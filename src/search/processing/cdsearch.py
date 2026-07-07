"""
CD-Search processing module.

This module contains the core RPS-BLAST/CD-Search logic:

- Run RPS-BLAST for each protein sequence
- Parse BLAST XML output
- Export detailed hit tables
- Export query-to-PSSM candidate tables
- Export default top1 PSSM proposal tables
- Extract selected domain alignment blocks
- Extract selected domain FASTA files

The CD-Search stage itself does not force the final top1 domain selection.
Instead, it writes a default top1 proposal table that can be inspected,
accepted, or manually edited before running the extraction stage.
"""

from __future__ import annotations

# Standard Library Imports
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

# Third-Party Library Imports
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_cdsearch(
    *,
    query_fasta: Path,
    cdd_rpsblast_db: Path,
    output_dir: Path,
    evalue: float,
    outfmt: int,
    num_threads: int,
    resume: bool,
) -> None:
    """
    Run local CD-Search using RPS-BLAST against a local CDD database.

    This stage performs sequence-level RPS-BLAST and writes selection
    artifacts only. It does not directly extract domain FASTA files.

    Parameters
    ----------
    query_fasta : pathlib.Path
        Input protein FASTA file.

    cdd_rpsblast_db : pathlib.Path
        Local CDD RPS-BLAST database basename.

        Example:
            data/external/cdd/rpsblast/Cdd

    output_dir : pathlib.Path
        CD-Search output root directory.

    evalue : float
        RPS-BLAST E-value threshold.

    outfmt : int
        RPS-BLAST output format. This pipeline expects XML format 5.

    num_threads : int
        Number of RPS-BLAST threads.

    resume : bool
        If True, skip sequence-level XML files that already exist.

    Returns
    -------
    None
    """
    start_time = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║                 [ CD-Search Alignment Stage Started ]                ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    print(f"📂 Input FASTA      : {query_fasta}")
    print(f"📁 Output Directory : {output_dir}")
    print(f"🧬 CDD Database     : {cdd_rpsblast_db}\n")

    output_dir.mkdir(parents=True, exist_ok=True)

    seq_records = load_query_records(query_fasta)
    seq_list = list(seq_records.items())
    total = len(seq_list)

    if total == 0:
        raise ValueError(f"No FASTA records found in input file: {query_fasta}")

    print(f"🔬 Total sequences to process: {total}\n")

    all_results: List[Dict[str, Any]] = []
    metadata: List[Dict[str, str]] = []

    success = 0
    skipped = 0
    failed = 0
    nohit = 0

    for idx, (seq_id, seq_record) in enumerate(seq_list, start=1):
        progress = (idx / total) * 100

        print("─" * 70)
        print(f"▶️  [{idx:3d} / {total:<3d} | {progress:5.1f}% ]  {seq_id}")
        print("─" * 70)

        xml_file, status, note = run_rpsblast_for_record(
            seq_record=seq_record,
            cdd_rpsblast_db=cdd_rpsblast_db,
            output_dir=output_dir,
            evalue=evalue,
            outfmt=outfmt,
            num_threads=num_threads,
            resume=resume,
        )

        if status in {"success", "skipped"} and xml_file:
            parsed_hits = parse_rpsblast_xml(
                xml_path=xml_file,
                seq_id=seq_id,
            )

            if parsed_hits:
                all_results.extend(parsed_hits)

                if status == "success":
                    success += 1
                    print(f"   ✅ Alignment success for {seq_id}\n")
                else:
                    skipped += 1
                    print("   ⏩ Skipped existing XML and parsed results\n")
            else:
                nohit += 1
                print("   ⚠️  No domains parsed from XML\n")

        elif status == "no_hit":
            nohit += 1
            print("   ⚠️  No conserved domain detected\n")

        else:
            failed += 1
            print("   ❌ RPS-BLAST error encountered\n")

        metadata.append(
            {
                "seq_id": seq_id,
                "alignment": "yes" if status in {"success", "skipped"} else "no",
                "note": note,
            }
        )

    output_paths = prepare_cdsearch_output_paths(output_dir)

    if all_results:
        write_cdsearch_result_tables(
            all_results=all_results,
            output_paths=output_paths,
        )
    else:
        print("⚠️  No CD-Search hits were collected. Result tables were not written.\n")

    pd.DataFrame(metadata).to_csv(
        output_paths["metadata"],
        sep="\t",
        index=False,
    )

    elapsed = time.time() - start_time

    print_cdsearch_summary(
        output_paths=output_paths,
        total=total,
        success=success,
        skipped=skipped,
        nohit=nohit,
        failed=failed,
        elapsed=elapsed,
    )


def load_query_records(query_fasta: Path) -> Dict[str, SeqRecord]:
    """
    Load query FASTA records and normalize sequence IDs.

    The original implementation used seq_record.id.replace("|", "_")
    when writing XML files. This helper applies the same normalization at
    FASTA loading time so that query_id values remain consistent across:

    - XML file names
    - parsed CD-Search tables
    - downstream domain extraction
    """
    seq_records: Dict[str, SeqRecord] = {}

    for rec in SeqIO.parse(str(query_fasta), "fasta"):
        safe_id = rec.id.replace("|", "_")

        rec.id = safe_id
        rec.name = safe_id
        rec.description = safe_id

        seq_records[safe_id] = rec

    return seq_records


def run_rpsblast_for_record(
    *,
    seq_record: SeqRecord,
    cdd_rpsblast_db: Path,
    output_dir: Path,
    evalue: float,
    outfmt: int,
    num_threads: int,
    resume: bool,
) -> Tuple[Path | None, str, str]:
    """
    Run RPS-BLAST for a single protein sequence.

    This function preserves the original behavior:

    - Create a temporary FASTA file for each sequence
    - Write XML output to output_dir/intermediate/<seq_id>.xml
    - Remove the temporary FASTA file after execution
    - Keep XML files for resume/reuse
    """
    seq_id = seq_record.id

    intermediate_dir = output_dir / "intermediate"
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    fasta_tmp = intermediate_dir / f"{seq_id}.fasta"
    xml_out = intermediate_dir / f"{seq_id}.xml"

    if resume and xml_out.exists():
        print(f"⏩ Skip {seq_id} (already completed)")
        return xml_out, "skipped", ""

    SeqIO.write(seq_record, str(fasta_tmp), "fasta")

    cmd = [
        "rpsblast",
        "-query",
        str(fasta_tmp),
        "-db",
        str(cdd_rpsblast_db),
        "-outfmt",
        str(outfmt),
        "-evalue",
        str(evalue),
        "-num_threads",
        str(num_threads),
        "-out",
        str(xml_out),
    ]

    try:
        subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )

        if not xml_out.exists() or xml_out.stat().st_size == 0:
            if xml_out.exists():
                xml_out.unlink()
            return None, "no_hit", "no domain alignment"

        return xml_out, "success", ""

    except subprocess.CalledProcessError as exc:
        error_log = output_dir / "cdsearch_error.log"
        with error_log.open("a", encoding="utf-8") as log:
            log.write(f"{seq_id}\t{exc}\n")
            if exc.stderr:
                log.write(f"{exc.stderr}\n")

        return None, "error", str(exc)

    finally:
        if fasta_tmp.exists():
            try:
                fasta_tmp.unlink()
            except OSError:
                pass


def parse_rpsblast_xml(
    *,
    xml_path: Path,
    seq_id: str,
) -> List[Dict[str, Any]]:
    """
    Parse RPS-BLAST XML output and extract detailed alignment information.

    Output columns are kept compatible with the original implementation.
    """
    results: List[Dict[str, Any]] = []

    try:
        with xml_path.open("r", encoding="utf-8") as handle:
            blast_records = NCBIXML.parse(handle)

            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        qseq = hsp.query
                        hseq = hsp.sbjct
                        midline = hsp.match

                        alignment_map = "".join(
                            "M" if char == "|" else "-" for char in midline
                        )

                        results.append(
                            {
                                "query_id": seq_id,
                                "PSSM_ID": alignment.hit_id.replace("gnl|CDD|", ""),
                                "title": alignment.hit_def,
                                "bitscore": hsp.bits,
                                "evalue": hsp.expect,
                                "pident": (
                                    sum(c1 == c2 for c1, c2 in zip(qseq, hseq))
                                    / len(qseq)
                                )
                                * 100,
                                "qstart": hsp.query_start,
                                "qend": hsp.query_end,
                                "sstart": hsp.sbjct_start,
                                "send": hsp.sbjct_end,
                                "qseq": qseq,
                                "hseq": hseq,
                                "midline": midline,
                                "alignment_map": alignment_map,
                            }
                        )

        return results

    except Exception as exc:  # pylint: disable=broad-exception-caught
        print(f"⚠️ Error parsing XML for {seq_id}: {exc}")
        return []


def prepare_cdsearch_output_paths(output_dir: Path) -> Dict[str, Path]:
    """
    Prepare and return all CD-Search output paths.

    The CD-Search stage writes search and selection artifacts only.
    Alignment blocks and domain FASTA files are generated by the
    cdsearch_extract stage after the selected PSSM table is confirmed.

    The default top1 table is intentionally detailed. It contains enough
    alignment-level information to uniquely identify the selected hit during
    downstream extraction.
    """
    return {
        "detailed": output_dir / "cdsearch_all_hits_detailed.tsv",
        "query_pssm_candidates": output_dir / "cdsearch_query_pssm_candidates.tsv",
        "default_top1": output_dir / "cdsearch_default_top1.tsv",
        "metadata": output_dir / "cdsearch_metadata.tsv",
    }


def write_cdsearch_result_tables(
    *,
    all_results: List[Dict[str, Any]],
    output_paths: Dict[str, Path],
) -> None:
    """
    Write CD-Search result tables.

    This function writes:

    - cdsearch_all_hits_detailed.tsv
    - cdsearch_query_pssm_candidates.tsv
    - cdsearch_default_top1.tsv

    The default top1 table is detailed rather than minimal. It contains the
    selected PSSM ID together with alignment coordinates and scores so that
    downstream extraction can uniquely recover the intended hit.
    """
    df = pd.DataFrame(all_results)

    df.to_csv(
        output_paths["detailed"],
        sep="\t",
        index=False,
    )

    candidates = build_query_pssm_candidates(df)
    candidates.to_csv(
        output_paths["query_pssm_candidates"],
        sep="\t",
        index=False,
    )

    default_top1 = select_default_top1_detailed(df)
    default_top1.to_csv(
        output_paths["default_top1"],
        sep="\t",
        index=False,
    )


def build_query_pssm_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a human-readable query-to-PSSM candidate table.

    This table contains all detected PSSM hits for each query protein,
    ranked using the default rule:

    - highest bitscore first
    - lowest evalue second

    Users can inspect this table when deciding whether the default top1
    PSSM assignment is biologically appropriate.
    """
    candidate_cols = [
        "query_id",
        "PSSM_ID",
        "rank",
        "bitscore",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "sstart",
        "send",
        "title",
    ]

    ranked = df.sort_values(
        by=["query_id", "bitscore", "evalue"],
        ascending=[True, False, True],
    ).copy()

    ranked["rank"] = ranked.groupby("query_id").cumcount() + 1

    return ranked[candidate_cols]  # type: ignore


def select_default_top1_detailed(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select the default top1 PSSM hit for each query.

    The output table is detailed and is intended to be used directly by
    the cdsearch_extract stage. Each selected hit keeps alignment-level
    information so that downstream extraction can uniquely recover the
    intended record from the all-hits table.

    Default sorting rule:
    - highest bitscore first
    - lowest evalue second

    The rank column follows the same ranking rule used by
    cdsearch_query_pssm_candidates.tsv.
    """
    ranked = df.sort_values(
        by=["query_id", "bitscore", "evalue"],
        ascending=[True, False, True],
    ).copy()

    ranked["rank"] = ranked.groupby("query_id").cumcount() + 1

    top1 = ranked[ranked["rank"] == 1].copy()

    output_cols = [
        "query_id",
        "PSSM_ID",
        "rank",
        "bitscore",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "sstart",
        "send",
        "title",
    ]

    return top1[output_cols].copy()  # type: ignore


def extract_selected_cdsearch_domains(
    *,
    query_fasta: Path,
    cdsearch_all_hits: Path,
    selected_top1: Path,
    output_dir: Path,
) -> None:
    """
    Extract selected CD-Search domain alignments and FASTA files.

    Parameters
    ----------
    query_fasta : pathlib.Path
        Input protein FASTA file used in the CD-Search stage.

    cdsearch_all_hits : pathlib.Path
        Full detailed hit table generated by the CD-Search stage.

    selected_top1 : pathlib.Path
        Selected top1 table. Must contain:

        - query_id
        - selected_PSSM_ID

        This can be the default proposal generated by CD-Search, or a
        manually edited table.

    output_dir : pathlib.Path
        Output directory for extracted domain artifacts.

    Returns
    -------
    None
    """
    start_time = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║              [ CD-Search Domain Extraction Started ]                 ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    print(f"📂 Query FASTA      : {query_fasta}")
    print(f"📄 All hits table   : {cdsearch_all_hits}")
    print(f"⭐ Selected top1     : {selected_top1}")
    print(f"📁 Output Directory : {output_dir}\n")

    output_dir.mkdir(parents=True, exist_ok=True)

    seq_records = load_query_records(query_fasta)

    if not seq_records:
        raise ValueError(f"No FASTA records found in input file: {query_fasta}")

    all_hits_df = pd.read_csv(cdsearch_all_hits, sep="\t")
    selected_df = pd.read_csv(selected_top1, sep="\t")

    selected_hits = build_selected_domain_hits(
        all_hits_df=all_hits_df,
        selected_df=selected_df,
        output_dir=output_dir,
    )

    output_paths = prepare_cdsearch_extract_output_paths(output_dir)

    selected_hits.to_csv(
        output_paths["selected_hits"],
        sep="\t",
        index=False,
    )

    write_alignment_blocks(
        top_hits=selected_hits,
        output_dir=output_paths["alignment_blocks"],
    )

    write_domain_fastas(
        top_hits=selected_hits,
        seq_records=seq_records,
        query_dir=output_paths["query_fasta"],
        query_with_gap_dir=output_paths["query_with_gap_fasta"],
        hseq_with_gap_dir=output_paths["hseq_with_gap_fasta"],
    )

    elapsed = time.time() - start_time

    print_cdsearch_extract_summary(
        output_paths=output_paths,
        total_selected=len(selected_df),
        total_extracted=len(selected_hits),
        elapsed=elapsed,
    )


def build_selected_domain_hits(
    *,
    all_hits_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    """
    Build a detailed selected-domain table from all hits and selected top1.

    The selected top1 table must contain detailed hit-level columns:

    - query_id
    - PSSM_ID
    - bitscore
    - evalue
    - pident
    - qstart
    - qend
    - sstart
    - send

    The rank column may be present in selected_df, but it is not used for
    matching because cdsearch_all_hits_detailed.tsv does not contain rank.

    This detailed matching avoids ambiguous extraction when RPS-BLAST reports
    multiple HSPs for the same query_id / PSSM_ID pair.
    """
    required_selected_cols = {
        "query_id",
        "PSSM_ID",
        "bitscore",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "sstart",
        "send",
    }

    required_hit_cols = {
        "query_id",
        "PSSM_ID",
        "bitscore",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "sstart",
        "send",
    }

    missing_selected_cols = required_selected_cols - set(selected_df.columns)
    if missing_selected_cols:
        raise ValueError(
            "Selected top1 table is missing required detailed columns: "
            f"{sorted(missing_selected_cols)}. "
            "Use results/cdsearch/cdsearch_default_top1.tsv generated by the "
            "updated cdsearch stage, or provide a manually edited selected table "
            "with the same detailed columns."
        )

    missing_hit_cols = required_hit_cols - set(all_hits_df.columns)
    if missing_hit_cols:
        raise ValueError(
            "CD-Search all-hits table is missing required columns: "
            f"{sorted(missing_hit_cols)}"
        )

    selected_df = selected_df.copy()
    all_hits_df = all_hits_df.copy()

    # Normalize string columns
    selected_df["query_id"] = selected_df["query_id"].astype(str)
    selected_df["PSSM_ID"] = selected_df["PSSM_ID"].astype(str)

    all_hits_df["query_id"] = all_hits_df["query_id"].astype(str)
    all_hits_df["PSSM_ID"] = all_hits_df["PSSM_ID"].astype(str)

    # Normalize numeric columns
    numeric_cols = [
        "bitscore",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "sstart",
        "send",
    ]

    for col in numeric_cols:
        selected_df[col] = pd.to_numeric(selected_df[col], errors="coerce")
        all_hits_df[col] = pd.to_numeric(all_hits_df[col], errors="coerce")

    # Integer coordinate columns should match exactly.
    coord_cols = [
        "qstart",
        "qend",
        "sstart",
        "send",
    ]

    for col in coord_cols:
        selected_df[col] = selected_df[col].astype("Int64")
        all_hits_df[col] = all_hits_df[col].astype("Int64")

    # Prevent multiple selected rows for the same query.
    duplicated_selected = selected_df[
        selected_df.duplicated(subset=["query_id"], keep=False)
    ]

    if not duplicated_selected.empty:
        duplicated_path = output_dir / "duplicated_selected_top1.tsv"
        duplicated_selected.to_csv(
            duplicated_path,
            sep="\t",
            index=False,
        )

        raise ValueError(
            "Selected top1 table contains duplicated query_id values. "
            f"See: {duplicated_path}"
        )

    merge_keys = [
        "query_id",
        "PSSM_ID",
        "bitscore",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "sstart",
        "send",
    ]

    selected_key_df = selected_df[merge_keys].copy()

    selected_hits = all_hits_df.merge(
        selected_key_df,
        on=merge_keys,
        how="inner",
    )

    # Check whether every selected row was found in all_hits_df.
    matched_keys = selected_hits[merge_keys].drop_duplicates()

    missing_selection = selected_key_df.merge(
        matched_keys,
        on=merge_keys,
        how="left",
        indicator=True,
    )

    missing_selection = missing_selection[missing_selection["_merge"] == "left_only"]

    if not missing_selection.empty:
        missing_path = output_dir / "missing_selected_top1.tsv"
        missing_selection.drop(columns=["_merge"]).to_csv(
            missing_path,
            sep="\t",
            index=False,
        )

        raise ValueError(
            "Some selected detailed hits were not found in the all-hits table. "
            "This usually means selected_top1 was edited inconsistently or was "
            "generated from a different cdsearch_all_hits_detailed.tsv. "
            f"See: {missing_path}"
        )

    # Check whether detailed matching still produced duplicated query hits.
    duplicated_selected_hits = selected_hits[
        selected_hits.duplicated(subset=["query_id"], keep=False)
    ]

    if not duplicated_selected_hits.empty:
        duplicated_path = output_dir / "duplicated_selected_domain_hits.tsv"
        duplicated_selected_hits.to_csv(
            duplicated_path,
            sep="\t",
            index=False,
        )

        raise ValueError(
            "Selected domain extraction produced duplicated query_id values "
            "even after detailed matching. "
            f"See: {duplicated_path}"
        )

    return selected_hits


def prepare_cdsearch_extract_output_paths(output_dir: Path) -> Dict[str, Path]:
    """
    Prepare output paths for selected CD-Search domain extraction.
    """
    align_dir = output_dir / "alignment_blocks"
    fasta_dir = output_dir / "domains_fasta"
    query_dir = fasta_dir / "query"
    query_with_gap_dir = fasta_dir / "query_with_gap"
    hseq_with_gap_dir = fasta_dir / "hseq_with_gap"

    align_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)
    query_dir.mkdir(parents=True, exist_ok=True)
    query_with_gap_dir.mkdir(parents=True, exist_ok=True)
    hseq_with_gap_dir.mkdir(parents=True, exist_ok=True)

    return {
        "selected_hits": output_dir / "selected_domain_hits_detailed.tsv",
        "alignment_blocks": align_dir,
        "domains_fasta": fasta_dir,
        "query_fasta": query_dir,
        "query_with_gap_fasta": query_with_gap_dir,
        "hseq_with_gap_fasta": hseq_with_gap_dir,
    }


def write_alignment_blocks(
    *,
    top_hits: pd.DataFrame,
    output_dir: Path,
) -> None:
    """
    Write plain-text alignment blocks for selected domain hits.
    """
    for _, row in top_hits.iterrows():
        qid = row["query_id"]
        pssm_id = row["PSSM_ID"]
        qseq = row["qseq"]
        hseq = row["hseq"]
        midline = row["midline"]
        qstart = int(row["qstart"])
        qend = int(row["qend"])
        sstart = int(row["sstart"])
        send = int(row["send"])

        align_txt = output_dir / f"{qid}_{pssm_id}.txt"

        with align_txt.open("w", encoding="utf-8") as handle:
            handle.write(f">{qid} vs PSSM{pssm_id}\n")
            handle.write(f"Query  {qstart:4d}  {qseq}  {qend}\n")
            handle.write(f"        {midline}\n")
            handle.write(f"Sbjct  {sstart:4d}  {hseq}  {send}\n")


def write_domain_fastas(
    *,
    top_hits: pd.DataFrame,
    seq_records: Dict[str, SeqRecord],
    query_dir: Path,
    query_with_gap_dir: Path,
    hseq_with_gap_dir: Path,
) -> None:
    """
    Write domain FASTA files from selected domain hits.

    Output folders and file naming are kept compatible with the original
    implementation:

    - domains_fasta/query/
    - domains_fasta/query_with_gap/
    - domains_fasta/hseq_with_gap/
    """
    for _, row in top_hits.iterrows():
        qid = row["query_id"]
        pssm_id = row["PSSM_ID"]
        qseq = row["qseq"]
        hseq = row["hseq"]
        qstart = int(row["qstart"])
        qend = int(row["qend"])

        if qid not in seq_records:
            raise KeyError(
                f"Query ID '{qid}' was found in selected hits but not in FASTA records."
            )

        query_fragment = seq_records[qid].seq[qstart - 1 : qend]

        SeqIO.write(
            SeqRecord(
                query_fragment,
                id=f"{qid}_{pssm_id}_query_fragment",
                description=f"domain_fragment {qstart}-{qend}",
            ),
            str(query_dir / f"{qid}_{pssm_id}_query_fragment.fasta"),
            "fasta",
        )

        SeqIO.write(
            SeqRecord(
                seq=Seq(qseq),
                id=f"{qid}_{pssm_id}_query_with_gap",
                description=f"aligned_query_with_gap {qstart}-{qend}",
            ),
            str(query_with_gap_dir / f"{qid}_{pssm_id}_query_with_gap.fasta"),
            "fasta",
        )

        SeqIO.write(
            SeqRecord(
                seq=Seq(hseq),
                id=f"{qid}_{pssm_id}_hseq_with_gap",
                description=(
                    f"aligned_hseq_with_gap "
                    f"sstart={row['sstart']} send={row['send']}"
                ),
            ),
            str(hseq_with_gap_dir / f"{qid}_{pssm_id}_hseq_with_gap.fasta"),
            "fasta",
        )


def print_cdsearch_summary(
    *,
    output_paths: Dict[str, Path],
    total: int,
    success: int,
    skipped: int,
    nohit: int,
    failed: int,
    elapsed: float,
) -> None:
    """
    Print CD-Search summary in a format close to the original implementation.
    """
    print("╔" + "═" * 70 + "╗")
    print("║" + " " * 24 + "CD-Search Alignment Summary" + " " * 19 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("\n📄 Output Files")
    print("─" * 70)
    print(f"✅ Detailed hits          : {output_paths['detailed']}")
    print(f"📋 Query-PSSM candidates  : {output_paths['query_pssm_candidates']}")
    print(f"⭐ Default top1 table     : {output_paths['default_top1']}")
    print(f"🧾 Metadata summary       : {output_paths['metadata']}\n")

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total sequences : {total}")
    print(f"Successful hits : {success}")
    print(f"No hits         : {nohit}")
    print(f"Skipped         : {skipped}")
    print(f"Failed          : {failed}\n")

    print(f"⏱️  Elapsed time : {elapsed:.2f} seconds")
    print("✅  CD-Search stage completed successfully!\n")

    print("Next step")
    print("─" * 70)
    print(
        "Inspect the detailed default top1 proposal or manually edit a selected "
        "top1 table, then run the cdsearch_extract stage to generate alignment "
        "blocks and domain FASTA files.\n"
    )


def print_cdsearch_extract_summary(
    *,
    output_paths: Dict[str, Path],
    total_selected: int,
    total_extracted: int,
    elapsed: float,
) -> None:
    """
    Print CD-Search extraction summary.
    """
    print("╔" + "═" * 70 + "╗")
    print("║" + " " * 20 + "CD-Search Extraction Summary" + " " * 22 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("\n📄 Output Files")
    print("─" * 70)
    print(f"✅ Selected hits    : {output_paths['selected_hits']}")
    print(f"📁 Alignment blocks : {output_paths['alignment_blocks']}")
    print(f"📁 Domain FASTAs    : {output_paths['domains_fasta']}\n")

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Selected records : {total_selected}")
    print(f"Extracted hits   : {total_extracted}\n")

    print(f"⏱️  Elapsed time : {elapsed:.2f} seconds")
    print("✅  CD-Search domain extraction completed successfully!\n")
