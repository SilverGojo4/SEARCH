"""
Scorecons conservation integration processing module.

This module implements the Scorecons integration stage:

- Parse external Scorecons output files
- Load reconstructed full-length PSSM tables
- Project alignment-level Scorecons conservation scores back to
  full-length protein residue positions
- Export integrated PSSM + Scorecons tables

The output structure and file names are designed to remain compatible
with the original PSSM/src/preprocess/run_scorecons_integrate.py implementation.
"""

from __future__ import annotations

# Standard Library Imports
import time
from pathlib import Path
from typing import Dict, List

# Third-Party Imports
import numpy as np
import pandas as pd
from Bio import SeqIO

DEFAULT_AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")
DEFAULT_PSSM_NUMERIC_COLS = DEFAULT_AA_ORDER + [
    "Po",
    "Hy",
    "Ch",
    "|Hy-Ch|",
    "|Hy-Po|",
]


def parse_scorecons_txt(path: Path) -> pd.DataFrame:
    """
    Parse Scorecons server output file.

    Expected format
    ---------------
    - Header lines, first 9 lines, are ignored.
    - Data lines contain:
        <score> # <two residues>

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:
        - aln_index
        - conservation
        - res_1
        - res_2
    """
    if not path.exists():
        raise FileNotFoundError(f"Scorecons file not found: {path}")

    records = []

    with path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    data_lines = lines[9:]

    aln_index = 1

    for line in data_lines:
        line = line.strip()

        if not line:
            continue

        try:
            left, right = line.split("#", 1)
            score = float(left.strip())
            residues = right.strip()

            records.append(
                {
                    "aln_index": aln_index,
                    "conservation": score,
                    "res_1": residues[0] if len(residues) > 0 else None,
                    "res_2": residues[1] if len(residues) > 1 else None,
                }
            )
            aln_index += 1

        except Exception:  # pylint: disable=broad-exception-caught
            continue

    if not records:
        raise ValueError(f"No valid Scorecons scores found: {path}")

    return pd.DataFrame(records)


def load_alignment_table(path: Path) -> pd.DataFrame:
    """
    Load selected CD-Search domain alignment table.
    """
    if not path.exists():
        raise FileNotFoundError(f"Selected CD-Search domain table not found: {path}")

    return pd.read_csv(path, sep="\t")


def load_original_sequences(fasta_path: Path) -> Dict[str, str]:
    """
    Load original sequences from FASTA.

    Returns
    -------
    dict
        Mapping from sequence ID to sequence string.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {fasta_path}")

    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}


def restore_numeric_types(
    df: pd.DataFrame,
    cols: List[str],
) -> pd.DataFrame:
    """
    Restore numeric dtypes.

    Rules
    -----
    - If a column contains only integer-like values, or NA, cast to Int64.
    - Otherwise, keep float.

    Notes
    -----
    Pandas can upcast integer columns to float when NA exists.
    This helper restores stable dtypes to avoid integer-to-float drift.
    """
    for col in cols:
        if col not in df.columns:
            continue

        series = pd.to_numeric(df[col], errors="coerce")

        if series.isna().all():
            continue

        non_na = series.dropna()

        if np.isclose(non_na, non_na.round()).all():
            df[col] = series.round().astype("Int64")
        else:
            df[col] = series.astype(float)

    return df


def load_pssm_reconstruct_table(
    *,
    reconstruct_dir: Path,
    query_id: str,
    numeric_cols: List[str],
) -> pd.DataFrame:
    """
    Load reconstructed full-length PSSM table.

    Expected path:
        <reconstruct_dir>/<query_id>.tsv
    """
    file_path = reconstruct_dir / f"{query_id}.tsv"

    if not file_path.exists():
        raise FileNotFoundError(f"PSSM reconstruct file not found: {file_path}")

    df = pd.read_csv(file_path, sep="\t")
    df = restore_numeric_types(df, numeric_cols)

    return df


def resolve_scorecons_path(
    *,
    scorecons_dir: Path,
    query_id: str,
) -> Path:
    """
    Resolve Scorecons txt file path.

    Expected naming:
        <scorecons_dir>/<query_id>.txt
    """
    file_path = scorecons_dir / f"{query_id}.txt"

    if not file_path.exists():
        raise FileNotFoundError(f"Scorecons file not found: {file_path}")

    return file_path


def integrate_scorecons(
    *,
    aln_row: pd.Series,
    scorecons_df: pd.DataFrame,
    base_table: pd.DataFrame,
) -> pd.DataFrame:
    """
    Integrate alignment-based Scorecons scores into a full-length table.

    Rules
    -----
    - Query gap positions are not projectable.
    - Only query residues receive scores.
    - Unaligned positions remain NA.

    Notes
    -----
    This uses qseq + qstart mapping from CD-Search and projects Scorecons
    aln_index back to full-length residue positions.
    """
    scorecons_df = scorecons_df.set_index("aln_index")

    qseq = str(aln_row["qseq"])
    orig_pos = int(aln_row["qstart"])
    aln_index = 1

    out = base_table.copy()

    title = str(aln_row["title"]).split(",")[0].strip()
    col_name = f"Conservation ({title})"

    if col_name in out.columns:
        out = out.drop(columns=[col_name])

    out.insert(3, col_name, pd.NA)

    for q_char in qseq:
        if q_char != "-":
            if aln_index in scorecons_df.index:
                score = scorecons_df.loc[aln_index, "conservation"]
                out.loc[out["Position"] == orig_pos, col_name] = score

            orig_pos += 1

        aln_index += 1

    return out


def run_scorecons_integrate(
    *,
    query_fasta: Path,
    cdsearch_table: Path,
    reconstruct_dir: Path,
    scorecons_dir: Path,
    output_dir: Path,
    aa_order: str = "GAILVMFWPCSTYNQHKRDE",
    resume: bool = True,
) -> None:
    """
    Integrate Scorecons conservation scores into reconstructed PSSM tables.

    Parameters
    ----------
    query_fasta : pathlib.Path
        Original full-length protein FASTA.
    cdsearch_table : pathlib.Path
        Selected CD-Search domain table generated by cdsearch_extract.
    reconstruct_dir : pathlib.Path
        Directory containing reconstructed full-length PSSM tables.
    scorecons_dir : pathlib.Path
        Directory containing external Scorecons output files.
    output_dir : pathlib.Path
        Directory to write integrated PSSM + Scorecons tables.
    aa_order : str, optional
        Amino-acid column order used by PSSM matrices.
    resume : bool, optional
        If True, skip integrated TSV files that already exist.

    Returns
    -------
    None
    """
    start_time = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║           [ Scorecons Integration Stage Started ]                   ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "scorecons_integrate_error.log"

    aa_order_list = list(aa_order)
    numeric_cols = aa_order_list + ["Po", "Hy", "Ch", "|Hy-Ch|", "|Hy-Po|"]

    print(f"📂 FASTA Path             : {query_fasta}")
    print(f"📄 Selected Domain Table  : {cdsearch_table}")
    print(f"📁 PSSM Reconstruct Dir   : {reconstruct_dir}")
    print(f"📁 Scorecons Dir          : {scorecons_dir}")
    print(f"📁 Integrated Output Dir  : {output_dir}")
    print(f"🐛 Error Log              : {log_path}\n")

    aln_table = load_alignment_table(cdsearch_table)
    original_seqs = load_original_sequences(query_fasta)

    total = len(aln_table)
    success = 0
    skipped = 0
    failed = 0

    for idx, row in aln_table.iterrows():
        query_id = row["query_id"]

        progress = (idx + 1) / total * 100

        print("─" * 70)
        print(f"▶️  [{idx + 1:3d}/{total:<3d} | {progress:5.1f}% ] {query_id}")
        print("─" * 70)

        out_path = output_dir / f"{query_id}.tsv"

        if resume and out_path.exists():
            skipped += 1
            print(f"   ⏩ Skip (already integrated): {out_path}\n")
            continue

        try:
            if query_id not in original_seqs:
                raise KeyError(f"Query sequence not found in FASTA: {query_id}")

            base_table = load_pssm_reconstruct_table(
                reconstruct_dir=reconstruct_dir,
                query_id=query_id,
                numeric_cols=numeric_cols,
            )
            scorecons_path = resolve_scorecons_path(
                scorecons_dir=scorecons_dir,
                query_id=query_id,
            )
            scorecons_df = parse_scorecons_txt(scorecons_path)

            merged = integrate_scorecons(
                aln_row=row,
                scorecons_df=scorecons_df,
                base_table=base_table,
            )

            merged = restore_numeric_types(merged, numeric_cols)
            merged.to_csv(out_path, sep="\t", index=False)

            print(f"   ✅ Scorecons integrated -> {out_path}\n")
            success += 1

        except Exception as exc:  # pylint: disable=broad-exception-caught
            failed += 1
            print(f"   ❌ Failed: {exc}\n")

            with log_path.open("a", encoding="utf-8") as log:
                log.write(f"[{query_id}] {exc}\n")

    elapsed = time.time() - start_time

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║                 Scorecons Integration Summary                       ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    print(f"Total proteins : {total}")
    print(f"Successful     : {success}")
    print(f"Skipped        : {skipped}")
    print(f"Failed         : {failed}")
    print(f"🐛 Error log   : {log_path}")
    print(f"⏱️ Elapsed time: {elapsed:.2f} seconds\n")

    print("✅  Scorecons integration completed!\n")
