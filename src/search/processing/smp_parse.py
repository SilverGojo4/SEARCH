"""
CDD .smp PSSM extraction processing module.

This module implements Branch B of the original PSSM workflow:

- Read CD-Search top-hit table
- Resolve corresponding CDD .smp files from CD-Search titles
- Parse CDD .smp ASN.1 profiles
- Extract domain-level PSSM matrices
- Export matrix TSV files, metadata summary, and error log

The output structure and file names are designed to remain compatible
with the original PSSM/src/preprocess/run_smp_parse.py implementation.
"""

from __future__ import annotations

# Standard Library Imports
import time
from pathlib import Path
from typing import List, Optional

# Third-Party Imports
import numpy as np
import pandas as pd

# App Imports
from search.processing.smp_parser import ParsedCDDProfile, parse_smp_file

AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")

HYDROPHOBIC = {"A", "I", "L", "V", "M", "F", "W", "P", "C"}
CHARGED = {"H", "D", "E", "K", "R"}
POLAR = {"S", "T", "Y", "N", "Q"}


def resolve_smp_filename_from_title(title: str) -> str:
    """
    Resolve .smp filename from CD-Search title field.

    Example
    -------
    title:
        "NF041153, organophos_OPH, organophosphate hydrolase OPH..."
    resolved:
        "NF041153.smp"
    """
    if not title or not isinstance(title, str):
        raise ValueError("Empty or invalid title field (cannot resolve .smp filename)")

    base = title.split(",")[0].strip()

    if not base:
        raise ValueError(f"Cannot resolve .smp filename from title: {title}")

    return f"{base}.smp"


def resolve_smp_path(smp_root_dir: Path | str, title: str) -> Path:
    """
    Resolve .smp file path from a title-derived filename.

    Common CDD extraction layouts include:
        <root>/<name>.smp
        <root>/smps/<name>.smp

    This preserves the original lookup behavior.
    """
    smp_root_dir = Path(smp_root_dir)
    smp_fname = resolve_smp_filename_from_title(title)

    path_direct = smp_root_dir / smp_fname
    path_nested = smp_root_dir / "smps" / smp_fname

    if path_direct.exists():
        return path_direct

    if path_nested.exists():
        return path_nested

    raise FileNotFoundError(f".smp file not found: {smp_fname} under: {smp_root_dir}")


def select_profile_block(
    profiles: List[ParsedCDDProfile],
    pssm_id: int,
) -> ParsedCDDProfile:
    """
    Select the correct PSSM block from a .smp file.

    Strategy
    --------
    - If only one block exists: use it
    - Otherwise:
        prefer block whose title contains "cd{pssm_id}"
        else fallback to the first block

    This preserves the original selection logic.
    """
    if len(profiles) == 1:
        return profiles[0]

    key = f"cd{pssm_id}"

    for profile in profiles:
        if profile.title and key in profile.title:
            return profile

    # Fallback is valid for many .smp files.
    return profiles[0]


def convert_pssm_to_dataframe(
    pssm_aa20: np.ndarray,
    query_seq: Optional[str] = None,
    aa_order: str = "GAILVMFWPCSTYNQHKRDE",
) -> pd.DataFrame:
    """
    Convert a 20 x L matrix into long-form TSV with AA columns and features.

    pssm_aa20 is assumed to be in parser AA_ORDER_20:

        ARNDCQEGHILKMFPSTWYV

    It is reordered into the global Branch B AA order:

        GAILVMFWPCSTYNQHKRDE

    Before calculating biochemical features, scaled PSSM values are rounded
    to integer scores. This makes Po, Hy, Ch, |Hy-Ch|, and |Po-Ch| consistent
    with integer-valued PSSM downstream processing.

    The biochemical feature calculation preserves the original positive-only rule.
    """
    aa_order_list = list(aa_order)

    parser_order = list("ARNDCQEGHILKMFPSTWYV")
    parser_index = {aa: idx for idx, aa in enumerate(parser_order)}

    missing = [aa for aa in aa_order_list if aa not in parser_index]

    if missing:
        raise ValueError(f"Missing AA in parser output: {missing}")

    length = pssm_aa20.shape[1]

    # Build reordered matrix with shape L x 20.
    reordered = np.zeros((length, 20), dtype=float)

    for col_idx, aa in enumerate(aa_order_list):
        reordered[:, col_idx] = pssm_aa20[parser_index[aa], :]

    # Round scaled PSSM values to integers before calculating derived features.
    #
    # We avoid np.round here because NumPy uses bankers rounding:
    #   np.round(2.5) -> 2.0
    #
    # The formula below performs conventional half-up rounding by magnitude:
    #   2.5  -> 3
    #   -2.5 -> -3
    #
    # NaN values are preserved if missing scores exist.
    # reordered = np.where(
    #     np.isnan(reordered),
    #     np.nan,
    #     np.sign(reordered) * np.floor(np.abs(reordered) + 0.5),
    # )

    df = pd.DataFrame(reordered, columns=aa_order_list)

    # Position and Residue columns.
    df.insert(0, "Position", list(range(1, length + 1)))

    if query_seq and len(query_seq) == length:
        df.insert(1, "Residue", list(query_seq))
    else:
        df.insert(1, "Residue", ["-"] * length)

    po_vals = []
    hy_vals = []
    ch_vals = []
    absch_vals = []
    abspoch_vals = []

    for _, row in df.iterrows():
        po = _sum_positive(row, POLAR, aa_order_list)
        hy = _sum_positive(row, HYDROPHOBIC, aa_order_list)
        ch = _sum_positive(row, CHARGED, aa_order_list)

        po_vals.append(po)
        hy_vals.append(hy)
        ch_vals.append(ch)
        absch_vals.append(abs(hy - ch))
        abspoch_vals.append(abs(po - ch))

    df["Po"] = po_vals
    df["Hy"] = hy_vals
    df["Ch"] = ch_vals
    df["|Hy-Ch|"] = absch_vals
    df["|Po-Ch|"] = abspoch_vals

    return df


def _sum_positive(
    row: pd.Series,
    aa_set: set[str],
    aa_order: List[str],
) -> float:
    """
    Sum only positive PSSM values for a given biochemical amino-acid set.

    This preserves the original Branch B feature calculation rule.
    """
    return sum(float(row[aa]) for aa in aa_order if aa in aa_set and row[aa] > 0)


def run_smp_parse(
    *,
    cdsearch_table: Path,
    smp_root_dir: Path,
    output_dir: Path,
    aa_order: str = "GAILVMFWPCSTYNQHKRDE",
    resume: bool = True,
) -> None:
    """
    Extract domain-level PSSM matrices directly from CDD .smp files.

    Parameters
    ----------
    cdsearch_table : pathlib.Path
        CD-Search top-hit table.

        Expected columns:
        - query_id
        - PSSM_ID
        - title

    smp_root_dir : pathlib.Path
        Root directory containing CDD .smp files.

    output_dir : pathlib.Path
        Output directory for domain-level matrix TSV files.

    aa_order : str, optional
        Amino-acid column order used in output matrices.

    resume : bool, optional
        If True, skip matrix files that already exist.

    Returns
    -------
    None
    """
    start_time = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║              [ CDD .smp PSSM Extraction Stage Started ]              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    output_dir.mkdir(parents=True, exist_ok=True)

    meta_path = output_dir / "smp_metadata.tsv"
    log_path = output_dir / "smp_parse_error.log"

    print(f"📄 CD-Search Table      : {cdsearch_table}")
    print(f"📂 SMP Root Dir         : {smp_root_dir}")
    print(f"📁 Output Matrix Dir    : {output_dir}")
    print(f"🧾 Metadata Summary     : {meta_path}")
    print(f"🐛 Error Log            : {log_path}\n")

    df_cd = pd.read_csv(cdsearch_table, sep="\t")

    required_cols = {"query_id", "PSSM_ID", "title"}

    if not required_cols.issubset(df_cd.columns):
        raise ValueError(f"CD-search table missing required columns: {required_cols}")

    total = len(df_cd)
    print(f"🔬 Total domain hits to process: {total}\n")

    summary_rows = []

    success = 0
    skipped = 0
    failed = 0

    for idx, row in df_cd.iterrows():
        query_id = row["query_id"]
        pssm_id = int(row["PSSM_ID"])
        title = str(row["title"])

        out_name = f"{query_id}_{pssm_id}_hseq_with_gap.tsv"
        out_path = output_dir / out_name

        progress = (idx + 1) / total * 100

        print("─" * 70)
        print(
            f"▶️  [{idx + 1:3d}/{total:<3d} | {progress:5.1f}% ] {query_id} (PSSM {pssm_id})"
        )
        print("─" * 70)

        if resume and out_path.exists():
            skipped += 1
            summary_rows.append(
                {
                    "query_id": query_id,
                    "PSSM_ID": pssm_id,
                    "status": "skipped",
                    "note": "already exists",
                }
            )
            print(f"   ⏩ Skip (already extracted): {out_name}\n")
            continue

        try:
            smp_path = resolve_smp_path(smp_root_dir, title)
            smp_fname = smp_path.name

            profiles = parse_smp_file(smp_path)
            selected = select_profile_block(profiles, pssm_id)

            matrix_df = convert_pssm_to_dataframe(
                selected.pssm_aa20,
                query_seq=selected.query_seq,
                aa_order=aa_order,
            )

            matrix_df.to_csv(out_path, sep="\t", index=False)

            success += 1
            summary_rows.append(
                {
                    "query_id": query_id,
                    "PSSM_ID": pssm_id,
                    "smp_file": smp_fname,
                    "status": "success",
                    "note": str(smp_path),
                }
            )

            print(f"   ✅ SMP parsed successfully -> {out_name}")
            print(f"   📄 SMP file used: {smp_fname}\n")

        except Exception as exc:  # pylint: disable=broad-exception-caught
            failed += 1
            summary_rows.append(
                {
                    "query_id": query_id,
                    "PSSM_ID": pssm_id,
                    "status": "error",
                    "note": str(exc),
                }
            )

            with log_path.open("a", encoding="utf-8") as log:
                log.write(f"[{query_id} | PSSM {pssm_id}] {exc}\n")

            print(f"   ❌ SMP parsing failed: {exc}\n")

    pd.DataFrame(summary_rows).to_csv(meta_path, sep="\t", index=False)

    elapsed = time.time() - start_time

    print("\n" + "╔" + "═" * 70 + "╗")
    print("║" + " " * 22 + "CDD SMP Extraction Summary" + " " * 22 + "║")
    print("╚" + "═" * 70 + "╝\n")

    print("📊 Summary Statistics")
    print("─" * 70)
    print(f"Total hits     : {total}")
    print(f"Successful     : {success}")
    print(f"Skipped        : {skipped}")
    print(f"Failed         : {failed}\n")

    print("📄 Output Files")
    print("─" * 70)
    print(f"📁 Matrix Dir         : {output_dir}")
    print(f"🧾 Metadata Summary   : {meta_path}")
    print(f"🐛 Error Log          : {log_path}")
    print(f"⏱️  Elapsed time      : {elapsed:.2f} seconds\n")

    print("✅  CDD .smp extraction stage completed successfully!\n")
