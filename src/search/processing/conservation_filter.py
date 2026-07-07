"""
Conservation-based PSSM feature masking processing module.

This module implements the final conservation filtering stage:

- Load integrated full-length PSSM tables
- Detect the Scorecons conservation column
- Validate the ConSurf evolutionary conservation column
- Mask PSSM numeric feature columns when conservation criteria fail
- Write filtered tables to a separate output directory

The behavior is designed to remain compatible with the original
PSSM/src/postprocess/run_conservation_filter.py implementation.

Important
---------
All residue positions are retained.
Only PSSM profile features are masked as NA if conservation criteria
are not satisfied.
"""

from __future__ import annotations

# Standard Library Imports
import time
from pathlib import Path
from typing import List, Tuple

# Third-Party Imports
import numpy as np
import pandas as pd

DEFAULT_AA_ORDER = list("GAILVMFWPCSTYNQHKRDE")
DEFAULT_FEATURE_COLS = ["Po", "Hy", "Ch", "|Hy-Ch|", "|Po-Ch|"]
DEFAULT_PSSM_NUMERIC_COLS = DEFAULT_AA_ORDER + DEFAULT_FEATURE_COLS


def detect_scorecons_column(df: pd.DataFrame) -> str:
    """
    Detect the Scorecons conservation column.

    Expected format
    ---------------
    Conservation (XXXX)

    Returns
    -------
    str
        Detected Scorecons column name.

    Raises
    ------
    ValueError
        If no Scorecons column is found, or if multiple candidate columns exist.
    """
    matches = [
        col
        for col in df.columns
        if col.startswith("Conservation (") and col.endswith(")")
    ]

    if len(matches) == 0:
        raise ValueError("No Scorecons column found (expected: Conservation (xxxx))")

    if len(matches) > 1:
        raise ValueError(
            f"Multiple Scorecons columns found: {matches}. "
            "Ambiguous integration result."
        )

    return matches[0]


def validate_required_columns(
    df: pd.DataFrame,
    pssm_numeric_cols: List[str],
) -> Tuple[str, List[str]]:
    """
    Validate required columns exist.

    Required columns
    ----------------
    - Position
    - Evolutionary conservation
    - Conservation (xxxx)
    - PSSM numeric columns

    Parameters
    ----------
    df : pandas.DataFrame
        Integrated PSSM table.

    pssm_numeric_cols : list[str]
        PSSM feature columns that will be masked.

    Returns
    -------
    tuple
        cons_col, pssm_cols
    """
    if "Position" not in df.columns:
        raise ValueError("Missing required column: Position")

    if "Evolutionary conservation" not in df.columns:
        raise ValueError("Missing required column: Evolutionary conservation")

    cons_col = detect_scorecons_column(df)

    missing = [col for col in pssm_numeric_cols if col not in df.columns]

    if missing:
        raise ValueError(f"Missing required PSSM feature columns: {missing}")

    return cons_col, pssm_numeric_cols


def restore_numeric_types(
    df: pd.DataFrame,
    cols: List[str],
) -> pd.DataFrame:
    """
    Restore numeric dtypes.

    Rules
    -----
    - If integer-like, cast to nullable Int64.
    - Otherwise, keep float.

    Notes
    -----
    Pandas may auto-upcast integer columns into float when NA exists.
    This helper restores stable dtypes after masking.
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


def mask_pssm_features(
    *,
    df: pd.DataFrame,
    ec_min: float,
    ec_max: float,
    cons_min: float,
    cons_max: float,
    pssm_numeric_cols: List[str],
) -> pd.DataFrame:
    """
    Mask PSSM features where conservation thresholds are not satisfied.

    Rules
    -----
    - All residue rows are retained.
    - Rows passing both Scorecons and ConSurf thresholds keep PSSM values.
    - Rows failing either threshold have PSSM numeric columns set to NA.
    """
    cons_col, pssm_cols = validate_required_columns(
        df=df,
        pssm_numeric_cols=pssm_numeric_cols,
    )

    evo = pd.to_numeric(df["Evolutionary conservation"], errors="coerce")
    cons = pd.to_numeric(df[cons_col], errors="coerce")

    valid_mask = (
        evo.notna()
        & cons.notna()
        & (evo >= ec_min)
        & (evo <= ec_max)
        & (cons >= cons_min)
        & (cons <= cons_max)
    )

    out = df.copy()
    out.loc[~valid_mask, pssm_cols] = pd.NA
    out = restore_numeric_types(out, pssm_cols)

    return out


def run_conservation_filter(
    *,
    integrated_dir: Path,
    filtered_output_dir: Path,
    ec_min: float,
    ec_max: float,
    cons_min: float,
    cons_max: float,
    aa_order: str = "GAILVMFWPCSTYNQHKRDE",
    resume: bool = True,
) -> None:
    """
    Mask PSSM profile features by conservation thresholds.

    Parameters
    ----------
    integrated_dir : pathlib.Path
        Directory containing integrated PSSM + conservation TSV files.

    filtered_output_dir : pathlib.Path
        Directory to write filtered / masked TSV files.

    ec_min : float
        Lower bound for ConSurf evolutionary conservation.

    ec_max : float
        Upper bound for ConSurf evolutionary conservation.

    cons_min : float
        Lower bound for Scorecons conservation.

    cons_max : float
        Upper bound for Scorecons conservation.

    aa_order : str, optional
        Amino-acid column order used by PSSM matrices.

    resume : bool, optional
        If True, skip filtered TSV files that already exist.

    Returns
    -------
    None
    """
    start_time = time.time()

    aa_order_list = list(aa_order)
    feature_cols = ["Po", "Hy", "Ch", "|Hy-Ch|", "|Po-Ch|"]
    pssm_numeric_cols = aa_order_list + feature_cols

    filtered_output_dir.mkdir(parents=True, exist_ok=True)

    log_path = filtered_output_dir / "conservation_filter_error.log"

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║     [ Conservation-Based PSSM Feature Masking Started ]              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    print(f"📥 Integrated Dir : {integrated_dir}")
    print(f"📤 Filtered Dir   : {filtered_output_dir}")
    print(f"🐛 Error Log      : {log_path}")
    print(f"📌 EC range       : {ec_min} ~ {ec_max}")
    print(f"📌 CONS range     : {cons_min} ~ {cons_max}\n")

    files = sorted(integrated_dir.glob("*.tsv"))
    files = [path for path in files if not path.name.endswith("_metadata.tsv")]

    total = len(files)

    if total == 0:
        print(f"⚠️ No TSV files found in: {integrated_dir}\n")
        return

    success = 0
    skipped = 0
    failed = 0

    for idx, path in enumerate(files):
        query_id = path.stem
        progress = (idx + 1) / total * 100
        out_path = filtered_output_dir / f"{query_id}.tsv"

        print("─" * 70)
        print(f"▶️  [{idx + 1:3d}/{total:<3d} | {progress:5.1f}% ] {query_id}")
        print("─" * 70)

        if resume and out_path.exists():
            skipped += 1
            print(f"   ⏩ Skip (already filtered): {out_path}\n")
            continue

        try:
            df = pd.read_csv(path, sep="\t")

            masked = mask_pssm_features(
                df=df,
                ec_min=ec_min,
                ec_max=ec_max,
                cons_min=cons_min,
                cons_max=cons_max,
                pssm_numeric_cols=pssm_numeric_cols,
            )

            masked.to_csv(out_path, sep="\t", index=False)

            print(f"   ✅ Filtered -> {out_path}\n")
            success += 1

        except Exception as exc:  # pylint: disable=broad-exception-caught
            failed += 1
            print(f"   ❌ Failed: {exc}\n")

            with log_path.open("a", encoding="utf-8") as log:
                log.write(f"[{query_id}] {exc}\n")

    elapsed = time.time() - start_time

    print("\n══════════════════════════════════════════════════════════════════════")
    print("Conservation Filtering Summary")
    print("──────────────────────────────────────────────────────────────────────")
    print(f"Total       : {total}")
    print(f"Successful  : {success}")
    print(f"Skipped     : {skipped}")
    print(f"Failed      : {failed}")
    print(f"🐛 Error log: {log_path}")
    print(f"Elapsed     : {elapsed:.2f} sec")
    print("══════════════════════════════════════════════════════════════════════\n")
