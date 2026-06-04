"""
ConSurf evolutionary conservation integration processing module.

This module implements the ConSurf integration stage:

- Load integrated full-length PSSM + Scorecons tables
- Parse external ConSurf grades files
- Build ConSurf sequence from parsed grades
- Map ConSurf sequence back to the original full-length query sequence
  using local window matching
- Add / overwrite the "Evolutionary conservation" column
- Update integrated TSV files in place

The behavior is designed to remain compatible with the original
PSSM/src/preprocess/run_consurf_integrate.py implementation.
"""

from __future__ import annotations

# Standard Library Imports
import time
from glob import glob
from pathlib import Path
from typing import Dict, List, Tuple

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


def load_base_integrated_table(
    path: Path,
    numeric_cols: List[str],
) -> pd.DataFrame:
    """
    Load integrated TSV table.

    Important
    ---------
    We restore numeric dtypes after loading. Otherwise integer PSSM matrices
    may be auto-cast to float when NA exists.
    """
    if not path.exists():
        raise FileNotFoundError(f"Integrated table not found: {path}")

    df = pd.read_csv(path, sep="\t")
    df = restore_numeric_types(df, numeric_cols)

    return df


def parse_consurf_grades(path: Path) -> pd.DataFrame:
    """
    Parse ConSurf grades file.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:
        - pos
        - aa
        - score
    """
    records = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 3:
                continue

            try:
                pos = int(parts[0])
                aa = parts[1]
            except Exception:  # pylint: disable=broad-exception-caught
                continue

            score = None

            for token in parts[2:]:
                try:
                    score = float(token)
                    break
                except Exception:  # pylint: disable=broad-exception-caught
                    continue

            if score is None:
                continue

            records.append(
                {
                    "pos": pos,
                    "aa": aa,
                    "score": score,
                }
            )

    if not records:
        raise ValueError(f"No valid ConSurf scores found: {path}")

    return pd.DataFrame(records)


def build_consurf_sequence(consurf_df: pd.DataFrame) -> str:
    """
    Build amino-acid sequence string from parsed ConSurf grades.
    """
    return "".join(consurf_df.sort_values("pos")["aa"].tolist())


def parse_mutation(query_id: str) -> Tuple[int | None, str | None, str | None]:
    """
    Parse mutation information from query_id if present.

    Expected format
    ---------------
    UniProt_A123B

    Returns
    -------
    tuple
        (mutation_position, from_aa, to_aa)
    """
    try:
        mutation = query_id.split("_")[1]
        from_aa = mutation[0]
        to_aa = mutation[-1]
        pos = int(mutation[1:-1])

        return pos, from_aa, to_aa

    except Exception:  # pylint: disable=broad-exception-caught
        return None, None, None


def restore_reference_sequence(
    seq: str,
    mut_pos: int | None,
    mut_from: str | None,
) -> str:
    """
    Restore original residue in a mutated sequence if mutation exists.

    ConSurf grades are often computed on the reference sequence.
    """
    if mut_pos is None:
        return seq

    if mut_from is None:
        return seq

    idx = mut_pos - 1

    if 0 <= idx < len(seq):
        return seq[:idx] + mut_from + seq[idx + 1 :]

    return seq


def find_best_local_mapping(
    ref_seq: str,
    consurf_seq: str,
    window: int = 25,
    max_mismatch: int = 5,
) -> Tuple[int, int, int]:
    """
    Find best local matching region between two sequences.

    Parameters
    ----------
    ref_seq : str
        Reference full-length sequence.

    consurf_seq : str
        Sequence reconstructed from ConSurf grades.

    window : int
        Local window size used for sequence matching.

    max_mismatch : int
        Maximum allowed mismatch count within the local window.

    Returns
    -------
    tuple
        ref_start, consurf_start, mismatch

        ref_start and consurf_start are 0-based.
    """
    best = None

    def scan(seq_a: str, seq_b: str, label: str) -> None:
        nonlocal best

        for i in range(len(seq_a) - window + 1):
            anchor = seq_a[i : i + window]

            for j in range(len(seq_b) - window + 1):
                compare = seq_b[j : j + window]
                mismatch = sum(a != b for a, b in zip(anchor, compare))

                if mismatch <= max_mismatch:
                    candidate = (i, j, mismatch, label)

                    if best is None or mismatch < best[2]:
                        best = candidate

    scan(ref_seq, consurf_seq, "ref_in_consurf")
    scan(consurf_seq, ref_seq, "consurf_in_ref")

    if best is None:
        raise ValueError("No reliable local mapping found")

    ref_i, consurf_i, mismatch, label = best

    if label == "ref_in_consurf":
        ref_start = ref_i
        consurf_start = consurf_i
    else:
        ref_start = consurf_i
        consurf_start = ref_i

    return ref_start, consurf_start, mismatch


def integrate_consurf(
    *,
    base_table: pd.DataFrame,
    consurf_df: pd.DataFrame,
    ref_start: int,
    consurf_start: int,
    overwrite_existing_column: bool = True,
) -> pd.DataFrame:
    """
    Integrate ConSurf scores into an integrated full-length table.

    Adds / overwrites the column:
        Evolutionary conservation
    """
    out = base_table.copy()
    col_name = "Evolutionary conservation"

    if col_name in out.columns and not overwrite_existing_column:
        raise ValueError(f"Column already exists and overwrite is disabled: {col_name}")

    # Original implementation always overwrites this column when rerun.
    out[col_name] = pd.NA

    for _, row in consurf_df.iterrows():
        consurf_pos = int(row["pos"]) - 1
        score = row["score"]

        ref_pos = ref_start + (consurf_pos - consurf_start) + 1

        if 1 <= ref_pos <= len(out):
            out.loc[out["Position"] == ref_pos, col_name] = score

    # Move to the 4th column, matching original behavior.
    cols = out.columns.tolist()

    if col_name in cols:
        cols.remove(col_name)

    cols.insert(3, col_name)
    out = out[cols]

    return out


def find_consurf_grades_file(
    *,
    consurf_dir: Path,
    query_id: str,
) -> Path:
    """
    Find ConSurf grades file for query_id.

    Search priority
    ---------------
    - query_id-based
    - UniProt-based

    Example patterns
    ----------------
    {query_id}_*_consurf_grades.txt
    {uniprot}_*_consurf_grades.txt
    {uniprot}_consurf_grades.txt
    """
    uniprot = query_id.split("_")[0]

    patterns = [
        f"{query_id}_*_consurf_grades.txt",
        f"{uniprot}_*_consurf_grades.txt",
        f"{uniprot}_consurf_grades.txt",
    ]

    for pattern in patterns:
        matches = glob(str(consurf_dir / pattern))

        if matches:
            return Path(sorted(matches)[0])

    raise FileNotFoundError(f"No ConSurf grades file found for {query_id}")


def run_consurf_integrate(
    *,
    query_fasta: Path,
    integrated_dir: Path,
    consurf_dir: Path,
    aa_order: str = "GAILVMFWPCSTYNQHKRDE",
    window: int = 25,
    max_mismatch: int = 5,
    overwrite_existing_column: bool = True,
) -> None:
    """
    Integrate ConSurf evolutionary conservation scores into integrated PSSM tables.

    Parameters
    ----------
    query_fasta : pathlib.Path
        Original full-length protein FASTA.

    integrated_dir : pathlib.Path
        Directory containing integrated PSSM + Scorecons tables.

    consurf_dir : pathlib.Path
        Directory containing external ConSurf grades files.

    aa_order : str, optional
        Amino-acid column order used by PSSM matrices.

    window : int, optional
        Local mapping window size.

    max_mismatch : int, optional
        Maximum allowed mismatches in local mapping window.

    overwrite_existing_column : bool, optional
        Whether to overwrite an existing "Evolutionary conservation" column.

    Returns
    -------
    None
    """
    start_time = time.time()

    print("\n╔══════════════════════════════════════════════════════════════════════╗")
    print("║          [ ConSurf Conservation Integration Stage Started ]          ║")
    print("╚══════════════════════════════════════════════════════════════════════╝\n")

    log_path = integrated_dir / "consurf_integrate_error.log"

    aa_order_list = list(aa_order)
    numeric_cols = aa_order_list + ["Po", "Hy", "Ch", "|Hy-Ch|", "|Hy-Po|"]

    print(f"📂 FASTA Path           : {query_fasta}")
    print(f"📁 Integrated Dir       : {integrated_dir}")
    print(f"📁 ConSurf Dir          : {consurf_dir}")
    print(f"🐛 Error Log            : {log_path}\n")

    original_seqs = load_original_sequences(query_fasta)

    integrated_files = sorted(integrated_dir.glob("*.tsv"))
    integrated_files = [
        path for path in integrated_files if not path.name.endswith("_metadata.tsv")
    ]

    total = len(integrated_files)

    if total == 0:
        raise ValueError(f"No integrated TSV files found in: {integrated_dir}")

    success = 0
    failed = 0

    for idx, path in enumerate(integrated_files):
        query_id = path.stem
        progress = (idx + 1) / total * 100

        print("─" * 70)
        print(f"▶️  [{idx + 1:3d}/{total:<3d} | {progress:5.1f}% ]  {query_id}")
        print("─" * 70)

        try:
            if query_id not in original_seqs:
                raise KeyError(f"Query sequence not found in FASTA: {query_id}")

            base_table = load_base_integrated_table(
                path=path,
                numeric_cols=numeric_cols,
            )

            consurf_path = find_consurf_grades_file(
                consurf_dir=consurf_dir,
                query_id=query_id,
            )

            consurf_df = parse_consurf_grades(consurf_path)

            ref_seq = original_seqs[query_id]
            mut_pos, mut_from, _ = parse_mutation(query_id)
            restored_seq = restore_reference_sequence(ref_seq, mut_pos, mut_from)

            consurf_seq = build_consurf_sequence(consurf_df)

            ref_start, consurf_start, mismatch = find_best_local_mapping(
                ref_seq=restored_seq,
                consurf_seq=consurf_seq,
                window=window,
                max_mismatch=max_mismatch,
            )

            merged = integrate_consurf(
                base_table=base_table,
                consurf_df=consurf_df,
                ref_start=ref_start,
                consurf_start=consurf_start,
                overwrite_existing_column=overwrite_existing_column,
            )

            merged = restore_numeric_types(merged, numeric_cols)
            merged.to_csv(path, sep="\t", index=False)

            print(
                f"   ✅ mapped (ref_start={ref_start + 1}, "
                f"consurf_start={consurf_start + 1}, mismatch={mismatch})"
            )
            print(f"   🧬 integrated -> {path}\n")

            success += 1

        except Exception as exc:  # pylint: disable=broad-exception-caught
            failed += 1
            print(f"   ❌ Failed: {exc}\n")

            with log_path.open("a", encoding="utf-8") as log:
                log.write(f"[{query_id}] {exc}\n")

    elapsed = time.time() - start_time

    print("\n══════════════════════════════════════════════════════════════════════")
    print("ConSurf Integration Summary")
    print("──────────────────────────────────────────────────────────────────────")
    print(f"Total       : {total}")
    print(f"Successful  : {success}")
    print(f"Failed      : {failed}")
    print(f"🐛 Error log: {log_path}")
    print(f"Elapsed     : {elapsed:.2f} sec")
    print("══════════════════════════════════════════════════════════════════════\n")
