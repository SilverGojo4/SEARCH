"""
PSSM reconstruction pipeline orchestration module.

This module validates configuration sections, checks required input
resources, prepares output directories, and dispatches the actual
full-length PSSM reconstruction logic.

The actual reconstruction logic is implemented in:

    search.processing.pssm_reconstruct
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Any, Dict

# Third-Party Imports
import pandas as pd

# App Imports
from search.exceptions import PipelineError
from search.utils.fs_utils import remove_directory_tree

PIPELINE_NAME = "pssm_reconstruct"


def run_pssm_reconstruct_pipeline(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> None:
    """
    Run the full-length PSSM reconstruction pipeline.

    This function is the pipeline-level entry point called by:

        search run --pipeline pssm_reconstruct \
          --config configs/stages/pssm_reconstruct_smp.yaml

    Parameters
    ----------
    inputs : dict
        Input paths.

        Expected keys:
        - query_fasta
        - cdsearch_table
        - matrix_dir

    parameters : dict
        Tool parameters that affect scientific results.

        Expected keys:
        - aa_order

    outputs : dict
        Output locations.

        Expected keys:
        - reconstruct_output_dir

    execution : dict
        Execution behavior.

        Expected keys:
        - overwrite
        - resume

    Returns
    -------
    None
    """
    config = _normalize_pssm_reconstruct_config(
        inputs=inputs,
        parameters=parameters,
        outputs=outputs,
        execution=execution,
    )

    _validate_pssm_reconstruct_resources(config)
    _prepare_pssm_reconstruct_output_dir(config)
    _print_pssm_reconstruct_configuration(config)

    # Import here to keep CLI startup lightweight.
    from search.processing.pssm_reconstruct import (  # pylint: disable=import-outside-toplevel
        run_pssm_reconstruct,
    )

    run_pssm_reconstruct(
        query_fasta=config["query_fasta"],
        cdsearch_table=config["cdsearch_table"],
        matrix_dir=config["matrix_dir"],
        output_dir=config["reconstruct_output_dir"],
        aa_order=config["aa_order"],
        resume=config["resume"],
    )


def _normalize_pssm_reconstruct_config(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Normalize PSSM reconstruction configuration sections into a flat dictionary.

    Defaults preserve the original Branch B reconstruction behavior.
    """
    required_inputs = ["query_fasta", "cdsearch_table", "matrix_dir"]
    required_outputs = ["reconstruct_output_dir"]

    for key in required_inputs:
        if key not in inputs or inputs[key] in (None, ""):
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="validate_config",
                reason=f"Missing required input field: inputs.{key}",
                action=(
                    "Add the missing field to "
                    "configs/stages/pssm_reconstruct_smp.yaml under "
                    "the 'inputs' section."
                ),
                context={
                    "missing_key": f"inputs.{key}",
                },
            )

    for key in required_outputs:
        if key not in outputs or outputs[key] in (None, ""):
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="validate_config",
                reason=f"Missing required output field: outputs.{key}",
                action=(
                    "Add the missing field to "
                    "configs/stages/pssm_reconstruct_smp.yaml under "
                    "the 'outputs' section."
                ),
                context={
                    "missing_key": f"outputs.{key}",
                },
            )

    config = {
        # Inputs
        "query_fasta": Path(inputs["query_fasta"]),
        "cdsearch_table": Path(inputs["cdsearch_table"]),
        "matrix_dir": Path(inputs["matrix_dir"]),
        # Parameters
        # Keep original Branch B AA column order by default.
        "aa_order": parameters.get("aa_order", "GAILVMFWPCSTYNQHKRDE"),
        # Outputs
        "reconstruct_output_dir": Path(outputs["reconstruct_output_dir"]),
        # Execution behavior
        "overwrite": bool(execution.get("overwrite", False)),
        "resume": bool(execution.get("resume", True)),
    }

    _validate_pssm_reconstruct_parameters(config)

    return config


def _validate_pssm_reconstruct_parameters(config: Dict[str, Any]) -> None:
    """
    Validate PSSM reconstruction parameter values.
    """
    aa_order = config["aa_order"]

    if not isinstance(aa_order, str):
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="Invalid aa_order parameter.",
            action=(
                "Set parameters.aa_order to a 20-character amino-acid string, "
                "for example: GAILVMFWPCSTYNQHKRDE."
            ),
            context={
                "aa_order": aa_order,
            },
        )

    if len(aa_order) != 20:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="aa_order must contain exactly 20 amino-acid characters.",
            action="Use the original Branch B order: GAILVMFWPCSTYNQHKRDE.",
            context={
                "aa_order": aa_order,
                "length": len(aa_order),
            },
        )

    if len(set(aa_order)) != 20:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="aa_order contains duplicated amino-acid characters.",
            action="Use a unique 20-amino-acid order, preferably GAILVMFWPCSTYNQHKRDE.",
            context={
                "aa_order": aa_order,
            },
        )


def _validate_pssm_reconstruct_resources(config: Dict[str, Any]) -> None:
    """
    Validate FASTA, CD-Search table, matrix directory, and required table columns.
    """
    query_fasta = config["query_fasta"]
    cdsearch_table = config["cdsearch_table"]
    matrix_dir = config["matrix_dir"]

    if not query_fasta.is_file():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="Input FASTA file does not exist.",
            action=(
                "Check inputs.query_fasta in "
                "configs/stages/pssm_reconstruct_smp.yaml. "
                "This should be the same FASTA used in the cdsearch stage."
            ),
            context={
                "query_fasta": str(query_fasta),
            },
        )

    if not cdsearch_table.is_file():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="Selected CD-Search domain table does not exist.",
            action=(
                "Run the cdsearch and cdsearch_extract stages first, or check "
                "inputs.cdsearch_table in configs/stages/pssm_reconstruct.yaml."
            ),
            context={
                "cdsearch_table": str(cdsearch_table),
            },
        )

    if not matrix_dir.is_dir():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="PSSM matrix directory does not exist.",
            action=(
                "Run the smp_parse stage first, or check inputs.matrix_dir "
                "in configs/stages/pssm_reconstruct_smp.yaml."
            ),
            context={
                "matrix_dir": str(matrix_dir),
            },
        )

    _validate_cdsearch_table_columns(cdsearch_table)


def _validate_cdsearch_table_columns(cdsearch_table: Path) -> None:
    """
    Validate required columns in selected_domain_hits_detailed.tsv.

    The reconstruction stage needs both coordinate columns and aligned
    sequence columns to project domain-level matrices back to full-length
    protein coordinates.
    """
    required_cols = {
        "query_id",
        "PSSM_ID",
        "qstart",
        "qend",
        "qseq",
        "hseq",
    }

    try:
        df = pd.read_csv(cdsearch_table, sep="\t", nrows=5)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="Failed to read CD-Search top-hit table.",
            action="Check whether the file is a valid tab-separated table.",
            context={
                "cdsearch_table": str(cdsearch_table),
                "error": str(exc),
            },
        ) from exc

    missing_cols = required_cols.difference(df.columns)

    if missing_cols:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="CD-Search table is missing required columns.",
            action=(
                "Use the cdsearch_top_hits_detailed.tsv generated by the "
                "cdsearch stage without modifying its required columns."
            ),
            context={
                "cdsearch_table": str(cdsearch_table),
                "required_columns": sorted(required_cols),
                "missing_columns": sorted(missing_cols),
                "available_columns": list(df.columns),
            },
        )


def _prepare_pssm_reconstruct_output_dir(config: Dict[str, Any]) -> None:
    """
    Prepare the reconstruction output directory according to execution settings.
    """
    output_dir = config["reconstruct_output_dir"]
    overwrite = config["overwrite"]
    resume = config["resume"]

    if output_dir.exists():
        if overwrite:
            remove_directory_tree(str(output_dir), include_self=True)
            output_dir.mkdir(parents=True, exist_ok=True)
            return

        if not resume and any(output_dir.iterdir()):
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="prepare_outputs",
                reason="Output directory already exists and is not empty.",
                action=(
                    "Either set execution.overwrite: true to rebuild from scratch, "
                    "or set execution.resume: true to skip existing reconstructed TSV files."
                ),
                context={
                    "output_dir": str(output_dir),
                    "overwrite": overwrite,
                    "resume": resume,
                },
            )

        output_dir.mkdir(parents=True, exist_ok=True)
        return

    output_dir.mkdir(parents=True, exist_ok=True)


def _print_pssm_reconstruct_configuration(config: Dict[str, Any]) -> None:
    """
    Print a simple PSSM reconstruction configuration summary.
    """
    print()
    print("[PSSM Reconstruction Configuration]")
    print(f"Input FASTA         : {config['query_fasta']}")
    print(f"Selected domain tbl : {config['cdsearch_table']}")
    print(f"PSSM matrix dir     : {config['matrix_dir']}")
    print(f"Reconstruct out dir : {config['reconstruct_output_dir']}")
    print()
    print("[Parameters]")
    print(f"aa_order            : {config['aa_order']}")
    print()
    print("[Execution]")
    print(f"overwrite           : {config['overwrite']}")
    print(f"resume              : {config['resume']}")
    print()
