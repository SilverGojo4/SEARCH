"""
SMP parse pipeline orchestration module.

This module validates configuration sections, checks required input
resources, prepares output directories, and dispatches the actual CDD
.smp parsing logic.

The actual .smp parsing and matrix generation logic is implemented in:

    search.processing.smp_parse
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

PIPELINE_NAME = "smp_parse"


def run_smp_parse_pipeline(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> None:
    """
    Run the CDD .smp parsing pipeline.

    This function is the pipeline-level entry point called by:

        search run --pipeline smp_parse --config configs/stages/smp_parse.yaml

    Parameters
    ----------
    inputs : dict
        Input paths and external resources.

        Expected keys:
        - cdsearch_table
        - smp_root_dir

    parameters : dict
        Tool parameters that affect scientific results.

        Expected keys:
        - aa_order

    outputs : dict
        Output locations.

        Expected keys:
        - matrix_output_dir

    execution : dict
        Execution behavior.

        Expected keys:
        - overwrite
        - resume

    Returns
    -------
    None
    """
    config = _normalize_smp_parse_config(
        inputs=inputs,
        parameters=parameters,
        outputs=outputs,
        execution=execution,
    )

    _validate_smp_parse_resources(config)
    _prepare_smp_parse_output_dir(config)
    _print_smp_parse_configuration(config)

    # Import here to keep CLI startup lightweight.
    from search.processing.smp_parse import (
        run_smp_parse,  # pylint: disable=import-outside-toplevel
    )

    run_smp_parse(
        cdsearch_table=config["cdsearch_table"],
        smp_root_dir=config["smp_root_dir"],
        output_dir=config["matrix_output_dir"],
        aa_order=config["aa_order"],
        resume=config["resume"],
    )


def _normalize_smp_parse_config(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Normalize SMP parse configuration sections into a flat dictionary.

    Defaults preserve the original Branch B behavior.
    """
    required_inputs = ["cdsearch_table", "smp_root_dir"]
    required_outputs = ["matrix_output_dir"]

    for key in required_inputs:
        if key not in inputs or inputs[key] in (None, ""):
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="validate_config",
                reason=f"Missing required input field: inputs.{key}",
                action=(
                    "Add the missing field to configs/stages/smp_parse.yaml "
                    "under the 'inputs' section."
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
                    "Add the missing field to configs/stages/smp_parse.yaml "
                    "under the 'outputs' section."
                ),
                context={
                    "missing_key": f"outputs.{key}",
                },
            )

    config = {
        # Inputs
        "cdsearch_table": Path(inputs["cdsearch_table"]),
        "smp_root_dir": Path(inputs["smp_root_dir"]),
        # Parameters
        # Keep original Branch B AA column order by default.
        "aa_order": parameters.get("aa_order", "GAILVMFWPCSTYNQHKRDE"),
        # Outputs
        "matrix_output_dir": Path(outputs["matrix_output_dir"]),
        # Execution behavior
        "overwrite": bool(execution.get("overwrite", False)),
        "resume": bool(execution.get("resume", True)),
    }

    _validate_smp_parse_parameters(config)

    return config


def _validate_smp_parse_parameters(config: Dict[str, Any]) -> None:
    """
    Validate SMP parse parameter values.
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


def _validate_smp_parse_resources(config: Dict[str, Any]) -> None:
    """
    Validate CD-Search table and SMP root directory.
    """
    cdsearch_table = config["cdsearch_table"]
    smp_root_dir = config["smp_root_dir"]

    if not cdsearch_table.is_file():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="CD-Search top-hit table does not exist.",
            action=(
                "Run the cdsearch stage first, or check inputs.cdsearch_table "
                "in configs/stages/smp_parse.yaml."
            ),
            context={
                "cdsearch_table": str(cdsearch_table),
            },
        )

    if not smp_root_dir.is_dir():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="SMP root directory does not exist.",
            action=(
                "Check inputs.smp_root_dir in configs/stages/smp_parse.yaml. "
                "For the current project structure, it is usually: "
                "data/external/cdd/profiles."
            ),
            context={
                "smp_root_dir": str(smp_root_dir),
            },
        )

    _validate_cdsearch_table_columns(cdsearch_table)


def _validate_cdsearch_table_columns(cdsearch_table: Path) -> None:
    """
    Validate required columns in cdsearch_top_hits_detailed.tsv.
    """
    required_cols = {"query_id", "PSSM_ID", "title"}

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


def _prepare_smp_parse_output_dir(config: Dict[str, Any]) -> None:
    """
    Prepare the SMP matrix output directory according to execution settings.
    """
    output_dir = config["matrix_output_dir"]
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
                    "or set execution.resume: true to skip existing matrix TSV files."
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


def _print_smp_parse_configuration(config: Dict[str, Any]) -> None:
    """
    Print a simple SMP parse configuration summary.
    """
    print()
    print("[SMP Parse Configuration]")
    print(f"CD-Search table  : {config['cdsearch_table']}")
    print(f"SMP root dir     : {config['smp_root_dir']}")
    print(f"Matrix output dir: {config['matrix_output_dir']}")
    print()
    print("[Parameters]")
    print(f"aa_order         : {config['aa_order']}")
    print()
    print("[Execution]")
    print(f"overwrite        : {config['overwrite']}")
    print(f"resume           : {config['resume']}")
    print()
