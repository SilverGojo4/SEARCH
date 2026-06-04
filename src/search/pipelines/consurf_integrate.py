"""
ConSurf integration pipeline orchestration module.

This module validates configuration sections, checks required input
resources, and dispatches the actual ConSurf conservation integration logic.

The actual ConSurf parsing, local sequence mapping, and in-place table update
logic is implemented in:

    search.processing.consurf_integrate
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Any, Dict

# App Imports
from search.exceptions import PipelineError

PIPELINE_NAME = "consurf_integrate"


def run_consurf_integrate_pipeline(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> None:
    """
    Run the ConSurf evolutionary conservation integration pipeline.

    This function is the pipeline-level entry point called by:

        search run --pipeline consurf_integrate \\
          --config configs/stages/consurf_integrate.yaml

    Parameters
    ----------
    inputs : dict
        Input paths.

        Expected keys:
        - query_fasta
        - integrated_dir
        - consurf_dir

    parameters : dict
        Tool parameters that affect mapping behavior.

        Expected keys:
        - aa_order
        - window
        - max_mismatch

    outputs : dict
        Output locations.

        Expected keys:
        - integrated_dir

        Notes
        -----
        ConSurf integration updates integrated tables in place.
        Therefore outputs.integrated_dir is expected to be the same as
        inputs.integrated_dir.

    execution : dict
        Execution behavior.

        Expected keys:
        - overwrite_existing_column

    Returns
    -------
    None
    """
    config = _normalize_consurf_integrate_config(
        inputs=inputs,
        parameters=parameters,
        outputs=outputs,
        execution=execution,
    )

    _validate_consurf_integrate_resources(config)
    _print_consurf_integrate_configuration(config)

    # Import here to keep CLI startup lightweight.
    from search.processing.consurf_integrate import (  # pylint: disable=import-outside-toplevel
        run_consurf_integrate,
    )

    run_consurf_integrate(
        query_fasta=config["query_fasta"],
        integrated_dir=config["integrated_dir"],
        consurf_dir=config["consurf_dir"],
        aa_order=config["aa_order"],
        window=config["window"],
        max_mismatch=config["max_mismatch"],
        overwrite_existing_column=config["overwrite_existing_column"],
    )


def _normalize_consurf_integrate_config(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Normalize ConSurf integration configuration sections into a flat dictionary.

    Defaults preserve the original ConSurf integration behavior.
    """
    required_inputs = [
        "query_fasta",
        "integrated_dir",
        "consurf_dir",
    ]
    required_outputs = ["integrated_dir"]

    for key in required_inputs:
        if key not in inputs or inputs[key] in (None, ""):
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="validate_config",
                reason=f"Missing required input field: inputs.{key}",
                action=(
                    "Add the missing field to "
                    "configs/stages/consurf_integrate.yaml under "
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
                    "configs/stages/consurf_integrate.yaml under "
                    "the 'outputs' section."
                ),
                context={
                    "missing_key": f"outputs.{key}",
                },
            )

    config = {
        # Inputs
        "query_fasta": Path(inputs["query_fasta"]),
        "integrated_dir": Path(inputs["integrated_dir"]),
        "consurf_dir": Path(inputs["consurf_dir"]),
        # Outputs
        # This should normally be identical to inputs.integrated_dir because
        # ConSurf integration updates files in place.
        "output_integrated_dir": Path(outputs["integrated_dir"]),
        # Parameters
        # Keep original Branch B AA column order by default.
        "aa_order": parameters.get("aa_order", "GAILVMFWPCSTYNQHKRDE"),
        # Original local mapping defaults from run_consurf_integrate.py.
        "window": parameters.get("window", 25),
        "max_mismatch": parameters.get("max_mismatch", 5),
        # Execution behavior
        # Original implementation always overwrites the Evolutionary conservation column.
        "overwrite_existing_column": bool(
            execution.get("overwrite_existing_column", True)
        ),
    }

    _validate_consurf_integrate_parameters(config)

    return config


def _validate_consurf_integrate_parameters(config: Dict[str, Any]) -> None:
    """
    Validate ConSurf integration parameter values.
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

    try:
        config["window"] = int(config["window"])
    except (TypeError, ValueError) as exc:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="Invalid window parameter.",
            action="Set parameters.window to a positive integer, for example 25.",
            context={
                "window": config.get("window"),
            },
        ) from exc

    if config["window"] < 1:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="window must be greater than or equal to 1.",
            action="Set parameters.window to a positive integer, for example 25.",
            context={
                "window": config["window"],
            },
        )

    try:
        config["max_mismatch"] = int(config["max_mismatch"])
    except (TypeError, ValueError) as exc:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="Invalid max_mismatch parameter.",
            action="Set parameters.max_mismatch to a non-negative integer, for example 5.",
            context={
                "max_mismatch": config.get("max_mismatch"),
            },
        ) from exc

    if config["max_mismatch"] < 0:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="max_mismatch must be greater than or equal to 0.",
            action="Set parameters.max_mismatch to a non-negative integer.",
            context={
                "max_mismatch": config["max_mismatch"],
            },
        )


def _validate_consurf_integrate_resources(config: Dict[str, Any]) -> None:
    """
    Validate FASTA, integrated directory, ConSurf directory, and integrated TSV files.
    """
    query_fasta = config["query_fasta"]
    integrated_dir = config["integrated_dir"]
    output_integrated_dir = config["output_integrated_dir"]
    consurf_dir = config["consurf_dir"]

    if not query_fasta.is_file():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="Input FASTA file does not exist.",
            action=(
                "Check inputs.query_fasta in "
                "configs/stages/consurf_integrate.yaml. "
                "This should be the same FASTA used in the previous stages."
            ),
            context={
                "query_fasta": str(query_fasta),
            },
        )

    if not integrated_dir.is_dir():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="Integrated directory does not exist.",
            action=(
                "Run the scorecons_integrate stage first, or check "
                "inputs.integrated_dir in configs/stages/consurf_integrate.yaml."
            ),
            context={
                "integrated_dir": str(integrated_dir),
            },
        )

    if not consurf_dir.is_dir():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="ConSurf directory does not exist.",
            action=(
                "Prepare external ConSurf grades files first, or check "
                "inputs.consurf_dir in configs/stages/consurf_integrate.yaml. "
                "For the current project structure, it is usually: "
                "data/external/conservation/consurf."
            ),
            context={
                "consurf_dir": str(consurf_dir),
            },
        )

    if integrated_dir != output_integrated_dir:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_outputs",
            reason="ConSurf integration is configured with different input and output directories.",
            action=(
                "Set outputs.integrated_dir to the same path as inputs.integrated_dir. "
                "This stage updates integrated TSV files in place to preserve the "
                "original implementation behavior."
            ),
            context={
                "inputs.integrated_dir": str(integrated_dir),
                "outputs.integrated_dir": str(output_integrated_dir),
            },
        )

    integrated_files = sorted(integrated_dir.glob("*.tsv"))
    integrated_files = [
        path for path in integrated_files if not path.name.endswith("_metadata.tsv")
    ]

    if not integrated_files:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="No integrated TSV files were found.",
            action=(
                "Run the scorecons_integrate stage first and check that "
                "results/pssm/integrated contains per-query TSV files."
            ),
            context={
                "integrated_dir": str(integrated_dir),
            },
        )


def _print_consurf_integrate_configuration(config: Dict[str, Any]) -> None:
    """
    Print a simple ConSurf integration configuration summary.
    """
    print()
    print("[ConSurf Integration Configuration]")
    print(f"Input FASTA              : {config['query_fasta']}")
    print(f"Integrated dir           : {config['integrated_dir']}")
    print(f"ConSurf dir              : {config['consurf_dir']}")
    print(f"Output integrated dir    : {config['output_integrated_dir']}")
    print()
    print("[Parameters]")
    print(f"aa_order                 : {config['aa_order']}")
    print(f"window                   : {config['window']}")
    print(f"max_mismatch             : {config['max_mismatch']}")
    print()
    print("[Execution]")
    print(f"overwrite_existing_column: {config['overwrite_existing_column']}")
    print()
