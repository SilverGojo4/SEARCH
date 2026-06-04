"""
Conservation filtering pipeline orchestration module.

This module validates configuration sections, checks required input
resources, prepares output directories, and dispatches the actual
conservation-based PSSM feature masking logic.

The actual filtering logic is implemented in:

    search.processing.conservation_filter
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Any, Dict

# App Imports
from search.exceptions import PipelineError
from search.utils.fs_utils import remove_directory_tree

PIPELINE_NAME = "conservation_filter"


def run_conservation_filter_pipeline(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> None:
    """
    Run the conservation-based PSSM feature masking pipeline.

    This function is the pipeline-level entry point called by:

        search run --pipeline conservation_filter \\
          --config configs/stages/conservation_filter.yaml

    Parameters
    ----------
    inputs : dict
        Input paths.

        Expected keys:
        - integrated_dir

    parameters : dict
        Filtering parameters.

        Expected keys:
        - cons_min
        - cons_max
        - ec_min
        - ec_max
        - aa_order

    outputs : dict
        Output locations.

        Expected keys:
        - filtered_output_dir

    execution : dict
        Execution behavior.

        Expected keys:
        - overwrite
        - resume

    Returns
    -------
    None
    """
    config = _normalize_conservation_filter_config(
        inputs=inputs,
        parameters=parameters,
        outputs=outputs,
        execution=execution,
    )

    _validate_conservation_filter_resources(config)
    _prepare_conservation_filter_output_dir(config)
    _print_conservation_filter_configuration(config)

    # Import here to keep CLI startup lightweight.
    from search.processing.conservation_filter import (  # pylint: disable=import-outside-toplevel
        run_conservation_filter,
    )

    run_conservation_filter(
        integrated_dir=config["integrated_dir"],
        filtered_output_dir=config["filtered_output_dir"],
        cons_min=config["cons_min"],
        cons_max=config["cons_max"],
        ec_min=config["ec_min"],
        ec_max=config["ec_max"],
        aa_order=config["aa_order"],
        resume=config["resume"],
    )


def _normalize_conservation_filter_config(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Normalize conservation filter configuration sections into a flat dictionary.

    Defaults preserve the original Branch B conservation filtering behavior.
    """
    required_inputs = ["integrated_dir"]
    required_outputs = ["filtered_output_dir"]

    for key in required_inputs:
        if key not in inputs or inputs[key] in (None, ""):
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="validate_config",
                reason=f"Missing required input field: inputs.{key}",
                action=(
                    "Add the missing field to "
                    "configs/stages/conservation_filter.yaml under "
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
                    "configs/stages/conservation_filter.yaml under "
                    "the 'outputs' section."
                ),
                context={
                    "missing_key": f"outputs.{key}",
                },
            )

    config = {
        # Inputs
        "integrated_dir": Path(inputs["integrated_dir"]),
        # Parameters
        "cons_min": parameters.get("cons_min", 0.15),
        "cons_max": parameters.get("cons_max", 1.0),
        "ec_min": parameters.get("ec_min", -1.5),
        "ec_max": parameters.get("ec_max", 1.5),
        # Keep original Branch B AA column order by default.
        "aa_order": parameters.get("aa_order", "GAILVMFWPCSTYNQHKRDE"),
        # Outputs
        "filtered_output_dir": Path(outputs["filtered_output_dir"]),
        # Execution behavior
        "overwrite": bool(execution.get("overwrite", False)),
        "resume": bool(execution.get("resume", True)),
    }

    _validate_conservation_filter_parameters(config)

    return config


def _validate_conservation_filter_parameters(config: Dict[str, Any]) -> None:
    """
    Validate conservation filter parameter values.
    """
    for key in ["cons_min", "cons_max", "ec_min", "ec_max"]:
        try:
            config[key] = float(config[key])
        except (TypeError, ValueError) as exc:
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="validate_parameters",
                reason=f"Invalid numeric threshold parameter: {key}.",
                action=(
                    f"Set parameters.{key} to a numeric value in "
                    "configs/stages/conservation_filter.yaml."
                ),
                context={
                    key: config.get(key),
                },
            ) from exc

    if config["cons_min"] > config["cons_max"]:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="Invalid Scorecons conservation range.",
            action="Ensure parameters.cons_min is less than or equal to cons_max.",
            context={
                "cons_min": config["cons_min"],
                "cons_max": config["cons_max"],
            },
        )

    if config["ec_min"] > config["ec_max"]:
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_parameters",
            reason="Invalid evolutionary conservation range.",
            action="Ensure parameters.ec_min is less than or equal to ec_max.",
            context={
                "ec_min": config["ec_min"],
                "ec_max": config["ec_max"],
            },
        )

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


def _validate_conservation_filter_resources(config: Dict[str, Any]) -> None:
    """
    Validate integrated directory and integrated TSV files.
    """
    integrated_dir = config["integrated_dir"]

    if not integrated_dir.is_dir():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="Integrated directory does not exist.",
            action=(
                "Run scorecons_integrate and consurf_integrate first, or check "
                "inputs.integrated_dir in configs/stages/conservation_filter.yaml."
            ),
            context={
                "integrated_dir": str(integrated_dir),
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
                "Check that the conservation integration stages produced "
                "per-query TSV files in inputs.integrated_dir."
            ),
            context={
                "integrated_dir": str(integrated_dir),
            },
        )


def _prepare_conservation_filter_output_dir(config: Dict[str, Any]) -> None:
    """
    Prepare the filtered output directory according to execution settings.
    """
    output_dir = config["filtered_output_dir"]
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
                    "or set execution.resume: true to skip existing filtered TSV files."
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


def _print_conservation_filter_configuration(config: Dict[str, Any]) -> None:
    """
    Print a simple conservation filter configuration summary.
    """
    print()
    print("[Conservation Filter Configuration]")
    print(f"Integrated dir     : {config['integrated_dir']}")
    print(f"Filtered output dir: {config['filtered_output_dir']}")
    print()
    print("[Parameters]")
    print(f"cons_min           : {config['cons_min']}")
    print(f"cons_max           : {config['cons_max']}")
    print(f"ec_min             : {config['ec_min']}")
    print(f"ec_max             : {config['ec_max']}")
    print(f"aa_order           : {config['aa_order']}")
    print()
    print("[Execution]")
    print(f"overwrite          : {config['overwrite']}")
    print(f"resume             : {config['resume']}")
    print()
