"""
CD-Search extraction pipeline orchestration module.

This module validates configuration sections, checks required input
resources, prepares output directories, and dispatches the selected-domain
extraction logic.

The CD-Search extraction stage reads:

- Query FASTA file
- CD-Search all-hits detailed table
- Selected top1 PSSM table

and generates:

- Selected domain hits detailed table
- Alignment blocks
- Domain FASTA files

This stage is intentionally separated from the CD-Search alignment stage so
that users can inspect, accept, or manually edit the selected top1 PSSM table
before extracting downstream domain artifacts.

The actual selected-domain extraction logic is implemented in:

    search.processing.cdsearch
"""

from __future__ import annotations

# Standard Library Imports
from pathlib import Path
from typing import Any, Dict

# App Imports
from search.exceptions import PipelineError
from search.utils.fs_utils import remove_directory_tree

PIPELINE_NAME = "cdsearch_extract"


def run_cdsearch_extract_pipeline(
    *,
    inputs: Dict[str, Any],
    parameters: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> None:
    """
    Run the CD-Search extraction pipeline.

    This function is the pipeline-level entry point called by:

        search run --pipeline cdsearch_extract --config configs/stages/cdsearch_extract.yaml

    This stage extracts selected domain alignment blocks and domain FASTA
    files based on a selected top1 PSSM table.

    Parameters
    ----------
    inputs : dict
        Input paths and external resources.

        Expected keys:
        - query_fasta
        - cdsearch_all_hits
        - selected_top1

    parameters : dict
        Tool parameters that affect scientific results.

        This stage does not currently require scientific parameters.

    outputs : dict
        Output locations.

        Expected keys:
        - extract_root

    execution : dict
        Execution behavior.

        Expected keys:
        - overwrite

    Returns
    -------
    None
    """
    del parameters  # This stage currently has no scientific parameters.

    config = _normalize_cdsearch_extract_config(
        inputs=inputs,
        outputs=outputs,
        execution=execution,
    )

    _validate_cdsearch_extract_resources(config)
    _prepare_cdsearch_extract_output_dir(config)
    _print_cdsearch_extract_configuration(config)

    # Import here to keep CLI startup lightweight and to avoid importing
    # Biopython/pandas until the pipeline is actually executed.
    from search.processing.cdsearch import (
        extract_selected_cdsearch_domains,  # pylint: disable=import-outside-toplevel
    )

    extract_selected_cdsearch_domains(
        query_fasta=config["query_fasta"],
        cdsearch_all_hits=config["cdsearch_all_hits"],
        selected_top1=config["selected_top1"],
        output_dir=config["extract_root"],
    )


def _normalize_cdsearch_extract_config(
    *,
    inputs: Dict[str, Any],
    outputs: Dict[str, Any],
    execution: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Normalize CD-Search extraction configuration sections into a flat dictionary.
    """
    required_inputs = [
        "query_fasta",
        "cdsearch_all_hits",
        "selected_top1",
    ]
    required_outputs = ["extract_root"]

    for key in required_inputs:
        if key not in inputs or inputs[key] in (None, ""):
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="validate_config",
                reason=f"Missing required input field: inputs.{key}",
                action=(
                    "Add the missing field to configs/stages/cdsearch_extract.yaml "
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
                    "Add the missing field to configs/stages/cdsearch_extract.yaml "
                    "under the 'outputs' section."
                ),
                context={
                    "missing_key": f"outputs.{key}",
                },
            )

    config = {
        # Inputs
        "query_fasta": Path(inputs["query_fasta"]),
        "cdsearch_all_hits": Path(inputs["cdsearch_all_hits"]),
        "selected_top1": Path(inputs["selected_top1"]),
        # Outputs
        "extract_root": Path(outputs["extract_root"]),
        # Execution behavior
        "overwrite": bool(execution.get("overwrite", False)),
    }

    return config


def _validate_cdsearch_extract_resources(config: Dict[str, Any]) -> None:
    """
    Validate query FASTA, CD-Search all-hits table, and selected top1 table.
    """
    query_fasta = config["query_fasta"]
    cdsearch_all_hits = config["cdsearch_all_hits"]
    selected_top1 = config["selected_top1"]

    if not query_fasta.is_file():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="Input FASTA file does not exist.",
            action=(
                "Check inputs.query_fasta in " "configs/stages/cdsearch_extract.yaml."
            ),
            context={
                "query_fasta": str(query_fasta),
            },
        )

    if not cdsearch_all_hits.is_file():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="CD-Search all-hits table does not exist.",
            action=(
                "Run the cdsearch stage first, or check inputs.cdsearch_all_hits "
                "in configs/stages/cdsearch_extract.yaml."
            ),
            context={
                "cdsearch_all_hits": str(cdsearch_all_hits),
            },
        )

    if not selected_top1.is_file():
        raise PipelineError(
            pipeline=PIPELINE_NAME,
            stage="validate_inputs",
            reason="Selected top1 PSSM table does not exist.",
            action=(
                "Use results/cdsearch/cdsearch_default_top1.tsv, or provide a "
                "manually edited selected top1 table in "
                "configs/stages/cdsearch_extract.yaml."
            ),
            context={
                "selected_top1": str(selected_top1),
            },
        )


def _prepare_cdsearch_extract_output_dir(config: Dict[str, Any]) -> None:
    """
    Prepare the CD-Search extraction output directory according to execution settings.

    If overwrite is true, the existing extraction output directory is removed
    and rebuilt from scratch.

    If overwrite is false, this function refuses to write into a non-empty
    output directory to avoid mixing old and new extracted artifacts.
    """
    output_dir = config["extract_root"]
    overwrite = config["overwrite"]

    if output_dir.exists():
        if overwrite:
            remove_directory_tree(str(output_dir), include_self=True)
            output_dir.mkdir(parents=True, exist_ok=True)
            return

        if any(output_dir.iterdir()):
            raise PipelineError(
                pipeline=PIPELINE_NAME,
                stage="prepare_outputs",
                reason="Output directory already exists and is not empty.",
                action=(
                    "Either set execution.overwrite: true to rebuild from scratch, "
                    "or choose a new outputs.extract_root path."
                ),
                context={
                    "output_dir": str(output_dir),
                    "overwrite": overwrite,
                },
            )

        output_dir.mkdir(parents=True, exist_ok=True)
        return

    output_dir.mkdir(parents=True, exist_ok=True)


def _print_cdsearch_extract_configuration(config: Dict[str, Any]) -> None:
    """
    Print a simple CD-Search extraction configuration summary.
    """
    print()
    print("[CD-Search Extract Configuration]")
    print(f"Query FASTA       : {config['query_fasta']}")
    print(f"All hits table    : {config['cdsearch_all_hits']}")
    print(f"Selected top1     : {config['selected_top1']}")
    print(f"Output directory  : {config['extract_root']}")
    print()
    print("[Execution]")
    print(f"overwrite         : {config['overwrite']}")
    print()
    print("[Expected Outputs]")
    print(
        "selected hits     : "
        f"{config['extract_root'] / 'selected_domain_hits_detailed.tsv'}"
    )
    print("alignment blocks  : " f"{config['extract_root'] / 'alignment_blocks'}")
    print("domain FASTAs     : " f"{config['extract_root'] / 'domains_fasta'}")
    print()
