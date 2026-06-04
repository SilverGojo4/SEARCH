# pylint: disable=missing-function-docstring
"""
SEARCH command-line interface.

This module provides the command-line interface for running SEARCH
pipeline stages.

The CLI supports:

- Listing registered pipelines
- Validating YAML configuration files
- Running registered pipeline stages

Pipeline implementations are registered through the PIPELINES registry.
Each registry entry should use the format:

    "module.path:function_name"

Example:

    "search.pipelines.cdsearch:run_cdsearch_pipeline"
"""

from __future__ import annotations

# Standard Library Imports
import argparse
import importlib
from pathlib import Path
from typing import Any, Callable, Dict

# App Imports
from search.exceptions import PipelineError
from search.utils.config_utils import load_yaml

# Pipeline registry
#
# Each registered pipeline maps a user-facing pipeline identifier to
# a callable entry point.
#
# Format:
#   "module.path:function_name"
#
# Example:
#   "search.pipelines.cdsearch:run_cdsearch_pipeline"
PIPELINES: Dict[str, str] = {
    "cdsearch": "search.pipelines.cdsearch:run_cdsearch_pipeline",
    "cdsearch_extract": (
        "search.pipelines.cdsearch_extract:run_cdsearch_extract_pipeline"
    ),
    "smp_parse": "search.pipelines.smp_parse:run_smp_parse_pipeline",
    "pssm_reconstruct": (
        "search.pipelines.pssm_reconstruct:run_pssm_reconstruct_pipeline"
    ),
}


def main() -> None:
    """
    Entry point for the SEARCH command-line interface.
    """
    parser = argparse.ArgumentParser(
        prog="search",
        description="SEARCH pipeline runner",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    # search list
    subparsers.add_parser(
        "list",
        help="List registered pipelines",
    )

    # search validate-config
    validate_parser = subparsers.add_parser(
        "validate-config",
        help="Validate that a YAML config file can be loaded",
    )
    validate_parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to a YAML config file",
    )

    # search run
    run_parser = subparsers.add_parser(
        "run",
        help="Run a registered pipeline",
    )
    run_parser.add_argument(
        "--pipeline",
        required=True,
        help="Pipeline identifier. Example: cdsearch",
    )
    run_parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to pipeline YAML config file",
    )

    args = parser.parse_args()

    if args.command == "list":
        _list_pipelines()
        return

    if args.command == "validate-config":
        _validate_config(config_path=args.config)
        return

    if args.command == "run":
        _run_pipeline(
            pipeline_name=args.pipeline,
            config_path=args.config,
        )
        return


def _list_pipelines() -> None:
    """
    Print registered pipeline names.
    """
    print()
    print("=" * 80)
    print("[SEARCH Registered Pipelines]")
    print("=" * 80)

    if not PIPELINES:
        print("No pipelines are registered yet.")
        print()
        print("Add pipeline implementations under:")
        print("  src/search/pipelines/")
        print()
        print("Then register them in:")
        print("  src/search/cli.py")

    else:
        for name in sorted(PIPELINES):
            print(f"- {name}")

    print("=" * 80)
    print()


def _validate_config(
    *,
    config_path: Path,
) -> None:
    """
    Validate that a YAML config file can be loaded.
    """
    try:
        config = load_yaml(config_path)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        _print_pipeline_error(
            PipelineError(
                pipeline="execution",
                stage="load_config",
                reason="Failed to load YAML configuration file.",
                action="Check whether the config path exists and the YAML syntax is valid.",
                context={
                    "config_path": str(config_path),
                    "error": str(exc),
                },
            )
        )
        return

    print()
    print("=" * 80)
    print("[SEARCH Config Validation]")
    print("=" * 80)
    print(f"Config : {config_path}")
    print("Status : valid")

    sections = ["inputs", "parameters", "outputs", "execution"]
    found_sections = [section for section in sections if section in config]

    if found_sections:
        print()
        print("Detected sections:")
        for section in found_sections:
            print(f"- {section}")

    print("=" * 80)
    print()


def _run_pipeline(
    *,
    pipeline_name: str,
    config_path: Path,
) -> None:
    """
    Run a registered pipeline.
    """
    if pipeline_name not in PIPELINES:
        _print_pipeline_error(
            PipelineError(
                pipeline=pipeline_name,
                stage="resolve_pipeline",
                reason="Pipeline is not registered.",
                action=(
                    "Run 'search list' to see available pipelines, or register "
                    "the pipeline in PIPELINES inside src/search/cli.py."
                ),
                context={
                    "requested_pipeline": pipeline_name,
                    "config_path": str(config_path),
                    "registered_pipelines": list(PIPELINES.keys()),
                },
            )
        )
        return

    try:
        config = load_yaml(config_path)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        _print_pipeline_error(
            PipelineError(
                pipeline=pipeline_name,
                stage="load_config",
                reason="Failed to load YAML configuration file.",
                action="Check whether the config path exists and the YAML syntax is valid.",
                context={
                    "config_path": str(config_path),
                    "error": str(exc),
                },
            )
        )
        return

    try:
        pipeline_func = _load_pipeline_entrypoint(
            pipeline_name=pipeline_name,
            entrypoint=PIPELINES[pipeline_name],
        )
    except PipelineError as exc:
        _print_pipeline_error(exc)
        return

    try:
        pipeline_func(
            inputs=config.get("inputs", {}),
            parameters=config.get("parameters", {}),
            outputs=config.get("outputs", {}),
            execution=config.get("execution", {}),
        )
    except PipelineError as exc:
        _print_pipeline_error(exc)
        return
    except Exception as exc:  # pylint: disable=broad-exception-caught
        _print_pipeline_error(
            PipelineError(
                pipeline=pipeline_name,
                stage="run_pipeline",
                reason="Unexpected error occurred while running the pipeline.",
                action=(
                    "Inspect the error message below. If this is a code-level "
                    "bug, rerun with the same config after fixing the pipeline."
                ),
                context={
                    "config_path": str(config_path),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                },
            )
        )
        return


def _load_pipeline_entrypoint(
    *,
    pipeline_name: str,
    entrypoint: str,
) -> Callable[..., Any]:
    """
    Load a pipeline entry point from a registry string.

    The registry string must use the format:

        module.path:function_name
    """
    if ":" not in entrypoint:
        raise PipelineError(
            pipeline=pipeline_name,
            stage="load_pipeline",
            reason="Invalid pipeline entrypoint format.",
            action=(
                "Use the format 'module.path:function_name' in the PIPELINES "
                "registry."
            ),
            context={
                "entrypoint": entrypoint,
            },
        )

    module_path, function_name = entrypoint.split(":", maxsplit=1)

    try:
        module = importlib.import_module(module_path)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        raise PipelineError(
            pipeline=pipeline_name,
            stage="load_pipeline",
            reason="Failed to import pipeline module.",
            action="Check whether the pipeline module path is correct.",
            context={
                "module_path": module_path,
                "function_name": function_name,
                "error": str(exc),
            },
        ) from exc

    try:
        pipeline_func = getattr(module, function_name)
    except AttributeError as exc:
        raise PipelineError(
            pipeline=pipeline_name,
            stage="load_pipeline",
            reason="Pipeline function was not found in the module.",
            action="Check whether the function name in PIPELINES is correct.",
            context={
                "module_path": module_path,
                "function_name": function_name,
            },
        ) from exc

    if not callable(pipeline_func):
        raise PipelineError(
            pipeline=pipeline_name,
            stage="load_pipeline",
            reason="Pipeline entrypoint is not callable.",
            action="Check that the registered entrypoint points to a function.",
            context={
                "module_path": module_path,
                "function_name": function_name,
            },
        )

    return pipeline_func


def _print_pipeline_error(exc: PipelineError) -> None:
    """
    Render pipeline errors in a human-readable CLI format.
    """
    print()
    print("=" * 80)
    print("[PIPELINE ERROR]")
    print("=" * 80)
    print(f"Pipeline : {exc.pipeline}")
    print(f"Stage    : {exc.stage}")

    print()
    print("Reason:")
    print(f"  {exc.reason}")

    if exc.action:
        print()
        print("Action:")
        print(f"  {exc.action}")

    if exc.context:
        print()
        print("Context:")
        for key, value in exc.context.items():
            print(f"  {key}: {value}")

    print("=" * 80)
    print()


if __name__ == "__main__":
    main()
