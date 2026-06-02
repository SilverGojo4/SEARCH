# pylint: disable=missing-function-docstring
"""
SEARCH command-line interface.

This module provides the initial CLI scaffold for the SEARCH project.

The full pipeline registry will be added together with the actual
pipeline implementations in a later feature commit.
"""

from __future__ import annotations

# Standard Library Imports
import argparse
from pathlib import Path
from typing import Dict

# App Imports
from search.exceptions import PipelineError
from search.utils.config_utils import load_yaml

# Pipeline registry
#
# The registry is intentionally empty in the initial project scaffold.
# Real pipeline entries should be added together with their implementation
# modules, for example:
#
# PIPELINES = {
#     "cdsearch": "search.pipelines.cdsearch:run_cdsearch_pipeline",
# }
PIPELINES: Dict[str, str] = {}


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

    In the initial scaffold, this function only checks whether the
    requested pipeline has been registered. Actual pipeline dispatching
    should be added together with real pipeline implementations.
    """
    if pipeline_name not in PIPELINES:
        _print_pipeline_error(
            PipelineError(
                pipeline=pipeline_name,
                stage="resolve_pipeline",
                reason="Pipeline is not registered.",
                action=(
                    "This initial scaffold does not include real pipeline "
                    "implementations yet. Add the pipeline module under "
                    "'src/search/pipelines/' and register it in 'PIPELINES'."
                ),
                context={
                    "requested_pipeline": pipeline_name,
                    "config_path": str(config_path),
                    "registered_pipelines": list(PIPELINES.keys()),
                },
            )
        )
        return

    # Real dispatch logic will be added in the feature commit.
    _print_pipeline_error(
        PipelineError(
            pipeline=pipeline_name,
            stage="run_pipeline",
            reason="Pipeline dispatch is not implemented in the scaffold CLI.",
            action=(
                "Add importlib-based dispatch when the real pipeline "
                "implementations are committed."
            ),
            context={
                "config_path": str(config_path),
            },
        )
    )


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
