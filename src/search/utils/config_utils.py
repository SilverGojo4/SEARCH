# pylint: disable=
"""
Configuration Utility Module.

This module provides helper functions for loading and writing
configuration and manifest files (e.g., YAML, JSON) used by pipelines,
datasets, and CLI execution.
"""

from __future__ import annotations

# Standard Library Imports
import json
from pathlib import Path
from typing import Any, Dict

# Third-Party Library Imports
import yaml

# App Imports
from search.utils.fs_utils import is_file


def load_yaml(path: Path | str) -> Dict[str, Any]:
    """
    Load a YAML configuration or manifest file.

    Parameters
    ----------
    path : pathlib.Path or str
        Path to the YAML file.

    Returns
    -------
    dict
        Parsed YAML content.
    """
    path = Path(path)

    if not is_file(str(path)):
        raise FileNotFoundError(f"YAML file not found: '{path}'")

    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    if data is None:
        raise ValueError(f"YAML file is empty: '{path}'")

    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping: '{path}'")

    return data


def write_yaml(
    path: Path | str,
    data: Dict[str, Any],
    *,
    overwrite: bool = True,
) -> None:
    """
    Write a dictionary to a YAML file.

    Parameters
    ----------
    path : pathlib.Path or str
        Output YAML path.
    data : dict
        Data to serialize.
    overwrite : bool, optional
        Whether to overwrite an existing file.

    Returns
    -------
    None
    """
    path = Path(path)

    if path.exists() and not overwrite:
        raise FileExistsError(f"YAML file already exists: '{path}'")

    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            data,
            f,
            sort_keys=False,
            allow_unicode=True,
        )


def load_json(path: Path | str) -> Dict[str, Any]:
    """
    Load a JSON configuration file.

    Parameters
    ----------
    path : pathlib.Path or str
        Path to the JSON file.

    Returns
    -------
    dict
        Parsed JSON content.
    """

    path = Path(path)

    if not is_file(str(path)):
        raise FileNotFoundError(f"JSON file not found: '{path}'")

    with path.open("r", encoding="utf-8") as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON format: '{path}'") from exc

    if data is None:
        raise ValueError(f"JSON file is empty: '{path}'")

    if not isinstance(data, dict):
        raise ValueError(f"JSON root must be a mapping: '{path}'")

    return data
