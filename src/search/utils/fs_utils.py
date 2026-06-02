# pylint: disable=
"""
Filesystem Utility Module.

This module provides low-level filesystem and OS-related helper functions
that are independent of domain logic, data schemas, or pipeline semantics.
"""

from __future__ import annotations

# Standard Library Imports
import gzip
import os
import shutil
from typing import List, Optional


# File System Predicates
def is_file(file_path: str) -> bool:
    """
    Determine whether the given path exists and is a regular file.

    Parameters
    ----------
    file_path : str
        Absolute or relative path.

    Returns
    -------
    bool
        True if the path exists and is a file.
    """
    return os.path.isfile(file_path)


def is_directory(dir_path: str) -> bool:
    """
    Determine whether the given path exists and is a directory.

    Parameters
    ----------
    dir_path : str
        Absolute or relative path.

    Returns
    -------
    bool
        True if the path exists and is a directory.
    """
    return os.path.isdir(dir_path)


# Directory Manipulation Utilities
def remove_directory_tree(
    dir_path: str,
    include_self: bool = True,
    extensions: Optional[List[str]] = None,
) -> None:
    """
    Remove a directory tree or its contents.

    This function can either:
    - Delete the entire directory tree, or
    - Delete only its contents while preserving the directory itself.

    Parameters
    ----------
    dir_path : str
        Target directory path.
    include_self : bool, optional
        If True, remove the directory itself.
        If False, remove only its contents.
    extensions : list of str, optional
        If provided, only files with matching extensions are removed.

    Returns
    -------
    None
    """
    # Validate path existence
    if not os.path.exists(dir_path):
        raise FileNotFoundError(f"Directory not found: '{dir_path}'")

    if not is_directory(dir_path):
        raise NotADirectoryError(f"Path is not a directory: '{dir_path}'")

    # Delete all contents
    if include_self:
        shutil.rmtree(dir_path)
        return

    # Only delete contents inside (keep the folder)
    for entry in os.listdir(dir_path):
        entry_path = os.path.join(dir_path, entry)

        if os.path.isfile(entry_path):
            if extensions and not any(entry_path.endswith(ext) for ext in extensions):
                continue
            os.remove(entry_path)
        elif os.path.isdir(entry_path):
            shutil.rmtree(entry_path)


# Compression Utilities
def decompress_gzip_files(
    input_dir: str,
    recursive: bool = True,
) -> None:
    """
    Decompress all gzip-compressed (.gz) files in a directory.

    Notes
    -----
    - Decompressed files are written alongside the originals.
    - Original '.gz' files are deleted after successful extraction.
    - Safe to call multiple times (no-op if no '.gz' files exist).

    Parameters
    ----------
    input_dir : str
        Directory containing gzip files.
    recursive : bool, optional
        Whether to search subdirectories recursively.

    Returns
    -------
    None
    """
    # Validate input path
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Directory not found: '{input_dir}'")

    if not is_directory(input_dir):
        raise NotADirectoryError(f"Path is not a directory: '{input_dir}'")

    # Find .gz files
    gz_files: List[str] = []
    if recursive:
        for root, _, files in os.walk(input_dir):
            gz_files.extend(os.path.join(root, f) for f in files if f.endswith(".gz"))
    else:
        gz_files = [
            os.path.join(input_dir, f)
            for f in os.listdir(input_dir)
            if f.endswith(".gz")
        ]

    # If no compressed files — do nothing
    if not gz_files:
        return

    # Decompress all files
    for gz_path in gz_files:
        output_path = gz_path[:-3]
        try:
            with gzip.open(gz_path, "rb") as f_in, open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_path)
        except Exception as e:
            raise RuntimeError(f"Failed to decompress '{gz_path}': {e}") from e
