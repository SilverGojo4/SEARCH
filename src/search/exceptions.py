# pylint: disable=
"""
Custom exception classes used across SEARCH pipelines.

These exceptions are designed to carry structured pipeline
error information for CLI reporting and debugging.
"""

from __future__ import annotations


class PipelineError(Exception):
    """
    Base class for all pipeline-related errors.
    """

    # pylint: disable=too-many-arguments
    def __init__(
        self,
        *,
        pipeline: str,
        stage: str,
        reason: str,
        action: str,
        context: dict | None = None,
    ):
        self.pipeline = pipeline
        self.stage = stage
        self.reason = reason
        self.action = action
        self.context = context or {}

        super().__init__(reason)
