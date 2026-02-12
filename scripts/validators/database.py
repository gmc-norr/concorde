"""Database integrity validation.

This module provides validators for checking database integrity,
such as detecting duplicate runs with the same parameters.
"""

from typing import Any

from sqlalchemy.orm import Session

from .base import ValidationError


def check_duplicate_run(
    session: Session,
    params: dict[str, Any]
) -> Any | None:
    """Check if run with same parameters already exists in database.

    Args:
        session: SQLAlchemy database session
        params: Dictionary with run parameters:
            - sample: Sample name
            - caller: Variant caller
            - pipeline_version: Pipeline version
            - decomposition_mode: Decomposition mode
            - comparison_tool: Comparison tool (optional)

    Returns:
        Existing Run object if found, None otherwise

    Example:
        >>> existing_run = check_duplicate_run(session, {
        ...     "sample": "NA12878",
        ...     "caller": "GATK",
        ...     "pipeline_version": "v1.0",
        ...     "decomposition_mode": "decomposed"
        ... })
        >>> if existing_run:
        ...     print(f"Duplicate run found: {existing_run.id}")
    """
    # Import here to avoid circular dependency
    try:
        from models.run import Run
    except ImportError:
        # Models not available (e.g., in testing without backend)
        return None

    # Extract parameters
    sample = params.get("sample")
    caller = params.get("caller")
    pipeline_version = params.get("pipeline_version")
    decomposition_mode = params.get("decomposition_mode")
    comparison_tool = params.get("comparison_tool")

    if not all([sample, caller, pipeline_version]):
        return None

    # Query for existing run with same key parameters
    query = session.query(Run).filter(
        Run.sample == sample,
        Run.caller == caller,
        Run.pipeline_version == pipeline_version
    )

    # Add optional filters
    if decomposition_mode:
        query = query.filter(Run.decomposition_mode == decomposition_mode)

    if comparison_tool:
        query = query.filter(Run.comparison_tool == comparison_tool)

    existing_run = query.first()

    return existing_run


def validate_no_duplicate_run(
    session: Session,
    params: dict[str, Any],
    allow_duplicates: bool = False
) -> None:
    """Validate that no run with same parameters already exists.

    Args:
        session: SQLAlchemy database session
        params: Run parameters dictionary
        allow_duplicates: If False, raise error if duplicate found

    Raises:
        ValidationError: If duplicate run exists and not allowed

    Example:
        >>> validate_no_duplicate_run(session, {
        ...     "sample": "NA12878",
        ...     "caller": "GATK",
        ...     "pipeline_version": "v1.0",
        ...     "decomposition_mode": "decomposed"
        ... })
    """
    existing_run = check_duplicate_run(session, params)

    if existing_run and not allow_duplicates:
        raise ValidationError(
            f"Duplicate run already exists in database:\n"
            f"  Run ID: {existing_run.id}\n"
            f"  Sample: {existing_run.sample}\n"
            f"  Caller: {existing_run.caller}\n"
            f"  Pipeline Version: {existing_run.pipeline_version}\n"
            f"  Decomposition Mode: {existing_run.decomposition_mode}\n"
            f"  Date: {existing_run.date}\n\n"
            f"To proceed with a new run, either:\n"
            f"  1. Delete the existing run from the database\n"
            f"  2. Use different parameters (caller, version, or mode)\n"
            f"  3. Set allow_duplicates=True if duplicate runs are intentional"
        )
