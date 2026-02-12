"""Data quality validation for genomic data.

This module provides validators for DNA sequences, probabilities,
percentages, and QC metrics to ensure data quality.
"""

from typing import Any

from .base import ValidationError


def validate_dna_sequence(seq: str, field_name: str = "sequence") -> None:
    """Validate DNA sequence contains only valid bases (ACGTN).

    Args:
        seq: DNA sequence string
        field_name: Field name for error messages

    Raises:
        ValidationError: If sequence contains invalid characters

    Example:
        >>> validate_dna_sequence("ACGT")  # Pass
        >>> validate_dna_sequence("ACGX")  # Fail
    """
    if not seq:
        raise ValidationError(
            f"{field_name} cannot be empty"
        )

    valid_bases = set("ACGTNacgtn")
    invalid_chars = set(seq) - valid_bases

    if invalid_chars:
        raise ValidationError(
            f"{field_name} contains invalid DNA characters: {sorted(invalid_chars)}. "
            f"Valid characters are: A, C, G, T, N (case-insensitive). "
            f"Sequence: {seq[:50]}{'...' if len(seq) > 50 else ''}"
        )


def validate_probability(value: float, field_name: str) -> None:
    """Validate value is a valid probability (between 0.0 and 1.0).

    Args:
        value: Numeric value to validate
        field_name: Field name for error messages

    Raises:
        ValidationError: If value is not in [0.0, 1.0]

    Example:
        >>> validate_probability(0.95, "AF")  # Pass
        >>> validate_probability(1.5, "AF")   # Fail
    """
    try:
        value_float = float(value)
    except (ValueError, TypeError):
        raise ValidationError(
            f"{field_name} must be a number, got: {value} ({type(value).__name__})"
        )

    if not (0.0 <= value_float <= 1.0):
        raise ValidationError(
            f"{field_name} must be between 0.0 and 1.0 (inclusive), "
            f"got: {value_float}"
        )


def validate_percentage(value: float, field_name: str) -> None:
    """Validate value is a valid percentage (between 0 and 100).

    Args:
        value: Numeric value to validate
        field_name: Field name for error messages

    Raises:
        ValidationError: If value is not in [0, 100]

    Example:
        >>> validate_percentage(95.5, "mapping_rate")  # Pass
        >>> validate_percentage(150, "mapping_rate")   # Fail
    """
    try:
        value_float = float(value)
    except (ValueError, TypeError):
        raise ValidationError(
            f"{field_name} must be a number, got: {value} ({type(value).__name__})"
        )

    if not (0.0 <= value_float <= 100.0):
        raise ValidationError(
            f"{field_name} must be between 0 and 100 (inclusive), "
            f"got: {value_float}"
        )


def validate_coverage_metrics(qc_summary: dict[str, Any]) -> None:
    """Validate QC coverage metrics are reasonable.

    Checks:
    - Coverage values are non-negative
    - Percentages are in valid range [0, 100]
    - Mean coverage >= median coverage (generally true)

    Args:
        qc_summary: Dictionary of QC summary metrics

    Raises:
        ValidationError: If coverage metrics are unreasonable

    Example:
        >>> validate_coverage_metrics({
        ...     "mean_coverage": 30.5,
        ...     "median_coverage": 32.0,
        ...     "pct_bases_10x": 99.1,
        ...     "pct_bases_20x": 98.5,
        ...     "mapping_rate": 99.8
        ... })
    """
    # Validate coverage values are non-negative
    coverage_fields = ["mean_coverage", "median_coverage"]
    for field in coverage_fields:
        if field in qc_summary and qc_summary[field] is not None:
            value = qc_summary[field]
            try:
                value_float = float(value)
            except (ValueError, TypeError):
                raise ValidationError(
                    f"QC metric '{field}' must be a number, "
                    f"got: {value} ({type(value).__name__})"
                )

            if value_float < 0:
                raise ValidationError(
                    f"QC metric '{field}' must be non-negative, "
                    f"got: {value_float}"
                )

    # Validate percentage fields are in [0, 100]
    percentage_fields = [
        "pct_bases_10x", "pct_bases_20x", "pct_bases_30x",
        "mapping_rate", "duplicate_rate"
    ]
    for field in percentage_fields:
        if field in qc_summary and qc_summary[field] is not None:
            value = qc_summary[field]
            try:
                validate_percentage(value, f"QC metric '{field}'")
            except ValidationError:
                raise

    # Validate mean >= median for coverage (with tolerance for rounding)
    if ("mean_coverage" in qc_summary and "median_coverage" in qc_summary
            and qc_summary["mean_coverage"] is not None
            and qc_summary["median_coverage"] is not None):

        mean_cov = float(qc_summary["mean_coverage"])
        median_cov = float(qc_summary["median_coverage"])

        # Allow small tolerance for rounding differences
        if mean_cov < median_cov - 0.1:
            raise ValidationError(
                f"QC metrics: mean_coverage ({mean_cov}) is less than "
                f"median_coverage ({median_cov}). This is unusual and may "
                f"indicate an error in QC calculation."
            )

    # Validate insert size is positive
    if "median_insert_size" in qc_summary and qc_summary["median_insert_size"] is not None:
        value = qc_summary["median_insert_size"]
        try:
            value_float = float(value)
        except (ValueError, TypeError):
            raise ValidationError(
                f"QC metric 'median_insert_size' must be a number, "
                f"got: {value} ({type(value).__name__})"
            )

        if value_float <= 0:
            raise ValidationError(
                f"QC metric 'median_insert_size' must be positive, "
                f"got: {value_float}"
            )

    # Validate contamination is percentage
    if "contamination_pct" in qc_summary and qc_summary["contamination_pct"] is not None:
        value = qc_summary["contamination_pct"]
        try:
            validate_percentage(value, "QC metric 'contamination_pct'")
        except ValidationError:
            raise
