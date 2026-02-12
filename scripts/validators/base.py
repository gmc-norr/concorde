"""Base validation utilities for pipeline scripts.

This module provides core validation functions to ensure data integrity
before processing, preventing silent failures and improving error messages.
"""

from pathlib import Path

import pandas as pd


class ValidationError(Exception):
    """Custom exception for data validation failures."""

    pass


def validate_required_columns(
    df: pd.DataFrame, required: list[str], file_name: str
) -> None:
    """Validate that DataFrame has all required columns.

    Args:
        df: DataFrame to validate
        required: List of required column names
        file_name: Name of file being validated (for error messages)

    Raises:
        ValidationError: If required columns are missing
    """
    missing = set(required) - set(df.columns)
    if missing:
        raise ValidationError(
            f"{file_name} missing required columns: {sorted(missing)}"
        )


def validate_file_exists(file_path: str | Path, file_description: str) -> Path:
    """Validate that a file exists and is readable.

    Args:
        file_path: Path to file
        file_description: Description of file for error messages

    Returns:
        Resolved Path object

    Raises:
        ValidationError: If file doesn't exist or isn't readable
    """
    path = Path(file_path)
    if not path.exists():
        raise ValidationError(f"{file_description} not found: {path}")
    if not path.is_file():
        raise ValidationError(f"{file_description} is not a file: {path}")
    return path


def validate_enum_values(
    df: pd.DataFrame, column: str, valid_values: set[str], file_name: str
) -> None:
    """Validate that a column contains only valid enum values.

    Args:
        df: DataFrame to validate
        column: Column name to check
        valid_values: Set of valid values
        file_name: Name of file being validated

    Raises:
        ValidationError: If invalid values are found
    """
    if column not in df.columns:
        return  # Column validation should be done separately

    actual_values = set(df[column].dropna().unique())
    invalid = actual_values - valid_values

    if invalid:
        raise ValidationError(
            f"{file_name} column '{column}' contains invalid values: {sorted(invalid)}. "
            f"Valid values are: {sorted(valid_values)}"
        )


def validate_non_empty(df: pd.DataFrame, file_name: str) -> None:
    """Validate that DataFrame is not empty.

    Args:
        df: DataFrame to validate
        file_name: Name of file being validated

    Raises:
        ValidationError: If DataFrame is empty
    """
    if len(df) == 0:
        raise ValidationError(f"{file_name} contains no data rows")


def validate_positive_numeric(
    df: pd.DataFrame, columns: list[str], file_name: str
) -> None:
    """Validate that numeric columns contain non-negative values.

    Args:
        df: DataFrame to validate
        columns: List of column names that should be non-negative
        file_name: Name of file being validated

    Raises:
        ValidationError: If negative values are found
    """
    for col in columns:
        if col not in df.columns:
            continue

        numeric_vals = pd.to_numeric(df[col], errors="coerce").dropna()
        if (numeric_vals < 0).any():
            min_val = numeric_vals.min()
            raise ValidationError(
                f"{file_name} column '{col}' contains negative values (min: {min_val})"
            )


def validate_variant_coordinates(
    df: pd.DataFrame, file_name: str
) -> None:
    """Validate variant coordinate columns.

    Args:
        df: DataFrame to validate
        file_name: Name of file being validated

    Raises:
        ValidationError: If coordinate validation fails
    """
    # Check chromosome format
    if "chrom" in df.columns:
        invalid_chroms = df[df["chrom"].isna()]
        if len(invalid_chroms) > 0:
            raise ValidationError(
                f"{file_name} contains {len(invalid_chroms)} rows with missing chromosome"
            )

    # Check position is positive integer
    if "pos" in df.columns:
        try:
            positions = pd.to_numeric(df["pos"], errors="coerce")
            if positions.isna().any():
                raise ValidationError(
                    f"{file_name} contains non-numeric position values"
                )
            if (positions <= 0).any():
                raise ValidationError(
                    f"{file_name} contains non-positive position values"
                )
        except Exception as e:
            raise ValidationError(f"{file_name} position validation failed: {e}")


def validate_no_duplicate_variants(
    df: pd.DataFrame, file_name: str
) -> None:
    """Check for duplicate variants at the same genomic position.

    Args:
        df: DataFrame to validate (must have chrom, pos, ref, alt columns)
        file_name: Name of file being validated

    Raises:
        ValidationError: If duplicate variants are found

    Example:
        >>> df = pd.DataFrame({
        ...     "chrom": ["chr1", "chr1", "chr2"],
        ...     "pos": [100, 100, 200],
        ...     "ref": ["A", "A", "G"],
        ...     "alt": ["T", "T", "C"]
        ... })
        >>> validate_no_duplicate_variants(df, "test.tsv")
        ValidationError: test.tsv contains 1 duplicate variants
    """
    required_cols = ["chrom", "pos", "ref", "alt"]
    missing = set(required_cols) - set(df.columns)
    if missing:
        # Skip validation if required columns don't exist
        return

    duplicates = df[df.duplicated(subset=required_cols, keep=False)]
    if len(duplicates) > 0:
        # Count unique duplicate groups
        n_duplicate_groups = len(df[df.duplicated(subset=required_cols, keep="first")])
        raise ValidationError(
            f"{file_name} contains {n_duplicate_groups} duplicate variant(s). "
            f"First duplicate: {duplicates.iloc[0]['chrom']}:{duplicates.iloc[0]['pos']} "
            f"{duplicates.iloc[0]['ref']}>{duplicates.iloc[0]['alt']}"
        )


def validate_no_nulls(
    df: pd.DataFrame, columns: list[str], file_name: str
) -> None:
    """Validate that specified columns have no null/NA values.

    Args:
        df: DataFrame to validate
        columns: List of column names that should not have nulls
        file_name: Name of file being validated

    Raises:
        ValidationError: If null values are found

    Example:
        >>> df = pd.DataFrame({"required_col": [1, None, 3]})
        >>> validate_no_nulls(df, ["required_col"], "test.tsv")
        ValidationError: test.tsv column 'required_col' has 1 null values
    """
    for col in columns:
        if col not in df.columns:
            continue

        null_count = df[col].isna().sum()
        if null_count > 0:
            raise ValidationError(
                f"{file_name} column '{col}' has {null_count} null value(s). "
                f"This column must not contain missing values."
            )


def validate_tsv_schema(
    file_path: str | Path,
    required_columns: list[str],
    enum_columns: dict[str, set[str]] | None = None,
    positive_columns: list[str] | None = None,
) -> pd.DataFrame:
    """Comprehensive TSV validation with schema checking.

    Args:
        file_path: Path to TSV file
        required_columns: List of required column names
        enum_columns: Dict mapping column names to valid value sets
        positive_columns: List of columns that must be non-negative

    Returns:
        Validated DataFrame

    Raises:
        ValidationError: If any validation fails

    Example:
        >>> df = validate_tsv_schema(
        ...     "variants.tsv",
        ...     required_columns=["chrom", "pos", "ref", "alt", "type"],
        ...     enum_columns={"type": {"SNP", "INDEL"}},
        ...     positive_columns=["pos", "dp"]
        ... )
    """
    file_name = Path(file_path).name

    # Validate file exists
    path = validate_file_exists(file_path, f"TSV file {file_name}")

    # Load TSV
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:
        raise ValidationError(f"Failed to read {file_name}: {e}")

    # Validate non-empty
    validate_non_empty(df, file_name)

    # Validate required columns
    validate_required_columns(df, required_columns, file_name)

    # Validate enum values
    if enum_columns:
        for col, valid_vals in enum_columns.items():
            if col in df.columns:
                validate_enum_values(df, col, valid_vals, file_name)

    # Validate positive numeric columns
    if positive_columns:
        validate_positive_numeric(df, positive_columns, file_name)

    return df
