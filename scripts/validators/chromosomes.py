"""Chromosome normalization and validation utilities.

This module provides utilities for normalizing chromosome names and
validating chromosome format consistency.
"""

from typing import Literal

import pandas as pd

from .base import ValidationError

# Standard human chromosomes
STANDARD_CHROMOSOMES_WITH_CHR = {
    f"chr{i}" for i in range(1, 23)
} | {"chrX", "chrY", "chrM", "chrMT"}

STANDARD_CHROMOSOMES_WITHOUT_CHR = {
    str(i) for i in range(1, 23)
} | {"X", "Y", "M", "MT"}

ALL_STANDARD_CHROMOSOMES = STANDARD_CHROMOSOMES_WITH_CHR | STANDARD_CHROMOSOMES_WITHOUT_CHR

# Thresholds for chr prefix format detection
CHR_FORMAT_MAJORITY_THRESHOLD = 0.7
CHR_FORMAT_MINORITY_THRESHOLD = 0.3


def normalize_chromosome(chrom: str) -> str:
    """Normalize chromosome name to standard format with 'chr' prefix.

    Normalization rules:
    - Add 'chr' prefix if missing: '1' -> 'chr1'
    - Convert MT variations: 'MT' -> 'chrM', 'chrMT' -> 'chrM'
    - Preserve non-standard chromosomes as-is with 'chr' prefix
    - Case-insensitive for standard chromosomes

    Args:
        chrom: Chromosome name to normalize

    Returns:
        Normalized chromosome name with 'chr' prefix

    Example:
        >>> normalize_chromosome("1")
        'chr1'
        >>> normalize_chromosome("chr1")
        'chr1'
        >>> normalize_chromosome("MT")
        'chrM'
        >>> normalize_chromosome("chrMT")
        'chrM'
        >>> normalize_chromosome("X")
        'chrX'
    """
    if not chrom:
        return chrom

    # Handle MT/M variations
    if chrom.upper() in ("MT", "CHRMT"):
        return "chrM"
    if chrom.upper() == "M" or chrom.upper() == "CHRM":
        return "chrM"

    # Add chr prefix if not present
    if not chrom.startswith("chr"):
        return f"chr{chrom}"

    return chrom


def get_chromosome_format(chromosomes: set[str]) -> Literal["with_chr", "without_chr", "mixed"]:
    """Detect chromosome naming format from a set of chromosome names.

    Args:
        chromosomes: Set of chromosome names

    Returns:
        "with_chr" if majority use chr prefix,
        "without_chr" if majority don't use chr prefix,
        "mixed" if roughly equal split

    Example:
        >>> get_chromosome_format({"chr1", "chr2", "chrX"})
        'with_chr'
        >>> get_chromosome_format({"1", "2", "X"})
        'without_chr'
        >>> get_chromosome_format({"chr1", "2", "X"})
        'mixed'
    """
    if not chromosomes:
        return "without_chr"

    with_chr = sum(1 for c in chromosomes if str(c).startswith("chr"))
    total = len(chromosomes)

    # Majority have chr prefix
    if with_chr > total * CHR_FORMAT_MAJORITY_THRESHOLD:
        return "with_chr"
    # Majority don't have chr prefix
    elif with_chr < total * CHR_FORMAT_MINORITY_THRESHOLD:
        return "without_chr"
    # Mixed format
    else:
        return "mixed"


def validate_chromosome_format(
    df: pd.DataFrame,
    expected_format: Literal["with_chr", "without_chr"],
    file_name: str
) -> None:
    """Validate chromosome format in DataFrame matches expected format.

    Args:
        df: DataFrame with 'chrom' column
        expected_format: Expected chromosome format
        file_name: File name for error messages

    Raises:
        ValidationError: If chromosome format doesn't match expected

    Example:
        >>> df = pd.DataFrame({"chrom": ["chr1", "chr2", "chrX"]})
        >>> validate_chromosome_format(df, "with_chr", "test.tsv")  # Pass
        >>> validate_chromosome_format(df, "without_chr", "test.tsv")  # Fail
    """
    if "chrom" not in df.columns:
        return

    chromosomes = set(df["chrom"].dropna().unique())
    actual_format = get_chromosome_format(chromosomes)

    if actual_format == "mixed":
        # Get examples of each format
        with_chr = [c for c in chromosomes if str(c).startswith("chr")][:3]
        without_chr = [c for c in chromosomes if not str(c).startswith("chr")][:3]

        raise ValidationError(
            f"{file_name} contains mixed chromosome formats:\n"
            f"  With 'chr' prefix: {with_chr}\n"
            f"  Without 'chr' prefix: {without_chr}\n"
            f"All chromosomes must use the same format."
        )

    if actual_format != expected_format:
        examples = sorted(list(chromosomes))[:5]
        raise ValidationError(
            f"{file_name} has {actual_format} chromosome format "
            f"(e.g., {examples}) but expected {expected_format} format. "
            f"Chromosomes must match the format used in the reference genome."
        )


def validate_standard_chromosomes(
    chromosomes: set[str],
    file_name: str,
    allow_nonstandard: bool = True
) -> None:
    """Validate chromosomes are standard human chromosomes.

    Args:
        chromosomes: Set of chromosome names to validate
        file_name: File name for error messages
        allow_nonstandard: If True, allow non-standard chromosomes (e.g., contigs)

    Raises:
        ValidationError: If non-standard chromosomes found and not allowed

    Example:
        >>> validate_standard_chromosomes({"chr1", "chr2", "chrX"}, "test.vcf")
        >>> validate_standard_chromosomes({"chr1", "chr99"}, "test.vcf", allow_nonstandard=False)
        ValidationError: ...
    """
    if not chromosomes:
        return

    # Normalize for comparison
    normalized = {normalize_chromosome(c) for c in chromosomes}

    # Check against standard chromosomes
    nonstandard = normalized - STANDARD_CHROMOSOMES_WITH_CHR

    if nonstandard and not allow_nonstandard:
        # Filter out common alternative contigs that are usually OK
        common_alts = {f"chr{i}_" for i in range(1, 23)} | {
            "chrUn", "chrEBV", "chrX_", "chrY_"
        }
        truly_nonstandard = {
            c for c in nonstandard
            if not any(c.startswith(alt) for alt in common_alts)
        }

        if truly_nonstandard:
            raise ValidationError(
                f"{file_name} contains non-standard chromosomes: "
                f"{sorted(truly_nonstandard)[:10]}. "
                f"Expected standard human chromosomes (1-22, X, Y, M)."
            )
