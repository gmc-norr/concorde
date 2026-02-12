"""Shared utility functions for pipeline scripts.

This module provides common helper functions used across multiple
pipeline scripts to avoid code duplication and ensure consistency.
"""

from typing import Any


def safe_float(val: Any, default: float | None = None) -> float | None:
    """Convert value to float, returning default for invalid values.

    Args:
        val: Value to convert to float
        default: Default value to return if conversion fails

    Returns:
        Float value or default if conversion fails

    Examples:
        >>> safe_float("3.14")
        3.14
        >>> safe_float("invalid")
        None
        >>> safe_float("invalid", default=0.0)
        0.0
    """
    try:
        result = float(val)
        # Check for NaN and infinity, return default instead
        import math
        if math.isnan(result) or math.isinf(result):
            return default
        return result
    except (ValueError, TypeError):
        return default


def safe_int(val: Any, default: int | None = None) -> int | None:
    """Convert value to int, returning default for invalid values.

    Handles float strings by truncating to int (e.g., "3.14" -> 3).

    Args:
        val: Value to convert to int
        default: Default value to return if conversion fails

    Returns:
        Integer value or default if conversion fails

    Examples:
        >>> safe_int("42")
        42
        >>> safe_int("3.14")
        3
        >>> safe_int("invalid")
        None
        >>> safe_int("invalid", default=0)
        0
    """
    try:
        # Try direct int conversion first
        return int(val)
    except (ValueError, TypeError):
        # If that fails, try converting to float first (handles "3.14" case)
        try:
            return int(float(val))
        except (ValueError, TypeError):
            return default


def safe_str(val: Any, default: str | None = None) -> str | None:
    """Convert value to string, returning default for None/NaN values.

    Strips whitespace and returns default if whitespace-only.

    Args:
        val: Value to convert to string
        default: Default value to return if val is None/NaN/whitespace-only

    Returns:
        String value (stripped) or default if val is None/NaN/whitespace-only

    Examples:
        >>> safe_str("text")
        'text'
        >>> safe_str("  hello  ")
        'hello'
        >>> safe_str("")
        ''
        >>> safe_str(None)
        None
        >>> safe_str("   ")
        None
        >>> safe_str(None, default="N/A")
        'N/A'
    """
    if val is None:
        return default
    # Handle pandas NaN
    try:
        import math

        if math.isnan(float(val)):
            return default
    except (ValueError, TypeError):
        pass

    # Convert to string
    str_val = str(val)

    # If already empty, return as-is
    if str_val == "":
        return ""

    # Strip whitespace
    result = str_val.strip()

    # Return default if empty after stripping (was whitespace-only)
    if not result:
        return default

    return result


def compute_classification_metrics(
    tp: int, fp: int, fn: int
) -> dict[str, float | None]:
    """Compute precision, recall, and F1 from TP/FP/FN counts.

    Single source of truth for these metrics. Returns None for undefined
    metrics (e.g. precision when TP+FP=0) to distinguish "not computable"
    from "zero".

    Args:
        tp: True positive count
        fp: False positive count
        fn: False negative count

    Returns:
        Dict with precision, recall, f1 (each float or None)
    """
    precision = tp / (tp + fp) if (tp + fp) > 0 else None
    recall = tp / (tp + fn) if (tp + fn) > 0 else None
    f1: float | None
    if precision is not None and recall is not None and (precision + recall) > 0:
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = None
    return {"precision": precision, "recall": recall, "f1": f1}


def format_variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Format variant coordinates into a standardized key.

    Args:
        chrom: Chromosome name
        pos: Position
        ref: Reference allele
        alt: Alternate allele

    Returns:
        Formatted variant key string

    Examples:
        >>> format_variant_key("chr1", 12345, "A", "T")
        'chr1:12345:A>T'
    """
    return f"{chrom}:{pos}:{ref}>{alt}"


def parse_variant_key(key: str) -> tuple[str, int, str, str]:
    """Parse a variant key string into components.

    Args:
        key: Variant key in format "chrom:pos:ref>alt"

    Returns:
        Tuple of (chrom, pos, ref, alt)

    Raises:
        ValueError: If key format is invalid

    Examples:
        >>> parse_variant_key("chr1:12345:A>T")
        ('chr1', 12345, 'A', 'T')
    """
    try:
        chrom_pos, alleles = key.split(":")[:2], key.split(":")[2]
        chrom, pos_str = chrom_pos
        ref, alt = alleles.split(">")
        return chrom, int(pos_str), ref, alt
    except (ValueError, IndexError) as e:
        raise ValueError(f"Invalid variant key format: {key}") from e
