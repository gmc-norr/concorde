"""VCF-specific utility functions for pipeline scripts.

This module provides helper functions for working with VCF files using pysam.
"""

import math
from pathlib import Path
from typing import Any

import pysam


def safe_format(sample: Any, field: str, default: Any = None) -> Any:
    """Safely extract a FORMAT field from a pysam sample record.

    Args:
        sample: Pysam sample object
        field: FORMAT field name to extract
        default: Default value if field is missing or invalid

    Returns:
        Field value or default if extraction fails

    Examples:
        >>> safe_format(sample, "DP")
        30
        >>> safe_format(sample, "MISSING_FIELD", default=0)
        0
    """
    try:
        val = sample[field]
        if isinstance(val, tuple):
            return val[0] if len(val) == 1 else val
        if isinstance(val, float) and math.isnan(val):
            return default
        return val
    except (KeyError, TypeError):
        return default


def determine_variant_type(ref: str, alt: str) -> str:
    """Determine variant type from REF and ALT alleles.

    Args:
        ref: Reference allele
        alt: Alternate allele

    Returns:
        "SNP", "INDEL", or "COMPLEX"
    """
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    elif len(ref) == len(alt):
        return "COMPLEX"
    else:
        return "INDEL"


def build_query_annotations(query_vcf_path: str | Path) -> dict[tuple, dict]:
    """Index annotations from the original query VCF by variant coordinates.

    Args:
        query_vcf_path: Path to query VCF file

    Returns:
        Dictionary mapping (chrom, pos, ref, alt) tuples to annotation dicts

    Note:
        Annotation dicts contain keys: dp, gq, af, mq, qual, filter
    """
    annotations = {}
    with pysam.VariantFile(str(query_vcf_path)) as qvcf:
        sample_name = qvcf.header.samples[0] if len(qvcf.header.samples) > 0 else None
        for rec in qvcf:
            for alt in rec.alts or []:
                key = (rec.chrom, rec.pos, rec.ref, alt)
                sample_data = rec.samples[sample_name] if sample_name else None
                annotations[key] = {
                    "dp": safe_format(sample_data, "DP") if sample_data else None,
                    "gq": safe_format(sample_data, "GQ") if sample_data else None,
                    "af": safe_format(sample_data, "AF") if sample_data else None,
                    "mq": rec.info.get("MQ", None),
                    "qual": rec.qual,
                    "filter": ";".join(rec.filter.keys()) if rec.filter else "PASS",
                }
    return annotations
