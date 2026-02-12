"""Cross-file consistency validation.

This module provides validators to ensure consistency across multiple
input files (VCFs, reference, etc.) before pipeline execution.
"""

from pathlib import Path

import pysam

from .base import ValidationError
from .input_files import validate_vcf_format, validate_vcf_chromosomes


def validate_chromosome_consistency(
    query_vcf: str | Path,
    truth_vcf: str | Path,
    reference: str | Path
) -> dict[str, set[str]]:
    """Ensure chromosome naming is consistent across files.

    Args:
        query_vcf: Path to query VCF file
        truth_vcf: Path to truth VCF file
        reference: Path to reference FASTA file

    Returns:
        Dictionary with chromosome sets from each file:
        {"query": {...}, "truth": {...}, "reference": {...}}

    Raises:
        ValidationError: If major chromosome inconsistencies found

    Example:
        >>> chroms = validate_chromosome_consistency(
        ...     "query.vcf.gz",
        ...     "truth.vcf.gz",
        ...     "ref.fa"
        ... )
    """
    # Get chromosomes from each file
    query_chroms = validate_vcf_chromosomes(query_vcf, quick_check=True)
    truth_chroms = validate_vcf_chromosomes(truth_vcf, quick_check=True)

    # Get chromosomes from reference FASTA
    try:
        with pysam.FastaFile(str(reference)) as fasta:
            ref_chroms = set(fasta.references)
    except Exception as e:
        raise ValidationError(
            f"Reference FASTA {Path(reference).name} could not be opened: {e}"
        )

    # Detect chromosome format for each file
    def _has_chr_prefix(chroms: set[str]) -> bool:
        """Check if majority of chromosomes use 'chr' prefix."""
        with_chr = sum(1 for c in chroms if str(c).startswith("chr"))
        return with_chr > len(chroms) / 2

    query_format = "with_chr" if _has_chr_prefix(query_chroms) else "without_chr"
    truth_format = "with_chr" if _has_chr_prefix(truth_chroms) else "without_chr"
    ref_format = "with_chr" if _has_chr_prefix(ref_chroms) else "without_chr"

    # Check for format inconsistencies
    issues = []

    if query_format != truth_format:
        issues.append(
            f"Query VCF uses {query_format} format (e.g., {sorted(list(query_chroms))[:3]}) "
            f"but truth VCF uses {truth_format} format (e.g., {sorted(list(truth_chroms))[:3]})"
        )

    if query_format != ref_format:
        issues.append(
            f"Query VCF uses {query_format} format "
            f"but reference FASTA uses {ref_format} format (e.g., {sorted(list(ref_chroms))[:3]})"
        )

    if truth_format != ref_format:
        issues.append(
            f"Truth VCF uses {truth_format} format "
            f"but reference FASTA uses {ref_format} format"
        )

    # Check for major chromosome mismatches
    # Normalize chromosome names for comparison (strip 'chr' prefix)
    def _normalize_chroms(chroms: set[str]) -> set[str]:
        """Normalize chromosome names by removing 'chr' prefix."""
        return {c.replace("chr", "") if c.startswith("chr") else c for c in chroms}

    query_norm = _normalize_chroms(query_chroms)
    truth_norm = _normalize_chroms(truth_chroms)
    ref_norm = _normalize_chroms(ref_chroms)

    # Find chromosomes present in VCFs but missing from reference
    query_missing = query_norm - ref_norm
    truth_missing = truth_norm - ref_norm

    if query_missing:
        issues.append(
            f"Query VCF contains chromosomes not in reference: {sorted(query_missing)[:10]}"
        )

    if truth_missing:
        issues.append(
            f"Truth VCF contains chromosomes not in reference: {sorted(truth_missing)[:10]}"
        )

    if issues:
        raise ValidationError(
            "Chromosome naming inconsistency detected:\n" +
            "\n".join(f"  - {issue}" for issue in issues) +
            "\n\nRecommendation: Ensure all files use the same chromosome naming convention "
            "(either all with 'chr' prefix or all without)."
        )

    return {
        "query": query_chroms,
        "truth": truth_chroms,
        "reference": ref_chroms
    }


def validate_sample_consistency(
    vcf_path: str | Path,
    expected_sample: str
) -> None:
    """Verify VCF header contains expected sample name.

    Args:
        vcf_path: Path to VCF file
        expected_sample: Expected sample name from config

    Raises:
        ValidationError: If sample name doesn't match

    Example:
        >>> validate_sample_consistency("sample.vcf.gz", "NA12878")
    """
    path = Path(vcf_path)
    vcf = validate_vcf_format(path)

    try:
        # Get sample names from VCF header
        samples = list(vcf.header.samples) if vcf.header.samples else []

        if not samples:
            raise ValidationError(
                f"VCF file {path.name} has no samples in header. "
                f"Expected sample: '{expected_sample}'"
            )

        if expected_sample not in samples:
            raise ValidationError(
                f"VCF file {path.name} does not contain expected sample '{expected_sample}'. "
                f"Found samples: {samples}"
            )

    finally:
        vcf.close()
