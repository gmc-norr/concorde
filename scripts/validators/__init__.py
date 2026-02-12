"""Data validation utilities for Concorde pipeline.

This package provides comprehensive validation functions to ensure data integrity
before and during pipeline processing, preventing silent failures and improving
error messages.

Modules:
    base: Core validators (columns, files, enums, numeric, duplicates, nulls)
    input_files: VCF, FASTA, BED format validation
    config: Configuration structure validation
    data_quality: DNA sequences, probabilities, percentages, QC metrics
    cross_file: Chromosome consistency, sample matching
    chromosomes: Chromosome normalization utilities
    database: Database integrity checks

Example:
    >>> from validators import ValidationError, validate_vcf_indexed
    >>> try:
    ...     validate_vcf_indexed("sample.vcf.gz")
    ... except ValidationError as e:
    ...     print(f"Validation failed: {e}")
"""

# Import all validators for backward compatibility
from .base import (
    ValidationError,
    validate_required_columns,
    validate_file_exists,
    validate_enum_values,
    validate_non_empty,
    validate_positive_numeric,
    validate_variant_coordinates,
    validate_no_duplicate_variants,
    validate_no_nulls,
    validate_tsv_schema,
)

from .input_files import (
    validate_vcf_indexed,
    validate_vcf_format,
    validate_vcf_chromosomes,
    validate_reference_indexed,
    validate_bed_format,
)

from .config import (
    validate_comparison_tool,
    validate_decomposition_modes,
    validate_gene_sets_config,
    validate_vep_config,
    validate_config,
)

from .data_quality import (
    validate_dna_sequence,
    validate_probability,
    validate_percentage,
    validate_coverage_metrics,
)

from .cross_file import (
    validate_chromosome_consistency,
    validate_sample_consistency,
)

from .chromosomes import (
    normalize_chromosome,
    get_chromosome_format,
    validate_chromosome_format,
    validate_standard_chromosomes,
    STANDARD_CHROMOSOMES_WITH_CHR,
    STANDARD_CHROMOSOMES_WITHOUT_CHR,
    ALL_STANDARD_CHROMOSOMES,
)

from .database import (
    check_duplicate_run,
    validate_no_duplicate_run,
)

# Export all validators
__all__ = [
    # Exception
    "ValidationError",
    # Base validators
    "validate_required_columns",
    "validate_file_exists",
    "validate_enum_values",
    "validate_non_empty",
    "validate_positive_numeric",
    "validate_variant_coordinates",
    "validate_no_duplicate_variants",
    "validate_no_nulls",
    "validate_tsv_schema",
    # Input file validators
    "validate_vcf_indexed",
    "validate_vcf_format",
    "validate_vcf_chromosomes",
    "validate_reference_indexed",
    "validate_bed_format",
    # Config validators
    "validate_comparison_tool",
    "validate_decomposition_modes",
    "validate_gene_sets_config",
    "validate_vep_config",
    "validate_config",
    # Data quality validators
    "validate_dna_sequence",
    "validate_probability",
    "validate_percentage",
    "validate_coverage_metrics",
    # Cross-file validators
    "validate_chromosome_consistency",
    "validate_sample_consistency",
    # Chromosome utilities
    "normalize_chromosome",
    "get_chromosome_format",
    "validate_chromosome_format",
    "validate_standard_chromosomes",
    "STANDARD_CHROMOSOMES_WITH_CHR",
    "STANDARD_CHROMOSOMES_WITHOUT_CHR",
    "ALL_STANDARD_CHROMOSOMES",
    # Database validators
    "check_duplicate_run",
    "validate_no_duplicate_run",
]
