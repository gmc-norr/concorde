# Validator Usage Guide

Quick reference for using Concorde validators in your scripts.

## Basic Usage

```python
from validators import ValidationError, validate_vcf_indexed

try:
    validate_vcf_indexed("sample.vcf.gz")
    print("✓ VCF is properly indexed")
except ValidationError as e:
    print(f"✗ Validation failed: {e}")
```

## Available Validators

### Input File Validators

```python
from validators import (
    validate_vcf_indexed,
    validate_vcf_format,
    validate_vcf_chromosomes,
    validate_reference_indexed,
    validate_bed_format,
)

# Check VCF is bgzipped with tabix index
validate_vcf_indexed("sample.vcf.gz")

# Open and validate VCF format
vcf = validate_vcf_format("sample.vcf.gz")
vcf.close()

# Get chromosomes from VCF
chroms = validate_vcf_chromosomes("sample.vcf.gz")
print(f"Found chromosomes: {chroms}")

# Validate against expected chromosomes
validate_vcf_chromosomes("sample.vcf.gz", expected_chroms={"chr1", "chr2"})

# Check FASTA has .fai index
validate_reference_indexed("GRCh38.fa")

# Validate BED file format
validate_bed_format("regions.bed")
```

### Configuration Validators

```python
from validators import validate_config, validate_comparison_tool

# Validate entire config
config = {
    "sample": "NA12878",
    "caller": "GATK",
    "pipeline_version": "v1.0",
    "comparison_tool": "happy",
    "reference": "/path/to/ref.fa",
    "truth_vcf": "/path/to/truth.vcf.gz",
    "query_vcfs": ["/path/to/query.vcf.gz"],
    "decomposition_modes": ["decomposed"],
    "database": "data/concorde.db",
    "output_dir": "results"
}
validate_config(config)

# Validate specific fields
validate_comparison_tool("happy")  # Must be "happy" or "rtg"
```

### Data Quality Validators

```python
from validators import (
    validate_dna_sequence,
    validate_probability,
    validate_percentage,
    validate_coverage_metrics,
)

# Validate DNA sequence (ACGTN only)
validate_dna_sequence("ACGTACGT", "ref_allele")

# Validate probability (0.0 - 1.0)
validate_probability(0.95, "allele_frequency")

# Validate percentage (0 - 100)
validate_percentage(98.5, "mapping_rate")

# Validate QC metrics
qc_metrics = {
    "mean_coverage": 30.5,
    "median_coverage": 28.0,
    "pct_bases_10x": 99.1,
    "mapping_rate": 99.8
}
validate_coverage_metrics(qc_metrics)
```

### DataFrame Validators

```python
from validators import (
    validate_required_columns,
    validate_enum_values,
    validate_positive_numeric,
    validate_no_duplicate_variants,
    validate_no_nulls,
    validate_tsv_schema,
)
import pandas as pd

df = pd.DataFrame({
    "chrom": ["chr1", "chr2"],
    "pos": [100, 200],
    "ref": ["A", "G"],
    "alt": ["T", "C"],
    "type": ["SNP", "SNP"]
})

# Check required columns
validate_required_columns(df, ["chrom", "pos", "ref", "alt"], "variants.tsv")

# Validate enum values
validate_enum_values(df, "type", {"SNP", "INDEL"}, "variants.tsv")

# Validate numeric columns are non-negative
validate_positive_numeric(df, ["pos"], "variants.tsv")

# Check for duplicate variants
validate_no_duplicate_variants(df, "variants.tsv")

# Ensure required fields have no nulls
validate_no_nulls(df, ["chrom", "pos", "ref", "alt"], "variants.tsv")

# Comprehensive TSV validation
df = validate_tsv_schema(
    "variants.tsv",
    required_columns=["chrom", "pos", "ref", "alt", "type"],
    enum_columns={"type": {"SNP", "INDEL"}},
    positive_columns=["pos"]
)
```

### Cross-File Validators

```python
from validators import (
    validate_chromosome_consistency,
    validate_sample_consistency,
)

# Check chromosome naming consistency
chroms = validate_chromosome_consistency(
    "query.vcf.gz",
    "truth.vcf.gz",
    "reference.fa"
)
print(f"Query chroms: {chroms['query']}")
print(f"Truth chroms: {chroms['truth']}")
print(f"Reference chroms: {chroms['reference']}")

# Verify VCF sample matches expected
validate_sample_consistency("sample.vcf.gz", "NA12878")
```

### Chromosome Utilities

```python
from validators import (
    normalize_chromosome,
    get_chromosome_format,
    validate_chromosome_format,
)

# Normalize chromosome name
print(normalize_chromosome("1"))      # "chr1"
print(normalize_chromosome("chr1"))   # "chr1"
print(normalize_chromosome("MT"))     # "chrM"

# Detect chromosome format
chroms = {"chr1", "chr2", "chrX"}
fmt = get_chromosome_format(chroms)
print(fmt)  # "with_chr"

# Validate DataFrame chromosome format
df = pd.DataFrame({"chrom": ["chr1", "chr2"]})
validate_chromosome_format(df, "with_chr", "variants.tsv")
```

## Error Handling

### ValidationError Exception

All validators raise `ValidationError` on failure:

```python
from validators import ValidationError, validate_vcf_indexed

try:
    validate_vcf_indexed("missing.vcf.gz")
except ValidationError as e:
    # Error message includes:
    # - File/context that failed
    # - Specific problem
    # - Expected values or fix instructions
    print(f"Validation error: {e}")
    # Example output:
    # VCF file missing.vcf.gz not found: /path/to/missing.vcf.gz
```

### Collect-All-Errors Pattern

For validating multiple inputs, collect errors instead of failing fast:

```python
from validators import ValidationError

class ValidationContext:
    def __init__(self):
        self.errors = []

    def validate(self, func, *args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ValidationError as e:
            self.errors.append(e)
            return None

    def raise_if_errors(self):
        if self.errors:
            msg = "\n".join(f"  - {e}" for e in self.errors)
            raise ValidationError(f"Validation failed:\n{msg}")

# Usage
ctx = ValidationContext()
ctx.validate(validate_vcf_indexed, "file1.vcf.gz")
ctx.validate(validate_vcf_indexed, "file2.vcf.gz")
ctx.validate(validate_reference_indexed, "ref.fa")
ctx.raise_if_errors()  # Raises with all errors at once
```

## Best Practices

### 1. Validate Early

Validate inputs at the start of your script:

```python
def main():
    # Validate first
    validate_vcf_indexed(args.vcf)
    validate_reference_indexed(args.reference)

    # Then process
    process_vcf(args.vcf)
```

### 2. Use TSV Schema Validation

For TSV files, use `validate_tsv_schema` for comprehensive checking:

```python
# Instead of multiple separate checks
df = pd.read_csv("variants.tsv", sep="\t")
validate_required_columns(df, ["chrom", "pos"], "variants.tsv")
validate_enum_values(df, "type", {"SNP"}, "variants.tsv")
# ... etc

# Use single comprehensive validation
df = validate_tsv_schema(
    "variants.tsv",
    required_columns=["chrom", "pos", "type"],
    enum_columns={"type": {"SNP", "INDEL"}},
    positive_columns=["pos"]
)
```

### 3. Provide Context in Error Messages

Always pass descriptive file names to validators:

```python
# Good
validate_required_columns(df, ["col1"], "query_variants.tsv")

# Bad
validate_required_columns(df, ["col1"], "file.tsv")
```

### 4. Handle Both Paths and Strings

Validators accept both `str` and `Path` objects:

```python
from pathlib import Path

# Both work
validate_vcf_indexed("sample.vcf.gz")
validate_vcf_indexed(Path("sample.vcf.gz"))
```

## Examples by Use Case

### Validating Pipeline Inputs

```python
from validators import (
    validate_config,
    validate_vcf_indexed,
    validate_reference_indexed,
    validate_chromosome_consistency,
)

def validate_pipeline_inputs(config_path):
    # Load config
    import yaml
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Validate config structure
    validate_config(config)

    # Validate input files
    validate_vcf_indexed(config["truth_vcf"])
    validate_reference_indexed(config["reference"])

    for query_vcf in config["query_vcfs"]:
        validate_vcf_indexed(query_vcf)

    # Validate consistency
    for query_vcf in config["query_vcfs"]:
        validate_chromosome_consistency(
            query_vcf,
            config["truth_vcf"],
            config["reference"]
        )
```

### Validating Parsed Data

```python
from validators import validate_tsv_schema, validate_no_duplicates

def load_variants(file_path):
    # Validate and load in one step
    df = validate_tsv_schema(
        file_path,
        required_columns=["chrom", "pos", "ref", "alt", "type"],
        enum_columns={"type": {"SNP", "INDEL", "COMPLEX"}},
        positive_columns=["pos", "dp", "gq"]
    )

    # Additional checks
    validate_no_duplicate_variants(df, file_path)

    return df
```

### Custom Validation Functions

Build on existing validators:

```python
from validators import ValidationError, validate_dna_sequence

def validate_variant_alleles(ref, alt, field_prefix=""):
    """Validate reference and alternate alleles."""
    validate_dna_sequence(ref, f"{field_prefix}ref")
    validate_dna_sequence(alt, f"{field_prefix}alt")

    # Custom validation
    if ref == alt:
        raise ValidationError(
            f"{field_prefix}ref and alt alleles cannot be identical: {ref}"
        )

    if len(ref) == 0 or len(alt) == 0:
        raise ValidationError(
            f"{field_prefix}alleles cannot be empty"
        )
```

## See Also

- [VALIDATION_IMPLEMENTATION_SUMMARY.md](archive/VALIDATION_IMPLEMENTATION_SUMMARY.md) - Complete implementation details
- [Getting Started](getting-started.md) - Pipeline usage with validation
- [test_validators/](../tests/validators/) - Test examples for all validators
