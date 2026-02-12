# Testing Strategy for Concorde Pipeline

**Date**: 2026-02-06
**Current Coverage**: 133 tests, ~45-50% estimated coverage
**Goal**: Achieve 80-85% coverage with 300+ tests

## Overview

This document outlines the testing strategy for the Concorde pipeline, with specific guidance for testing Snakemake scripts and other challenging components.

## Current Status

### ✅ Well Tested (9 files, 133 tests)
- **Validators** (95% coverage): base, chromosomes, config, data_quality, input_files, cross_file
- **Utils** (90% coverage): utils.py, vcf_utils.py
- **Test Infrastructure**: Excellent fixtures for VCF/FASTA/BED creation

### ❌ Untested (12 files, 0 tests)
- **Parsers** (5 files): All Snakemake scripts with zero tests
- **Ingestion** (2 files): Database interaction scripts
- **Analysis** (2 files): Gene set intersection, LLM analysis
- **Models** (7 files): SQLAlchemy database models
- **Tools** (1 file): query_llm_analysis.py

## Testing Challenges

### Challenge 1: Snakemake Scripts

**Problem**: Snakemake scripts receive the `snakemake` object at runtime via injection. This object contains:
- `snakemake.input` - Input file paths
- `snakemake.output` - Output file paths
- `snakemake.params` - Parameters
- `snakemake.log` - Log file path
- `snakemake.config` - Configuration dictionary

**Example**:
```python
# In parse_happy_vcf.py
vcf_path = snakemake.input.vcf
output_path = snakemake.output.tsv
sample = snakemake.params.sample
```

**Solutions**:

#### Option A: Mock the Snakemake Object (Recommended for Unit Tests)
```python
# tests/parsers/test_parse_happy_vcf.py
from unittest.mock import Mock
import sys

def test_parse_happy_vcf_basic(tmp_path, temp_vcf_file):
    """Test parsing of hap.py VCF output."""
    # Create test VCF
    vcf_path = temp_vcf_file([
        {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T",
         "info": {"BD": "TP", "BK": ".", "type": "SNP"}}
    ])

    output_path = tmp_path / "output.tsv"
    log_path = tmp_path / "test.log"

    # Mock snakemake object
    mock_snakemake = Mock()
    mock_snakemake.input.vcf = str(vcf_path)
    mock_snakemake.output.tsv = str(output_path)
    mock_snakemake.params.run_id = 1
    mock_snakemake.params.sample = "NA12878"
    mock_snakemake.log = [str(log_path)]

    # Inject mock into module namespace
    sys.modules['__main__'].snakemake = mock_snakemake

    # Import and run parser
    import scripts.parsers.parse_happy_vcf as parser
    # Parser runs on import due to Snakemake pattern

    # Verify output
    import pandas as pd
    df = pd.read_csv(output_path, sep='\t')
    assert len(df) == 1
    assert df.iloc[0]['chrom'] == 'chr1'
    assert df.iloc[0]['type'] == 'SNP'
```

**Challenges with Option A**:
- Snakemake scripts execute on import (not ideal for testing)
- Need to mock the entire snakemake object structure
- Difficult to isolate specific functions

#### Option B: Refactor for Testability (Recommended for New Code)
```python
# In parse_happy_vcf.py - refactored version
def parse_happy_vcf_file(vcf_path: str, run_id: int, sample: str) -> pd.DataFrame:
    """Parse hap.py VCF and return DataFrame.

    This function is testable without Snakemake.
    """
    variants = []
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            variant = {
                "run_id": run_id,
                "sample": sample,
                "chrom": record.chrom,
                "pos": record.pos,
                "ref": record.ref,
                "alt": record.alts[0],
                "type": record.info.get("type", "."),
                "bd": record.info.get("BD", "."),
            }
            variants.append(variant)
    return pd.DataFrame(variants)


# Snakemake wrapper (minimal, hard to break)
if __name__ == "__main__" or "snakemake" in globals():
    df = parse_happy_vcf_file(
        vcf_path=snakemake.input.vcf,
        run_id=snakemake.params.run_id,
        sample=snakemake.params.sample
    )
    df.to_csv(snakemake.output.tsv, sep='\t', index=False)
```

```python
# tests/parsers/test_parse_happy_vcf.py - much simpler!
def test_parse_happy_vcf_basic(temp_vcf_file):
    """Test parsing of hap.py VCF output."""
    vcf_path = temp_vcf_file([
        {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}
    ])

    from scripts.parsers.parse_happy_vcf import parse_happy_vcf_file
    df = parse_happy_vcf_file(str(vcf_path), run_id=1, sample="NA12878")

    assert len(df) == 1
    assert df.iloc[0]['chrom'] == 'chr1'
```

**Benefits of Option B**:
- Pure functions easy to test
- No mocking required
- Better separation of concerns
- Snakemake wrapper is trivial (unlikely to break)

#### Option C: Integration Tests with Real Snakemake
```bash
# Run actual Snakemake workflow with test data
snakemake --cores 1 --configfile tests/data/test_config.yaml -n
```

**Use for**: End-to-end validation, not unit tests

### Challenge 2: Database Interactions

**Problem**: SQLAlchemy sessions, transactions, and database state

**Solutions**:

#### Use In-Memory SQLite for Tests
```python
# tests/ingestion/conftest.py
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from models.base import Base

@pytest.fixture
def test_db():
    """Create in-memory test database."""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    yield session
    session.close()

# tests/ingestion/test_ingest_helpers.py
def test_tsv_bulk_inserter(test_db, temp_tsv_file):
    """Test TSVBulkInserter with metrics data."""
    from scripts.ingestion.ingest_helpers import TSVBulkInserter
    from models import Metric

    tsv_path = temp_tsv_file(
        headers=["run_id", "type", "filter", "recall"],
        rows=[["1", "SNP", "PASS", "0.95"]]
    )

    inserter = TSVBulkInserter(test_db, run_id=1)
    count = inserter.insert_from_tsv(
        tsv_path=tsv_path,
        model_class=Metric,
        required_columns=["run_id", "type", "filter", "recall"],
        row_mapper=lambda row, run_id: {
            "run_id": run_id,
            "type": row["type"],
            "filter": row["filter"],
            "recall": float(row["recall"])
        },
        entity_name="metrics"
    )

    assert count == 1
    test_db.commit()

    # Verify data in database
    metrics = test_db.query(Metric).all()
    assert len(metrics) == 1
    assert metrics[0].type == "SNP"
    assert metrics[0].recall == 0.95
```

### Challenge 3: LLM Analysis Scripts

**Problem**: Dependency on external Ollama service

**Solutions**:

#### Mock the Ollama Client
```python
# tests/analysis/test_llm_qc_analyzer.py
from unittest.mock import Mock, patch

def test_llm_analysis_integration(test_db, sample_qc_data):
    """Test LLM analysis with mocked Ollama."""
    mock_response = {
        "response": "Overall Quality: Good\n\nCoverage is excellent at 35X..."
    }

    with patch('ollama.chat') as mock_chat:
        mock_chat.return_value = mock_response

        from scripts.analysis.llm_qc_analyzer import analyze_qc_metrics
        result = analyze_qc_metrics(sample_qc_data, model="llama3.2:3b")

        assert "Overall Quality: Good" in result
        mock_chat.assert_called_once()
```

## Implementation Roadmap

### Phase 1: Refactor Parsers for Testability (3-4 days)

1. **Extract pure functions from Snakemake scripts**
   - parse_happy_vcf.py → extract parse_happy_vcf_file()
   - parse_happy_metrics.py → extract parse_happy_metrics_file()
   - parse_nfcore_qc.py → extract parse_multiqc_json()
   - parse_rtg_vcf.py → extract parse_rtg_vcf_file()
   - parse_rtg_metrics.py → extract parse_rtg_metrics_file()

2. **Keep Snakemake wrappers minimal**
   ```python
   if __name__ == "__main__" or "snakemake" in globals():
       df = parse_function(snakemake.input.file, snakemake.params.X)
       df.to_csv(snakemake.output.tsv, sep='\t', index=False)
   ```

3. **Create comprehensive parser tests**
   - Valid inputs
   - Missing columns
   - Invalid values (NaN, Inf, negative)
   - Type errors
   - Empty files
   - Malformed data

### Phase 2: Add Ingestion Tests (2 days)

1. **Test TSVBulkInserter abstraction**
   - Valid data insertion
   - Row mapper functions
   - Error handling
   - Transaction rollback

2. **Test database constraints**
   - Foreign keys
   - Unique constraints
   - Not null constraints
   - Cascade deletes

### Phase 3: Expand Validator Edge Cases (2 days)

**Already added (11 tests)**:
- ✅ Infinity values in safe_float
- ✅ Extreme float values
- ✅ Large integers
- ✅ IUPAC ambiguity codes
- ✅ Probability infinity handling
- ✅ Extreme percentage boundaries
- ✅ Empty/partial coverage metrics
- ✅ Mitochondrial chromosome variants
- ✅ Leading zeros in chromosomes
- ✅ Non-standard contigs

**Still needed**:
- .csi index support
- Multi-allelic variants
- Empty VCF files
- Symlink handling
- File permission errors
- Mixed line endings

### Phase 4: Integration Tests (2 days)

1. **End-to-end pipeline test**
   ```python
   def test_full_pipeline_small_sample(tmp_path):
       """Test complete pipeline with small sample data."""
       # Setup test data (10 variants)
       # Run Snakemake workflow
       # Verify database contents
       # Check all expected files created
   ```

2. **Validation workflow test**
   ```python
   def test_validation_collects_all_errors():
       """Test ValidationContext error collection."""
       # Create config with multiple errors
       # Run validate_inputs.py
       # Verify all errors collected and reported
   ```

### Phase 5: Model Tests (1 day)

Test SQLAlchemy models:
- Relationships (one-to-many, many-to-many)
- Cascade deletes
- Constraints
- Default values

## Testing Best Practices

### 1. Test File Organization

Mirror source structure:
```
tests/
├── parsers/
│   ├── conftest.py           # Parser-specific fixtures
│   ├── test_parse_happy_vcf.py
│   ├── test_parse_happy_metrics.py
│   └── ...
├── ingestion/
│   ├── conftest.py           # Database fixtures
│   ├── test_ingest_helpers.py
│   └── test_ingest.py
├── analysis/
│   ├── test_intersect_gene_sets.py
│   └── test_llm_qc_analyzer.py
├── integration/
│   ├── test_end_to_end.py
│   └── test_validation_workflow.py
└── validators/               # Already well-tested
```

### 2. Fixture Design

**Good**: Factory fixtures that create test data on demand
```python
@pytest.fixture
def temp_vcf_file(tmp_path):
    def _create(variants):
        # Create VCF with specified variants
        return vcf_path
    return _create
```

**Bad**: Fixed test data that doesn't adapt to test needs

### 3. Test Naming

**Good**: `test_parse_happy_vcf_with_missing_column_raises_error()`
**Bad**: `test_parser()`, `test_case1()`

### 4. Error Message Validation

Always validate error messages:
```python
with pytest.raises(ValidationError, match="column 'recall' is required"):
    validate_tsv_schema(tsv_path, required_columns=["recall"])
```

### 5. Parametrize Similar Tests

```python
@pytest.mark.parametrize("invalid_value,error_msg", [
    (float('inf'), "infinity"),
    (float('-inf'), "infinity"),
    (float('nan'), "NaN"),
    (150, "between 0 and 100"),
])
def test_validate_percentage_invalid(invalid_value, error_msg):
    with pytest.raises(ValidationError, match=error_msg):
        validate_percentage(invalid_value, "test")
```

## Coverage Goals

| Component | Current | Target | Priority |
|-----------|---------|--------|----------|
| Validators | 95% | 98% | Low |
| Utils | 90% | 95% | Low |
| Parsers | 0% | 85% | **HIGH** |
| Ingestion | 10% | 80% | **HIGH** |
| Analysis | 0% | 70% | Medium |
| Models | 0% | 75% | Medium |
| Integration | 0% | 60% | Medium |
| **Overall** | **~45%** | **80-85%** | |

## Success Metrics

### Quantitative
- ✅ 133 → 300+ tests (2.3x increase)
- ✅ 45% → 80% coverage (1.8x increase)
- ✅ 0 parser tests → 100 parser tests
- ✅ All critical paths tested

### Qualitative
- ✅ Tests catch bugs before production
- ✅ Refactored code is more testable
- ✅ Clear testing patterns for future development
- ✅ Fast test execution (<10 seconds for full suite)
- ✅ Reliable (no flaky tests)

## Conclusion

The testing strategy prioritizes:
1. **Refactoring for testability** - Extract pure functions from Snakemake scripts
2. **Critical path coverage** - Test parsers and ingestion first (highest risk)
3. **Practical approaches** - Use in-memory databases, mock external services
4. **Incremental progress** - Can be implemented in phases

**Next Steps**:
1. Start with Phase 1 (refactor parsers for testability)
2. Add comprehensive parser tests
3. Measure coverage improvement
4. Continue through phases 2-5

**Estimated Total Effort**: 10-12 days to reach 80% coverage with ~300 tests
