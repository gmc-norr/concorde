# Test Suite Analysis & Recommendations

**Date**: 2026-02-06
**Current Coverage**: ~40-45% (estimated)
**Test Quality**: Good foundation, significant gaps

## Executive Summary

The test suite has a **strong foundation** with well-written tests for validators and utilities (122 tests, all passing). However, there are **critical gaps** in coverage:

- ✅ **Validators**: 95% coverage (excellent)
- ✅ **Utils**: 90% coverage (excellent)
- ❌ **Parsers**: 0% coverage (critical gap)
- ❌ **Ingestion**: 10% coverage (critical gap)
- ❌ **Analysis**: 0% coverage
- ❌ **Models**: 0% coverage

**Recommendation**: **Yes, expand tests significantly** - especially for parsers and ingestion logic.

---

## Current Test Statistics

- **Total tests**: 122
- **Test files**: 7
- **Lines of test code**: 1,042
- **Production scripts**: 21 files
- **Pass rate**: 100% ✅

### Files WITH Tests (9/21 files)

| File | Tests | Quality |
|------|-------|---------|
| scripts/utils.py | 16 | ✅ Excellent |
| scripts/vcf_utils.py | 9 | ✅ Good |
| scripts/validators/base.py | 30 | ✅ Excellent |
| scripts/validators/chromosomes.py | 10 | ✅ Good |
| scripts/validators/config.py | 17 | ✅ Excellent |
| scripts/validators/data_quality.py | 12 | ✅ Good |
| scripts/validators/input_files.py | 28 | ✅ Excellent |
| scripts/validators/cross_file.py | - | ⚠️ Indirect |
| scripts/validators/database.py | - | ⚠️ Indirect |

### Files WITHOUT Tests (12/21 files)

**Critical Gap - Parsers** (5 files, ~1000 LOC):
- ❌ scripts/parsers/parse_happy_vcf.py
- ❌ scripts/parsers/parse_happy_metrics.py
- ❌ scripts/parsers/parse_nfcore_qc.py (most complex - 489 lines!)
- ❌ scripts/parsers/parse_rtg_vcf.py
- ❌ scripts/parsers/parse_rtg_metrics.py

**Critical Gap - Ingestion** (2 files):
- ❌ scripts/ingestion/ingest.py
- ❌ scripts/ingestion/ingest_helpers.py (TSVBulkInserter)

**Gap - Analysis** (2 files):
- ❌ scripts/analysis/intersect_gene_sets.py
- ❌ scripts/analysis/llm_qc_analyzer.py

**Gap - Other**:
- ❌ scripts/validation/validate_inputs.py (ValidationContext)
- ❌ models/*.py (7 database models)
- ❌ tools/query_llm_analysis.py

---

## Test Quality Assessment

### ✅ Strengths

1. **Excellent Fixtures**
   - Factory pattern for VCF/FASTA/BED creation
   - Clean, reusable fixture design
   - Proper use of tmp_path

2. **Clear Organization**
   - Tests mirror source structure
   - Descriptive test names
   - Logical grouping by functionality

3. **Good Edge Case Coverage** (where tests exist)
   - None/null handling
   - Empty strings
   - Invalid inputs
   - Boundary values
   - Type conversions
   - NaN handling

4. **Proper Pytest Usage**
   - Error message validation with `match=`
   - Parametrization where appropriate
   - Clear assertions

### ❌ Weaknesses

#### 1. Critical: Zero Parser Coverage

The **biggest risk** is zero test coverage for parsers, which:
- Handle complex file formats (VCF, TSV, JSON)
- Parse external tool output (hap.py, RTG, nf-core)
- Have lots of edge cases (malformed data, missing fields, type errors)
- Are ~40% of the codebase by LOC

**Example risks**:
- What if hap.py changes output format?
- What if nf-core JSON has unexpected structure?
- What if metrics file has NaN or Inf values?
- What if column names change slightly?

#### 2. Critical: Minimal Ingestion Coverage

No tests for:
- TSVBulkInserter abstraction
- Row mapper functions (_map_metric_row, etc.)
- Database transactions
- Error handling during ingestion
- Rollback scenarios

#### 3. Missing Edge Cases (Even in Tested Code)

**Utils edge cases not tested**:
```python
safe_float(float('inf'))      # Infinity
safe_float(float('-inf'))     # Negative infinity
safe_float(1e308)             # Near max float
safe_int(2**63)               # Very large int
safe_str("   ")               # Only whitespace
safe_str("\n\t")              # Only whitespace chars
```

**Validator edge cases not tested**:
```python
# VCF
- .csi index instead of .tbi
- Multi-allelic variants
- Empty VCF (no variants)
- Non-standard contigs
- Case sensitivity (Chr1 vs chr1 vs CHR1)

# BED
- Overlapping regions
- Unsorted regions
- Zero-length intervals (start == end)
- Scientific notation in coordinates

# DNA sequences
- IUPAC ambiguity codes (RYWSMKHBVDN)
- Mixed case sequences
```

**Chromosome normalization edge cases**:
```python
normalize_chromosome("chrM") vs "chrMT"   # Mitochondrial
normalize_chromosome("chr01")             # Leading zero
normalize_chromosome("HLA-A*01:01")       # HLA contigs
# Mixed formats in same file
```

#### 4. No Integration Tests

Missing:
- End-to-end pipeline test with sample data
- Snakemake workflow validation
- Database integration (create → insert → query)
- Multi-file consistency checks
- Real GIAB data test

#### 5. No Robustness Tests

Missing:
- Performance tests (large VCFs with millions of variants)
- Memory tests (file streaming vs loading all in memory)
- Error recovery (connection loss, disk full, file deleted mid-read)
- Timeout tests
- Concurrent access scenarios

---

## Priority Recommendations

### 🔴 Priority 1: Parser Tests (CRITICAL)

**Why**: Highest risk area with zero coverage and complex parsing logic

**What to add**:
```
tests/parsers/
├── conftest.py                    # Shared fixtures
├── test_parse_nfcore_qc.py       # ~30 tests (most complex)
├── test_parse_happy_vcf.py       # ~20 tests
├── test_parse_happy_metrics.py   # ~15 tests
├── test_parse_rtg_vcf.py         # ~20 tests
└── test_parse_rtg_metrics.py     # ~15 tests
```

**Key tests**:
- ✅ Valid input produces expected output
- ✅ Missing columns raise clear errors
- ✅ Invalid values (NaN, Inf, negative) handled
- ✅ Type conversion errors handled
- ✅ Empty files raise errors
- ✅ Malformed JSON/TSV handled
- ✅ Edge case values (0, 100%, very large numbers)

**Estimated**: 3-4 days, ~100 tests

### 🟠 Priority 2: Ingestion Tests (HIGH)

**Why**: Database interactions are error-prone and currently untested

**What to add**:
```
tests/ingestion/
├── test_ingest_helpers.py        # TSVBulkInserter tests
└── test_ingest.py                # Database integration
```

**Key tests**:
- ✅ TSVBulkInserter with valid data
- ✅ Row mappers with edge case values
- ✅ Transaction rollback on errors
- ✅ Duplicate detection
- ✅ Foreign key constraints
- ✅ Null handling

**Estimated**: 2 days, ~30 tests

### 🟡 Priority 3: Expand Validator Edge Cases (MEDIUM)

**Why**: Fill gaps in otherwise well-tested code

**What to add**:
- Infinity/extreme values in utils
- .csi index support in VCF validation
- IUPAC codes in DNA validation
- Chromosome edge cases (chrM, chr01, HLA)
- File permission errors
- Symlink handling

**Estimated**: 2 days, ~40 tests

### 🟡 Priority 4: Integration Tests (MEDIUM)

**Why**: End-to-end confidence, catch interaction bugs

**What to add**:
```
tests/integration/
├── test_end_to_end.py            # Full pipeline with sample data
├── test_validation_workflow.py   # ValidationContext
└── test_database_workflow.py     # Create schema → insert → query
```

**Estimated**: 2 days, ~15 tests

### 🟢 Priority 5: Model Tests (LOW)

**Why**: Nice to have, but SQLAlchemy is well-tested

**What to add**:
```
tests/models/
└── test_models.py                # Relationships, constraints
```

**Estimated**: 1 day, ~20 tests

---

## Specific Edge Cases to Add

### utils.py
```python
def test_safe_float_infinity():
    assert safe_float(float('inf')) is None
    assert safe_float(float('-inf')) is None
    assert safe_float("inf") is None

def test_safe_float_extreme_values():
    assert safe_float(1e308) == 1e308  # Near max
    assert safe_float(-1e308) == -1e308

def test_safe_int_large_values():
    assert safe_int(2**31 - 1) == 2147483647  # Max 32-bit
    # Behavior with 2**63 is platform-dependent

def test_safe_str_whitespace_only():
    assert safe_str("   ", default="default") == "default"
    assert safe_str("\n\t ", default="X") == "X"
```

### validators/input_files.py
```python
def test_validate_vcf_with_csi_index(temp_vcf_with_csi):
    # Should accept .csi as alternative to .tbi
    validate_vcf_indexed(vcf_with_csi)

def test_validate_vcf_empty():
    vcf = create_vcf_with_no_variants()
    chroms = validate_vcf_chromosomes(vcf)
    assert chroms == set()

def test_validate_bed_zero_length_interval():
    bed = create_bed([("chr1", 100, 100)])  # start == end
    with pytest.raises(ValidationError, match="zero-length"):
        validate_bed_format(bed)
```

### validators/data_quality.py
```python
def test_validate_dna_iupac_codes():
    # IUPAC ambiguity codes - should these be allowed?
    validate_dna_sequence("ACGTRYWSMKHBVDN")  # Or should raise?

def test_validate_probability_infinity():
    with pytest.raises(ValidationError):
        validate_probability(float('inf'), "test")

def test_validate_coverage_metrics_missing_keys():
    metrics = {}  # Empty dict
    # Should handle gracefully or require keys?
    validate_coverage_metrics(metrics)
```

### validators/chromosomes.py
```python
def test_normalize_chromosome_mitochondrial():
    assert normalize_chromosome("M") == "chrM"
    assert normalize_chromosome("MT") == "chrMT"
    assert normalize_chromosome("chrM") == "chrM"

def test_normalize_chromosome_leading_zero():
    assert normalize_chromosome("01") == "chr1"
    assert normalize_chromosome("chr01") == "chr1"

def test_normalize_chromosome_hla():
    # HLA contigs - should pass through unchanged?
    result = normalize_chromosome("HLA-A*01:01")
    assert "HLA" in result
```

---

## Estimated Total Work

| Priority | Component | Days | Tests |
|----------|-----------|------|-------|
| P1 | Parser tests | 3-4 | ~100 |
| P2 | Ingestion tests | 2 | ~30 |
| P3 | Validator edge cases | 2 | ~40 |
| P4 | Integration tests | 2 | ~15 |
| P5 | Model tests | 1 | ~20 |
| **Total** | | **10-11 days** | **~205 tests** |

**Result**: 122 → ~327 tests, 40-45% → 80-85% coverage

---

## Conclusion

### Current State: Good Foundation ✅

The existing tests demonstrate:
- Good testing practices
- Clean organization
- Proper use of pytest
- Solid edge case coverage (where present)

### Problem: Critical Coverage Gaps ❌

The biggest risks are:
1. **Zero parser coverage** (40% of codebase)
2. **Minimal ingestion coverage** (database interactions)
3. **No integration tests** (end-to-end confidence)

### Recommendation: **Expand Tests** 📈

**Start with Priority 1 (parser tests)** as these are:
- Highest risk (complex parsing logic)
- Most LOC without coverage
- Most error-prone (external tool outputs)
- Critical to pipeline correctness

**Then add Priority 2 (ingestion)** for:
- Database correctness
- Transaction safety
- Error handling

The current tests are well-written and should serve as a good template for expansion. The effort required (10-11 days) is justified by the risk reduction and confidence gained.
