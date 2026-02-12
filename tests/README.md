# Pipeline Tests

This directory contains unit tests for pipeline scripts.

## Test Files

- **test_utils.py** - Tests for utility functions (safe_float, safe_int, safe_str)
- **test_validators.py** - Tests for validation functions (TSV schema validation, enum validation, etc.)
- **test_vcf_utils.py** - Tests for VCF-specific utilities (safe_format, query annotation building)

## Running Tests

### Using pixi (recommended)

After installing pytest dependencies:
```bash
pixi run test-pipeline
```

### Using pytest directly

```bash
python -m pytest pipeline/tests/ -v
```

### Running specific test files

```bash
python -m pytest pipeline/tests/test_utils.py -v
python -m pytest pipeline/tests/test_validators.py -v
python -m pytest pipeline/tests/test_vcf_utils.py -v
```

### Running specific test classes

```bash
python -m pytest pipeline/tests/test_utils.py::TestSafeFloat -v
```

## Test Coverage

The test suite covers:

### utils.py (test_utils.py)
- ✅ `safe_float()` - 5 test cases
  - Valid float conversion
  - None handling with defaults
  - Invalid strings return defaults
  - NaN detection
  - Empty string handling

- ✅ `safe_int()` - 6 test cases
  - Valid int conversion
  - None handling with defaults
  - Invalid strings return defaults
  - Empty string handling
  - Float string truncation

- ✅ `safe_str()` - 6 test cases
  - Valid string conversion
  - Number to string conversion
  - None handling with defaults
  - Boolean to string conversion
  - Whitespace stripping
  - Empty string after strip

### validators.py (test_validators.py)
- ✅ `validate_required_columns()` - 3 test cases
- ✅ `validate_file_exists()` - 3 test cases
- ✅ `validate_enum_values()` - 4 test cases
- ✅ `validate_non_empty()` - 3 test cases
- ✅ `validate_positive_numeric()` - 5 test cases
- ✅ `validate_tsv_schema()` - 7 test cases (comprehensive integration tests)

### vcf_utils.py (test_vcf_utils.py)
- ✅ `safe_format()` - 10 test cases
  - Simple value extraction
  - Missing fields with defaults
  - Tuple unpacking (single and multiple values)
  - NaN detection
  - None value handling
  - String, float, and complex type handling

## Total Test Coverage

- **Total test files**: 3
- **Total test classes**: 12
- **Total test cases**: ~45

All tests use pytest fixtures and follow best practices:
- Descriptive test names
- Comprehensive edge case coverage
- Proper use of pytest.raises for error testing
- Temporary file handling for file I/O tests
- Mock objects for external dependencies
