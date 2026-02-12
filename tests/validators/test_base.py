"""Tests for base validator functions."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from validators import (
    ValidationError,
    validate_enum_values,
    validate_file_exists,
    validate_non_empty,
    validate_positive_numeric,
    validate_required_columns,
    validate_tsv_schema,
    validate_no_duplicate_variants,
    validate_no_nulls,
)


class TestValidateRequiredColumns:
    """Tests for validate_required_columns function."""

    def test_all_columns_present(self):
        """Test validation passes when all columns present."""
        df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
        validate_required_columns(df, ["a", "b"], "test.tsv")  # Should not raise

    def test_missing_columns_raises(self):
        """Test validation fails when columns missing."""
        df = pd.DataFrame({"a": [1], "b": [2]})
        with pytest.raises(ValidationError, match="missing required columns"):
            validate_required_columns(df, ["a", "b", "c"], "test.tsv")

    def test_error_message_contains_file_name(self):
        """Test error message includes file name."""
        df = pd.DataFrame({"a": [1]})
        with pytest.raises(ValidationError, match="test.tsv"):
            validate_required_columns(df, ["a", "b"], "test.tsv")


class TestValidateFileExists:
    """Tests for validate_file_exists function."""

    def test_existing_file_returns_path(self):
        """Test validation passes for existing file."""
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            result = validate_file_exists(tmp_path, "test file")
            assert result == tmp_path
        finally:
            tmp_path.unlink()

    def test_missing_file_raises(self):
        """Test validation fails for missing file."""
        with pytest.raises(ValidationError, match="not found"):
            validate_file_exists("/nonexistent/file.txt", "test file")

    def test_directory_raises(self):
        """Test validation fails for directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValidationError, match="is not a file"):
                validate_file_exists(tmpdir, "test file")


class TestValidateEnumValues:
    """Tests for validate_enum_values function."""

    def test_valid_enum_values(self):
        """Test validation passes for valid enum values."""
        df = pd.DataFrame({"type": ["SNP", "INDEL", "SNP"]})
        validate_enum_values(
            df, "type", {"SNP", "INDEL", "COMPLEX"}, "test.tsv"
        )  # Should not raise

    def test_invalid_enum_values_raises(self):
        """Test validation fails for invalid enum values."""
        df = pd.DataFrame({"type": ["SNP", "INVALID", "INDEL"]})
        with pytest.raises(ValidationError, match="invalid values"):
            validate_enum_values(df, "type", {"SNP", "INDEL"}, "test.tsv")

    def test_missing_column_is_ignored(self):
        """Test validation passes if column doesn't exist."""
        df = pd.DataFrame({"other": [1, 2, 3]})
        validate_enum_values(
            df, "type", {"SNP", "INDEL"}, "test.tsv"
        )  # Should not raise

    def test_null_values_are_ignored(self):
        """Test that null/NA values are ignored."""
        df = pd.DataFrame({"type": ["SNP", None, "INDEL", pd.NA]})
        validate_enum_values(
            df, "type", {"SNP", "INDEL"}, "test.tsv"
        )  # Should not raise


class TestValidateNonEmpty:
    """Tests for validate_non_empty function."""

    def test_non_empty_dataframe(self):
        """Test validation passes for non-empty dataframe."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        validate_non_empty(df, "test.tsv")  # Should not raise

    def test_empty_dataframe_raises(self):
        """Test validation fails for empty dataframe."""
        df = pd.DataFrame()
        with pytest.raises(ValidationError, match="contains no data rows"):
            validate_non_empty(df, "test.tsv")

    def test_dataframe_with_columns_but_no_rows_raises(self):
        """Test validation fails for dataframe with columns but no data."""
        df = pd.DataFrame(columns=["a", "b", "c"])
        with pytest.raises(ValidationError, match="contains no data rows"):
            validate_non_empty(df, "test.tsv")


class TestValidatePositiveNumeric:
    """Tests for validate_positive_numeric function."""

    def test_all_positive_values(self):
        """Test validation passes for all positive values."""
        df = pd.DataFrame({"pos": [1, 10, 100], "dp": [30, 40, 50]})
        validate_positive_numeric(df, ["pos", "dp"], "test.tsv")  # Should not raise

    def test_zero_is_allowed(self):
        """Test that zero values are allowed."""
        df = pd.DataFrame({"value": [0, 1, 2]})
        validate_positive_numeric(df, ["value"], "test.tsv")  # Should not raise

    def test_negative_values_raise(self):
        """Test validation fails for negative values."""
        df = pd.DataFrame({"value": [1, -5, 10]})
        with pytest.raises(ValidationError, match="negative values"):
            validate_positive_numeric(df, ["value"], "test.tsv")

    def test_missing_column_is_ignored(self):
        """Test validation passes if column doesn't exist."""
        df = pd.DataFrame({"other": [1, 2, 3]})
        validate_positive_numeric(df, ["missing"], "test.tsv")  # Should not raise

    def test_null_values_are_ignored(self):
        """Test that null/NA values are ignored."""
        df = pd.DataFrame({"value": [1, None, 3, pd.NA]})
        validate_positive_numeric(df, ["value"], "test.tsv")  # Should not raise


class TestValidateNoDuplicateVariants:
    """Tests for validate_no_duplicate_variants function."""

    def test_no_duplicates_passes(self):
        """Test validation passes when no duplicates exist."""
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2"],
            "pos": [100, 200, 100],
            "ref": ["A", "A", "G"],
            "alt": ["T", "G", "C"]
        })
        validate_no_duplicate_variants(df, "test.tsv")  # Should not raise

    def test_duplicate_variants_raises(self):
        """Test validation fails when duplicate variants exist."""
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2"],
            "pos": [100, 100, 200],
            "ref": ["A", "A", "G"],
            "alt": ["T", "T", "C"]
        })
        with pytest.raises(ValidationError, match="duplicate variant"):
            validate_no_duplicate_variants(df, "test.tsv")

    def test_multiple_duplicates_raises(self):
        """Test validation fails with multiple duplicate groups."""
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "pos": [100, 100, 200, 200],
            "ref": ["A", "A", "G", "G"],
            "alt": ["T", "T", "C", "C"]
        })
        with pytest.raises(ValidationError, match="2 duplicate variant"):
            validate_no_duplicate_variants(df, "test.tsv")

    def test_missing_columns_skips_validation(self):
        """Test validation is skipped if required columns missing."""
        df = pd.DataFrame({"other": [1, 2, 3]})
        validate_no_duplicate_variants(df, "test.tsv")  # Should not raise


class TestValidateNoNulls:
    """Tests for validate_no_nulls function."""

    def test_no_nulls_passes(self):
        """Test validation passes when no nulls exist."""
        df = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
        validate_no_nulls(df, ["col1", "col2"], "test.tsv")  # Should not raise

    def test_null_values_raise(self):
        """Test validation fails when null values exist."""
        df = pd.DataFrame({"col1": [1, None, 3]})
        with pytest.raises(ValidationError, match="null value"):
            validate_no_nulls(df, ["col1"], "test.tsv")

    def test_na_values_raise(self):
        """Test validation fails when NA values exist."""
        df = pd.DataFrame({"col1": [1, pd.NA, 3]})
        with pytest.raises(ValidationError, match="null value"):
            validate_no_nulls(df, ["col1"], "test.tsv")

    def test_multiple_null_columns(self):
        """Test validation reports multiple nulls."""
        df = pd.DataFrame({"col1": [1, None, 3, None]})
        with pytest.raises(ValidationError, match="2 null value"):
            validate_no_nulls(df, ["col1"], "test.tsv")

    def test_missing_column_is_ignored(self):
        """Test validation skips missing columns."""
        df = pd.DataFrame({"col1": [1, 2, 3]})
        validate_no_nulls(df, ["col1", "missing_col"], "test.tsv")  # Should not raise


class TestValidateTSVSchema:
    """Tests for validate_tsv_schema function."""

    def test_valid_tsv_file(self):
        """Test validation of valid TSV file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
            tmp.write("chrom\tpos\ttype\n")
            tmp.write("chr1\t100\tSNP\n")
            tmp.write("chr2\t200\tINDEL\n")
            tmp_path = Path(tmp.name)

        try:
            df = validate_tsv_schema(
                tmp_path,
                required_columns=["chrom", "pos", "type"],
                enum_columns={"type": {"SNP", "INDEL"}},
                positive_columns=["pos"],
            )
            assert len(df) == 2
            assert list(df.columns) == ["chrom", "pos", "type"]
        finally:
            tmp_path.unlink()

    def test_missing_required_column_raises(self):
        """Test validation fails for missing required columns."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
            tmp.write("chrom\tpos\n")
            tmp.write("chr1\t100\n")
            tmp_path = Path(tmp.name)

        try:
            with pytest.raises(ValidationError, match="missing required columns"):
                validate_tsv_schema(tmp_path, required_columns=["chrom", "pos", "type"])
        finally:
            tmp_path.unlink()

    def test_invalid_enum_value_raises(self):
        """Test validation fails for invalid enum values."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
            tmp.write("type\n")
            tmp.write("INVALID\n")
            tmp_path = Path(tmp.name)

        try:
            with pytest.raises(ValidationError, match="invalid values"):
                validate_tsv_schema(
                    tmp_path,
                    required_columns=["type"],
                    enum_columns={"type": {"SNP", "INDEL"}},
                )
        finally:
            tmp_path.unlink()

    def test_negative_value_raises(self):
        """Test validation fails for negative values in positive columns."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
            tmp.write("pos\n")
            tmp.write("-100\n")
            tmp_path = Path(tmp.name)

        try:
            with pytest.raises(ValidationError, match="negative values"):
                validate_tsv_schema(
                    tmp_path, required_columns=["pos"], positive_columns=["pos"]
                )
        finally:
            tmp_path.unlink()

    def test_empty_file_raises(self):
        """Test validation fails for empty file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
            tmp.write("col1\tcol2\n")  # Only header, no data
            tmp_path = Path(tmp.name)

        try:
            with pytest.raises(ValidationError, match="contains no data rows"):
                validate_tsv_schema(tmp_path, required_columns=["col1"])
        finally:
            tmp_path.unlink()

    def test_missing_file_raises(self):
        """Test validation fails for missing file."""
        with pytest.raises(ValidationError, match="not found"):
            validate_tsv_schema("/nonexistent/file.tsv", required_columns=["col1"])


class TestValidateRequiredColumnsEdgeCases:
    """Edge case tests for validate_required_columns."""

    def test_unicode_column_names(self):
        """Test handling of unicode characters in column names."""
        df = pd.DataFrame({"col_α": [1, 2], "col_β": [3, 4]})
        validate_required_columns(df, ["col_α", "col_β"], "test.tsv")

    def test_column_names_with_spaces(self):
        """Test handling of column names with spaces."""
        df = pd.DataFrame({"column one": [1], "column two": [2]})
        validate_required_columns(df, ["column one", "column two"], "test.tsv")

    def test_case_sensitive_columns(self):
        """Test that column matching is case-sensitive."""
        df = pd.DataFrame({"Column": [1, 2]})
        with pytest.raises(ValidationError, match="missing required columns"):
            validate_required_columns(df, ["column"], "test.tsv")


class TestValidateEnumValuesEdgeCases:
    """Edge case tests for validate_enum_values."""

    def test_numeric_enum_values(self):
        """Test enum validation with numeric values."""
        df = pd.DataFrame({"type": [1, 2, 1, 2]})
        validate_enum_values(df, "type", {1, 2}, "test.tsv")

    def test_mixed_type_enum_values(self):
        """Test enum validation with mixed types (should fail)."""
        df = pd.DataFrame({"type": ["SNP", 1, "INDEL"]})
        # Mixed types in enum values
        with pytest.raises(ValidationError):
            validate_enum_values(df, "type", {"SNP", "INDEL"}, "test.tsv")

    def test_whitespace_in_enum_values(self):
        """Test enum values with leading/trailing whitespace."""
        df = pd.DataFrame({"type": [" SNP ", "INDEL "]})
        # Whitespace matters in enum matching
        with pytest.raises(ValidationError, match="invalid values"):
            validate_enum_values(df, "type", {"SNP", "INDEL"}, "test.tsv")


class TestValidatePositiveNumericEdgeCases:
    """Edge case tests for validate_positive_numeric."""

    def test_very_small_positive_value(self):
        """Test with very small positive values."""
        df = pd.DataFrame({"coverage": [0.0001, 0.00001, 1e-10]})
        validate_positive_numeric(df, ["coverage"], "test.tsv")

    def test_infinity_values(self):
        """Test that infinity values are allowed (they're technically positive).

        Note: The validator only checks for negative values, not infinity.
        If infinity should be rejected, the validator needs to be updated.
        """
        df = pd.DataFrame({"coverage": [10.0, float('inf'), 20.0]})
        # Positive infinity passes (not < 0)
        validate_positive_numeric(df, ["coverage"], "test.tsv")

        # Negative infinity should fail
        df_neg = pd.DataFrame({"coverage": [10.0, float('-inf'), 20.0]})
        with pytest.raises(ValidationError):
            validate_positive_numeric(df_neg, ["coverage"], "test.tsv")

    def test_mixed_numeric_types(self):
        """Test with mixed int and float values."""
        df = pd.DataFrame({"value": [10, 20.5, 30, 40.1]})
        validate_positive_numeric(df, ["value"], "test.tsv")
