"""Tests for pipeline/scripts/utils.py utility functions."""

from __future__ import annotations

import math
import sys
from pathlib import Path


# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent / "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from utils import safe_float, safe_int, safe_str


class TestSafeFloat:
    """Tests for safe_float function."""

    def test_valid_float(self):
        """Test conversion of valid float values."""
        assert safe_float(3.14) == 3.14
        assert safe_float("2.5") == 2.5
        assert safe_float(42) == 42.0

    def test_none_returns_default(self):
        """Test that None returns default value."""
        assert safe_float(None) is None
        assert safe_float(None, default=0.0) == 0.0

    def test_invalid_string_returns_default(self):
        """Test that invalid strings return default."""
        assert safe_float("not a number") is None
        assert safe_float("not a number", default=-1.0) == -1.0

    def test_nan_returns_default(self):
        """Test that NaN returns default value."""
        assert safe_float(float("nan")) is None
        assert safe_float(math.nan, default=0.0) == 0.0

    def test_empty_string_returns_default(self):
        """Test that empty string returns default."""
        assert safe_float("") is None
        assert safe_float("", default=999.0) == 999.0

    def test_infinity_returns_default(self):
        """Test that infinity values return default."""
        assert safe_float(float('inf')) is None
        assert safe_float(float('-inf')) is None
        assert safe_float("inf", default=0.0) == 0.0
        assert safe_float("-inf", default=0.0) == 0.0

    def test_extreme_values(self):
        """Test handling of extreme float values."""
        very_large = 1e308  # Near max float
        assert safe_float(very_large) == very_large
        assert safe_float(-1e308) == -1e308


class TestSafeInt:
    """Tests for safe_int function."""

    def test_valid_int(self):
        """Test conversion of valid integer values."""
        assert safe_int(42) == 42
        assert safe_int("123") == 123
        assert safe_int(3.14) == 3

    def test_none_returns_default(self):
        """Test that None returns default value."""
        assert safe_int(None) is None
        assert safe_int(None, default=0) == 0

    def test_invalid_string_returns_default(self):
        """Test that invalid strings return default."""
        assert safe_int("not a number") is None
        assert safe_int("not a number", default=-1) == -1

    def test_empty_string_returns_default(self):
        """Test that empty string returns default."""
        assert safe_int("") is None
        assert safe_int("", default=999) == 999

    def test_float_string_truncates(self):
        """Test that float strings are converted to int."""
        assert safe_int("3.14") == 3
        assert safe_int("9.99") == 9

    def test_very_large_int(self):
        """Test handling of very large integers."""
        large_int = 2**31 - 1  # Max 32-bit signed int
        assert safe_int(large_int) == large_int
        assert safe_int(str(large_int)) == large_int


class TestSafeStr:
    """Tests for safe_str function."""

    def test_valid_string(self):
        """Test conversion of valid string values."""
        assert safe_str("hello") == "hello"
        assert safe_str("") == ""

    def test_number_to_string(self):
        """Test conversion of numbers to strings."""
        assert safe_str(123) == "123"
        assert safe_str(3.14) == "3.14"

    def test_none_returns_default(self):
        """Test that None returns default value."""
        assert safe_str(None) is None
        assert safe_str(None, default="N/A") == "N/A"

    def test_boolean_to_string(self):
        """Test conversion of booleans to strings."""
        assert safe_str(True) == "True"
        assert safe_str(False) == "False"

    def test_strips_whitespace(self):
        """Test that strings are stripped of whitespace."""
        assert safe_str("  hello  ") == "hello"
        assert safe_str("\tworld\n") == "world"

    def test_empty_after_strip_returns_default(self):
        """Test that whitespace-only strings return default."""
        assert safe_str("   ") is None
        assert safe_str("\t\n", default="empty") == "empty"
