"""Shared test fixtures for pipeline script tests."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add scripts directory to path so all tests can import from it
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)


@pytest.fixture
def temp_tsv_file(tmp_path):
    """Fixture to create temporary TSV files for testing."""

    def _create_tsv(content: str, filename: str = "test.tsv") -> Path:
        """Create a TSV file with given content.

        Args:
            content: TSV content (including headers)
            filename: Name of the file

        Returns:
            Path to the created file
        """
        file_path = tmp_path / filename
        file_path.write_text(content)
        return file_path

    return _create_tsv
