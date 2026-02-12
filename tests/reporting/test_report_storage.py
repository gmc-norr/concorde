"""Tests for report storage utilities (Spec S16.4)."""

from __future__ import annotations

import sys
from datetime import UTC, datetime
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from reporting.report_storage import (  # noqa: E402
    compute_file_checksum,
    generate_report_name,
    save_report,
)


class TestGenerateReportName:
    def test_default_json(self):
        name = generate_report_name(
            "NA12878", "germline", "initial_validation",
            date=datetime(2026, 2, 7, 10, 30, 0, tzinfo=UTC),
        )
        assert name == "NA12878_germline_20260207_103000_initial_validation.json"

    def test_html_extension(self):
        name = generate_report_name(
            "NA12878", "germline", "initial_validation",
            extension="html",
            date=datetime(2026, 2, 7, 10, 30, 0, tzinfo=UTC),
        )
        assert name.endswith(".html")

    def test_pdf_extension(self):
        name = generate_report_name(
            "NA12878", "somatic", "continuous_verification",
            extension="pdf",
            date=datetime(2026, 2, 7, 10, 30, 0, tzinfo=UTC),
        )
        assert "somatic" in name
        assert name.endswith(".pdf")

    def test_sanitizes_sample_name(self):
        name = generate_report_name(
            "sample/with spaces", "germline", "initial_validation",
            date=datetime(2026, 2, 7, 10, 30, 0, tzinfo=UTC),
        )
        assert " " not in name
        assert "/" not in name

    def test_default_date(self):
        name = generate_report_name("NA12878", "germline", "initial_validation")
        assert "NA12878" in name
        assert "germline" in name


class TestComputeFileChecksum:
    def test_checksum(self, tmp_path):
        f = tmp_path / "report.json"
        f.write_text('{"test": true}')
        checksum = compute_file_checksum(str(f))
        assert len(checksum) == 64

    def test_deterministic(self, tmp_path):
        f = tmp_path / "report.json"
        f.write_text("content")
        c1 = compute_file_checksum(str(f))
        c2 = compute_file_checksum(str(f))
        assert c1 == c2


class TestSaveReport:
    def test_save_json(self, tmp_path):
        path = save_report('{"test": true}', str(tmp_path), "report.json")
        assert Path(path).exists()
        assert Path(path).read_text() == '{"test": true}'

    def test_save_creates_directory(self, tmp_path):
        out_dir = tmp_path / "reports" / "sub"
        path = save_report("content", str(out_dir), "report.html")
        assert Path(path).exists()
        assert out_dir.exists()

    def test_save_returns_path(self, tmp_path):
        path = save_report("content", str(tmp_path), "test.json")
        assert path == str(tmp_path / "test.json")
