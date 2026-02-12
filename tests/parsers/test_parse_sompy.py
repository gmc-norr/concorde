"""Tests for som.py parsers (VCF and metrics)."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from parsers.parse_sompy_metrics import parse_sompy_stats  # noqa: E402
from vcf_utils import determine_variant_type  # noqa: E402


class TestDetermineVariantType:
    """Tests for variant type classification."""

    def test_snp(self):
        assert determine_variant_type("A", "T") == "SNP"

    def test_indel_insertion(self):
        assert determine_variant_type("A", "AT") == "INDEL"

    def test_indel_deletion(self):
        assert determine_variant_type("AT", "A") == "INDEL"

    def test_complex(self):
        assert determine_variant_type("AC", "TG") == "COMPLEX"


class TestParseSompyStats:
    """Tests for som.py stats CSV parsing."""

    def test_parse_valid_csv(self, tmp_path):
        csv_path = tmp_path / "stats.csv"
        csv_path.write_text(
            "type,subtype,filter,tp,fp,fn,recall,precision\n"
            "SNP,*,ALL,100,5,10,0.909,0.952\n"
            "INDEL,*,ALL,50,3,8,0.862,0.943\n"
        )

        records = parse_sompy_stats(str(csv_path))

        assert len(records) == 2
        assert records[0]["variant_type"] == "SNP"
        assert records[0]["tp_count"] == 100
        assert records[0]["fp_count"] == 5
        assert records[0]["fn_count"] == 10
        assert records[0]["recall"] == 0.909
        assert records[0]["precision"] == 0.952
        assert records[0]["stratification"] is None

    def test_parse_pass_filter(self, tmp_path):
        csv_path = tmp_path / "stats.csv"
        csv_path.write_text(
            "type,subtype,filter,tp,fp,fn,recall,precision\n"
            "SNP,*,PASS,90,2,5,0.947,0.978\n"
        )

        records = parse_sompy_stats(str(csv_path))

        assert len(records) == 1
        assert records[0]["stratification"] == "filter=PASS"

    def test_parse_subset_filter(self, tmp_path):
        csv_path = tmp_path / "stats.csv"
        csv_path.write_text(
            "type,subtype,filter,tp,fp,fn,recall,precision\n"
            "SNP,ti,PASS,60,1,3,0.952,0.984\n"
        )

        records = parse_sompy_stats(str(csv_path))

        assert records[0]["stratification"] == "filter=PASS;subset=ti"

    def test_f1_computed(self, tmp_path):
        csv_path = tmp_path / "stats.csv"
        csv_path.write_text(
            "type,subtype,filter,tp,fp,fn,recall,precision\n"
            "SNP,*,ALL,100,0,0,1.0,1.0\n"
        )

        records = parse_sompy_stats(str(csv_path))
        assert records[0]["f1"] == pytest.approx(1.0)

    def test_f1_zero_when_no_precision_recall(self, tmp_path):
        csv_path = tmp_path / "stats.csv"
        csv_path.write_text(
            "type,subtype,filter,tp,fp,fn,recall,precision\n"
            "SNP,*,ALL,0,0,0,,\n"
        )

        records = parse_sompy_stats(str(csv_path))
        assert records[0]["f1"] == 0.0

    def test_parse_missing_file(self, tmp_path):
        records = parse_sompy_stats(str(tmp_path / "nonexistent.csv"))
        assert records == []

    def test_parse_empty_csv(self, tmp_path):
        csv_path = tmp_path / "stats.csv"
        csv_path.write_text("type,subtype,filter,tp,fp,fn,recall,precision\n")

        records = parse_sompy_stats(str(csv_path))
        assert records == []
