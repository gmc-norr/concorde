"""Tests for VEP annotation merge during ingestion."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from ingestion.ingest_helpers import merge_vep_annotations  # noqa: E402


def _make_variant(chrom, pos, ref, alt):
    """Create a mock Variant object for testing."""
    v = MagicMock()
    v.chrom = chrom
    v.pos = pos
    v.ref = ref
    v.alt = alt
    v.consequence = None
    v.impact = None
    v.gene_symbol = None
    v.gene_id = None
    v.transcript_id = None
    v.hgvsc = None
    v.hgvsp = None
    v.exon = None
    v.protein_position = None
    v.sift_prediction = None
    v.polyphen_prediction = None
    return v


class TestMergeVepAnnotations:
    """Tests for VEP annotation merge function."""

    def test_merge_basic_annotations(self, tmp_path):
        tsv = tmp_path / "vep.tsv"
        tsv.write_text(
            "chrom\tpos\tref\talt\tconsequence\timpact\tgene_symbol\tgene_id\ttranscript_id\thgvsc\thgvsp\texon\tprotein_position\tsift_prediction\tpolyphen_prediction\n"
            "chr17\t7577538\tC\tT\tmissense_variant\tMODERATE\tTP53\tENSG00000141510\tENST00000269305\tc.743G>A\tp.Arg248Gln\t7/11\t248\tdeleterious(0.0)\tprobably_damaging(1.0)\n"
        )

        v1 = _make_variant("chr17", 7577538, "C", "T")
        v2 = _make_variant("chr1", 100, "A", "G")

        count = merge_vep_annotations(str(tsv), [v1, v2])

        assert count == 1
        assert v1.consequence == "missense_variant"
        assert v1.impact == "MODERATE"
        assert v1.gene_symbol == "TP53"
        assert v1.protein_position == 248
        assert v2.consequence is None  # Not annotated

    def test_merge_empty_file(self, tmp_path):
        tsv = tmp_path / "vep.tsv"
        tsv.write_text(
            "chrom\tpos\tref\talt\tconsequence\timpact\tgene_symbol\tgene_id\ttranscript_id\thgvsc\thgvsp\texon\tprotein_position\tsift_prediction\tpolyphen_prediction\n"
        )

        v1 = _make_variant("chr1", 100, "A", "T")
        count = merge_vep_annotations(str(tsv), [v1])

        assert count == 0
        assert v1.consequence is None

    def test_merge_missing_file(self, tmp_path):
        count = merge_vep_annotations(str(tmp_path / "nonexistent.tsv"), [])
        assert count == 0

    def test_merge_no_matching_variants(self, tmp_path):
        tsv = tmp_path / "vep.tsv"
        tsv.write_text(
            "chrom\tpos\tref\talt\tconsequence\timpact\tgene_symbol\tgene_id\ttranscript_id\thgvsc\thgvsp\texon\tprotein_position\tsift_prediction\tpolyphen_prediction\n"
            "chrX\t1000\tA\tT\tintron_variant\tMODIFIER\tGENE1\tENSG1\tENST1\t\t\t\t\t\t\n"
        )

        v1 = _make_variant("chr1", 100, "A", "T")
        count = merge_vep_annotations(str(tsv), [v1])

        assert count == 0
        assert v1.consequence is None

    def test_merge_multiple_variants(self, tmp_path):
        tsv = tmp_path / "vep.tsv"
        tsv.write_text(
            "chrom\tpos\tref\talt\tconsequence\timpact\tgene_symbol\tgene_id\ttranscript_id\thgvsc\thgvsp\texon\tprotein_position\tsift_prediction\tpolyphen_prediction\n"
            "chr1\t100\tA\tT\tstop_gained\tHIGH\tBRCA1\tENSG1\tENST1\tc.1A>T\tp.Lys1*\t1/10\t1\tdeleterious\tprobably_damaging\n"
            "chr2\t200\tG\tC\tsynonymous_variant\tLOW\tTP53\tENSG2\tENST2\tc.2G>C\t\t\t\tbenign\tbenign\n"
        )

        v1 = _make_variant("chr1", 100, "A", "T")
        v2 = _make_variant("chr2", 200, "G", "C")
        v3 = _make_variant("chr3", 300, "C", "G")  # No annotation

        count = merge_vep_annotations(str(tsv), [v1, v2, v3])

        assert count == 2
        assert v1.consequence == "stop_gained"
        assert v1.impact == "HIGH"
        assert v2.consequence == "synonymous_variant"
        assert v2.impact == "LOW"
        assert v3.consequence is None
