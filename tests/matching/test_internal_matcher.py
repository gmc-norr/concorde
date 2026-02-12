"""Tests for internal somatic variant matcher."""

from __future__ import annotations

import sys
from pathlib import Path

import pysam
import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from matching.internal_matcher import (  # noqa: E402
    _extract_variants_from_vcf,
    match_variants,
)
from vcf_utils import determine_variant_type  # noqa: E402


def _create_vcf(tmp_path, filename, variants, sample="SAMPLE"):
    """Create a VCF file with given variants for testing.

    Args:
        tmp_path: pytest tmp_path fixture
        filename: Name for the VCF file
        variants: List of dicts with chrom, pos, ref, alt, and optional af, dp
        sample: Sample name

    Returns:
        Path to bgzipped, indexed VCF
    """
    vcf_path = tmp_path / f"{filename}.vcf"

    header = pysam.VariantHeader()
    header.add_sample(sample)
    header.add_line("##fileformat=VCFv4.2")
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
    header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
    header.add_line('##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping Quality">')

    chromosomes = sorted(set(v["chrom"] for v in variants))
    for chrom in chromosomes:
        header.add_line(f"##contig=<ID={chrom}>")

    with pysam.VariantFile(str(vcf_path), "w", header=header) as vcf:
        for variant in variants:
            rec = vcf.new_record()
            rec.chrom = variant["chrom"]
            rec.pos = variant["pos"]
            rec.ref = variant["ref"]
            rec.alts = (variant["alt"],)
            rec.samples[sample]["GT"] = (0, 1)
            if "af" in variant:
                rec.samples[sample]["AF"] = variant["af"]
            if "dp" in variant:
                rec.samples[sample]["DP"] = variant["dp"]
            vcf.write(rec)

    bgz_path = tmp_path / f"{filename}.vcf.gz"
    pysam.tabix_compress(str(vcf_path), str(bgz_path), force=True)
    pysam.tabix_index(str(bgz_path), preset="vcf", force=True)
    vcf_path.unlink()

    return bgz_path


class TestDetermineVariantType:
    """Tests for variant type determination."""

    def test_snp(self):
        assert determine_variant_type("A", "T") == "SNP"

    def test_insertion(self):
        assert determine_variant_type("A", "AT") == "INDEL"

    def test_deletion(self):
        assert determine_variant_type("AT", "A") == "INDEL"

    def test_complex(self):
        assert determine_variant_type("AC", "TG") == "COMPLEX"


class TestExtractVariants:
    """Tests for VCF variant extraction."""

    def test_extract_basic_variants(self, tmp_path):
        vcf_path = _create_vcf(
            tmp_path,
            "query",
            [
                {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.35, "dp": 100},
                {"chrom": "chr1", "pos": 200, "ref": "G", "alt": "C", "af": 0.50, "dp": 80},
            ],
        )

        variants = _extract_variants_from_vcf(vcf_path)

        assert len(variants) == 2
        key1 = ("chr1", 100, "A", "T")
        assert key1 in variants
        assert variants[key1]["af"] == pytest.approx(0.35)
        assert variants[key1]["dp"] == 100

    def test_extract_empty_vcf(self, tmp_path):
        vcf_path = _create_vcf(tmp_path, "empty", [])
        variants = _extract_variants_from_vcf(vcf_path)
        assert len(variants) == 0


class TestMatchVariants:
    """Tests for variant matching logic."""

    def test_exact_match_tp(self, tmp_path):
        truth = _create_vcf(
            tmp_path,
            "truth",
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.30, "dp": 100}],
        )
        query = _create_vcf(
            tmp_path,
            "query",
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.32, "dp": 120}],
        )

        records, metrics = match_variants(truth, query, vaf_tolerance=0.10)

        tp_records = [r for r in records if r["classification"] == "TP"]
        assert len(tp_records) == 1
        assert tp_records[0]["chrom"] == "chr1"
        assert tp_records[0]["pos"] == 100

    def test_vaf_discordant(self, tmp_path):
        truth = _create_vcf(
            tmp_path,
            "truth",
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.50, "dp": 100}],
        )
        query = _create_vcf(
            tmp_path,
            "query",
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.10, "dp": 100}],
        )

        records, metrics = match_variants(truth, query, vaf_tolerance=0.10)

        disc = [r for r in records if r["classification"] == "TP_VAF_DISCORDANT"]
        assert len(disc) == 1

    def test_false_positive(self, tmp_path):
        truth = _create_vcf(tmp_path, "truth", [])
        query = _create_vcf(
            tmp_path,
            "query",
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.30, "dp": 100}],
        )

        records, metrics = match_variants(truth, query)

        fp = [r for r in records if r["classification"] == "FP"]
        assert len(fp) == 1

    def test_false_negative(self, tmp_path):
        truth = _create_vcf(
            tmp_path,
            "truth",
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.30, "dp": 100}],
        )
        query = _create_vcf(tmp_path, "query", [])

        records, metrics = match_variants(truth, query)

        fn = [r for r in records if r["classification"] == "FN"]
        assert len(fn) == 1

    def test_partial_match(self, tmp_path):
        truth = _create_vcf(
            tmp_path,
            "truth",
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.30, "dp": 100}],
        )
        query = _create_vcf(
            tmp_path,
            "query",
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G", "af": 0.30, "dp": 100}],
        )

        records, metrics = match_variants(truth, query, track_partial_matches=True)

        partial = [r for r in records if r["classification"] == "PARTIAL"]
        assert len(partial) == 1

    def test_metrics_computed(self, tmp_path):
        truth = _create_vcf(
            tmp_path,
            "truth",
            [
                {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.30, "dp": 100},
                {"chrom": "chr1", "pos": 200, "ref": "G", "alt": "C", "af": 0.40, "dp": 80},
                {"chrom": "chr1", "pos": 300, "ref": "C", "alt": "G", "af": 0.50, "dp": 90},
            ],
        )
        query = _create_vcf(
            tmp_path,
            "query",
            [
                {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.32, "dp": 120},
                {"chrom": "chr1", "pos": 200, "ref": "G", "alt": "C", "af": 0.38, "dp": 70},
                {"chrom": "chr1", "pos": 400, "ref": "A", "alt": "G", "af": 0.20, "dp": 50},
            ],
        )

        records, metrics = match_variants(truth, query, vaf_tolerance=0.10)

        all_metric = [m for m in metrics if m["variant_type"] == "ALL"][0]
        assert all_metric["tp_count"] == 2
        assert all_metric["fp_count"] == 1
        assert all_metric["fn_count"] == 1
        assert all_metric["precision"] == pytest.approx(2 / 3)
        assert all_metric["recall"] == pytest.approx(2 / 3)

    def test_multiple_variant_types(self, tmp_path):
        truth = _create_vcf(
            tmp_path,
            "truth",
            [
                {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.30, "dp": 100},
                {"chrom": "chr1", "pos": 200, "ref": "AT", "alt": "A", "af": 0.40, "dp": 80},
            ],
        )
        query = _create_vcf(
            tmp_path,
            "query",
            [
                {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "af": 0.32, "dp": 120},
                {"chrom": "chr1", "pos": 200, "ref": "AT", "alt": "A", "af": 0.42, "dp": 70},
            ],
        )

        records, metrics = match_variants(truth, query, vaf_tolerance=0.10)

        snp_metric = [m for m in metrics if m["variant_type"] == "SNP"]
        indel_metric = [m for m in metrics if m["variant_type"] == "INDEL"]
        assert len(snp_metric) == 1
        assert len(indel_metric) == 1
        assert snp_metric[0]["tp_count"] == 1
        assert indel_metric[0]["tp_count"] == 1

    def test_indel_size_computed(self, tmp_path):
        truth = _create_vcf(
            tmp_path,
            "truth",
            [{"chrom": "chr1", "pos": 100, "ref": "ATG", "alt": "A", "af": 0.30, "dp": 100}],
        )
        query = _create_vcf(
            tmp_path,
            "query",
            [{"chrom": "chr1", "pos": 100, "ref": "ATG", "alt": "A", "af": 0.30, "dp": 100}],
        )

        records, metrics = match_variants(truth, query)

        assert records[0]["indel_size"] == 2
        assert records[0]["type"] == "INDEL"
