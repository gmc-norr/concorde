"""Tests for extended metrics: GT concordance, Ti/Tv, Het/Hom, CIs, MCC."""

from __future__ import annotations

import sys
from pathlib import Path
import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.extended_metrics import (  # noqa: E402
    compute_confidence_intervals,
    compute_genotype_concordance,
    compute_het_hom_ratio,
    compute_mcc,
    compute_titv_ratio,
    wilson_score_interval,
)


class _MockVariant:
    """Simple mock variant that returns None for unset attributes."""

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __getattr__(self, name):
        return None


def _make_variant(**kwargs):
    """Create a mock variant with given attributes."""
    return _MockVariant(**kwargs)


# ── Genotype Concordance ──


class TestGenotypeConcordance:
    def test_perfect_concordance(self):
        variants = [
            _make_variant(classification="TP", truth_gt="0/1", query_gt="0/1"),
            _make_variant(classification="TP", truth_gt="1/1", query_gt="1/1"),
        ]
        result = compute_genotype_concordance(variants)
        assert result["gt_concordant_count"] == 2
        assert result["gt_discordant_count"] == 0
        assert result["genotype_concordance"] == 1.0

    def test_partial_concordance(self):
        variants = [
            _make_variant(classification="TP", truth_gt="0/1", query_gt="0/1"),
            _make_variant(classification="TP", truth_gt="1/1", query_gt="0/1"),
            _make_variant(classification="TP", truth_gt="0/1", query_gt="0/1"),
        ]
        result = compute_genotype_concordance(variants)
        assert result["gt_concordant_count"] == 2
        assert result["gt_discordant_count"] == 1
        assert result["genotype_concordance"] == pytest.approx(2 / 3)

    def test_no_gt_data(self):
        variants = [
            _make_variant(classification="TP", truth_gt=None, query_gt=None),
            _make_variant(classification="TP"),
        ]
        result = compute_genotype_concordance(variants)
        assert result["gt_concordant_count"] == 0
        assert result["gt_discordant_count"] == 0
        assert result["genotype_concordance"] is None

    def test_fp_excluded(self):
        variants = [
            _make_variant(classification="FP", truth_gt="0/1", query_gt="0/1"),
            _make_variant(classification="FN", truth_gt="0/1", query_gt="0/1"),
        ]
        result = compute_genotype_concordance(variants)
        assert result["gt_concordant_count"] == 0
        assert result["genotype_concordance"] is None

    def test_tp_vaf_discordant_included(self):
        variants = [
            _make_variant(classification="TP_VAF_DISCORDANT", truth_gt="0/1", query_gt="0/1"),
        ]
        result = compute_genotype_concordance(variants)
        assert result["gt_concordant_count"] == 1
        assert result["genotype_concordance"] == 1.0

    def test_empty_variants(self):
        result = compute_genotype_concordance([])
        assert result["genotype_concordance"] is None

    def test_mixed_with_missing_gt(self):
        variants = [
            _make_variant(classification="TP", truth_gt="0/1", query_gt="0/1"),
            _make_variant(classification="TP", truth_gt="0/1", query_gt=None),
            _make_variant(classification="TP", truth_gt=None, query_gt="0/1"),
        ]
        result = compute_genotype_concordance(variants)
        assert result["gt_concordant_count"] == 1
        assert result["gt_discordant_count"] == 0
        assert result["genotype_concordance"] == 1.0


# ── Ti/Tv Ratio ──


class TestTiTvRatio:
    def test_typical_ratio(self):
        # 4 transitions: A>G, G>A, C>T, T>C; 2 transversions: A>C, G>T
        variants = [
            _make_variant(type="SNP", ref="A", alt="G"),
            _make_variant(type="SNP", ref="G", alt="A"),
            _make_variant(type="SNP", ref="C", alt="T"),
            _make_variant(type="SNP", ref="T", alt="C"),
            _make_variant(type="SNP", ref="A", alt="C"),
            _make_variant(type="SNP", ref="G", alt="T"),
        ]
        result = compute_titv_ratio(variants)
        assert result["transitions_count"] == 4
        assert result["transversions_count"] == 2
        assert result["titv_ratio"] == pytest.approx(2.0)

    def test_all_transitions(self):
        variants = [
            _make_variant(type="SNP", ref="A", alt="G"),
            _make_variant(type="SNP", ref="C", alt="T"),
        ]
        result = compute_titv_ratio(variants)
        assert result["transitions_count"] == 2
        assert result["transversions_count"] == 0
        assert result["titv_ratio"] is None  # Division by zero

    def test_all_transversions(self):
        variants = [
            _make_variant(type="SNP", ref="A", alt="T"),
            _make_variant(type="SNP", ref="C", alt="G"),
        ]
        result = compute_titv_ratio(variants)
        assert result["transitions_count"] == 0
        assert result["transversions_count"] == 2
        assert result["titv_ratio"] == pytest.approx(0.0)

    def test_no_snps(self):
        variants = [
            _make_variant(type="INDEL", ref="A", alt="AT"),
            _make_variant(type="COMPLEX", ref="AC", alt="GT"),
        ]
        result = compute_titv_ratio(variants)
        assert result["transitions_count"] == 0
        assert result["transversions_count"] == 0
        assert result["titv_ratio"] is None

    def test_empty_variants(self):
        result = compute_titv_ratio([])
        assert result["titv_ratio"] is None

    def test_non_snp_filtered(self):
        variants = [
            _make_variant(type="SNP", ref="A", alt="G"),
            _make_variant(type="INDEL", ref="A", alt="G"),  # INDEL skipped
        ]
        result = compute_titv_ratio(variants)
        assert result["transitions_count"] == 1
        assert result["transversions_count"] == 0


# ── Het/Hom Ratio ──


class TestHetHomRatio:
    def test_typical_ratio(self):
        variants = [
            _make_variant(zygosity="HET"),
            _make_variant(zygosity="HET"),
            _make_variant(zygosity="HET"),
            _make_variant(zygosity="HOM_ALT"),
        ]
        result = compute_het_hom_ratio(variants)
        assert result["het_count"] == 3
        assert result["hom_count"] == 1
        assert result["het_hom_ratio"] == pytest.approx(3.0)

    def test_all_het(self):
        variants = [
            _make_variant(zygosity="HET"),
            _make_variant(zygosity="HET"),
        ]
        result = compute_het_hom_ratio(variants)
        assert result["het_count"] == 2
        assert result["hom_count"] == 0
        assert result["het_hom_ratio"] is None

    def test_no_zygosity(self):
        variants = [
            _make_variant(zygosity=None),
            _make_variant(),
        ]
        result = compute_het_hom_ratio(variants)
        assert result["het_count"] == 0
        assert result["hom_count"] == 0
        assert result["het_hom_ratio"] is None

    def test_empty_variants(self):
        result = compute_het_hom_ratio([])
        assert result["het_hom_ratio"] is None

    def test_equal_het_hom(self):
        variants = [
            _make_variant(zygosity="HET"),
            _make_variant(zygosity="HOM_ALT"),
        ]
        result = compute_het_hom_ratio(variants)
        assert result["het_hom_ratio"] == pytest.approx(1.0)


# ── Wilson Score Interval ──


class TestWilsonScoreInterval:
    def test_perfect_score(self):
        lower, upper = wilson_score_interval(100, 100)
        assert lower > 0.95
        assert upper == pytest.approx(1.0)

    def test_zero_score(self):
        lower, upper = wilson_score_interval(0, 100)
        assert lower == pytest.approx(0.0, abs=1e-10)
        assert upper < 0.05

    def test_half_score(self):
        lower, upper = wilson_score_interval(50, 100)
        assert lower < 0.5
        assert upper > 0.5
        assert lower > 0.35
        assert upper < 0.65

    def test_small_n(self):
        lower, upper = wilson_score_interval(1, 2)
        assert lower < 0.5
        assert upper > 0.5

    def test_large_n(self):
        lower, upper = wilson_score_interval(950, 1000)
        assert upper - lower < 0.03  # Narrow interval for large N

    def test_zero_total(self):
        lower, upper = wilson_score_interval(0, 0)
        assert lower == 0.0
        assert upper == 0.0

    def test_bounds_within_zero_one(self):
        for n in [1, 5, 10, 50, 100]:
            for s in range(n + 1):
                lower, upper = wilson_score_interval(s, n)
                assert 0.0 <= lower <= upper <= 1.0


# ── Confidence Intervals ──


class TestConfidenceIntervals:
    def test_basic_ci(self):
        result = compute_confidence_intervals(tp=90, fp=10, fn=5)
        assert 0.0 <= result["precision_ci"][0] <= result["precision_ci"][1] <= 1.0
        assert 0.0 <= result["recall_ci"][0] <= result["recall_ci"][1] <= 1.0
        assert 0.0 <= result["f1_ci"][0] <= result["f1_ci"][1] <= 1.0

    def test_zero_counts(self):
        result = compute_confidence_intervals(tp=0, fp=0, fn=0)
        assert result["precision_ci"] == (0.0, 0.0)
        assert result["recall_ci"] == (0.0, 0.0)
        assert result["f1_ci"] == (0.0, 0.0)

    def test_perfect_precision(self):
        result = compute_confidence_intervals(tp=100, fp=0, fn=10)
        assert result["precision_ci"][0] > 0.95
        assert result["precision_ci"][1] == pytest.approx(1.0)

    def test_perfect_recall(self):
        result = compute_confidence_intervals(tp=100, fp=10, fn=0)
        assert result["recall_ci"][0] > 0.95
        assert result["recall_ci"][1] == pytest.approx(1.0)

    def test_f1_between_precision_and_recall(self):
        result = compute_confidence_intervals(tp=80, fp=20, fn=10)
        # F1 point estimate should be between precision and recall
        prec = 80 / 100
        rec = 80 / 90
        f1_point = 2 * prec * rec / (prec + rec)
        # CI midpoint should be close to point estimate
        f1_mid = (result["f1_ci"][0] + result["f1_ci"][1]) / 2
        assert abs(f1_mid - f1_point) < 0.1


# ── MCC ──


class TestMCC:
    def test_perfect_classification(self):
        # All TP, no errors
        mcc = compute_mcc(tp=100, fp=0, fn=0, tn=100)
        assert mcc == pytest.approx(1.0)

    def test_zero_denominator(self):
        mcc = compute_mcc(tp=0, fp=0, fn=0, tn=0)
        assert mcc is None

    def test_balanced_errors(self):
        mcc = compute_mcc(tp=50, fp=50, fn=50, tn=50)
        assert mcc == pytest.approx(0.0)

    def test_all_fp(self):
        mcc = compute_mcc(tp=0, fp=100, fn=0, tn=0)
        assert mcc is None  # denominator has (TP+FP)*(TP+FN) with TP+FN=0

    def test_typical_variant_calling(self):
        # TN=0 is typical for variant calling (we don't count true negatives)
        mcc = compute_mcc(tp=90, fp=10, fn=5, tn=0)
        assert mcc is not None
        # With TN=0, numerator = 0*0 - 10*5 = -50
        # MCC should be negative (not ideal but reflects TN=0 scenario)
        assert mcc < 0

    def test_negative_mcc(self):
        mcc = compute_mcc(tp=10, fp=90, fn=90, tn=10)
        assert mcc is not None
        assert mcc < 0

    def test_partial_denominator_zero(self):
        # TP+FN = 0 → denominator = 0
        mcc = compute_mcc(tp=0, fp=5, fn=0, tn=10)
        assert mcc is None
