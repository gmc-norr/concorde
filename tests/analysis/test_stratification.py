"""Tests for the stratification engine."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.stratification import (  # noqa: E402
    StratificationEngine,
    _assign_functional_impact,
    _assign_indel_size,
    _assign_vaf_bin,
    _assign_variant_class,
    _assign_zygosity,
    _compute_stratum_metrics,
)


def _make_variant(
    vtype="SNP",
    classification="TP",
    impact=None,
    zygosity=None,
    tumor_af=None,
    af=None,
    indel_size=None,
    gene_sets=None,
):
    """Create a mock variant for testing."""
    v = MagicMock()
    v.type = vtype
    v.classification = classification
    v.impact = impact
    v.zygosity = zygosity
    v.tumor_af = tumor_af
    v.af = af
    v.indel_size = indel_size
    v.gene_sets = gene_sets or []
    return v


class TestAssignVariantClass:
    def test_snp(self):
        assert _assign_variant_class(_make_variant(vtype="SNP")) == "SNP"

    def test_indel(self):
        assert _assign_variant_class(_make_variant(vtype="INDEL")) == "INDEL"

    def test_unknown(self):
        assert _assign_variant_class(_make_variant(vtype="OTHER")) is None


class TestAssignIndelSize:
    def test_1bp(self):
        v = _make_variant(vtype="INDEL", indel_size=1)
        assert _assign_indel_size(v) == "1bp"

    def test_3bp(self):
        v = _make_variant(vtype="INDEL", indel_size=3)
        assert _assign_indel_size(v) == "2-5bp"

    def test_10bp(self):
        v = _make_variant(vtype="INDEL", indel_size=10)
        assert _assign_indel_size(v) == "6-15bp"

    def test_25bp(self):
        v = _make_variant(vtype="INDEL", indel_size=25)
        assert _assign_indel_size(v) == "16-50bp"

    def test_100bp(self):
        v = _make_variant(vtype="INDEL", indel_size=100)
        assert _assign_indel_size(v) == ">50bp"

    def test_snp_returns_none(self):
        v = _make_variant(vtype="SNP", indel_size=None)
        assert _assign_indel_size(v) is None

    def test_no_size_returns_none(self):
        v = _make_variant(vtype="INDEL", indel_size=None)
        assert _assign_indel_size(v) is None


class TestAssignFunctionalImpact:
    def test_high(self):
        assert _assign_functional_impact(_make_variant(impact="HIGH")) == "HIGH"

    def test_moderate(self):
        assert _assign_functional_impact(_make_variant(impact="MODERATE")) == "MODERATE"

    def test_none(self):
        assert _assign_functional_impact(_make_variant(impact=None)) is None


class TestAssignZygosity:
    def test_het(self):
        assert _assign_zygosity(_make_variant(zygosity="HET")) == "HET"

    def test_hom_alt(self):
        assert _assign_zygosity(_make_variant(zygosity="HOM_ALT")) == "HOM_ALT"

    def test_none(self):
        assert _assign_zygosity(_make_variant(zygosity=None)) is None


class TestAssignVafBin:
    def test_low_vaf(self):
        v = _make_variant(tumor_af=0.03)
        assert _assign_vaf_bin(v) == "<0.05"

    def test_mid_vaf(self):
        v = _make_variant(tumor_af=0.15)
        assert _assign_vaf_bin(v) == "0.10-0.20"

    def test_high_vaf(self):
        v = _make_variant(tumor_af=0.60)
        assert _assign_vaf_bin(v) == ">0.50"

    def test_no_af(self):
        v = _make_variant(tumor_af=None, af=None)
        assert _assign_vaf_bin(v) is None

    def test_falls_back_to_af(self):
        v = _make_variant(tumor_af=None, af=0.25)
        assert _assign_vaf_bin(v) == "0.20-0.50"


class TestComputeStratumMetrics:
    def test_basic_metrics(self):
        variants = [
            _make_variant(classification="TP"),
            _make_variant(classification="TP"),
            _make_variant(classification="FP"),
            _make_variant(classification="FN"),
        ]
        metrics = _compute_stratum_metrics(variants)

        assert metrics["tp_count"] == 2
        assert metrics["fp_count"] == 1
        assert metrics["fn_count"] == 1
        assert metrics["precision"] == pytest.approx(2 / 3)
        assert metrics["recall"] == pytest.approx(2 / 3)
        assert metrics["total_variants"] == 4

    def test_perfect_precision_recall(self):
        variants = [
            _make_variant(classification="TP"),
            _make_variant(classification="TP"),
        ]
        metrics = _compute_stratum_metrics(variants)
        assert metrics["precision"] == 1.0
        assert metrics["recall"] == 1.0
        assert metrics["f1"] == 1.0

    def test_empty_returns_none_metrics(self):
        metrics = _compute_stratum_metrics([])
        assert metrics["precision"] is None
        assert metrics["recall"] is None
        assert metrics["f1"] is None

    def test_vaf_discordant_counts_as_tp(self):
        variants = [
            _make_variant(classification="TP_VAF_DISCORDANT"),
            _make_variant(classification="FN"),
        ]
        metrics = _compute_stratum_metrics(variants)
        assert metrics["tp_count"] == 1
        assert metrics["fn_count"] == 1
        assert metrics["recall"] == pytest.approx(0.5)


class TestStratificationEngine:
    def test_stratify_by_variant_class(self):
        engine = StratificationEngine(mode="germline")
        variants = [
            _make_variant(vtype="SNP"),
            _make_variant(vtype="SNP"),
            _make_variant(vtype="INDEL"),
        ]

        strata = engine.stratify(variants)

        assert len(strata[("variant_class", "SNP")]) == 2
        assert len(strata[("variant_class", "INDEL")]) == 1

    def test_stratify_skips_somatic_dimensions_in_germline(self):
        config = {
            "dimensions": {
                "variant_class": {"enabled": True},
                "vaf_bins": {"enabled": True, "mode": "somatic"},
            }
        }
        engine = StratificationEngine(config=config, mode="germline")
        variants = [_make_variant(vtype="SNP", tumor_af=0.30)]

        strata = engine.stratify(variants)

        assert ("variant_class", "SNP") in strata
        assert not any(k[0] == "vaf_bins" for k in strata)

    def test_stratify_includes_somatic_dimensions(self):
        config = {
            "dimensions": {
                "variant_class": {"enabled": True},
                "vaf_bins": {"enabled": True, "mode": "somatic"},
            }
        }
        engine = StratificationEngine(config=config, mode="somatic")
        variants = [_make_variant(vtype="SNP", tumor_af=0.30)]

        strata = engine.stratify(variants)

        assert ("variant_class", "SNP") in strata
        assert ("vaf_bins", "0.20-0.50") in strata

    def test_disabled_dimension_skipped(self):
        config = {
            "dimensions": {
                "variant_class": {"enabled": False},
                "functional_impact": {"enabled": True},
            }
        }
        engine = StratificationEngine(config=config, mode="germline")
        variants = [_make_variant(vtype="SNP", impact="HIGH")]

        strata = engine.stratify(variants)

        assert not any(k[0] == "variant_class" for k in strata)
        assert ("functional_impact", "HIGH") in strata

    def test_compute_metrics_with_low_confidence(self):
        engine = StratificationEngine(min_variants=5)
        strata = {
            ("variant_class", "SNP"): [_make_variant(classification="TP")] * 10,
            ("variant_class", "INDEL"): [_make_variant(classification="TP")] * 3,
        }

        metrics = engine.compute_metrics(strata)

        snp_m = [m for m in metrics if m["stratum"] == "SNP"][0]
        indel_m = [m for m in metrics if m["stratum"] == "INDEL"][0]

        assert snp_m["low_confidence"] is False
        assert indel_m["low_confidence"] is True

    def test_gene_panel_off_panel(self):
        engine = StratificationEngine(mode="germline")
        variants = [_make_variant(vtype="SNP")]  # No gene_sets

        strata = engine.stratify(variants)

        assert ("gene_panel", "off_panel") in strata
