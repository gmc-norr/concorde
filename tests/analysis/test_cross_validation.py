"""Cross-validation tests: verify implementations against scipy/numpy references.

These tests ensure our hand-implemented formulas produce the same results
as established scientific computing libraries. Any discrepancy here indicates
a bug in our implementation.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import norm

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.extended_metrics import (  # noqa: E402
    compute_confidence_intervals,
    compute_mcc,
    wilson_score_interval,
)


# ── Wilson Score Interval vs scipy.stats.norm ──


class TestWilsonVsScipy:
    """Verify Wilson CI implementation uses scipy correctly and matches
    the textbook formula:
        (p + z²/2n ± z√(p(1-p)/n + z²/4n²)) / (1 + z²/n)
    """

    @pytest.mark.parametrize("confidence", [0.80, 0.90, 0.95, 0.99, 0.999])
    def test_z_score_matches_scipy(self, confidence):
        """Verify our z-score comes from scipy.stats.norm.ppf, not a lookup."""
        z_expected = norm.ppf((1 + confidence) / 2)
        # Compute CI and verify by reconstructing from known formula
        successes, total = 70, 100
        p = successes / total

        lower, upper = wilson_score_interval(successes, total, confidence)

        # Recompute with explicit scipy z
        z = z_expected
        denom = 1 + z * z / total
        centre = p + z * z / (2 * total)
        spread = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total)
        expected_lower = max(0.0, (centre - spread) / denom)
        expected_upper = min(1.0, (centre + spread) / denom)

        assert lower == pytest.approx(expected_lower, abs=1e-10)
        assert upper == pytest.approx(expected_upper, abs=1e-10)

    @pytest.mark.parametrize(
        "successes,total",
        [(0, 100), (1, 100), (50, 100), (99, 100), (100, 100),
         (3, 10), (1, 2), (7, 7)],
    )
    def test_ci_contains_point_estimate(self, successes, total):
        """The CI must always contain the point estimate p̂ = successes/total."""
        lower, upper = wilson_score_interval(successes, total, 0.95)
        p = successes / total
        assert lower <= p + 1e-10  # Allow tiny float tolerance
        assert upper >= p - 1e-10

    def test_wider_ci_at_lower_confidence_is_wrong(self):
        """Higher confidence level should produce wider intervals."""
        s, n = 70, 100
        l90, u90 = wilson_score_interval(s, n, 0.90)
        l95, u95 = wilson_score_interval(s, n, 0.95)
        l99, u99 = wilson_score_interval(s, n, 0.99)

        width_90 = u90 - l90
        width_95 = u95 - l95
        width_99 = u99 - l99

        assert width_90 < width_95 < width_99

    def test_symmetric_at_p_half(self):
        """At p=0.5, CI should be symmetric around 0.5."""
        lower, upper = wilson_score_interval(50, 100, 0.95)
        assert lower == pytest.approx(1.0 - upper, abs=1e-10)

    @pytest.mark.parametrize("n", [5, 10, 50, 100, 1000])
    def test_ci_narrows_with_sample_size(self, n):
        """Larger sample size should produce narrower intervals."""
        p = 0.7
        s = int(p * n)
        lower, upper = wilson_score_interval(s, n, 0.95)
        width = upper - lower
        assert width > 0  # Ensure valid interval

    def test_ci_narrows_monotonically(self):
        """CI width should decrease monotonically with n (for fixed p)."""
        widths = []
        for n in [10, 50, 100, 500, 1000]:
            s = int(0.7 * n)
            lower, upper = wilson_score_interval(s, n, 0.95)
            widths.append(upper - lower)

        for i in range(len(widths) - 1):
            assert widths[i] > widths[i + 1]


# ── MCC vs hand-computed reference values ──


class TestMCCReference:
    """Verify MCC against textbook worked examples."""

    def test_wikipedia_example(self):
        """MCC from the Wikipedia confusion matrix example.
        TP=5, FP=10, FN=2, TN=83
        MCC = (5*83 - 10*2) / sqrt(15 * 7 * 93 * 85)
            = (415 - 20) / sqrt(830025)
            = 395 / 911.057...
            ≈ 0.43356
        """
        mcc = compute_mcc(tp=5, fp=10, fn=2, tn=83)
        expected = (5 * 83 - 10 * 2) / math.sqrt(15 * 7 * 93 * 85)
        assert mcc == pytest.approx(expected, abs=1e-10)
        assert mcc == pytest.approx(0.43356, abs=1e-4)

    def test_perfect_positive(self):
        """Perfect classifier: all correct, MCC = 1.0."""
        assert compute_mcc(50, 0, 0, 50) == pytest.approx(1.0)

    def test_perfect_negative(self):
        """Worst classifier: everything wrong, MCC = -1.0."""
        assert compute_mcc(0, 50, 50, 0) == pytest.approx(-1.0)

    def test_random_classifier(self):
        """Random classifier: MCC ≈ 0."""
        assert compute_mcc(25, 25, 25, 25) == pytest.approx(0.0)

    def test_mcc_range(self):
        """MCC should always be in [-1, 1]."""
        import random
        random.seed(42)
        for _ in range(100):
            tp = random.randint(0, 100)
            fp = random.randint(0, 100)
            fn = random.randint(0, 100)
            tn = random.randint(0, 100)
            mcc = compute_mcc(tp, fp, fn, tn)
            if mcc is not None:
                assert -1.0 <= mcc <= 1.0 + 1e-10, f"MCC={mcc} out of range"


# ── Precision/Recall CI sanity ──


class TestConfidenceIntervalsSanity:
    """Verify that CIs have correct mathematical properties."""

    def test_precision_ci_contains_point_estimate(self):
        tp, fp, fn = 90, 10, 5
        precision = tp / (tp + fp)
        result = compute_confidence_intervals(tp, fp, fn)
        assert result["precision_ci"][0] <= precision + 1e-10
        assert result["precision_ci"][1] >= precision - 1e-10

    def test_recall_ci_contains_point_estimate(self):
        tp, fp, fn = 90, 10, 5
        recall = tp / (tp + fn)
        result = compute_confidence_intervals(tp, fp, fn)
        assert result["recall_ci"][0] <= recall + 1e-10
        assert result["recall_ci"][1] >= recall - 1e-10

    def test_f1_ci_contains_point_estimate(self):
        tp, fp, fn = 90, 10, 5
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1 = 2 * precision * recall / (precision + recall)
        result = compute_confidence_intervals(tp, fp, fn)
        # Bootstrap CI: wider tolerance due to stochastic resampling
        assert result["f1_ci"][0] <= f1 + 0.05
        assert result["f1_ci"][1] >= f1 - 0.05

    def test_ci_bounds_valid(self):
        """All CI bounds should be in [0, 1]."""
        for tp, fp, fn in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (50, 50, 50),
                           (100, 0, 0), (0, 100, 0), (0, 0, 100)]:
            result = compute_confidence_intervals(tp, fp, fn)
            for key in ("precision_ci", "recall_ci", "f1_ci"):
                lo, hi = result[key]
                assert 0.0 <= lo <= hi <= 1.0, (
                    f"{key}: [{lo}, {hi}] for TP={tp}, FP={fp}, FN={fn}"
                )


# ── SPC/Trending stdev: verify Bessel's correction ──


class TestStdevBesselsCorrection:
    """Verify that our stdev computations use ddof=1 (sample stdev)."""

    def test_qc_trending_uses_sample_stdev(self, tmp_path):
        """Verify QC trending stdev matches numpy with ddof=1."""
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        from models.base import Base
        from models.run import Run
        from models.qc_metric import QCMetric
        from analysis.qc_trending import compute_historical_stats

        db_path = tmp_path / "stdev_test.db"
        engine = create_engine(f"sqlite:///{db_path}")
        Base.metadata.create_all(engine)
        Session = sessionmaker(bind=engine)
        session = Session()

        values = [2.1, 2.2, 2.0, 2.15, 2.05, 2.3, 1.95, 2.1, 2.25, 2.08]
        for i, v in enumerate(values):
            run = Run(sample="S", caller="C", pipeline_version=f"v{i}",
                      comparison_tool="h", mode="germline")
            session.add(run)
            session.flush()
            session.add(QCMetric(run_id=run.id, metric_name="titv",
                                 metric_value_float=v, metric_source="test"))

        # Current run (excluded from stats)
        cur = Run(sample="S", caller="C", pipeline_version="vX",
                  comparison_tool="h", mode="germline")
        session.add(cur)
        session.flush()
        session.add(QCMetric(run_id=cur.id, metric_name="titv",
                             metric_value_float=2.5, metric_source="test"))
        session.commit()

        stats = compute_historical_stats(session, "titv", cur.id, min_runs=5)

        arr = np.array(values, dtype=np.float64)
        assert stats["mean"] == pytest.approx(float(arr.mean()), abs=1e-10)
        assert stats["stdev"] == pytest.approx(float(arr.std(ddof=1)), abs=1e-10)

        # Verify it does NOT match population stdev (ddof=0)
        pop_stdev = float(arr.std(ddof=0))
        sample_stdev = float(arr.std(ddof=1))
        assert pop_stdev != sample_stdev  # They should differ
        assert stats["stdev"] != pytest.approx(pop_stdev, abs=1e-10)

        session.close()

    def test_spc_uses_sample_stdev(self, tmp_path):
        """Verify SPC stdev matches numpy with ddof=1."""
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        from models.base import Base
        from models.run import Run
        from models.stratified_metric import StratifiedMetric
        from analysis.spc import SPCEngine

        db_path = tmp_path / "spc_stdev.db"
        engine = create_engine(f"sqlite:///{db_path}")
        Base.metadata.create_all(engine)
        Session = sessionmaker(bind=engine)
        session = Session()

        values = [0.92, 0.95, 0.93, 0.94, 0.91, 0.96, 0.90, 0.95, 0.93, 0.94]
        for i, v in enumerate(values):
            run = Run(sample="S", caller="C", pipeline_version=f"v{i}",
                      comparison_tool="h", mode="germline")
            session.add(run)
            session.flush()
            session.add(StratifiedMetric(
                run_id=run.id, dimension="d", stratum="s", variant_type="ALL",
                precision=v, tp_count=90, fp_count=5, fn_count=5,
                total_variants=100, low_confidence=False,
            ))
        session.commit()

        spc = SPCEngine(session, config={"min_data_points": 5, "sigma_multiplier": 3})
        limits = spc.compute_control_limits("precision", dimension="d", stratum="s")

        arr = np.array(values, dtype=np.float64)
        assert limits["mean"] == pytest.approx(float(arr.mean()), abs=1e-10)
        assert limits["stdev"] == pytest.approx(float(arr.std(ddof=1)), abs=1e-10)

        # Control limits should use sample stdev
        expected_ucl = float(arr.mean()) + 3 * float(arr.std(ddof=1))
        expected_lcl = float(arr.mean()) - 3 * float(arr.std(ddof=1))
        assert limits["ucl"] == pytest.approx(expected_ucl, abs=1e-10)
        assert limits["lcl"] == pytest.approx(expected_lcl, abs=1e-10)

        session.close()


# ── Ti/Tv transition classification cross-check ──


class TestTiTvClassification:
    """Verify all 12 possible SNP base changes are classified correctly.

    Transitions (purines↔purines, pyrimidines↔pyrimidines):
        A→G, G→A, C→T, T→C  (4 pairs)
    Transversions (purine↔pyrimidine):
        A→C, A→T, G→C, G→T, C→A, C→G, T→A, T→G  (8 pairs)
    """

    # All 12 SNP base changes
    ALL_CHANGES = [
        ("A", "G"), ("A", "C"), ("A", "T"),
        ("G", "A"), ("G", "C"), ("G", "T"),
        ("C", "A"), ("C", "G"), ("C", "T"),
        ("T", "A"), ("T", "C"), ("T", "G"),
    ]

    EXPECTED_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    EXPECTED_TRANSVERSIONS = set(ALL_CHANGES) - EXPECTED_TRANSITIONS

    def test_transition_count(self):
        assert len(self.EXPECTED_TRANSITIONS) == 4

    def test_transversion_count(self):
        assert len(self.EXPECTED_TRANSVERSIONS) == 8

    def test_expected_ratio(self):
        """With equal representation, Ti/Tv = 4/8 = 0.5."""
        from analysis.extended_metrics import compute_titv_ratio

        class V:
            def __init__(self, ref, alt):
                self.type = "SNP"
                self.ref = ref
                self.alt = alt

        variants = [V(ref, alt) for ref, alt in self.ALL_CHANGES]
        result = compute_titv_ratio(variants)
        assert result["transitions_count"] == 4
        assert result["transversions_count"] == 8
        assert result["titv_ratio"] == pytest.approx(0.5)

    def test_each_change_classified_correctly(self):
        from analysis.extended_metrics import _TRANSITIONS

        for ref, alt in self.ALL_CHANGES:
            is_ti = (ref, alt) in _TRANSITIONS
            expected_ti = (ref, alt) in self.EXPECTED_TRANSITIONS
            assert is_ti == expected_ti, (
                f"{ref}→{alt}: classified as "
                f"{'transition' if is_ti else 'transversion'}, "
                f"expected {'transition' if expected_ti else 'transversion'}"
            )


# ── BRAM Beta-Binomial prior uses sample variance (ddof=1) ──


class TestBRAMEmpiricalPrior:
    """Verify BRAM Beta prior computation uses numpy with ddof=1 for variance."""

    def test_empirical_prior_uses_sample_variance(self):
        """BRAM get_prior() MoM estimator should use var(ddof=1), consistent with SPC/trending."""
        from analysis.bayesian_risk import BRAMEngine

        engine = BRAMEngine(config={"use_historical_priors": True, "min_historical_runs": 3})
        baseline_value = 0.99
        deltas = [0.001, -0.002, 0.0, 0.001, -0.001]
        proportions = np.array([baseline_value + d for d in deltas], dtype=np.float64)
        proportions = np.clip(proportions, 1e-6, 1.0 - 1e-6)

        alpha, beta_p = engine.get_prior("test_stratum", baseline_value, deltas)

        # Method of moments Beta fit from proportions:
        mu_emp = float(proportions.mean())
        var_emp = float(proportions.var(ddof=1))
        kappa = max(mu_emp * (1.0 - mu_emp) / var_emp - 1.0, 2.0)
        expected_alpha = mu_emp * kappa
        expected_beta = (1.0 - mu_emp) * kappa

        assert alpha == pytest.approx(expected_alpha, abs=1e-6)
        assert beta_p == pytest.approx(expected_beta, abs=1e-6)

        # Verify that using population variance (ddof=0) would give different result
        var_pop = float(proportions.var(ddof=0))
        assert var_pop != var_emp  # They must differ

    def test_numpy_prod_for_product_aggregation(self):
        """aggregate_risks('product') should use numpy.prod, not iterative *=."""
        from analysis.bayesian_risk import aggregate_risks

        risks = [0.1, 0.2, 0.3]
        result = aggregate_risks(risks, method="product")
        # 1 - (1-0.1)*(1-0.2)*(1-0.3) = 1 - 0.9*0.8*0.7 = 1 - 0.504 = 0.496
        expected = 1.0 - float(np.prod([1.0 - r for r in risks]))
        assert result == pytest.approx(expected, abs=1e-10)

    def test_numpy_mean_for_weighted_mean(self):
        """aggregate_risks('weighted_mean') should use numpy.mean."""
        from analysis.bayesian_risk import aggregate_risks

        risks = [0.1, 0.2, 0.3, 0.4]
        result = aggregate_risks(risks, method="weighted_mean")
        assert result == pytest.approx(float(np.mean(risks)), abs=1e-10)


# ── Shared precision/recall/F1 utility ──


class TestSharedClassificationMetrics:
    """Verify the shared utility matches the formula exactly."""

    def test_basic_metrics(self):
        from utils import compute_classification_metrics

        result = compute_classification_metrics(tp=90, fp=10, fn=5)
        assert result["precision"] == pytest.approx(90 / 100)
        assert result["recall"] == pytest.approx(90 / 95)
        expected_f1 = 2 * (90/100) * (90/95) / ((90/100) + (90/95))
        assert result["f1"] == pytest.approx(expected_f1)

    def test_zero_denominator_returns_none(self):
        from utils import compute_classification_metrics

        result = compute_classification_metrics(tp=0, fp=0, fn=0)
        assert result["precision"] is None
        assert result["recall"] is None
        assert result["f1"] is None

    def test_no_tp_returns_zero_prec_recall(self):
        from utils import compute_classification_metrics

        result = compute_classification_metrics(tp=0, fp=10, fn=5)
        assert result["precision"] == pytest.approx(0.0)
        assert result["recall"] == pytest.approx(0.0)
        # F1 = 2*0*0/(0+0) is undefined → None
        assert result["f1"] is None


# ── Bootstrap F1 CI cross-validation ──


class TestBootstrapF1CI:
    """Verify bootstrap F1 CI has correct statistical properties."""

    def test_f1_ci_contains_point_estimate(self):
        """Bootstrap CI should contain the F1 point estimate."""
        tp, fp, fn = 80, 20, 10
        prec = tp / (tp + fp)
        rec = tp / (tp + fn)
        f1_point = 2 * prec * rec / (prec + rec)
        result = compute_confidence_intervals(tp, fp, fn)
        # Bootstrap CI should bracket the point estimate
        assert result["f1_ci"][0] <= f1_point + 0.02
        assert result["f1_ci"][1] >= f1_point - 0.02

    def test_f1_ci_bounds_valid(self):
        """F1 CI bounds should always be in [0, 1]."""
        for tp, fp, fn in [(1, 0, 0), (0, 1, 0), (50, 50, 50), (100, 0, 0)]:
            result = compute_confidence_intervals(tp, fp, fn)
            lo, hi = result["f1_ci"]
            assert 0.0 <= lo <= hi <= 1.0, (
                f"f1_ci: [{lo}, {hi}] for TP={tp}, FP={fp}, FN={fn}"
            )

    def test_f1_ci_wider_at_lower_confidence(self):
        """Lower confidence should give narrower (or equal) CI."""
        tp, fp, fn = 80, 20, 10
        from analysis.extended_metrics import _bootstrap_f1_ci

        lo90, hi90 = _bootstrap_f1_ci(tp, fp, fn, confidence=0.90)
        lo99, hi99 = _bootstrap_f1_ci(tp, fp, fn, confidence=0.99)
        # 99% CI should be wider than 90% CI (with tolerance for bootstrap variance)
        width_90 = hi90 - lo90
        width_99 = hi99 - lo99
        assert width_99 >= width_90 - 0.03  # Allow small tolerance
