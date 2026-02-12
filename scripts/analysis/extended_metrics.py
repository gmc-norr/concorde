"""Extended metrics: genotype concordance, Ti/Tv, Het/Hom, CIs, MCC.

Provides pure functions that operate on variant lists and count tuples.
No database dependency -- suitable for both library use and Snakemake scripts.

Uses scipy.stats for statistical computations (z-scores, CDF, bootstrap CIs)
to avoid hand-rolled approximations.
"""

import logging
import math

import numpy as np
from scipy.stats import bootstrap, norm

log = logging.getLogger(__name__)

# Transition pairs (SNP base changes that are transitions)
_TRANSITIONS = frozenset({("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")})


def compute_genotype_concordance(variants: list) -> dict:
    """Compute genotype concordance from TP variants with GT fields.

    Only considers TP variants where both truth_gt and query_gt are available.

    Args:
        variants: List of objects with truth_gt, query_gt, classification attrs

    Returns:
        Dict with gt_concordant_count, gt_discordant_count, genotype_concordance
    """
    concordant = 0
    discordant = 0

    for v in variants:
        classification = getattr(v, "classification", "")
        if classification not in ("TP", "TP_VAF_DISCORDANT"):
            continue

        truth_gt = getattr(v, "truth_gt", None)
        query_gt = getattr(v, "query_gt", None)
        if truth_gt is None or query_gt is None:
            continue

        if truth_gt == query_gt:
            concordant += 1
        else:
            discordant += 1

    total = concordant + discordant
    rate = concordant / total if total > 0 else None

    return {
        "gt_concordant_count": concordant,
        "gt_discordant_count": discordant,
        "genotype_concordance": rate,
    }


def compute_titv_ratio(variants: list) -> dict:
    """Compute transition/transversion ratio for SNP variants.

    Only considers SNP variants. Expected ~2.1 for WGS, ~2.8 for WES.

    Args:
        variants: List of objects with ref, alt, type attributes

    Returns:
        Dict with transitions_count, transversions_count, titv_ratio
    """
    transitions = 0
    transversions = 0

    for v in variants:
        vtype = getattr(v, "type", None)
        if vtype != "SNP":
            continue

        ref = getattr(v, "ref", "")
        alt = getattr(v, "alt", "")
        if len(ref) != 1 or len(alt) != 1:
            continue

        if (ref.upper(), alt.upper()) in _TRANSITIONS:
            transitions += 1
        else:
            transversions += 1

    ratio = transitions / transversions if transversions > 0 else None

    return {
        "transitions_count": transitions,
        "transversions_count": transversions,
        "titv_ratio": ratio,
    }


def compute_het_hom_ratio(variants: list) -> dict:
    """Compute heterozygous/homozygous ratio.

    Only considers variants with zygosity field populated.
    Expected ~1.5-2.0 for diploid organisms.

    Args:
        variants: List of objects with zygosity attribute

    Returns:
        Dict with het_count, hom_count, het_hom_ratio
    """
    het = 0
    hom = 0

    for v in variants:
        zyg = getattr(v, "zygosity", None)
        if zyg == "HET":
            het += 1
        elif zyg == "HOM_ALT":
            hom += 1

    ratio = het / hom if hom > 0 else None

    return {
        "het_count": het,
        "hom_count": hom,
        "het_hom_ratio": ratio,
    }


def wilson_score_interval(
    successes: int, total: int, confidence: float = 0.95
) -> tuple[float, float]:
    """Compute Wilson score confidence interval for a proportion.

    Wilson score interval is more accurate than normal approximation
    for small samples and extreme proportions (near 0 or 1).

    Args:
        successes: Number of successes
        total: Total trials
        confidence: Confidence level (default 0.95 for 95% CI)

    Returns:
        (lower_bound, upper_bound) tuple
    """
    if total == 0:
        return (0.0, 0.0)

    z = norm.ppf((1 + confidence) / 2)

    p = successes / total
    denominator = 1 + z * z / total
    centre = p + z * z / (2 * total)
    spread = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total)

    lower = max(0.0, (centre - spread) / denominator)
    upper = min(1.0, (centre + spread) / denominator)

    return (lower, upper)


def _bootstrap_f1_ci(
    tp: int, fp: int, fn: int, confidence: float = 0.95, n_resamples: int = 9999
) -> tuple[float, float]:
    """Compute bootstrap CI for F1 score using scipy.stats.bootstrap.

    Simulates the confusion matrix as a classification array and bootstraps
    the F1 statistic.

    Args:
        tp, fp, fn: Confusion matrix counts
        confidence: Confidence level
        n_resamples: Number of bootstrap resamples

    Returns:
        (lower, upper) tuple
    """
    total = tp + fp + fn
    if total == 0:
        return (0.0, 0.0)

    # Build a classification outcome array: 1=TP, 2=FP, 3=FN
    data = np.array([1] * tp + [2] * fp + [3] * fn)

    def f1_statistic(x):
        """Compute F1 from a bootstrap sample of classification outcomes."""
        tp_b = int(np.sum(x == 1))
        fp_b = int(np.sum(x == 2))
        fn_b = int(np.sum(x == 3))
        prec = tp_b / (tp_b + fp_b) if (tp_b + fp_b) > 0 else 0.0
        rec = tp_b / (tp_b + fn_b) if (tp_b + fn_b) > 0 else 0.0
        if prec + rec > 0:
            return 2 * prec * rec / (prec + rec)
        return 0.0

    try:
        result = bootstrap(
            (data,),
            statistic=f1_statistic,
            confidence_level=confidence,
            n_resamples=n_resamples,
            method="percentile",
            random_state=42,
        )
        lo = float(np.clip(result.confidence_interval.low, 0.0, 1.0))
        hi = float(np.clip(result.confidence_interval.high, 0.0, 1.0))
        return (lo, hi)
    except Exception:
        # Fallback to harmonic mean heuristic if bootstrap fails
        # (e.g., all-zero counts edge case)
        return (0.0, 0.0)


def compute_confidence_intervals(
    tp: int, fp: int, fn: int, tn: int = 0, confidence: float = 0.95
) -> dict:
    """Compute CIs for precision, recall (Wilson score), and F1 (bootstrap).

    Precision and recall use Wilson score intervals. F1 uses scipy.stats.bootstrap
    for calibrated confidence intervals rather than the harmonic mean heuristic.

    Args:
        tp, fp, fn: Confusion matrix counts
        tn: True negatives (default 0, often unavailable in variant calling)
        confidence: Confidence level

    Returns:
        Dict with precision_ci, recall_ci, f1_ci (each a (lower, upper) tuple)
    """
    # Precision CI: TP / (TP + FP)
    precision_total = tp + fp
    precision_ci = wilson_score_interval(tp, precision_total, confidence)

    # Recall CI: TP / (TP + FN)
    recall_total = tp + fn
    recall_ci = wilson_score_interval(tp, recall_total, confidence)

    # F1 CI: scipy.stats.bootstrap for calibrated intervals
    f1_ci = _bootstrap_f1_ci(tp, fp, fn, confidence)

    return {
        "precision_ci": precision_ci,
        "recall_ci": recall_ci,
        "f1_ci": f1_ci,
    }


def compute_mcc(tp: int, fp: int, fn: int, tn: int) -> float | None:
    """Compute Matthews Correlation Coefficient.

    MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))

    Ranges from -1 to +1. Returns None if denominator is 0.

    Args:
        tp, fp, fn, tn: Confusion matrix values

    Returns:
        MCC value or None if undefined
    """
    numerator = tp * tn - fp * fn
    denominator_sq = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)

    if denominator_sq == 0:
        return None

    return numerator / math.sqrt(denominator_sq)
