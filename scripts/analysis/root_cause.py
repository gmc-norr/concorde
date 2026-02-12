"""Root cause analysis for biological differences (Spec S14).

Collects evidence (coverage, mapping quality, filter, quality score,
tool version) for variant transitions and scores explanatory power.
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)

# Evidence score definitions (S14.2)
STRONG = "STRONG"
MODERATE = "MODERATE"
WEAK = "WEAK"
NONE = "NONE"

# Thresholds for scoring
COVERAGE_DROP_STRONG = 0.50  # >50% drop = STRONG
COVERAGE_DROP_MODERATE = 0.25  # >25% drop = MODERATE
MQ_DROP_STRONG = 15  # >15 point MQ drop = STRONG
MQ_DROP_MODERATE = 5  # >5 point MQ drop = MODERATE
QUAL_DROP_STRONG = 0.50  # >50% relative QUAL drop = STRONG
QUAL_DROP_MODERATE = 0.25  # >25% relative QUAL drop = MODERATE


def collect_evidence(
    transition: dict,
    baseline_variant: dict | None = None,
    verification_variant: dict | None = None,
    tool_changes: list[dict] | None = None,
    bam_available: bool = False,
) -> list[dict]:
    """Collect root cause evidence for a single variant transition.

    Args:
        transition: Transition dict with chrom, pos, ref, alt, transition
        baseline_variant: Baseline variant data (dict with dp, mq, qual, filter_status)
        verification_variant: Verification variant data
        tool_changes: List of tool version changes [{tool, old_version, new_version}]
        bam_available: Whether BAM files were available for analysis

    Returns:
        List of evidence dicts with category, values, and score
    """
    evidence = []

    base = baseline_variant or {}
    verify = verification_variant or {}

    # Coverage evidence
    evidence.append(_score_coverage(base, verify, bam_available))

    # Mapping quality evidence
    evidence.append(_score_mapping_quality(base, verify, bam_available))

    # Filter evidence
    evidence.append(_score_filter_change(base, verify))

    # Quality score evidence
    evidence.append(_score_quality_change(base, verify))

    # Tool version evidence
    evidence.append(_score_tool_versions(tool_changes))

    return [e for e in evidence if e is not None]


def _score_coverage(base: dict, verify: dict, bam_available: bool) -> dict:
    """Score coverage evidence."""
    base_dp = base.get("dp")
    verify_dp = verify.get("dp")

    if not bam_available and base_dp is None and verify_dp is None:
        return {
            "category": "coverage",
            "baseline_value": None,
            "verification_value": None,
            "change_description": "BAM not available, coverage not assessed",
            "change_magnitude": None,
            "score": NONE,
        }

    if base_dp is None or verify_dp is None:
        return {
            "category": "coverage",
            "baseline_value": str(base_dp),
            "verification_value": str(verify_dp),
            "change_description": "Coverage data incomplete",
            "change_magnitude": None,
            "score": NONE,
        }

    if base_dp == 0:
        relative_change = 1.0 if verify_dp > 0 else 0.0
    else:
        relative_change = abs(verify_dp - base_dp) / base_dp

    score = NONE
    if relative_change >= COVERAGE_DROP_STRONG:
        score = STRONG
    elif relative_change >= COVERAGE_DROP_MODERATE:
        score = MODERATE
    elif relative_change > 0:
        score = WEAK

    return {
        "category": "coverage",
        "baseline_value": str(base_dp),
        "verification_value": str(verify_dp),
        "change_description": f"Depth changed from {base_dp} to {verify_dp} ({relative_change:.1%})",
        "change_magnitude": relative_change,
        "score": score,
    }


def _score_mapping_quality(base: dict, verify: dict, bam_available: bool) -> dict:
    """Score mapping quality evidence."""
    base_mq = base.get("mq")
    verify_mq = verify.get("mq")

    if not bam_available and base_mq is None and verify_mq is None:
        return {
            "category": "mapping_quality",
            "baseline_value": None,
            "verification_value": None,
            "change_description": "BAM not available, MQ not assessed",
            "change_magnitude": None,
            "score": NONE,
        }

    if base_mq is None or verify_mq is None:
        return {
            "category": "mapping_quality",
            "baseline_value": str(base_mq),
            "verification_value": str(verify_mq),
            "change_description": "MQ data incomplete",
            "change_magnitude": None,
            "score": NONE,
        }

    mq_change = abs(verify_mq - base_mq)

    score = NONE
    if mq_change >= MQ_DROP_STRONG:
        score = STRONG
    elif mq_change >= MQ_DROP_MODERATE:
        score = MODERATE
    elif mq_change > 0:
        score = WEAK

    return {
        "category": "mapping_quality",
        "baseline_value": str(base_mq),
        "verification_value": str(verify_mq),
        "change_description": f"MQ changed from {base_mq} to {verify_mq} (Δ{mq_change:.1f})",
        "change_magnitude": mq_change,
        "score": score,
    }


def _score_filter_change(base: dict, verify: dict) -> dict:
    """Score filter evidence."""
    base_filter = base.get("filter_status", ".")
    verify_filter = verify.get("filter_status", ".")

    if base_filter == verify_filter:
        return {
            "category": "filter",
            "baseline_value": str(base_filter),
            "verification_value": str(verify_filter),
            "change_description": "No filter change",
            "change_magnitude": 0.0,
            "score": NONE,
        }

    # PASS→non-PASS or non-PASS→PASS is strong evidence
    pass_to_nonpass = base_filter == "PASS" and verify_filter != "PASS"
    nonpass_to_pass = base_filter != "PASS" and verify_filter == "PASS"

    score = STRONG if (pass_to_nonpass or nonpass_to_pass) else MODERATE

    return {
        "category": "filter",
        "baseline_value": str(base_filter),
        "verification_value": str(verify_filter),
        "change_description": f"Filter changed: {base_filter}→{verify_filter}",
        "change_magnitude": 1.0,
        "score": score,
    }


def _score_quality_change(base: dict, verify: dict) -> dict:
    """Score quality score evidence."""
    base_qual = base.get("qual")
    verify_qual = verify.get("qual")

    if base_qual is None or verify_qual is None:
        return {
            "category": "quality_score",
            "baseline_value": str(base_qual),
            "verification_value": str(verify_qual),
            "change_description": "QUAL data incomplete",
            "change_magnitude": None,
            "score": NONE,
        }

    if base_qual == 0:
        relative_change = 1.0 if verify_qual > 0 else 0.0
    else:
        relative_change = abs(verify_qual - base_qual) / base_qual

    score = NONE
    if relative_change >= QUAL_DROP_STRONG:
        score = STRONG
    elif relative_change >= QUAL_DROP_MODERATE:
        score = MODERATE
    elif relative_change > 0:
        score = WEAK

    return {
        "category": "quality_score",
        "baseline_value": str(base_qual),
        "verification_value": str(verify_qual),
        "change_description": f"QUAL changed from {base_qual} to {verify_qual} ({relative_change:.1%})",
        "change_magnitude": relative_change,
        "score": score,
    }


def _score_tool_versions(tool_changes: list[dict] | None) -> dict:
    """Score tool version evidence."""
    if not tool_changes:
        return {
            "category": "tool_version",
            "baseline_value": None,
            "verification_value": None,
            "change_description": "No tool version changes detected",
            "change_magnitude": 0.0,
            "score": NONE,
        }

    changed = [f"{t['tool']}: {t['old_version']}→{t['new_version']}" for t in tool_changes]

    return {
        "category": "tool_version",
        "baseline_value": ", ".join(t["old_version"] for t in tool_changes),
        "verification_value": ", ".join(t["new_version"] for t in tool_changes),
        "change_description": f"Tool changes: {'; '.join(changed)}",
        "change_magnitude": float(len(tool_changes)),
        "score": WEAK,
    }


def determine_likely_cause(evidence: list[dict]) -> str:
    """Determine the likely root cause from evidence items.

    Args:
        evidence: List of evidence dicts with category and score

    Returns:
        Likely cause string: coverage_drop, filter_change,
        mapping_quality_shift, quality_change, tool_update, unknown
    """
    cause_map = {
        "coverage": "coverage_drop",
        "mapping_quality": "mapping_quality_shift",
        "filter": "filter_change",
        "quality_score": "quality_change",
        "tool_version": "tool_update",
    }

    # Find the strongest evidence
    score_rank = {STRONG: 4, MODERATE: 3, WEAK: 2, NONE: 1}
    best_score = NONE
    best_category = None

    for e in evidence:
        s = e.get("score", NONE)
        if score_rank.get(s, 0) > score_rank.get(best_score, 0):
            best_score = s
            best_category = e.get("category")

    if best_score == NONE or best_category is None:
        return "unknown"

    return cause_map.get(best_category, "unknown")


def requires_review(evidence: list[dict]) -> bool:
    """Determine if a variant transition requires manual review.

    Returns True if no STRONG evidence was found (S14.3).
    """
    return not any(e.get("score") == STRONG for e in evidence)
