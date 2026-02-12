"""Difference classifier for continuous verification (Spec S13.3).

Classifies each variant transition as equivalent, numerical_drift,
or biological_difference based on envelope comparisons and variant
characteristics.
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)

FLOAT_TOLERANCE = 1e-9


def classify_transitions(
    transitions: list[dict],
    envelope_violations: list[dict],
    high_impact_variants: set[tuple] | None = None,
    regulated_panel_variants: set[tuple] | None = None,
) -> list[dict]:
    """Classify each variant transition per S13.3 rules.

    Args:
        transitions: Transition dicts from compute_variant_diff
        envelope_violations: Violations from BaselineManager.check_envelopes
        high_impact_variants: Set of (chrom, pos, ref, alt) keys for
            HIGH/MODERATE impact variants
        regulated_panel_variants: Set of (chrom, pos, ref, alt) keys for
            variants in regulated gene panels

    Returns:
        transitions list with added 'diff_class' field:
        - "equivalent": no actual difference (shouldn't appear)
        - "numerical_drift": within envelope, borderline variant
        - "biological_difference": exceeds envelope or affects important variants
    """
    high_impact = high_impact_variants or set()
    regulated = regulated_panel_variants or set()

    # Build set of violated strata for quick lookup
    violated_strata = set()
    for v in envelope_violations:
        violated_strata.add((v["dimension"], v["stratum"]))

    for t in transitions:
        key = (t["chrom"], t["pos"], t["ref"], t["alt"])
        strata = t.get("strata", [])

        # Check if any of the variant's strata had an envelope breach
        has_envelope_breach = any(s in violated_strata for s in strata)

        # Check if variant is high impact or in a regulated panel
        is_high_impact = key in high_impact
        is_regulated = key in regulated

        if has_envelope_breach or is_high_impact or is_regulated:
            t["diff_class"] = "biological_difference"
        else:
            t["diff_class"] = "numerical_drift"

    classified = _count_by_class(transitions)
    log.info(
        "Classified %d transitions: %d drift, %d biological",
        len(transitions),
        classified.get("numerical_drift", 0),
        classified.get("biological_difference", 0),
    )
    return transitions


def _count_by_class(transitions: list[dict]) -> dict[str, int]:
    """Count transitions by diff_class."""
    counts: dict[str, int] = {}
    for t in transitions:
        cls = t.get("diff_class", "unknown")
        counts[cls] = counts.get(cls, 0) + 1
    return counts


def determine_verdict(
    transitions: list[dict],
    envelope_violations: list[dict],
) -> str:
    """Determine overall verification verdict (S13.5).

    Args:
        transitions: Classified transition list
        envelope_violations: Envelope violations

    Returns:
        "PASS", "FAIL", or "REVIEW_REQUIRED"
    """
    if not transitions and not envelope_violations:
        return "PASS"

    bio_count = sum(1 for t in transitions if t.get("diff_class") == "biological_difference")

    if bio_count > 0 or len(envelope_violations) > 0:
        return "FAIL"

    if transitions:
        return "REVIEW_REQUIRED"

    return "PASS"
