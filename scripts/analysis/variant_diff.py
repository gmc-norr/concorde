"""Variant-level diff engine for continuous verification (Spec S13.2).

Compares two sets of variant classifications (baseline vs verification)
and produces variant-level transitions (gained, lost, changed).
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)


def compute_variant_diff(
    baseline_variants: list,
    verification_variants: list,
) -> list[dict]:
    """Compute variant-level diff between baseline and verification runs.

    Args:
        baseline_variants: List of baseline Variant objects (or dicts with
            chrom, pos, ref, alt, classification attributes)
        verification_variants: List of verification Variant objects

    Returns:
        List of transition dicts with:
        - chrom, pos, ref, alt
        - baseline_classification (or None if new)
        - verification_classification (or None if lost)
        - transition (e.g. "TP→FN", "None→TP")
    """
    def _variant_key(v):
        if isinstance(v, dict):
            return (v["chrom"], v["pos"], v["ref"], v["alt"])
        return (v.chrom, v.pos, v.ref, v.alt)

    def _classification(v):
        if isinstance(v, dict):
            return v.get("classification")
        return getattr(v, "classification", None)

    # Build maps
    baseline_map = {}
    for v in baseline_variants:
        key = _variant_key(v)
        baseline_map[key] = _classification(v)

    verification_map = {}
    for v in verification_variants:
        key = _variant_key(v)
        verification_map[key] = _classification(v)

    all_keys = set(baseline_map.keys()) | set(verification_map.keys())
    transitions = []

    for key in sorted(all_keys):
        base_cls = baseline_map.get(key)
        verify_cls = verification_map.get(key)

        if base_cls == verify_cls:
            continue  # No change

        chrom, pos, ref, alt = key
        transition_str = f"{base_cls or 'None'}→{verify_cls or 'None'}"

        transitions.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "baseline_classification": base_cls,
                "verification_classification": verify_cls,
                "transition": transition_str,
            }
        )

    log.info(
        "Computed variant diff: %d transitions from %d total variants",
        len(transitions),
        len(all_keys),
    )
    return transitions


def assign_strata_to_transitions(
    transitions: list[dict],
    verification_variants: list,
    strata: dict,
) -> list[dict]:
    """Annotate transitions with the strata each variant belongs to.

    Args:
        transitions: List of transition dicts from compute_variant_diff
        verification_variants: Verification variant objects
        strata: Strata dict from StratificationEngine.stratify()

    Returns:
        transitions list with added 'strata' field (list of (dim, stratum) tuples)
    """
    # Build variant key → strata lookup
    strata_lookup = {}
    for (dim, stratum_val), variants in strata.items():
        for v in variants:
            key = (v.chrom, v.pos, v.ref, v.alt) if not isinstance(v, dict) else (v["chrom"], v["pos"], v["ref"], v["alt"])
            strata_lookup.setdefault(key, []).append((dim, stratum_val))

    for t in transitions:
        key = (t["chrom"], t["pos"], t["ref"], t["alt"])
        t["strata"] = strata_lookup.get(key, [])

    return transitions
