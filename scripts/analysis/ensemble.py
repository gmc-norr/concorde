"""Multi-caller ensemble comparison engine.

Computes cross-caller concordance and applies ensemble filtering
methods (majority vote, intersection, union) across multiple caller runs.
"""

from __future__ import annotations

import logging
from collections import defaultdict

log = logging.getLogger(__name__)


class EnsembleEngine:
    """Engine for cross-caller comparison and ensemble filtering.

    Args:
        session: SQLAlchemy session
        config: Optional ensemble configuration dict
    """

    def __init__(self, session, config=None):
        self.session = session
        self.config = config or {}

    def find_runs_for_sample(self, sample: str) -> list:
        """Find all runs for a given sample.

        Args:
            sample: Sample identifier

        Returns:
            List of Run objects for the sample
        """
        from models.run import Run

        return (
            self.session.query(Run)
            .filter(Run.sample == sample)
            .order_by(Run.id)
            .all()
        )

    def compute_cross_caller_concordance(self, run_ids: list[int]) -> dict:
        """Compute concordance between callers for the same sample.

        Args:
            run_ids: List of Run IDs to compare

        Returns:
            Dict with callers, concordance_matrix, variants_by_caller_count
        """
        from models.run import Run
        from models.variant import Variant

        # Load variants per run
        runs = {}
        variant_sets = {}
        for run_id in run_ids:
            run = self.session.query(Run).filter(Run.id == run_id).first()
            if not run:
                continue
            runs[run_id] = run

            variants = (
                self.session.query(Variant)
                .filter(Variant.run_id == run_id)
                .all()
            )
            variant_sets[run_id] = {
                (v.chrom, v.pos, v.ref, v.alt): v.classification
                for v in variants
            }

        callers = [runs[rid].caller for rid in run_ids if rid in runs]

        # Compute pairwise concordance
        concordance_matrix = {}
        ordered_ids = [rid for rid in run_ids if rid in runs]
        for i, rid_a in enumerate(ordered_ids):
            for j, rid_b in enumerate(ordered_ids):
                if i >= j:
                    continue
                set_a = set(variant_sets[rid_a].keys())
                set_b = set(variant_sets[rid_b].keys())
                intersection = len(set_a & set_b)
                union = len(set_a | set_b)
                jaccard = intersection / union if union > 0 else 0.0
                concordance_matrix[(rid_a, rid_b)] = {
                    "shared": intersection,
                    "total": union,
                    "jaccard": jaccard,
                }

        # Count variants by caller support level
        all_variant_keys = set()
        for vs in variant_sets.values():
            all_variant_keys.update(vs.keys())

        by_count = defaultdict(int)
        for key in all_variant_keys:
            count = sum(1 for vs in variant_sets.values() if key in vs)
            by_count[count] += 1

        log.info(
            "Cross-caller concordance: %d runs, %d total variants",
            len(runs), len(all_variant_keys),
        )

        return {
            "callers": callers,
            "run_ids": ordered_ids,
            "concordance_matrix": {
                f"{k[0]}_vs_{k[1]}": v for k, v in concordance_matrix.items()
            },
            "variants_by_caller_count": dict(by_count),
            "total_unique_variants": len(all_variant_keys),
        }

    def apply_ensemble_filter(
        self, run_ids: list[int], method: str = "majority_vote"
    ) -> list[dict]:
        """Apply ensemble filtering to select consensus variants.

        Args:
            run_ids: List of Run IDs
            method: Filtering method - "majority_vote", "intersection", "union"

        Returns:
            List of variant dicts that pass the ensemble filter
        """
        from models.variant import Variant

        # Load variants per run
        variant_maps = {}
        for run_id in run_ids:
            variants = (
                self.session.query(Variant)
                .filter(Variant.run_id == run_id)
                .all()
            )
            variant_maps[run_id] = {
                (v.chrom, v.pos, v.ref, v.alt): v
                for v in variants
            }

        # Collect all variant keys with their support count
        all_keys = set()
        for vm in variant_maps.values():
            all_keys.update(vm.keys())

        n_callers = len(run_ids)
        threshold = n_callers // 2 + 1 if method == "majority_vote" else n_callers

        result = []
        for key in sorted(all_keys):
            support = sum(1 for vm in variant_maps.values() if key in vm)

            if method == "intersection" and support < n_callers:
                continue
            elif method == "majority_vote" and support < threshold:
                continue
            # "union" includes everything

            # Use the variant from the first run that has it
            variant = None
            for vm in variant_maps.values():
                if key in vm:
                    variant = vm[key]
                    break

            if variant:
                result.append({
                    "chrom": variant.chrom,
                    "pos": variant.pos,
                    "ref": variant.ref,
                    "alt": variant.alt,
                    "type": variant.type,
                    "classification": variant.classification,
                    "caller_support": support,
                    "total_callers": n_callers,
                })

        log.info(
            "Ensemble filter (%s): %d/%d variants passed",
            method, len(result), len(all_keys),
        )
        return result
