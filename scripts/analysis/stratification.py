"""Stratification engine for variant-level stratified metric computation.

Assigns variants to strata across multiple dimensions (variant_class,
indel_size, functional_impact, gene_panel, zygosity, vaf_bins) and
computes precision/recall/F1 independently for each stratum.

Can be used as a library or called via Snakemake script.
"""

import logging
from collections import defaultdict

log = logging.getLogger(__name__)

# Default stratification config matching config/stratifications.yaml
DEFAULT_INDEL_BINS = [
    {"name": "1bp", "min": 1, "max": 1},
    {"name": "2-5bp", "min": 2, "max": 5},
    {"name": "6-15bp", "min": 6, "max": 15},
    {"name": "16-50bp", "min": 16, "max": 50},
    {"name": ">50bp", "min": 51, "max": None},
]

DEFAULT_VAF_BINS = [
    {"name": "<0.05", "min": 0.0, "max": 0.05},
    {"name": "0.05-0.10", "min": 0.05, "max": 0.10},
    {"name": "0.10-0.20", "min": 0.10, "max": 0.20},
    {"name": "0.20-0.50", "min": 0.20, "max": 0.50},
    {"name": ">0.50", "min": 0.50, "max": 1.0},
]


def _assign_variant_class(variant):
    """Assign variant to variant_class stratum."""
    vtype = getattr(variant, "type", None)
    if vtype in ("SNP", "INDEL", "COMPLEX", "SV"):
        return vtype
    return None


def _assign_indel_size(variant, bins=None):
    """Assign variant to indel_size stratum based on indel_size field."""
    if getattr(variant, "type", None) != "INDEL":
        return None

    size = getattr(variant, "indel_size", None)
    if size is None:
        return None

    bins = bins or DEFAULT_INDEL_BINS
    for b in bins:
        b_min = b["min"]
        b_max = b["max"]
        if b_max is None:
            if size >= b_min:
                return b["name"]
        elif b_min <= size <= b_max:
            return b["name"]
    return None


def _assign_functional_impact(variant):
    """Assign variant to functional_impact stratum based on VEP impact."""
    impact = getattr(variant, "impact", None)
    if impact in ("HIGH", "MODERATE", "LOW", "MODIFIER"):
        return impact
    return None


def _assign_gene_panel(variant, gene_set_names=None):
    """Assign variant to gene_panel strata.

    Returns list of panel names the variant belongs to, or ["off_panel"].
    """
    gene_sets = getattr(variant, "gene_sets", [])
    if gene_sets:
        return [gs.name for gs in gene_sets]
    return ["off_panel"]


def _assign_zygosity(variant):
    """Assign variant to zygosity stratum (germline only)."""
    zyg = getattr(variant, "zygosity", None)
    if zyg in ("HET", "HOM_ALT"):
        return zyg
    return None


def _assign_vaf_bin(variant, bins=None):
    """Assign variant to VAF bin stratum (somatic only)."""
    af = getattr(variant, "tumor_af", None) or getattr(variant, "af", None)
    if af is None:
        return None

    bins = bins or DEFAULT_VAF_BINS
    for b in bins:
        b_min = b["min"]
        b_max = b["max"]
        if b_max is None:
            if af >= b_min:
                return b["name"]
        elif b_min <= af < b_max:
            return b["name"]
        elif b_max == 1.0 and af >= b_min:
            # Include 1.0 in the last bin
            return b["name"]
    return None


def _assign_low_complexity(variant):
    """Assign variant to low_complexity stratum based on region annotation."""
    val = getattr(variant, "in_low_complexity", None)
    if val is None:
        return None
    return "low_complexity" if val else "non_lcr"


def _assign_gc_content(variant, bins=None):
    """Assign variant to gc_content stratum bin."""
    gc = getattr(variant, "gc_content", None)
    if gc is None:
        return None
    if bins is None:
        return None
    for b in bins:
        b_min = b["min"]
        b_max = b["max"]
        if b_max is None:
            if gc >= b_min:
                return b["name"]
        elif b_min <= gc < b_max:
            return b["name"]
        elif b_max >= 1.0 and gc >= b_min:
            return b["name"]
    return None


def _assign_segdup(variant):
    """Assign variant to segdup stratum based on region annotation."""
    val = getattr(variant, "in_segdup", None)
    if val is None:
        return None
    return "in_segdup" if val else "non_segdup"


def _assign_mappability(variant, bins=None):
    """Assign variant to mappability stratum bin."""
    score = getattr(variant, "mappability_score", None)
    if score is None:
        return None
    if bins is None:
        return None
    for b in bins:
        b_min = b["min"]
        b_max = b["max"]
        if b_max is None:
            if score >= b_min:
                return b["name"]
        elif b_min <= score < b_max:
            return b["name"]
        elif b_max > 1.0 and score >= b_min:
            return b["name"]
    return None


def _assign_coverage_depth(variant, bins=None):
    """Assign variant to coverage_depth stratum based on DP field."""
    dp = getattr(variant, "dp", None)
    if dp is None:
        return None
    if bins is None:
        return None
    for b in bins:
        b_min = b["min"]
        b_max = b["max"]
        if b_max is None:
            if dp >= b_min:
                return b["name"]
        elif b_min <= dp < b_max:
            return b["name"]
    return None


def _compute_stratum_metrics(variants):
    """Compute TP/FP/FN counts, precision/recall/F1, CIs, MCC, and GT concordance."""
    from analysis.extended_metrics import (
        compute_confidence_intervals,
        compute_genotype_concordance,
        compute_mcc,
    )
    from utils import compute_classification_metrics

    tp = sum(1 for v in variants if getattr(v, "classification", "") in ("TP", "TP_VAF_DISCORDANT"))
    fp = sum(1 for v in variants if getattr(v, "classification", "") == "FP")
    fn = sum(1 for v in variants if getattr(v, "classification", "") == "FN")

    prf = compute_classification_metrics(tp, fp, fn)
    precision = prf["precision"]
    recall = prf["recall"]
    f1 = prf["f1"]

    ci_data = compute_confidence_intervals(tp, fp, fn)
    mcc = compute_mcc(tp, fp, fn, tn=0)
    gt_data = compute_genotype_concordance(variants)

    return {
        "tp_count": tp,
        "fp_count": fp,
        "fn_count": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "total_variants": len(variants),
        "precision_ci_lower": ci_data["precision_ci"][0],
        "precision_ci_upper": ci_data["precision_ci"][1],
        "recall_ci_lower": ci_data["recall_ci"][0],
        "recall_ci_upper": ci_data["recall_ci"][1],
        "f1_ci_lower": ci_data["f1_ci"][0],
        "f1_ci_upper": ci_data["f1_ci"][1],
        "mcc": mcc,
        "genotype_concordance": gt_data["genotype_concordance"],
    }


class StratificationEngine:
    """Engine for stratifying variants and computing per-stratum metrics.

    Args:
        config: Stratification configuration dictionary (from stratifications.yaml)
        mode: Execution mode ("germline" or "somatic")
        gene_set_names: List of configured gene set names
        min_variants: Minimum variants per stratum for high-confidence flag
    """

    def __init__(self, config=None, mode="germline", gene_set_names=None, min_variants=20):
        self.config = config or {}
        self.mode = mode
        self.gene_set_names = gene_set_names or []
        self.min_variants = self.config.get("min_variants_per_stratum", min_variants)
        self.dimensions = self.config.get("dimensions", {})

    def _is_dimension_enabled(self, dim_name):
        """Check if a dimension is enabled and applies to the current mode."""
        dim_config = self.dimensions.get(dim_name, {})
        if not dim_config.get("enabled", True):
            return False
        dim_mode = dim_config.get("mode")
        if dim_mode and dim_mode != self.mode:
            return False
        return True

    def stratify(self, variants):
        """Assign each variant to its applicable strata across all dimensions.

        Args:
            variants: List of Variant objects (or objects with matching attributes)

        Returns:
            Dictionary mapping (dimension, stratum) to list of variants
        """
        strata = defaultdict(list)

        for variant in variants:
            # variant_class
            if self._is_dimension_enabled("variant_class"):
                vc = _assign_variant_class(variant)
                if vc:
                    strata[("variant_class", vc)].append(variant)

            # indel_size
            if self._is_dimension_enabled("indel_size"):
                indel_bins = self.dimensions.get("indel_size", {}).get("bins", DEFAULT_INDEL_BINS)
                isize = _assign_indel_size(variant, indel_bins)
                if isize:
                    strata[("indel_size", isize)].append(variant)

            # functional_impact
            if self._is_dimension_enabled("functional_impact"):
                fi = _assign_functional_impact(variant)
                if fi:
                    strata[("functional_impact", fi)].append(variant)

            # gene_panel
            if self._is_dimension_enabled("gene_panel"):
                panels = _assign_gene_panel(variant, self.gene_set_names)
                for panel in panels:
                    strata[("gene_panel", panel)].append(variant)

            # zygosity (germline only)
            if self._is_dimension_enabled("zygosity"):
                zyg = _assign_zygosity(variant)
                if zyg:
                    strata[("zygosity", zyg)].append(variant)

            # vaf_bins (somatic only)
            if self._is_dimension_enabled("vaf_bins"):
                vaf_bins = self.dimensions.get("vaf_bins", {}).get("bins", DEFAULT_VAF_BINS)
                vb = _assign_vaf_bin(variant, vaf_bins)
                if vb:
                    strata[("vaf_bins", vb)].append(variant)

            # low_complexity (requires region annotation)
            if self._is_dimension_enabled("low_complexity"):
                lcr = _assign_low_complexity(variant)
                if lcr:
                    strata[("low_complexity", lcr)].append(variant)

            # gc_content (requires region annotation or reference)
            if self._is_dimension_enabled("gc_content"):
                gc_bins = self.dimensions.get("gc_content", {}).get("bins")
                gc = _assign_gc_content(variant, gc_bins)
                if gc:
                    strata[("gc_content", gc)].append(variant)

            # segdup (requires region annotation)
            if self._is_dimension_enabled("segdup"):
                sd = _assign_segdup(variant)
                if sd:
                    strata[("segdup", sd)].append(variant)

            # mappability (requires region annotation)
            if self._is_dimension_enabled("mappability"):
                mapp_bins = self.dimensions.get("mappability", {}).get("bins")
                mp = _assign_mappability(variant, mapp_bins)
                if mp:
                    strata[("mappability", mp)].append(variant)

            # coverage_depth (uses variant DP field)
            if self._is_dimension_enabled("coverage_depth"):
                cov_bins = self.dimensions.get("coverage_depth", {}).get("bins")
                cd = _assign_coverage_depth(variant, cov_bins)
                if cd:
                    strata[("coverage_depth", cd)].append(variant)

        log.info(
            "Stratified %d variants across %d strata",
            len(variants),
            len(strata),
        )
        return dict(strata)

    def compute_metrics(self, strata):
        """Compute metrics for each stratum.

        Args:
            strata: Dictionary from stratify() - (dimension, stratum) -> [variants]

        Returns:
            List of metric dictionaries ready for database insertion
        """
        results = []

        for (dimension, stratum), variants in sorted(strata.items()):
            metrics = _compute_stratum_metrics(variants)
            metrics["dimension"] = dimension
            metrics["stratum"] = stratum
            metrics["variant_type"] = "ALL"
            metrics["low_confidence"] = metrics["total_variants"] < self.min_variants
            results.append(metrics)

        log.info("Computed metrics for %d strata", len(results))
        return results
