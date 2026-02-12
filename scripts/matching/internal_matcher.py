"""Internal somatic variant matcher for position + allele matching.

Provides a fallback matching strategy when external tools (som.py, RTG)
are not suitable. Matches truth and query variants by genomic position
and alleles, with configurable VAF concordance checking.

Can be called as a Snakemake script or used as a library.
"""

import logging
from pathlib import Path

import pandas as pd
import pysam

from utils import safe_float, safe_int
from vcf_utils import determine_variant_type

log = logging.getLogger(__name__)


def _extract_variants_from_vcf(vcf_path, sample_name=None):
    """Extract variant records from a VCF file.

    Args:
        vcf_path: Path to VCF file
        sample_name: Sample name to extract FORMAT fields from

    Returns:
        Dictionary mapping (chrom, pos, ref, alt) to variant info
    """
    variants = {}

    with pysam.VariantFile(str(vcf_path)) as vcf:
        samples = list(vcf.header.samples)
        target_sample = sample_name if sample_name in samples else (samples[0] if samples else None)

        for rec in vcf:
            for alt in rec.alts or []:
                if alt == "*":
                    continue

                key = (rec.chrom, rec.pos, rec.ref, alt)
                sample_data = rec.samples[target_sample] if target_sample else None

                # Extract AF
                af = None
                if sample_data:
                    try:
                        af_val = sample_data.get("AF")
                        if af_val is not None:
                            af = float(af_val[0]) if isinstance(af_val, tuple) else float(af_val)
                    except (ValueError, TypeError, IndexError):
                        pass

                # Extract DP
                dp = None
                if sample_data:
                    try:
                        dp = int(sample_data.get("DP", 0) or 0)
                    except (ValueError, TypeError):
                        pass

                # Extract GQ
                gq = None
                if sample_data:
                    try:
                        gq = int(sample_data.get("GQ", 0) or 0)
                    except (ValueError, TypeError):
                        pass

                qual = safe_float(rec.qual)
                filt = ";".join(rec.filter.keys()) if rec.filter else None

                variants[key] = {
                    "chrom": rec.chrom,
                    "pos": rec.pos,
                    "ref": rec.ref,
                    "alt": alt,
                    "af": safe_float(af),
                    "dp": safe_int(dp),
                    "gq": safe_int(gq),
                    "qual": qual,
                    "filter": filt,
                    "mq": safe_float(rec.info.get("MQ", None)),
                }

    return variants


def _build_variant_record(chrom, pos, ref, alt, variant_type, classification, variant_data, indel_size, is_query=True):
    """Build a variant record dict for output.

    Args:
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternate allele
        variant_type: SNP, INDEL, or COMPLEX
        classification: TP, FP, FN, etc.
        variant_data: Dict with dp, mq, af, gq, qual, filter
        indel_size: Indel size for stratification (None for non-indels)
        is_query: If True, populate somatic fields from variant_data

    Returns:
        Variant record dict
    """
    record = {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "type": variant_type,
        "classification": classification,
        "dp": variant_data.get("dp"),
        "mq": variant_data.get("mq"),
        "af": variant_data.get("af"),
        "gq": variant_data.get("gq"),
        "qual": variant_data.get("qual"),
        "filter_status": variant_data.get("filter"),
        "tumor_dp": variant_data.get("dp") if is_query else None,
        "tumor_af": variant_data.get("af") if is_query else None,
        "normal_dp": None,
        "normal_af": None,
        "somatic_quality": variant_data.get("qual") if is_query else None,
        "caller_filter": variant_data.get("filter") if is_query else None,
        "indel_size": indel_size,
    }
    return record


def _compute_classification_counts(records):
    """Compute classification counts and derived metrics.

    Args:
        records: List of variant record dicts

    Returns:
        Tuple of (tp, fp, fn, sensitivity, precision, f1)
    """
    from utils import compute_classification_metrics

    tp = sum(1 for r in records if r["classification"] in ("TP", "TP_VAF_DISCORDANT"))
    fp = sum(1 for r in records if r["classification"] == "FP")
    fn = sum(1 for r in records if r["classification"] == "FN")
    prf = compute_classification_metrics(tp, fp, fn)
    # Return 0.0 for None to maintain backward compatibility with metric TSV output
    sensitivity = prf["recall"] if prf["recall"] is not None else 0.0
    precision = prf["precision"] if prf["precision"] is not None else 0.0
    f1 = prf["f1"] if prf["f1"] is not None else 0.0
    return tp, fp, fn, sensitivity, precision, f1


def match_variants(
    truth_vcf_path,
    query_vcf_path,
    vaf_tolerance=0.10,
    position_window=0,
    require_allele_match=True,
    track_partial_matches=True,
    truth_sample=None,
    query_sample=None,
):
    """Match truth and query somatic variants by position and allele.

    Args:
        truth_vcf_path: Path to truth VCF
        query_vcf_path: Path to query VCF
        vaf_tolerance: Maximum absolute VAF difference for concordance
        position_window: Position matching window in bp (0 = exact)
        require_allele_match: Whether to require exact REF+ALT match
        track_partial_matches: Whether to track partial matches
        truth_sample: Sample name in truth VCF
        query_sample: Sample name in query VCF

    Returns:
        Tuple of (variant_records, metric_records)
    """
    truth_variants = _extract_variants_from_vcf(truth_vcf_path, truth_sample)
    query_variants = _extract_variants_from_vcf(query_vcf_path, query_sample)

    log.info(
        "Internal matcher: %d truth variants, %d query variants",
        len(truth_variants),
        len(query_variants),
    )

    records = []
    matched_truth = set()
    matched_query = set()

    # Match query variants against truth
    for q_key, q_var in query_variants.items():
        q_chrom, q_pos, q_ref, q_alt = q_key

        best_match = None
        best_match_key = None

        # Search for matching truth variants
        for t_key, t_var in truth_variants.items():
            t_chrom, t_pos, t_ref, t_alt = t_key

            # Check chromosome
            if q_chrom != t_chrom:
                continue

            # Check position (with optional window)
            if abs(q_pos - t_pos) > position_window:
                continue

            # Check alleles
            if require_allele_match and (q_ref != t_ref or q_alt != t_alt):
                if track_partial_matches and q_pos == t_pos:
                    # Position match but allele mismatch
                    best_match = "PARTIAL"
                    best_match_key = t_key
                continue

            # Position + allele match found
            best_match_key = t_key

            # Check VAF concordance
            t_af = safe_float(t_var.get("af"))
            q_af = safe_float(q_var.get("af"))

            if t_af is not None and q_af is not None:
                if abs(t_af - q_af) <= vaf_tolerance:
                    best_match = "TP"
                else:
                    best_match = "TP_VAF_DISCORDANT"
            else:
                best_match = "TP"  # No VAF info, assume position+allele match is TP

            break  # Found exact match, stop searching

        variant_type = determine_variant_type(q_ref, q_alt)
        indel_size = abs(len(q_alt) - len(q_ref)) if variant_type == "INDEL" else None

        if best_match in ("TP", "TP_VAF_DISCORDANT"):
            classification = best_match
            matched_truth.add(best_match_key)
            matched_query.add(q_key)
        elif best_match == "PARTIAL" and track_partial_matches:
            classification = "PARTIAL"
            matched_query.add(q_key)
        else:
            classification = "FP"

        records.append(_build_variant_record(
            q_chrom, q_pos, q_ref, q_alt, variant_type,
            classification, q_var, indel_size, is_query=True,
        ))

    # Add FN records for unmatched truth variants
    for t_key, t_var in truth_variants.items():
        if t_key not in matched_truth:
            t_chrom, t_pos, t_ref, t_alt = t_key
            variant_type = determine_variant_type(t_ref, t_alt)
            indel_size = abs(len(t_alt) - len(t_ref)) if variant_type == "INDEL" else None

            records.append(_build_variant_record(
                t_chrom, t_pos, t_ref, t_alt, variant_type,
                "FN", t_var, indel_size, is_query=False,
            ))

    # Compute metrics
    tp, fp, fn, sensitivity, precision, f1 = _compute_classification_counts(records)
    vaf_disc = sum(1 for r in records if r["classification"] == "TP_VAF_DISCORDANT")
    partial = sum(1 for r in records if r["classification"] == "PARTIAL")

    metric_records = [
        {
            "variant_type": "ALL",
            "stratification": None,
            "precision": precision,
            "recall": sensitivity,
            "f1": f1,
            "tp_count": tp,
            "fp_count": fp,
            "fn_count": fn,
        }
    ]

    # Per variant type metrics
    for vtype in ("SNP", "INDEL", "COMPLEX"):
        vtype_records = [r for r in records if r["type"] == vtype]
        if not vtype_records:
            continue
        vt_tp, vt_fp, vt_fn, vt_sens, vt_prec, vt_f1 = _compute_classification_counts(vtype_records)

        metric_records.append(
            {
                "variant_type": vtype,
                "stratification": None,
                "precision": vt_prec,
                "recall": vt_sens,
                "f1": vt_f1,
                "tp_count": vt_tp,
                "fp_count": vt_fp,
                "fn_count": vt_fn,
            }
        )

    log.info(
        "Internal matcher results: TP=%d FP=%d FN=%d VAF_DISC=%d PARTIAL=%d",
        tp - vaf_disc,
        fp,
        fn,
        vaf_disc,
        partial,
    )

    return records, metric_records


def run_internal_matcher_snakemake():
    """Entry point for Snakemake script execution."""
    from typing import TYPE_CHECKING as _TC

    if _TC:
        from snakemake.script import Snakemake
        snakemake_var: Snakemake
    else:
        snakemake_var = snakemake  # type: ignore  # noqa: F821

    logging.basicConfig(
        filename=snakemake_var.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    truth_vcf = snakemake_var.input.truth_vcf
    query_vcf = snakemake_var.input.query_vcf

    vaf_tolerance = snakemake_var.params.get("vaf_tolerance", 0.10)
    position_window = snakemake_var.params.get("position_window", 0)
    require_allele_match = snakemake_var.params.get("require_allele_match", True)
    track_partial = snakemake_var.params.get("track_partial_matches", True)

    log.info("Internal matcher: truth=%s query=%s", truth_vcf, query_vcf)

    variant_records, metric_records = match_variants(
        truth_vcf,
        query_vcf,
        vaf_tolerance=vaf_tolerance,
        position_window=position_window,
        require_allele_match=require_allele_match,
        track_partial_matches=track_partial,
    )

    # Write variants TSV
    variants_df = pd.DataFrame(variant_records)
    variants_df.to_csv(snakemake_var.output.variants_tsv, sep="\t", index=False)
    log.info("Wrote %d variant records", len(variants_df))

    # Write metrics TSV
    metrics_df = pd.DataFrame(metric_records)
    metrics_df.to_csv(snakemake_var.output.metrics_tsv, sep="\t", index=False)
    log.info("Wrote %d metric records", len(metrics_df))
