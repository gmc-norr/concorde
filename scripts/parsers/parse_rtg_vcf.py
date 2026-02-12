"""Parse rtg vcfeval output VCFs into a variants TSV.

Called via Snakemake `script:` directive.
Reads the rtg vcfeval output VCFs (tp.vcf.gz, fp.vcf.gz, fn.vcf.gz) and the
original query VCF (for annotation fallback), and writes a tab-separated file
with one row per classified variant.

rtg vcfeval outputs separate VCF files for each classification:
- tp.vcf.gz: True positive variants from the query
- fp.vcf.gz: False positive variants from the query
- fn.vcf.gz: False negative variants from the baseline/truth
"""

import logging
import math
from typing import TYPE_CHECKING

import pandas as pd
import pysam

from vcf_utils import build_query_annotations, determine_variant_type, safe_format

if TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake
else:
    snakemake = snakemake  # type: ignore  # noqa: F821

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


def parse_rtg_vcf_file(vcf_path, classification, query_annots):
    """Parse a single rtg vcfeval output VCF and extract variant records."""
    records = []

    with pysam.VariantFile(str(vcf_path)) as vcf:
        sample_name = vcf.header.samples[0] if len(vcf.header.samples) > 0 else None

        for rec in vcf:
            for alt in rec.alts or []:
                if alt == "*":
                    continue  # Skip spanning deletions

                key = (rec.chrom, rec.pos, rec.ref, alt)
                sample_data = rec.samples[sample_name] if sample_name else None

                # Extract genotype (GT) from sample
                gt_raw = safe_format(sample_data, "GT") if sample_data else None
                gt_str = (
                    "/".join(str(a) for a in gt_raw)
                    if isinstance(gt_raw, tuple)
                    else str(gt_raw) if gt_raw is not None else None
                )

                # Extract annotations from rtg VCF
                dp = safe_format(sample_data, "DP") if sample_data else None
                gq = safe_format(sample_data, "GQ") if sample_data else None
                af = safe_format(sample_data, "AF") if sample_data else None
                mq = rec.info.get("MQ", None)
                qual = rec.qual
                filt = ";".join(rec.filter.keys()) if rec.filter else None

                # Fallback to original query VCF annotations
                orig = query_annots.get(key, {})
                if dp is None:
                    dp = orig.get("dp")
                if gq is None:
                    gq = orig.get("gq")
                if af is None:
                    af = orig.get("af")
                if mq is None:
                    mq = orig.get("mq")
                if qual is None:
                    qual = orig.get("qual")
                if filt is None:
                    filt = orig.get("filter")

                # Determine variant type
                variant_type = determine_variant_type(rec.ref, alt)

                # Determine GT fields based on classification
                # TP/FP files contain query genotypes; FN files contain truth genotypes
                if classification == "FN":
                    truth_gt = gt_str
                    query_gt = None
                else:
                    truth_gt = None
                    query_gt = gt_str

                # Compute indel_size for INDEL variants
                indel_size = (
                    abs(len(alt) - len(rec.ref))
                    if variant_type == "INDEL"
                    else None
                )

                records.append(
                    {
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": rec.ref,
                        "alt": alt,
                        "type": variant_type,
                        "classification": classification,
                        "dp": int(dp) if dp is not None else None,
                        "mq": float(mq) if mq is not None else None,
                        "af": float(af) if af is not None else None,
                        "gq": int(gq) if gq is not None else None,
                        "qual": float(qual)
                        if qual is not None and not math.isnan(qual)
                        else None,
                        "filter_status": str(filt) if filt else None,
                        "truth_gt": truth_gt,
                        "query_gt": query_gt,
                        "indel_size": indel_size,
                    }
                )

    return records


def parse_rtg_vcfeval(tp_vcf, fp_vcf, fn_vcf, query_vcf):
    """Parse all rtg vcfeval output VCFs and combine into unified records."""
    query_annots = build_query_annotations(query_vcf)

    records = []

    # Parse TP variants
    log.info("Parsing TP variants from: %s", tp_vcf)
    records.extend(parse_rtg_vcf_file(tp_vcf, "TP", query_annots))

    # Parse FP variants
    log.info("Parsing FP variants from: %s", fp_vcf)
    records.extend(parse_rtg_vcf_file(fp_vcf, "FP", query_annots))

    # Parse FN variants (from baseline/truth, may have limited annotations)
    log.info("Parsing FN variants from: %s", fn_vcf)
    records.extend(parse_rtg_vcf_file(fn_vcf, "FN", query_annots))

    return records


# --- Main execution via Snakemake ---
log.info("Parsing rtg vcfeval VCFs")
log.info("  TP VCF: %s", snakemake.input.tp_vcf)
log.info("  FP VCF: %s", snakemake.input.fp_vcf)
log.info("  FN VCF: %s", snakemake.input.fn_vcf)
log.info("  Query VCF for annotation fallback: %s", snakemake.input.query_vcf)

records = parse_rtg_vcfeval(
    snakemake.input.tp_vcf,
    snakemake.input.fp_vcf,
    snakemake.input.fn_vcf,
    snakemake.input.query_vcf,
)

df = pd.DataFrame(records)
df.to_csv(snakemake.output.tsv, sep="\t", index=False)

log.info("Wrote %d variant records to %s", len(df), snakemake.output.tsv)
