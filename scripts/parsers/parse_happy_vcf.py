"""Parse hap.py output VCF into a variants TSV.

Called via Snakemake `script:` directive.
Reads the hap.py annotated VCF and the original query VCF (for annotation
fallback), and writes a tab-separated file with one row per classified variant.
"""

import logging
import math
from typing import TYPE_CHECKING

import pandas as pd
import pysam

from vcf_utils import build_query_annotations, safe_format

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


def parse_happy_vcf(happy_vcf_path, query_vcf_path):
    """Parse hap.py output VCF and extract classified variant records."""
    query_annots = build_query_annotations(query_vcf_path)
    records = []

    with pysam.VariantFile(str(happy_vcf_path)) as vcf:
        for rec in vcf:
            truth_data = rec.samples["TRUTH"]
            query_data = rec.samples["QUERY"]

            truth_bd = safe_format(truth_data, "BD")
            query_bd = safe_format(query_data, "BD")
            truth_bvt = safe_format(truth_data, "BVT")
            query_bvt = safe_format(query_data, "BVT")

            # Determine classification
            if truth_bd == "FN":
                classification = "FN"
                variant_type = truth_bvt or "UNKNOWN"
            elif query_bd in ("TP", "FP"):
                classification = query_bd
                variant_type = query_bvt or "UNKNOWN"
            else:
                continue  # Not assessed

            for alt in rec.alts or []:
                if alt == "*":
                    continue  # Skip spanning deletions

                key = (rec.chrom, rec.pos, rec.ref, alt)

                # Extract genotype (GT) from both samples
                truth_gt_raw = safe_format(truth_data, "GT")
                query_gt_raw = safe_format(query_data, "GT")
                truth_gt_str = (
                    "/".join(str(a) for a in truth_gt_raw)
                    if isinstance(truth_gt_raw, tuple)
                    else str(truth_gt_raw) if truth_gt_raw is not None else None
                )
                query_gt_str = (
                    "/".join(str(a) for a in query_gt_raw)
                    if isinstance(query_gt_raw, tuple)
                    else str(query_gt_raw) if query_gt_raw is not None else None
                )

                # Extract annotations from hap.py VCF query sample
                dp = safe_format(query_data, "DP")
                gq = safe_format(query_data, "GQ")
                af = safe_format(query_data, "AF")
                mq = rec.info.get("MQ", None)
                qual = rec.qual
                filt = ";".join(rec.filter.keys()) if rec.filter else None

                # Fallback to original query VCF
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
                        "truth_gt": truth_gt_str,
                        "query_gt": query_gt_str,
                        "indel_size": indel_size,
                    }
                )

    return records


# --- Main execution via Snakemake ---
log.info("Parsing hap.py VCF: %s", snakemake.input.vcf)
log.info("Query VCF for annotation fallback: %s", snakemake.input.query_vcf)

records = parse_happy_vcf(snakemake.input.vcf, snakemake.input.query_vcf)

df = pd.DataFrame(records)
df.to_csv(snakemake.output.tsv, sep="\t", index=False)

log.info("Wrote %d variant records to %s", len(df), snakemake.output.tsv)
