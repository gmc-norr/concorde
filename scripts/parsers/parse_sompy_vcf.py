"""Parse som.py (hap.py somatic module) output VCF into a variants TSV.

Called via Snakemake `script:` directive.
Reads the som.py annotated VCF and the original query VCF (for annotation
fallback), and writes a tab-separated file with one row per classified variant.

Somatic-specific fields extracted: tumor_dp, tumor_af, normal_dp, normal_af,
somatic_quality, caller_filter.
"""

import logging
import math

import pandas as pd
import pysam

from utils import safe_float, safe_int
from vcf_utils import build_query_annotations, determine_variant_type, safe_format

log = logging.getLogger(__name__)


def parse_sompy_vcf(sompy_vcf_path, query_vcf_path, tumor_sample=None, normal_sample=None):
    """Parse som.py output VCF and extract classified somatic variant records.

    Args:
        sompy_vcf_path: Path to som.py output VCF
        query_vcf_path: Path to original query VCF for annotation fallback
        tumor_sample: Tumor sample name in the query VCF
        normal_sample: Normal sample name in the query VCF (optional)

    Returns:
        List of variant record dictionaries
    """
    query_annots = build_query_annotations(query_vcf_path)
    records = []

    with pysam.VariantFile(str(sompy_vcf_path)) as vcf:
        sample_names = list(vcf.header.samples)

        for rec in vcf:
            # som.py uses TRUTH and QUERY samples
            truth_data = rec.samples.get("TRUTH") if "TRUTH" in sample_names else None
            query_data = rec.samples.get("QUERY") if "QUERY" in sample_names else None

            if truth_data is None and query_data is None:
                continue

            # Determine classification from som.py BD tag
            truth_bd = safe_format(truth_data, "BD") if truth_data else None
            query_bd = safe_format(query_data, "BD") if query_data else None

            if truth_bd == "FN":
                classification = "FN"
            elif query_bd in ("TP", "FP"):
                classification = query_bd
            else:
                continue  # Not assessed

            for alt in rec.alts or []:
                if alt == "*":
                    continue  # Skip spanning deletions

                key = (rec.chrom, rec.pos, rec.ref, alt)
                variant_type = determine_variant_type(rec.ref, alt)

                # Extract query annotations
                dp = safe_format(query_data, "DP") if query_data else None
                af = safe_format(query_data, "AF") if query_data else None
                gq = safe_format(query_data, "GQ") if query_data else None
                mq = rec.info.get("MQ", None)
                qual = rec.qual
                filt = ";".join(rec.filter.keys()) if rec.filter else None

                # Fallback to original query VCF
                orig = query_annots.get(key, {})
                if dp is None:
                    dp = orig.get("dp")
                if af is None:
                    af = orig.get("af")
                if gq is None:
                    gq = orig.get("gq")
                if mq is None:
                    mq = orig.get("mq")
                if qual is None:
                    qual = orig.get("qual")
                if filt is None:
                    filt = orig.get("filter")

                # Extract somatic-specific fields from original query VCF
                tumor_dp_val = orig.get("dp")
                tumor_af_val = orig.get("af")

                # Compute indel size for stratification
                indel_size = abs(len(alt) - len(rec.ref)) if variant_type == "INDEL" else None

                records.append(
                    {
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": rec.ref,
                        "alt": alt,
                        "type": variant_type,
                        "classification": classification,
                        "dp": safe_int(dp),
                        "mq": safe_float(mq),
                        "af": safe_float(af),
                        "gq": safe_int(gq),
                        "qual": safe_float(qual) if qual is not None and not math.isnan(qual) else None,
                        "filter_status": str(filt) if filt else None,
                        "tumor_dp": safe_int(tumor_dp_val or dp),
                        "tumor_af": safe_float(tumor_af_val or af),
                        "normal_dp": None,
                        "normal_af": None,
                        "somatic_quality": safe_float(qual),
                        "caller_filter": str(filt) if filt else None,
                        "indel_size": indel_size,
                    }
                )

    return records


# --- Main execution via Snakemake ---
try:
    from typing import TYPE_CHECKING

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

    log.info("Parsing som.py VCF: %s", snakemake.input.vcf)
    log.info("Query VCF for annotation fallback: %s", snakemake.input.query_vcf)

    tumor_sample = snakemake.params.get("tumor_sample", None)
    normal_sample = snakemake.params.get("normal_sample", None)

    records = parse_sompy_vcf(
        snakemake.input.vcf,
        snakemake.input.query_vcf,
        tumor_sample=tumor_sample,
        normal_sample=normal_sample,
    )

    df = pd.DataFrame(records)
    df.to_csv(snakemake.output.tsv, sep="\t", index=False)

    log.info("Wrote %d somatic variant records to %s", len(df), snakemake.output.tsv)
except NameError:
    pass  # Not running via Snakemake (e.g., imported for testing)
