"""Parse VEP-annotated VCF and extract functional annotations.

Reads the CSQ (or ANN) INFO field from a VEP-annotated VCF and produces
a tab-separated annotations file keyed by (chrom, pos, ref, alt).

Can be called as a Snakemake script or used as a library.
"""

import logging

import pandas as pd
import pysam

from utils import safe_int

log = logging.getLogger(__name__)

# VEP CSQ fields we want to extract (order matches default VEP output)
VEP_FIELDS_OF_INTEREST = {
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature",
    "HGVSc",
    "HGVSp",
    "EXON",
    "Protein_position",
    "SIFT",
    "PolyPhen",
}


def _parse_csq_header(vcf_header):
    """Extract CSQ field names from VCF header.

    Args:
        vcf_header: pysam VariantHeader object

    Returns:
        List of CSQ field names, or None if CSQ not found
    """
    for rec in vcf_header.records:
        if rec.type == "INFO" and rec.get("ID") == "CSQ":
            desc = rec.get("Description", "")
            # Format: "Consequence annotations from Ensembl VEP. Format: A|B|C"
            if "Format: " in desc:
                format_str = desc.split("Format: ")[1].strip('"')
                return format_str.split("|")
    return None


def _parse_ann_header(vcf_header):
    """Extract ANN (SnpEff-style) field names from VCF header.

    Args:
        vcf_header: pysam VariantHeader object

    Returns:
        List of ANN field names, or None if ANN not found
    """
    for rec in vcf_header.records:
        if rec.type == "INFO" and rec.get("ID") == "ANN":
            desc = rec.get("Description", "")
            if "'" in desc:
                # ANN format: "Functional annotations: 'Allele | Annotation | ...'"
                format_str = desc.split("'")[1]
                return [f.strip() for f in format_str.split("|")]
    return None


def _select_transcript(annotations, selection_mode="canonical"):
    """Select the best transcript annotation from a list.

    Args:
        annotations: List of annotation dicts
        selection_mode: "canonical", "mane_select", or "most_severe"

    Returns:
        Best annotation dict, or None if empty
    """
    if not annotations:
        return None

    if selection_mode == "most_severe":
        # Sort by impact severity
        impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
        return min(
            annotations,
            key=lambda a: impact_order.get(a.get("impact", "MODIFIER"), 4),
        )

    if selection_mode == "mane_select":
        # Prefer MANE Select transcripts
        mane = [a for a in annotations if a.get("mane_select")]
        if mane:
            return mane[0]

    if selection_mode == "canonical":
        # Prefer canonical transcripts
        canonical = [a for a in annotations if a.get("canonical")]
        if canonical:
            return canonical[0]

    # Fallback: most severe impact
    impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    return min(
        annotations,
        key=lambda a: impact_order.get(a.get("impact", "MODIFIER"), 4),
    )


def parse_vep_vcf(vep_vcf_path, transcript_selection="canonical"):
    """Parse VEP-annotated VCF and extract annotations per variant.

    Args:
        vep_vcf_path: Path to VEP-annotated VCF
        transcript_selection: Transcript selection mode

    Returns:
        List of annotation record dictionaries
    """
    records = []

    with pysam.VariantFile(str(vep_vcf_path)) as vcf:
        # Determine annotation format (CSQ or ANN)
        csq_fields = _parse_csq_header(vcf.header)
        ann_fields = _parse_ann_header(vcf.header) if csq_fields is None else None

        if csq_fields is None and ann_fields is None:
            log.warning("No CSQ or ANN INFO field found in VCF header")
            return records

        field_names = csq_fields or ann_fields
        info_key = "CSQ" if csq_fields else "ANN"
        log.info("Using %s annotation format with %d fields", info_key, len(field_names))

        for rec in vcf:
            csq_values = rec.info.get(info_key)
            if csq_values is None:
                continue

            for alt in rec.alts or []:
                if alt == "*":
                    continue

                # Parse all transcript annotations for this alt allele
                alt_annotations = []
                for csq_str in csq_values:
                    fields = csq_str.split("|")
                    field_dict = dict(zip(field_names, fields))

                    # Check if this annotation is for this alt allele
                    allele = field_dict.get("Allele", "")
                    if allele and allele != alt and allele != "-":
                        continue

                    annotation = {
                        "consequence": field_dict.get("Consequence", ""),
                        "impact": field_dict.get("IMPACT", ""),
                        "gene_symbol": field_dict.get("SYMBOL", ""),
                        "gene_id": field_dict.get("Gene", ""),
                        "transcript_id": field_dict.get("Feature", ""),
                        "hgvsc": field_dict.get("HGVSc", "") or None,
                        "hgvsp": field_dict.get("HGVSp", "") or None,
                        "exon": field_dict.get("EXON", "") or None,
                        "protein_position": safe_int(
                            field_dict.get("Protein_position", "")
                        ),
                        "sift_prediction": field_dict.get("SIFT", "") or None,
                        "polyphen_prediction": field_dict.get("PolyPhen", "") or None,
                        "canonical": field_dict.get("CANONICAL", "") == "YES",
                        "mane_select": bool(field_dict.get("MANE_SELECT", "")),
                    }
                    alt_annotations.append(annotation)

                # Select best transcript
                best = _select_transcript(alt_annotations, transcript_selection)
                if best is None:
                    continue

                records.append(
                    {
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": rec.ref,
                        "alt": alt,
                        "consequence": best["consequence"],
                        "impact": best["impact"],
                        "gene_symbol": best["gene_symbol"],
                        "gene_id": best["gene_id"],
                        "transcript_id": best["transcript_id"],
                        "hgvsc": best["hgvsc"],
                        "hgvsp": best["hgvsp"],
                        "exon": best["exon"],
                        "protein_position": best["protein_position"],
                        "sift_prediction": best["sift_prediction"],
                        "polyphen_prediction": best["polyphen_prediction"],
                    }
                )

    log.info("Parsed %d VEP annotations", len(records))
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

    transcript_selection = snakemake.params.get("transcript_selection", "canonical")

    log.info("Parsing VEP annotations: %s", snakemake.input.vcf)
    records = parse_vep_vcf(snakemake.input.vcf, transcript_selection)

    df = pd.DataFrame(records)
    df.to_csv(snakemake.output.tsv, sep="\t", index=False)

    log.info("Wrote %d VEP annotation records to %s", len(df), snakemake.output.tsv)
except NameError:
    pass  # Not running via Snakemake (e.g., imported for testing)
