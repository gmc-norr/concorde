"""Tests for VEP annotation parser."""

from __future__ import annotations

import sys
from pathlib import Path

import pysam
import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from parsers.parse_vep_annotations import (  # noqa: E402
    _parse_csq_header,
    _select_transcript,
    parse_vep_vcf,
)


def _create_vep_vcf(tmp_path, variants, csq_format=None):
    """Create a VEP-annotated VCF for testing.

    Args:
        tmp_path: pytest tmp_path
        variants: List of dicts with chrom, pos, ref, alt, csq (CSQ string)
        csq_format: CSQ format fields (defaults to standard VEP fields)

    Returns:
        Path to bgzipped, indexed VCF
    """
    if csq_format is None:
        csq_format = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|EXON|Protein_position|SIFT|PolyPhen|CANONICAL|MANE_SELECT"

    vcf_path = tmp_path / "vep.vcf"

    header = pysam.VariantHeader()
    header.add_sample("SAMPLE")
    header.add_line("##fileformat=VCFv4.2")
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line(
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {csq_format}">'
    )

    chromosomes = sorted(set(v["chrom"] for v in variants))
    for chrom in chromosomes:
        header.add_line(f"##contig=<ID={chrom}>")

    with pysam.VariantFile(str(vcf_path), "w", header=header) as vcf:
        for variant in variants:
            rec = vcf.new_record()
            rec.chrom = variant["chrom"]
            rec.pos = variant["pos"]
            rec.ref = variant["ref"]
            rec.alts = (variant["alt"],)
            rec.samples["SAMPLE"]["GT"] = (0, 1)
            if "csq" in variant:
                rec.info["CSQ"] = variant["csq"] if isinstance(variant["csq"], tuple) else (variant["csq"],)
            vcf.write(rec)

    bgz_path = tmp_path / "vep.vcf.gz"
    pysam.tabix_compress(str(vcf_path), str(bgz_path), force=True)
    pysam.tabix_index(str(bgz_path), preset="vcf", force=True)
    vcf_path.unlink()

    return bgz_path


class TestParseCSQHeader:
    """Tests for CSQ header parsing."""

    def test_parse_csq_format(self, tmp_path):
        fmt = "Allele|Consequence|IMPACT|SYMBOL"
        vcf_path = _create_vep_vcf(tmp_path, [], csq_format=fmt)

        with pysam.VariantFile(str(vcf_path)) as vcf:
            fields = _parse_csq_header(vcf.header)

        assert fields == ["Allele", "Consequence", "IMPACT", "SYMBOL"]

    def test_no_csq_returns_none(self, tmp_path):
        vcf_path = tmp_path / "no_csq.vcf"
        header = pysam.VariantHeader()
        header.add_sample("SAMPLE")
        header.add_line("##fileformat=VCFv4.2")
        header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        header.add_line("##contig=<ID=chr1>")

        with pysam.VariantFile(str(vcf_path), "w", header=header) as vcf:
            pass

        bgz = tmp_path / "no_csq.vcf.gz"
        pysam.tabix_compress(str(vcf_path), str(bgz), force=True)
        pysam.tabix_index(str(bgz), preset="vcf", force=True)

        with pysam.VariantFile(str(bgz)) as vcf:
            fields = _parse_csq_header(vcf.header)

        assert fields is None


class TestSelectTranscript:
    """Tests for transcript selection logic."""

    def test_most_severe(self):
        annotations = [
            {"impact": "MODIFIER", "consequence": "intron_variant"},
            {"impact": "HIGH", "consequence": "stop_gained"},
            {"impact": "LOW", "consequence": "synonymous_variant"},
        ]
        best = _select_transcript(annotations, "most_severe")
        assert best["impact"] == "HIGH"

    def test_canonical_preferred(self):
        annotations = [
            {"impact": "MODERATE", "canonical": False},
            {"impact": "LOW", "canonical": True},
        ]
        best = _select_transcript(annotations, "canonical")
        assert best["canonical"] is True

    def test_mane_select_preferred(self):
        annotations = [
            {"impact": "MODERATE", "mane_select": False},
            {"impact": "LOW", "mane_select": True},
        ]
        best = _select_transcript(annotations, "mane_select")
        assert best["mane_select"] is True

    def test_fallback_to_most_severe_when_no_canonical(self):
        annotations = [
            {"impact": "MODIFIER", "canonical": False},
            {"impact": "HIGH", "canonical": False},
        ]
        best = _select_transcript(annotations, "canonical")
        assert best["impact"] == "HIGH"

    def test_empty_returns_none(self):
        assert _select_transcript([], "canonical") is None


class TestParseVepVcf:
    """Tests for VEP VCF parsing."""

    def test_parse_basic_annotation(self, tmp_path):
        csq = "T|missense_variant|MODERATE|TP53|ENSG00000141510|ENST00000269305|c.743G>A|p.Arg248Gln|7/11|248|deleterious(0.0)|probably_damaging(1.0)|YES|NM_000546"
        vcf_path = _create_vep_vcf(
            tmp_path,
            [{"chrom": "chr17", "pos": 7577538, "ref": "C", "alt": "T", "csq": csq}],
        )

        records = parse_vep_vcf(str(vcf_path))

        assert len(records) == 1
        r = records[0]
        assert r["chrom"] == "chr17"
        assert r["pos"] == 7577538
        assert r["consequence"] == "missense_variant"
        assert r["impact"] == "MODERATE"
        assert r["gene_symbol"] == "TP53"
        assert r["gene_id"] == "ENSG00000141510"
        assert r["transcript_id"] == "ENST00000269305"
        assert r["hgvsc"] == "c.743G>A"
        assert r["hgvsp"] == "p.Arg248Gln"
        assert r["exon"] == "7/11"
        assert r["protein_position"] == 248
        assert r["sift_prediction"] == "deleterious(0.0)"
        assert r["polyphen_prediction"] == "probably_damaging(1.0)"

    def test_parse_multiple_transcripts_selects_canonical(self, tmp_path):
        csq = (
            "T|intron_variant|MODIFIER|TP53|ENSG00000141510|ENST00000999999||||||||NO|",
            "T|missense_variant|MODERATE|TP53|ENSG00000141510|ENST00000269305|c.743G>A|p.Arg248Gln|7/11|248|deleterious(0.0)|probably_damaging(1.0)|YES|NM_000546",
        )
        vcf_path = _create_vep_vcf(
            tmp_path,
            [{"chrom": "chr17", "pos": 7577538, "ref": "C", "alt": "T", "csq": csq}],
        )

        records = parse_vep_vcf(str(vcf_path), transcript_selection="canonical")

        assert len(records) == 1
        assert records[0]["consequence"] == "missense_variant"
        assert records[0]["transcript_id"] == "ENST00000269305"

    def test_parse_empty_vcf(self, tmp_path):
        vcf_path = _create_vep_vcf(tmp_path, [])
        records = parse_vep_vcf(str(vcf_path))
        assert records == []

    def test_parse_variant_without_csq(self, tmp_path):
        vcf_path = _create_vep_vcf(
            tmp_path,
            [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}],
        )
        records = parse_vep_vcf(str(vcf_path))
        assert records == []

    def test_none_fields_handled(self, tmp_path):
        csq = "T|synonymous_variant|LOW|BRCA1|ENSG00000012048|ENST00000357654||||||benign(0.9)|YES|"
        vcf_path = _create_vep_vcf(
            tmp_path,
            [{"chrom": "chr17", "pos": 43094464, "ref": "G", "alt": "T", "csq": csq}],
        )

        records = parse_vep_vcf(str(vcf_path))

        assert len(records) == 1
        assert records[0]["hgvsc"] is None
        assert records[0]["hgvsp"] is None
        assert records[0]["exon"] is None
        assert records[0]["protein_position"] is None
