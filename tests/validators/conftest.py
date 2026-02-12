"""Shared test fixtures for validator tests."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import pysam
import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)


@pytest.fixture
def temp_vcf_file(tmp_path) -> Callable:
    """Factory fixture for creating temporary VCF files.

    Returns a function that creates a VCF file with specified variants.
    The VCF will be bgzipped and tab indexed.

    Example:
        >>> vcf_path = temp_vcf_file([
        ...     {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
        ...     {"chrom": "chr2", "pos": 200, "ref": "G", "alt": "C"}
        ... ], sample="NA12878")
    """
    def _create_vcf(variants: list[dict], sample: str = "SAMPLE") -> Path:
        vcf_path = tmp_path / "test.vcf"

        # Create VCF header
        header = pysam.VariantHeader()
        header.add_sample(sample)
        header.add_line('##fileformat=VCFv4.2')
        header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

        # Add contigs from variants
        chromosomes = sorted(set(v["chrom"] for v in variants))
        for chrom in chromosomes:
            header.add_line(f'##contig=<ID={chrom}>')

        # Write VCF
        with pysam.VariantFile(str(vcf_path), 'w', header=header) as vcf:
            for variant in variants:
                rec = vcf.new_record()
                rec.chrom = variant["chrom"]
                rec.pos = variant["pos"]
                rec.ref = variant["ref"]
                rec.alts = (variant["alt"],)
                rec.samples[sample]["GT"] = (0, 1)  # Heterozygous
                vcf.write(rec)

        # Bgzip and index
        bgz_path = tmp_path / "test.vcf.gz"
        pysam.tabix_compress(str(vcf_path), str(bgz_path), force=True)
        pysam.tabix_index(str(bgz_path), preset="vcf", force=True)

        vcf_path.unlink()  # Remove uncompressed VCF
        return bgz_path

    return _create_vcf


@pytest.fixture
def temp_fasta_file(tmp_path) -> Callable:
    """Factory fixture for creating temporary FASTA files with index.

    Returns a function that creates a FASTA file with specified sequences.

    Example:
        >>> fasta_path = temp_fasta_file({
        ...     "chr1": "ACGTACGTACGT",
        ...     "chr2": "TGCATGCATGCA"
        ... })
    """
    def _create_fasta(sequences: dict[str, str]) -> Path:
        fasta_path = tmp_path / "test.fa"

        # Write FASTA
        with open(fasta_path, 'w') as f:
            for chrom, seq in sequences.items():
                f.write(f">{chrom}\n")
                # Write sequence in lines of 60 characters
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")

        # Index with pysam
        pysam.faidx(str(fasta_path))

        return fasta_path

    return _create_fasta


@pytest.fixture
def temp_bed_file(tmp_path) -> Callable:
    """Factory fixture for creating temporary BED files.

    Returns a function that creates a BED file with specified regions.

    Example:
        >>> bed_path = temp_bed_file([
        ...     {"chrom": "chr1", "start": 100, "end": 200, "name": "GENE1"},
        ...     {"chrom": "chr2", "start": 300, "end": 400, "name": "GENE2"}
        ... ])
    """
    def _create_bed(regions: list[dict]) -> Path:
        bed_path = tmp_path / "test.bed"

        with open(bed_path, 'w') as f:
            for region in regions:
                chrom = region["chrom"]
                start = region["start"]
                end = region["end"]
                name = region.get("name", "")

                if name:
                    f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
                else:
                    f.write(f"{chrom}\t{start}\t{end}\n")

        return bed_path

    return _create_bed


@pytest.fixture
def sample_config(tmp_path) -> Callable:
    """Factory fixture for creating test configuration dictionaries.

    Returns a function that creates a config dict with sensible defaults.

    Example:
        >>> config = sample_config(sample="NA12878", caller="GATK")
    """
    def _create_config(**overrides) -> dict:
        # Create default paths
        ref_path = tmp_path / "ref.fa"
        truth_vcf = tmp_path / "truth.vcf.gz"
        query_vcf = tmp_path / "query.vcf.gz"

        config = {
            "sample": "NA12878",
            "caller": "GATK",
            "pipeline_version": "v1.0",
            "comparison_tool": "happy",
            "reference": str(ref_path),
            "truth_vcf": str(truth_vcf),
            "query_vcfs": [str(query_vcf)],
            "decomposition_modes": ["decomposed"],
            "database": "data/concorde.db",
            "output_dir": "results",
            "confident_regions": "",
            "gene_sets": [],
        }

        # Apply overrides
        config.update(overrides)
        return config

    return _create_config
