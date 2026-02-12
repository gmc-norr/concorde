"""Tests for input file validators (VCF, FASTA, BED)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from validators import (
    ValidationError,
    validate_vcf_indexed,
    validate_vcf_format,
    validate_vcf_chromosomes,
    validate_reference_indexed,
    validate_bed_format,
)


class TestValidateVCFIndexed:
    """Tests for validate_vcf_indexed function."""

    def test_valid_indexed_vcf_passes(self, temp_vcf_file):
        """Test validation passes for properly indexed VCF."""
        vcf_path = temp_vcf_file([
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}
        ])
        validate_vcf_indexed(vcf_path)  # Should not raise

    def test_missing_vcf_raises(self):
        """Test validation fails for missing VCF."""
        with pytest.raises(ValidationError, match="not found"):
            validate_vcf_indexed("/nonexistent/file.vcf.gz")

    def test_uncompressed_vcf_raises(self, tmp_path):
        """Test validation fails for uncompressed VCF."""
        vcf_path = tmp_path / "test.vcf"
        vcf_path.write_text("##fileformat=VCFv4.2\n")

        with pytest.raises(ValidationError, match="not bgzipped"):
            validate_vcf_indexed(vcf_path)

    def test_missing_index_raises(self, tmp_path):
        """Test validation fails when index is missing."""
        vcf_path = tmp_path / "test.vcf.gz"
        vcf_path.write_bytes(b"dummy content")  # Fake bgzipped file

        with pytest.raises(ValidationError, match="missing tabix index"):
            validate_vcf_indexed(vcf_path)


class TestValidateVCFFormat:
    """Tests for validate_vcf_format function."""

    def test_valid_vcf_returns_file_object(self, temp_vcf_file):
        """Test validation returns opened VCF file."""
        vcf_path = temp_vcf_file([
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}
        ])
        vcf = validate_vcf_format(vcf_path)
        assert vcf is not None
        assert hasattr(vcf, 'header')
        vcf.close()

    def test_invalid_vcf_raises(self, tmp_path):
        """Test validation fails for malformed VCF."""
        vcf_path = tmp_path / "invalid.vcf.gz"
        vcf_path.write_bytes(b"not a vcf file")

        with pytest.raises(ValidationError, match="cannot be opened"):
            validate_vcf_format(vcf_path)


class TestValidateVCFChromosomes:
    """Tests for validate_vcf_chromosomes function."""

    def test_returns_chromosome_set(self, temp_vcf_file):
        """Test validation returns set of chromosomes."""
        vcf_path = temp_vcf_file([
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
            {"chrom": "chr2", "pos": 200, "ref": "G", "alt": "C"},
        ])
        chroms = validate_vcf_chromosomes(vcf_path)
        assert chroms == {"chr1", "chr2"}

    def test_validates_against_expected_chromosomes(self, temp_vcf_file):
        """Test validation checks against expected chromosomes."""
        vcf_path = temp_vcf_file([
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}
        ])
        # Should pass - matches expected
        chroms = validate_vcf_chromosomes(vcf_path, expected_chroms={"chr1"})
        assert chroms == {"chr1"}

    def test_missing_expected_chromosomes_raises(self, temp_vcf_file):
        """Test validation fails when expected chromosomes missing."""
        vcf_path = temp_vcf_file([
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"}
        ])
        with pytest.raises(ValidationError, match="chromosome mismatch"):
            validate_vcf_chromosomes(vcf_path, expected_chroms={"chr1", "chr2"})

    def test_extra_chromosomes_raises(self, temp_vcf_file):
        """Test validation fails when extra chromosomes present."""
        vcf_path = temp_vcf_file([
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T"},
            {"chrom": "chr2", "pos": 200, "ref": "G", "alt": "C"},
        ])
        with pytest.raises(ValidationError, match="chromosome mismatch"):
            validate_vcf_chromosomes(vcf_path, expected_chroms={"chr1"})


class TestValidateReferenceIndexed:
    """Tests for validate_reference_indexed function."""

    def test_valid_indexed_fasta_passes(self, temp_fasta_file):
        """Test validation passes for indexed FASTA."""
        fasta_path = temp_fasta_file({"chr1": "ACGTACGT"})
        validate_reference_indexed(fasta_path)  # Should not raise

    def test_missing_fasta_raises(self):
        """Test validation fails for missing FASTA."""
        with pytest.raises(ValidationError, match="not found"):
            validate_reference_indexed("/nonexistent/ref.fa")

    def test_missing_index_raises(self, tmp_path):
        """Test validation fails when .fai index missing."""
        fasta_path = tmp_path / "test.fa"
        fasta_path.write_text(">chr1\nACGT\n")

        with pytest.raises(ValidationError, match="missing index"):
            validate_reference_indexed(fasta_path)

    def test_empty_index_raises(self, tmp_path):
        """Test validation fails when .fai index is empty."""
        fasta_path = tmp_path / "test.fa"
        fasta_path.write_text(">chr1\nACGT\n")

        fai_path = tmp_path / "test.fa.fai"
        fai_path.write_text("")  # Empty index

        with pytest.raises(ValidationError, match="empty"):
            validate_reference_indexed(fasta_path)


class TestValidateVCFEdgeCases:
    """Edge case tests for VCF validation."""

    def test_vcf_with_no_variants(self, tmp_path):
        """Test VCF file with header but no variants."""
        import pysam
        vcf_path = tmp_path / "empty.vcf"

        # Create VCF with header but no variants
        header = pysam.VariantHeader()
        header.add_sample("SAMPLE")
        header.add_line('##fileformat=VCFv4.2')
        header.add_line('##contig=<ID=chr1>')

        with pysam.VariantFile(str(vcf_path), 'w', header=header) as vcf:
            pass  # No variants written

        # Bgzip and index
        bgz_path = tmp_path / "empty.vcf.gz"
        pysam.tabix_compress(str(vcf_path), str(bgz_path), force=True)
        pysam.tabix_index(str(bgz_path), preset="vcf", force=True)

        # Should handle gracefully
        # Note: VCF header declares chr1 in ##contig line, so it's returned
        # even though there are no variant records
        chroms = validate_vcf_chromosomes(bgz_path)
        assert chroms == {'chr1'}  # Chromosomes from header, not variants

    def test_vcf_case_sensitivity(self, temp_vcf_file):
        """Test chromosome name case sensitivity."""
        vcf_path = temp_vcf_file([
            {"chrom": "Chr1", "pos": 100, "ref": "A", "alt": "T"},  # Mixed case
            {"chrom": "CHR2", "pos": 200, "ref": "G", "alt": "C"},  # Upper case
        ])

        chroms = validate_vcf_chromosomes(vcf_path)
        # Chromosomes are returned as-is from VCF
        assert "Chr1" in chroms
        assert "CHR2" in chroms


class TestValidateBEDFormat:
    """Tests for validate_bed_format function."""

    def test_valid_bed_passes(self, temp_bed_file):
        """Test validation passes for valid BED file."""
        bed_path = temp_bed_file([
            {"chrom": "chr1", "start": 100, "end": 200, "name": "GENE1"}
        ])
        validate_bed_format(bed_path)  # Should not raise

    def test_minimal_bed_passes(self, temp_bed_file):
        """Test validation passes for minimal 3-column BED."""
        bed_path = temp_bed_file([
            {"chrom": "chr1", "start": 100, "end": 200}
        ])
        validate_bed_format(bed_path)  # Should not raise

    def test_missing_bed_raises(self):
        """Test validation fails for missing BED file."""
        with pytest.raises(ValidationError, match="not found"):
            validate_bed_format("/nonexistent/file.bed")

    def test_too_few_columns_raises(self, tmp_path):
        """Test validation fails for BED with < 3 columns."""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text("chr1\t100\n")  # Only 2 columns

        with pytest.raises(ValidationError, match="2 columns"):
            validate_bed_format(bed_path)

    def test_non_numeric_coordinates_raises(self, tmp_path):
        """Test validation fails for non-numeric coordinates."""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text("chr1\tABC\t200\n")

        with pytest.raises(ValidationError, match="non-numeric"):
            validate_bed_format(bed_path)

    def test_start_greater_than_end_raises(self, tmp_path):
        """Test validation fails when start >= end."""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text("chr1\t200\t100\n")

        with pytest.raises(ValidationError, match="start >= end"):
            validate_bed_format(bed_path)

    def test_negative_start_raises(self, tmp_path):
        """Test validation fails for negative start coordinate."""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text("chr1\t-100\t200\n")

        with pytest.raises(ValidationError, match="negative start"):
            validate_bed_format(bed_path)

    def test_skips_comments_and_headers(self, tmp_path):
        """Test validation skips comment lines and track headers."""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text(
            "# This is a comment\n"
            "track name=test\n"
            "browser position chr1:100-200\n"
            "chr1\t100\t200\tGENE1\n"
        )
        validate_bed_format(bed_path)  # Should not raise


class TestValidateBEDEdgeCases:
    """Edge case tests for BED validation."""

    def test_zero_length_interval(self, tmp_path):
        """Test BED with zero-length interval (start == end).

        Note: While BED format allows start == end for insertion points,
        this validator enforces start < end for simplicity.
        """
        bed_path = tmp_path / "test.bed"
        bed_path.write_text("chr1\t100\t100\n")

        # Current validator rejects zero-length intervals
        with pytest.raises(ValidationError, match="start >= end"):
            validate_bed_format(bed_path)

    def test_very_large_coordinates(self, tmp_path):
        """Test BED with very large coordinate values."""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text("chr1\t1000000000\t1000000100\n")

        validate_bed_format(bed_path)  # Should handle large numbers

    def test_scientific_notation_coordinates(self, tmp_path):
        """Test BED with scientific notation (should fail)."""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text("chr1\t1e6\t2e6\n")

        with pytest.raises(ValidationError, match="non-numeric"):
            validate_bed_format(bed_path)

    def test_overlapping_regions(self, tmp_path):
        """Test BED with overlapping regions (valid in BED format)."""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text(
            "chr1\t100\t200\n"
            "chr1\t150\t250\n"  # Overlaps
        )

        # Overlapping regions are valid in BED files
        validate_bed_format(bed_path)
