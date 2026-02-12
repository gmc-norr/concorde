"""Validation for input files (VCF, FASTA, BED).

This module provides validators for genomic input file formats to ensure
they exist, are properly indexed, and can be opened correctly.
"""

from pathlib import Path

import pysam

from .base import ValidationError, validate_file_exists


def validate_vcf_indexed(vcf_path: str | Path) -> None:
    """Validate VCF is bgzipped and has tabix index.

    Args:
        vcf_path: Path to VCF file

    Raises:
        ValidationError: If VCF is not bgzipped or index missing

    Example:
        >>> validate_vcf_indexed("sample.vcf.gz")  # Pass
        >>> validate_vcf_indexed("sample.vcf")     # Fail - not bgzipped
    """
    path = Path(vcf_path)
    validate_file_exists(path, "VCF file")

    # Check for .gz extension
    if not str(path).endswith(".vcf.gz"):
        raise ValidationError(
            f"VCF file {path.name} is not bgzipped. "
            f"Expected .vcf.gz extension. Use bgzip to compress."
        )

    # Check for index (.tbi or .csi)
    tbi_path = Path(str(path) + ".tbi")
    csi_path = Path(str(path) + ".csi")

    if not tbi_path.exists() and not csi_path.exists():
        raise ValidationError(
            f"VCF file {path.name} is missing tabix index. "
            f"Expected {path.name}.tbi or {path.name}.csi. "
            f"Run: tabix -p vcf {path.name}"
        )


def validate_vcf_format(vcf_path: str | Path) -> pysam.VariantFile:
    """Validate VCF can be opened and has valid header.

    Args:
        vcf_path: Path to VCF file

    Returns:
        Opened pysam.VariantFile object (caller should close it)

    Raises:
        ValidationError: If VCF cannot be opened or header is invalid

    Example:
        >>> vcf = validate_vcf_format("sample.vcf.gz")
        >>> vcf.close()
    """
    path = Path(vcf_path)
    validate_file_exists(path, "VCF file")

    try:
        vcf = pysam.VariantFile(str(path))

        # Try to access header to ensure it's valid
        _ = vcf.header

        return vcf

    except Exception as e:
        raise ValidationError(
            f"VCF file {path.name} cannot be opened or has invalid header: {e}"
        )


def validate_vcf_chromosomes(
    vcf_path: str | Path,
    expected_chroms: set[str] | None = None,
    quick_check: bool = True
) -> set[str]:
    """Get chromosomes from VCF, optionally validate against expected set.

    Args:
        vcf_path: Path to VCF file
        expected_chroms: Optional set of expected chromosome names
        quick_check: If True, only check first 1000 variants (faster)

    Returns:
        Set of chromosome names found in VCF

    Raises:
        ValidationError: If chromosomes don't match expected set

    Example:
        >>> chroms = validate_vcf_chromosomes("sample.vcf.gz")
        >>> chroms = validate_vcf_chromosomes(
        ...     "sample.vcf.gz",
        ...     expected_chroms={"chr1", "chr2", "chrX"}
        ... )
    """
    path = Path(vcf_path)
    vcf = validate_vcf_format(path)

    try:
        # Get chromosomes from VCF header first (fast)
        header_chroms = set(vcf.header.contigs.keys()) if vcf.header.contigs else set()

        # Optionally sample actual records to verify
        record_chroms = set()
        max_records = 1000 if quick_check else None
        count = 0

        for rec in vcf:
            record_chroms.add(rec.chrom)
            count += 1
            if max_records and count >= max_records:
                break

        # Use header chromosomes if available, otherwise use sampled records
        found_chroms = header_chroms if header_chroms else record_chroms

        if not found_chroms:
            raise ValidationError(
                f"VCF file {path.name} contains no chromosomes in header or records"
            )

        # Validate against expected chromosomes if provided
        if expected_chroms:
            missing = expected_chroms - found_chroms
            extra = found_chroms - expected_chroms

            if missing or extra:
                msg_parts = [f"VCF file {path.name} chromosome mismatch:"]
                if missing:
                    msg_parts.append(f"  Missing: {sorted(missing)}")
                if extra:
                    msg_parts.append(f"  Extra: {sorted(extra)}")
                raise ValidationError("\n".join(msg_parts))

        return found_chroms

    finally:
        vcf.close()


def validate_reference_indexed(ref_path: str | Path) -> None:
    """Validate FASTA reference has .fai index.

    Args:
        ref_path: Path to FASTA file

    Raises:
        ValidationError: If FASTA is not indexed

    Example:
        >>> validate_reference_indexed("GRCh38.fa")
    """
    path = Path(ref_path)
    validate_file_exists(path, "Reference FASTA")

    # Check for .fai index
    fai_path = Path(str(path) + ".fai")

    if not fai_path.exists():
        raise ValidationError(
            f"Reference FASTA {path.name} is missing index. "
            f"Expected {path.name}.fai. "
            f"Run: samtools faidx {path.name}"
        )

    # Verify .fai file is not empty
    if fai_path.stat().st_size == 0:
        raise ValidationError(
            f"Reference FASTA index {fai_path.name} is empty. "
            f"Re-run: samtools faidx {path.name}"
        )


def validate_bed_format(bed_path: str | Path) -> None:
    """Validate BED file has proper format (at least 3 columns: chr, start, end).

    Args:
        bed_path: Path to BED file

    Raises:
        ValidationError: If BED format is invalid

    Example:
        >>> validate_bed_format("regions.bed")
    """
    path = Path(bed_path)
    validate_file_exists(path, "BED file")

    try:
        # Read first few lines to validate format
        with open(path, "r") as f:
            for i, line in enumerate(f):
                # Skip empty lines and comments
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                    continue

                # Check minimum 3 columns
                cols = line.split("\t")
                if len(cols) < 3:
                    raise ValidationError(
                        f"BED file {path.name} line {i+1} has {len(cols)} columns, "
                        f"expected at least 3 (chrom, start, end). Line: {line[:50]}"
                    )

                # Validate start and end are integers
                try:
                    start = int(cols[1])
                    end = int(cols[2])
                except ValueError:
                    raise ValidationError(
                        f"BED file {path.name} line {i+1} has non-numeric start/end coordinates. "
                        f"Start: '{cols[1]}', End: '{cols[2]}'"
                    )

                # Validate start < end
                if start >= end:
                    raise ValidationError(
                        f"BED file {path.name} line {i+1} has start >= end "
                        f"({start} >= {end}). BED coordinates must have start < end."
                    )

                # Validate start >= 0
                if start < 0:
                    raise ValidationError(
                        f"BED file {path.name} line {i+1} has negative start coordinate ({start}). "
                        f"BED coordinates must be >= 0."
                    )

                # Only check first data line (rest are assumed similar)
                break

    except ValidationError:
        raise
    except Exception as e:
        raise ValidationError(
            f"BED file {path.name} could not be read: {e}"
        )
