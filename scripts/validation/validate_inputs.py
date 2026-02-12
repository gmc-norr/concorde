"""Pre-flight validation script for Concorde pipeline.

This script validates all inputs before pipeline execution to catch errors early.
It uses a collect-all-errors strategy to report all issues at once.

Called from Snakemake via the validate_inputs rule.
"""

import logging
import sys
from pathlib import Path
from typing import Any, TYPE_CHECKING

import yaml

# Add validators to path
SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

if TYPE_CHECKING:
    from snakemake.script import Snakemake
    snakemake: Snakemake
else:
    snakemake = snakemake  # type: ignore  # noqa: F821

from validators import (
    ValidationError,
    validate_config,
    validate_vcf_indexed,
    validate_reference_indexed,
    validate_bed_format,
    validate_chromosome_consistency,
)


class ValidationContext:
    """Context for collecting validation errors without failing fast."""

    def __init__(self):
        self.errors: list[ValidationError] = []
        self.warnings: list[str] = []

    def validate(self, func, *args, **kwargs) -> Any:
        """Run validation function, collecting errors instead of raising.

        Returns:
            Result of validation function, or None if error occurred
        """
        try:
            return func(*args, **kwargs)
        except ValidationError as e:
            self.errors.append(e)
            return None

    def warn(self, message: str):
        """Add a warning message."""
        self.warnings.append(message)

    def raise_if_errors(self):
        """Raise combined error if any validations failed."""
        if self.errors:
            error_list = "\n".join(f"  {i+1}. {e}" for i, e in enumerate(self.errors))
            raise ValidationError(
                f"Validation failed with {len(self.errors)} error(s):\n{error_list}"
            )

    def has_errors(self) -> bool:
        """Check if any errors were collected."""
        return len(self.errors) > 0


def load_config(config_path: str | Path) -> dict:
    """Load configuration from YAML file.

    Args:
        config_path: Path to config.yaml

    Returns:
        Configuration dictionary

    Raises:
        ValidationError: If config cannot be loaded
    """
    path = Path(config_path)

    if not path.exists():
        raise ValidationError(f"Configuration file not found: {path}")

    try:
        with open(path, 'r') as f:
            config = yaml.safe_load(f)

        if not isinstance(config, dict):
            raise ValidationError(
                f"Configuration file {path} does not contain a valid YAML dictionary"
            )

        return config

    except yaml.YAMLError as e:
        raise ValidationError(f"Failed to parse configuration file {path}: {e}")
    except Exception as e:
        raise ValidationError(f"Failed to load configuration file {path}: {e}")


def validate_all_inputs(config_path: str | Path, log: logging.Logger) -> ValidationContext:
    """Validate all pipeline inputs using collect-all-errors strategy.

    Args:
        config_path: Path to configuration file
        log: Logger instance

    Returns:
        ValidationContext with collected errors and warnings
    """
    ctx = ValidationContext()

    # Load configuration
    log.info("Loading configuration from: %s", config_path)
    try:
        config = load_config(config_path)
    except ValidationError as e:
        ctx.errors.append(e)
        return ctx  # Can't proceed without config

    log.info("Configuration loaded successfully")

    # Validate configuration structure
    log.info("Validating configuration structure...")
    ctx.validate(validate_config, config)

    # Extract paths from config
    reference = config.get("reference", "")
    truth_vcf = config.get("truth_vcf", "")
    query_vcfs = config.get("query_vcfs", [])
    confident_regions = config.get("confident_regions", "")
    gene_sets = config.get("gene_sets", [])

    # Validate reference FASTA
    if reference:
        log.info("Validating reference FASTA: %s", reference)
        ctx.validate(validate_reference_indexed, reference)

    # Validate truth VCF
    if truth_vcf:
        log.info("Validating truth VCF: %s", truth_vcf)
        ctx.validate(validate_vcf_indexed, truth_vcf)

    # Validate query VCFs
    for query_vcf in query_vcfs:
        log.info("Validating query VCF: %s", query_vcf)
        ctx.validate(validate_vcf_indexed, query_vcf)

    # Validate confident regions BED (if provided)
    if confident_regions:
        log.info("Validating confident regions BED: %s", confident_regions)
        ctx.validate(validate_bed_format, confident_regions)

    # Validate gene set BED files
    for gene_set in gene_sets:
        if "bed" in gene_set:
            bed_path = gene_set["bed"]
            log.info("Validating gene set BED: %s (name=%s)",
                    bed_path, gene_set.get("name", "unknown"))
            ctx.validate(validate_bed_format, bed_path)

    # Validate chromosome consistency across files
    if reference and truth_vcf and query_vcfs and not ctx.has_errors():
        log.info("Validating chromosome consistency across files...")
        for query_vcf in query_vcfs:
            ctx.validate(
                validate_chromosome_consistency,
                query_vcf,
                truth_vcf,
                reference
            )

    # Report warnings
    if ctx.warnings:
        log.warning("Validation completed with %d warning(s):", len(ctx.warnings))
        for i, warning in enumerate(ctx.warnings, 1):
            log.warning("  %d. %s", i, warning)

    # Report summary
    if ctx.has_errors():
        log.error("Validation failed with %d error(s)", len(ctx.errors))
    else:
        log.info("✓ All input validations passed successfully")

    return ctx


def main(config_path: str, output_file: str, log_file: str):
    """Main entry point for validation script.

    Args:
        config_path: Path to configuration file
        output_file: Path to output sentinel file
        log_file: Path to log file
    """
    # Set up logging
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s"
    )
    log = logging.getLogger(__name__)

    # Also log to console
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    log.addHandler(console)

    log.info("=" * 60)
    log.info("Concorde Pipeline Pre-Flight Validation")
    log.info("=" * 60)

    try:
        # Run all validations
        ctx = validate_all_inputs(config_path, log)

        # Raise if errors occurred
        ctx.raise_if_errors()

        # Create output sentinel file
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        Path(output_file).touch()

        log.info("=" * 60)
        log.info("Validation completed successfully")
        log.info("Output: %s", output_file)
        log.info("=" * 60)

        return 0

    except ValidationError as e:
        log.error("=" * 60)
        log.error("VALIDATION FAILED")
        log.error("=" * 60)
        log.error(str(e))
        log.error("=" * 60)
        log.error("Please fix the errors above and re-run the pipeline")
        return 1

    except Exception:
        log.exception("Unexpected error during validation")
        return 1


if __name__ == "__main__":
    # When called from Snakemake, access snakemake object
    if "snakemake" in globals():
        config_path = snakemake.input.config
        output_file = snakemake.output.validation_done
        log_file = snakemake.log[0]
    else:
        # Standalone execution
        if len(sys.argv) < 2:
            print("Usage: python validate_inputs.py <config.yaml>")
            sys.exit(1)

        config_path = sys.argv[1]
        output_file = ".validation_done"
        log_file = "validation.log"

    sys.exit(main(config_path, output_file, log_file))
