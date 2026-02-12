"""Configuration validation for pipeline config.yaml.

This module provides validators to ensure config.yaml has required fields,
valid values, and proper structure before pipeline execution.
"""


from .base import ValidationError

VALID_MODES = {"germline", "somatic"}
VALID_GERMLINE_TOOLS = {"happy", "rtg"}
VALID_SOMATIC_TOOLS = {"sompy", "rtg", "internal"}
VALID_CALLING_MODES = {"tumor_only", "tumor_normal"}
VALID_GENOME_BUILDS = {"GRCh37", "GRCh38"}
VALID_TRANSCRIPT_SOURCES = {"ensembl", "refseq", "merged"}
VALID_TRANSCRIPT_SELECTIONS = {"canonical", "mane_select", "most_severe"}


def validate_mode(mode: str) -> None:
    """Validate execution mode is 'germline' or 'somatic'.

    Args:
        mode: Execution mode string

    Raises:
        ValidationError: If mode is not valid
    """
    if mode not in VALID_MODES:
        raise ValidationError(
            f"Invalid mode '{mode}'. "
            f"Valid options are: {sorted(VALID_MODES)}"
        )


def validate_comparison_tool(tool: str, mode: str = "germline") -> None:
    """Validate comparison tool for the given mode.

    Args:
        tool: Comparison tool name
        mode: Execution mode (germline or somatic)

    Raises:
        ValidationError: If tool is not valid for the mode
    """
    if mode == "somatic":
        valid_tools = VALID_SOMATIC_TOOLS
    else:
        valid_tools = VALID_GERMLINE_TOOLS
    if tool not in valid_tools:
        raise ValidationError(
            f"Invalid comparison_tool '{tool}' for mode '{mode}'. "
            f"Valid options are: {sorted(valid_tools)}"
        )


def validate_somatic_config(somatic: dict) -> None:
    """Validate somatic configuration block.

    Args:
        somatic: Somatic configuration dictionary

    Raises:
        ValidationError: If somatic config is invalid
    """
    if not isinstance(somatic, dict):
        raise ValidationError(
            f"Config field 'somatic' must be a dictionary, "
            f"got: {type(somatic).__name__}"
        )

    # tumor_sample is required
    if "tumor_sample" not in somatic or not somatic["tumor_sample"]:
        raise ValidationError(
            "Somatic config requires 'tumor_sample' to be a non-empty string."
        )

    # calling_mode is required
    calling_mode = somatic.get("calling_mode")
    if not calling_mode:
        raise ValidationError(
            "Somatic config requires 'calling_mode'. "
            f"Valid options are: {sorted(VALID_CALLING_MODES)}"
        )
    if calling_mode not in VALID_CALLING_MODES:
        raise ValidationError(
            f"Invalid somatic calling_mode '{calling_mode}'. "
            f"Valid options are: {sorted(VALID_CALLING_MODES)}"
        )

    # normal_sample required for tumor_normal mode
    if calling_mode == "tumor_normal":
        if "normal_sample" not in somatic or not somatic["normal_sample"]:
            raise ValidationError(
                "Somatic config requires 'normal_sample' when "
                "calling_mode is 'tumor_normal'."
            )

    # Validate numeric fields
    if "vaf_tolerance" in somatic:
        vaf_tol = somatic["vaf_tolerance"]
        if not isinstance(vaf_tol, (int, float)) or vaf_tol < 0 or vaf_tol > 1:
            raise ValidationError(
                f"Somatic 'vaf_tolerance' must be a number between 0 and 1, "
                f"got: {vaf_tol}"
            )

    if "min_vaf" in somatic:
        min_vaf = somatic["min_vaf"]
        if not isinstance(min_vaf, (int, float)) or min_vaf < 0 or min_vaf > 1:
            raise ValidationError(
                f"Somatic 'min_vaf' must be a number between 0 and 1, "
                f"got: {min_vaf}"
            )


def validate_truth_matching_config(truth_matching: dict, mode: str) -> None:
    """Validate truth matching configuration.

    Args:
        truth_matching: Truth matching configuration dictionary
        mode: Execution mode

    Raises:
        ValidationError: If truth matching config is invalid
    """
    if not isinstance(truth_matching, dict):
        raise ValidationError(
            f"Config field 'truth_matching' must be a dictionary, "
            f"got: {type(truth_matching).__name__}"
        )

    mode_config = truth_matching.get(mode)
    if mode_config is None:
        return  # No mode-specific matching config; use defaults

    if not isinstance(mode_config, dict):
        raise ValidationError(
            f"truth_matching.{mode} must be a dictionary, "
            f"got: {type(mode_config).__name__}"
        )

    if "tool" in mode_config:
        tool = mode_config["tool"]
        if mode == "somatic":
            valid = VALID_SOMATIC_TOOLS
        else:
            valid = VALID_GERMLINE_TOOLS
        if tool not in valid:
            raise ValidationError(
                f"Invalid truth_matching.{mode}.tool '{tool}'. "
                f"Valid options are: {sorted(valid)}"
            )


def validate_vep_config(vep: dict) -> None:
    """Validate VEP annotation configuration block.

    Args:
        vep: VEP configuration dictionary

    Raises:
        ValidationError: If VEP config is invalid
    """
    if not isinstance(vep, dict):
        raise ValidationError(
            f"Config field 'vep' must be a dictionary, "
            f"got: {type(vep).__name__}"
        )

    if not vep.get("enabled", False):
        return  # VEP disabled, no further validation needed

    # genome_build is required when enabled
    genome_build = vep.get("genome_build")
    if genome_build and genome_build not in VALID_GENOME_BUILDS:
        raise ValidationError(
            f"Invalid vep.genome_build '{genome_build}'. "
            f"Valid options are: {sorted(VALID_GENOME_BUILDS)}"
        )

    # transcript_source validation
    transcript_source = vep.get("transcript_source")
    if transcript_source and transcript_source not in VALID_TRANSCRIPT_SOURCES:
        raise ValidationError(
            f"Invalid vep.transcript_source '{transcript_source}'. "
            f"Valid options are: {sorted(VALID_TRANSCRIPT_SOURCES)}"
        )

    # transcript_selection validation
    transcript_selection = vep.get("transcript_selection")
    if transcript_selection and transcript_selection not in VALID_TRANSCRIPT_SELECTIONS:
        raise ValidationError(
            f"Invalid vep.transcript_selection '{transcript_selection}'. "
            f"Valid options are: {sorted(VALID_TRANSCRIPT_SELECTIONS)}"
        )

    # plugins must be a list if present
    plugins = vep.get("plugins")
    if plugins is not None and not isinstance(plugins, list):
        raise ValidationError(
            f"Config field 'vep.plugins' must be a list, "
            f"got: {type(plugins).__name__}"
        )


def validate_decomposition_modes(modes: list[str]) -> None:
    """Validate decomposition modes are valid.

    Args:
        modes: List of decomposition mode strings

    Raises:
        ValidationError: If any mode is invalid

    Example:
        >>> validate_decomposition_modes(["decomposed", "non-decomposed"])  # Pass
        >>> validate_decomposition_modes(["invalid"])  # Fail
    """
    valid_modes = {"decomposed", "non-decomposed"}
    invalid = set(modes) - valid_modes

    if invalid:
        raise ValidationError(
            f"Invalid decomposition_modes: {sorted(invalid)}. "
            f"Valid options are: {sorted(valid_modes)}"
        )

    if not modes:
        raise ValidationError(
            "decomposition_modes cannot be empty. "
            f"Valid options are: {sorted(valid_modes)}"
        )


def validate_gene_sets_config(gene_sets: list[dict]) -> None:
    """Validate gene sets configuration structure.

    Args:
        gene_sets: List of gene set dictionaries from config

    Raises:
        ValidationError: If gene set structure is invalid

    Example:
        >>> validate_gene_sets_config([{
        ...     "name": "ACMG59",
        ...     "description": "ACMG genes",
        ...     "version": "v3.1",
        ...     "bed": "path/to/acmg.bed"
        ... }])
    """
    required_fields = {"name", "description", "version", "bed"}

    for i, gene_set in enumerate(gene_sets):
        if not isinstance(gene_set, dict):
            raise ValidationError(
                f"gene_sets[{i}] is not a dictionary. "
                f"Each gene set must have: {sorted(required_fields)}"
            )

        missing = required_fields - set(gene_set.keys())
        if missing:
            raise ValidationError(
                f"gene_sets[{i}] (name='{gene_set.get('name', 'unknown')}') "
                f"missing required fields: {sorted(missing)}"
            )

        # Validate name is not empty
        if not gene_set["name"] or not isinstance(gene_set["name"], str):
            raise ValidationError(
                f"gene_sets[{i}] has invalid name: '{gene_set.get('name')}'. "
                f"Name must be a non-empty string."
            )


def validate_config(config: dict) -> None:
    """Comprehensive configuration validation.

    Validates:
    - Required fields are present
    - Comparison tool is valid
    - Decomposition modes are valid
    - Gene sets have proper structure
    - File paths are strings (existence checked separately)

    Args:
        config: Configuration dictionary from config.yaml

    Raises:
        ValidationError: If configuration is invalid

    Example:
        >>> validate_config({
        ...     "sample": "NA12878",
        ...     "caller": "GATK",
        ...     "pipeline_version": "v1.0",
        ...     "comparison_tool": "happy",
        ...     "reference": "/path/to/ref.fa",
        ...     "truth_vcf": "/path/to/truth.vcf.gz",
        ...     "query_vcfs": ["/path/to/query.vcf.gz"],
        ...     "decomposition_modes": ["decomposed"],
        ...     "database": "data/concorde.db",
        ...     "output_dir": "results"
        ... })
    """
    # Required top-level fields
    required_fields = {
        "sample",
        "caller",
        "pipeline_version",
        "comparison_tool",
        "reference",
        "truth_vcf",
        "query_vcfs",
        "decomposition_modes",
        "database",
        "output_dir",
    }

    missing = required_fields - set(config.keys())
    if missing:
        raise ValidationError(
            f"Config missing required fields: {sorted(missing)}"
        )

    # Validate sample name
    if not config["sample"] or not isinstance(config["sample"], str):
        raise ValidationError(
            f"Config field 'sample' must be a non-empty string, "
            f"got: {config.get('sample')}"
        )

    # Validate caller
    if not config["caller"] or not isinstance(config["caller"], str):
        raise ValidationError(
            f"Config field 'caller' must be a non-empty string, "
            f"got: {config.get('caller')}"
        )

    # Validate pipeline_version
    if not config["pipeline_version"] or not isinstance(config["pipeline_version"], str):
        raise ValidationError(
            f"Config field 'pipeline_version' must be a non-empty string, "
            f"got: {config.get('pipeline_version')}"
        )

    # Validate mode (default to germline for backward compatibility)
    mode = config.get("mode", "germline")
    validate_mode(mode)

    # Validate comparison tool (mode-aware)
    validate_comparison_tool(config["comparison_tool"], mode)

    # Validate somatic config if in somatic mode
    if mode == "somatic":
        if "somatic" not in config:
            raise ValidationError(
                "Config must include 'somatic' block when mode is 'somatic'."
            )
        validate_somatic_config(config["somatic"])

    # Validate truth_matching config (optional)
    if "truth_matching" in config and config["truth_matching"]:
        validate_truth_matching_config(config["truth_matching"], mode)

    # Validate decomposition modes
    if not isinstance(config["decomposition_modes"], list):
        raise ValidationError(
            f"Config field 'decomposition_modes' must be a list, "
            f"got: {type(config['decomposition_modes']).__name__}"
        )
    validate_decomposition_modes(config["decomposition_modes"])

    # Validate query_vcfs is a list
    if not isinstance(config["query_vcfs"], list):
        raise ValidationError(
            f"Config field 'query_vcfs' must be a list, "
            f"got: {type(config['query_vcfs']).__name__}"
        )

    if not config["query_vcfs"]:
        raise ValidationError(
            "Config field 'query_vcfs' cannot be empty. "
            "At least one query VCF must be specified."
        )

    # Validate file path fields are strings
    path_fields = ["reference", "truth_vcf", "database", "output_dir"]
    for field in path_fields:
        if not isinstance(config[field], str):
            raise ValidationError(
                f"Config field '{field}' must be a string, "
                f"got: {type(config[field]).__name__}"
            )

    # Validate confident_regions (optional)
    if "confident_regions" in config and config["confident_regions"]:
        if not isinstance(config["confident_regions"], str):
            raise ValidationError(
                f"Config field 'confident_regions' must be a string, "
                f"got: {type(config['confident_regions']).__name__}"
            )

    # Validate gene_sets (optional)
    if "gene_sets" in config and config["gene_sets"]:
        if not isinstance(config["gene_sets"], list):
            raise ValidationError(
                f"Config field 'gene_sets' must be a list, "
                f"got: {type(config['gene_sets']).__name__}"
            )
        validate_gene_sets_config(config["gene_sets"])

    # Validate VEP config (optional)
    if "vep" in config and config["vep"]:
        validate_vep_config(config["vep"])

    # Validate nfcore_qc_dir (optional)
    if "nfcore_qc_dir" in config and config["nfcore_qc_dir"]:
        if not isinstance(config["nfcore_qc_dir"], str):
            raise ValidationError(
                f"Config field 'nfcore_qc_dir' must be a string, "
                f"got: {type(config['nfcore_qc_dir']).__name__}"
            )

