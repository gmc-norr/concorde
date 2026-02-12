"""Constants and enumerations for the Concorde pipeline.

Centralizes magic strings into named constants for type safety and IDE support.
"""


class ExecutionMode:
    """Pipeline execution modes."""

    GERMLINE = "germline"
    SOMATIC = "somatic"

    ALL = {GERMLINE, SOMATIC}


class SomaticCallingMode:
    """Somatic calling modes."""

    TUMOR_ONLY = "tumor_only"
    TUMOR_NORMAL = "tumor_normal"

    ALL = {TUMOR_ONLY, TUMOR_NORMAL}


class ComparisonTool:
    """Variant comparison tools."""

    HAPPY = "happy"
    RTG = "rtg"
    SOMPY = "sompy"
    INTERNAL = "internal"

    GERMLINE_TOOLS = {HAPPY, RTG}
    SOMATIC_TOOLS = {SOMPY, RTG, INTERNAL}
    ALL = GERMLINE_TOOLS | SOMATIC_TOOLS


class ClassificationCategory:
    """Variant classification categories."""

    TP = "TP"
    FP = "FP"
    FN = "FN"
    TP_GT_MISMATCH = "TP_GT_MISMATCH"
    TP_VAF_DISCORDANT = "TP_VAF_DISCORDANT"
    PARTIAL = "PARTIAL"

    GERMLINE = {TP, FP, FN, TP_GT_MISMATCH}
    SOMATIC = {TP, FP, FN, TP_VAF_DISCORDANT, PARTIAL}
    ALL = GERMLINE | SOMATIC


class VariantType:
    """Variant type classifications."""

    SNP = "SNP"
    INDEL = "INDEL"
    COMPLEX = "COMPLEX"
    SV = "SV"

    ALL = {SNP, INDEL, COMPLEX, SV}


class DiffClassification:
    """Difference classification for continuous verification."""

    EQUIVALENT = "equivalent"
    NUMERICAL_DRIFT = "numerical_drift"
    BIOLOGICAL_DIFFERENCE = "biological_difference"

    ALL = {EQUIVALENT, NUMERICAL_DRIFT, BIOLOGICAL_DIFFERENCE}


class EvidenceScore:
    """Root cause evidence scoring."""

    STRONG = "STRONG"
    MODERATE = "MODERATE"
    WEAK = "WEAK"
    NONE = "NONE"

    ALL = {STRONG, MODERATE, WEAK, NONE}


class AcceptanceDecision:
    """Acceptance criteria outcomes."""

    PASS = "PASS"
    FAIL = "FAIL"
    CONDITIONAL_PASS = "CONDITIONAL_PASS"
    REVIEW_REQUIRED = "REVIEW_REQUIRED"

    ALL = {PASS, FAIL, CONDITIONAL_PASS, REVIEW_REQUIRED}


class RiskTier:
    """Risk-based tiering for acceptance criteria."""

    TIER_1 = 1  # Critical (regulated panels, HIGH impact)
    TIER_2 = 2  # Standard (coding variants)
    TIER_3 = 3  # Informational (non-coding, off-panel)

    ALL = {TIER_1, TIER_2, TIER_3}


class VEPImpact:
    """VEP functional impact tiers."""

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MODIFIER = "MODIFIER"

    ALL = {HIGH, MODERATE, LOW, MODIFIER}
