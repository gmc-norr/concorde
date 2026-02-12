# Acceptance Criteria

The acceptance engine evaluates per-stratum metrics against tiered thresholds to produce a deterministic pass/fail decision. When BRAM is enabled, a two-stage decision combines deterministic and probabilistic assessments.

## Three-Tier Risk Model

Strata are assigned to risk tiers based on clinical significance. Higher tiers enforce stricter thresholds, reflecting the proportionally greater consequence of missed variants in clinically actionable regions.

| Tier | Label | Assignment | Threshold Strictness |
|------|-------|------------|---------------------|
| Tier 1 | Critical | Regulated gene panels (e.g., ACMG SF v3.1), HIGH impact variants | Strictest |
| Tier 2 | Standard | MODERATE and LOW impact coding variants | Standard |
| Tier 3 | Informational | MODIFIER variants, off-panel regions | Relaxed |

## Tier Assignment

A stratum is assigned to a tier based on the following priority:

1. If the stratum's gene panel is listed in a tier's `panels`, it receives that tier
2. If the stratum's functional impact is listed in a tier's `impacts`, it receives that tier
3. Otherwise, the stratum defaults to Tier 3

## Decision Outcomes

| Outcome | Condition |
|---------|-----------|
| `PASS` | All strata meet their tier thresholds |
| `FAIL` | One or more Tier 1 or Tier 2 strata violate thresholds |
| `CONDITIONAL_PASS` | Only Tier 3 violations detected |
| `REVIEW_REQUIRED` | BRAM escalated a deterministic PASS or CONDITIONAL_PASS |

## Threshold Configuration

Thresholds can be set at three levels of specificity (most specific wins):

1. **Global** - Baseline thresholds for all strata
2. **Per-stratum** - Override for specific variant types (e.g., SNP, INDEL)
3. **Per-panel** - Override for specific gene panels (e.g., ACMG59)

```yaml
global:
  min_sensitivity: 0.99
  min_precision: 0.999

per_stratum:
  SNP:
    min_sensitivity: 0.995
  INDEL:
    min_sensitivity: 0.98

per_panel:
  ACMG59:
    min_sensitivity: 0.999
    min_precision: 0.9999
```

## Violations

When a stratum fails its thresholds, a violation record is created containing:

- The dimension and stratum name
- The assigned tier
- Which metric(s) violated the threshold
- The observed vs. required values

Violations are included in the acceptance section of the validation report.

## Two-Stage Decision (with BRAM)

When BRAM is enabled, the acceptance engine runs a two-stage evaluation:

1. **Stage 1 (Deterministic)**: Evaluate per-stratum metrics against thresholds as described above
2. **Stage 2 (Probabilistic)**: BRAM computes the posterior probability of clinically meaningful degradation for each stratum

The final verdict combines both stages. BRAM can escalate a PASS or CONDITIONAL_PASS to REVIEW_REQUIRED, but never overrides a deterministic FAIL.

| Deterministic | BRAM | Final Verdict |
|--------------|------|---------------|
| PASS | BRAM_PASS | PASS |
| PASS | BRAM_FLAG | REVIEW_REQUIRED |
| CONDITIONAL_PASS | BRAM_PASS | CONDITIONAL_PASS |
| CONDITIONAL_PASS | BRAM_FLAG | REVIEW_REQUIRED |
| REVIEW_REQUIRED | any | REVIEW_REQUIRED |
| FAIL | any | FAIL |

See [Bayesian Risk Assessment](bram.md) for the full BRAM methodology.

## Low-Count Strata

Strata with fewer than `low_count_threshold` (default: 20) variants receive special treatment based on the `low_count_policy`:

| Policy | Behavior |
|--------|----------|
| `warn` | Include with a low-confidence flag in the report (default) |
| `exclude` | Exclude from acceptance evaluation entirely |
| `fail` | Treat as an automatic failure |

## Configuration Reference

See `config/acceptance_criteria.yaml` for the full configuration. Key fields:

```yaml
acceptance_criteria:
  mode: "germline"
  global:
    min_sensitivity: 0.99
    min_precision: 0.999
  per_stratum: { ... }
  per_panel: { ... }
  tiers:
    tier_1:
      label: "Critical"
      panels: ["ACMG59"]
      impacts: ["HIGH"]
    tier_2:
      label: "Standard"
      impacts: ["MODERATE", "LOW"]
    tier_3:
      label: "Informational"
      impacts: ["MODIFIER"]
  low_count_threshold: 20
  low_count_policy: "warn"
```
