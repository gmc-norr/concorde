# Reporting

Concorde generates validation reports in three formats: JSON (machine-readable), HTML (human-readable), and PDF (print-ready). All formats derive from the same underlying data and include stratified metrics with confidence intervals, acceptance results, BRAM posteriors, and verification data.

## JSON Report (S16.1)

Machine-readable structured output. This is the primary report format and the source of truth for all other formats.

### Schema

```json
{
  "report_type": "initial_validation",
  "report_version": "1.0",
  "generated_at": "2026-02-12T12:00:00+00:00",
  "mode": "germline",
  "run_metadata": {
    "sample": "NA12878",
    "pipeline_version": "v4.5.0",
    "caller": "GATK-HC",
    "comparison_tool": "happy"
  },
  "input_checksums": { ... },
  "metrics": {
    "global": {
      "sensitivity": 0.995,
      "precision": 0.999,
      "f1": 0.997
    },
    "per_stratum": [
      {
        "dimension": "variant_class",
        "stratum": "SNP",
        "tp_count": 4850,
        "fp_count": 5,
        "fn_count": 20,
        "precision": 0.999,
        "recall": 0.996,
        "f1": 0.997,
        "precision_ci_lower": 0.997,
        "precision_ci_upper": 0.999,
        "recall_ci_lower": 0.994,
        "recall_ci_upper": 0.998,
        "f1_ci_lower": 0.996,
        "f1_ci_upper": 0.998,
        "mcc": 0.985,
        "genotype_concordance": 0.998,
        "low_confidence": false
      }
    ],
    "qc_metrics": {
      "titv_ratio": 2.08,
      "het_hom_ratio": 1.65
    }
  },
  "acceptance": {
    "overall_result": "PASS",
    "per_stratum_results": [ ... ],
    "violations": [ ... ]
  },
  "bayesian_risk": {
    "enabled": true,
    "verdict": "BRAM_PASS",
    "aggregate_risk_score": 0.12,
    "mean_risk_score": 0.08,
    "per_stratum": [ ... ]
  },
  "qc_trending": [
    {
      "metric_name": "sensitivity",
      "current_value": 0.995,
      "historical_mean": 0.996,
      "historical_stdev": 0.001,
      "z_score": -1.0,
      "flagged": false
    }
  ],
  "variant_diff": [ ... ],
  "root_cause_evidence": [ ... ]
}
```

### Report Types

| Type | Description |
|------|-------------|
| `initial_validation` | First-time validation against truth set |
| `continuous_verification` | Re-verification against a baseline (includes `variant_diff`, `root_cause_evidence`, and BRAM results) |

### Extended Metric Fields

Per-stratum metric records include:

| Field | Description |
|-------|-------------|
| `precision_ci_lower`, `precision_ci_upper` | Wilson score 95% CI for precision |
| `recall_ci_lower`, `recall_ci_upper` | Wilson score 95% CI for recall |
| `f1_ci_lower`, `f1_ci_upper` | Bootstrap percentile 95% CI for F1 |
| `mcc` | Matthews Correlation Coefficient (-1 to +1) |
| `genotype_concordance` | Fraction of TPs with matching truth/query GT |

## HTML Report (S16.2)

Human-readable report rendered from Jinja2 templates. Contains the following sections:

1. **Header** - Overall result banner (color-coded: green/red/amber), run metadata
2. **Executive Summary** - Global sensitivity, precision, F1; Ti/Tv ratio; Het/Hom ratio; violation count
3. **Per-Stratum Metrics** - Table with TP/FP/FN counts, precision, recall, F1, CI bounds, MCC, and genotype concordance per stratum. Low-confidence strata are italicized.
4. **Acceptance Results** - Per-stratum pass/fail with tier assignment and threshold comparison
5. **Bayesian Risk Assessment** - BRAM verdict banner, configuration summary, per-stratum risk table with flagged rows highlighted
6. **QC Trending** - Table of metric names with current value, historical mean, historical stdev, z-score, and flagged status (when trending is enabled)
7. **Variant Differences** - Variant-level transitions between baseline and verification with drift classification column (first 100 shown; verification reports only)
8. **Traceability** - Input file SHA-256 checksums

### Template

The HTML template is at `scripts/reporting/templates/report.html.j2`. Custom Jinja2 filters:

| Filter | Purpose | Example |
|--------|---------|---------|
| `format_pct` | Float to percentage | `0.995` -> `99.50%` |
| `format_float` | Float to N decimals | `0.12345` -> `0.1235` |
| `status_color` | Status to CSS class | `PASS` -> `green` |

## PDF Report (S16.3)

Print-ready PDF generated from the HTML report using weasyprint. Preserves all sections with proper page breaks and styling.

## Report Storage

Reports follow a naming convention:

```
{sample}_{caller}_{version}_{mode}_{timestamp}.{ext}
```

Example: `NA12878_GATK-HC_v4.5.0_germline_20260212T120000.json`

Each report includes a SHA-256 checksum computed from the serialized JSON for integrity verification.

## Output Location

```
results/
  {sample}_{caller}_{version}_{mode}/
    report.json
    report.html
    report.pdf
```
