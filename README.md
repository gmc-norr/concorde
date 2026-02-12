# Concorde Pipeline

> **Status: Under Active Development** — This project is not yet ready for production use. APIs, configuration formats, and database schemas may change without notice. Use at your own risk.

Concorde is a computational framework for the systematic validation and continuous verification of next-generation sequencing (NGS) variant calling pipelines. It implements a multi-stage analytical workflow encompassing truth-set concordance, stratified performance metrics, deterministic acceptance criteria, and Bayesian risk assessment to provide statistically rigorous evaluation of variant calling accuracy.

## Motivation

Clinical and research variant calling pipelines require quantitative evidence that sensitivity, precision, and genotype concordance meet defined performance thresholds. Simple global metrics (e.g., overall sensitivity) mask clinically important variation across variant classes, genomic contexts, and allele frequency ranges. Concorde addresses this by decomposing performance into fine-grained strata and applying both frequentist confidence intervals and Bayesian posterior inference to each stratum independently.

## Analytical Capabilities

### Variant Concordance

Concorde supports four concordance engines for truth-vs-query comparison:

| Engine | Application | Method |
|--------|------------|--------|
| hap.py | Germline | Haplotype-aware comparison with decomposition |
| RTG vcfeval | Germline / Somatic | Diploid-aware global alignment of variant representations |
| som.py | Somatic | Somatic-aware comparison with allele frequency tracking |
| Internal matcher | Somatic | Position + allele matching with configurable VAF concordance |

Each engine classifies variants as TP, FP, or FN. The internal matcher additionally tracks TP_VAF_DISCORDANT (correct position/allele but discordant allele frequency) and PARTIAL (position match, allele mismatch) classifications.

### Stratified Performance Metrics

Variants are assigned to strata across **11 independent dimensions**, and TP/FP/FN counts, precision, recall, and F1 are computed per stratum:

| Dimension | Mode | Strata | Data Source |
|-----------|------|--------|-------------|
| `variant_class` | Both | SNP, INDEL, COMPLEX | Variant type |
| `indel_size` | Both | 1bp, 2-5bp, 6-15bp, 16-50bp, >50bp | Allele length difference |
| `functional_impact` | Both | HIGH, MODERATE, LOW, MODIFIER | VEP annotation |
| `gene_panel` | Both | Configured gene sets + off_panel | BED intersection |
| `zygosity` | Germline | HET, HOM_ALT | Genotype field |
| `vaf_bins` | Somatic | <0.05, 0.05-0.10, 0.10-0.20, 0.20-0.50, >0.50 | Tumor allele frequency |
| `low_complexity` | Both | low_complexity, non_lcr | BED annotation |
| `gc_content` | Both | 0-25%, 25-30%, 30-55%, 55-65%, 65-100% | Reference sequence |
| `segdup` | Both | in_segdup, non_segdup | BED annotation |
| `mappability` | Both | low, medium, high | BED annotation |
| `coverage_depth` | Both | <10x, 10-30x, 30-50x, >50x | Variant DP field |

Strata with fewer than a configurable minimum number of variants (default: 20) are flagged as low-confidence.

### Extended Metrics

Beyond precision, recall, and F1, the following metrics are computed per stratum:

- **Wilson score confidence intervals** for precision and recall, using `scipy.stats.norm` for quantile computation. Wilson intervals are preferred over normal approximation for their accuracy at extreme proportions and small sample sizes.
- **Bootstrap confidence intervals** for F1, using `scipy.stats.bootstrap` with 9,999 percentile resamples. F1 is a nonlinear function of precision and recall; bootstrap resampling provides calibrated intervals without distributional assumptions.
- **Matthews Correlation Coefficient (MCC)**: (TP\*TN - FP\*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN)). MCC is a balanced measure that accounts for all four confusion matrix quadrants.
- **Genotype concordance**: fraction of true-positive variants where truth and query genotypes match.
- **Ti/Tv ratio**: transition-to-transversion ratio for SNPs. Expected values are approximately 2.0-2.1 for whole-genome and 2.8-3.0 for exome sequencing.
- **Het/Hom ratio**: heterozygous-to-homozygous-alternate ratio. Expected approximately 1.5-2.0 for diploid organisms.

### Acceptance Criteria

A three-tier risk-based acceptance model evaluates per-stratum metrics against configurable thresholds:

| Tier | Scope | Strictness |
|------|-------|------------|
| Tier 1 (Critical) | Regulated gene panels (e.g., ACMG SF v3.1), HIGH-impact variants | Strictest |
| Tier 2 (Standard) | MODERATE and LOW impact coding variants | Standard |
| Tier 3 (Informational) | MODIFIER variants, off-panel regions | Relaxed |

Decision outcomes: PASS, FAIL, CONDITIONAL_PASS (only Tier 3 violations), or REVIEW_REQUIRED (BRAM escalation).

### Bayesian Risk Assessment Module (BRAM)

BRAM augments deterministic acceptance with calibrated uncertainty estimates using a Beta-Binomial conjugate Bayesian model. Sensitivity and precision are binomial proportions, and the Beta distribution is the natural conjugate prior. For each stratum:

1. **Observation counts**: successes and trials from raw TP/FP/FN counts (or reconstructed from proportions when counts are unavailable)
2. **Beta prior**: centered at the baseline value with concentration parameter kappa: alpha = baseline \* kappa, beta = (1 - baseline) \* kappa. Prior source priority: per-stratum override > empirical (method-of-moments Beta fit from historical data, using `numpy.var(ddof=1)`) > default (kappa=20)
3. **Conjugate update**: alpha\_post = alpha\_prior + successes, beta\_post = beta\_prior + (trials - successes)
4. **Tail probability**: P(p\_true < baseline - threshold | data) via `scipy.stats.beta.cdf`
5. **Per-tier degradation thresholds**: Tier 1 = 0.005, Tier 2 = 0.01, Tier 3 = 0.02 (configurable)
6. **Aggregation**: per-stratum tail probabilities combined via max, weighted mean, or joint probability (1 - prod(1 - risk\_i), computed with `numpy.prod`)

BRAM produces a two-stage decision: it can escalate a deterministic PASS to REVIEW_REQUIRED but never overrides a deterministic FAIL.

### Continuous Verification

For regression testing across pipeline versions:

- **Baseline management**: locked, signed baselines with statistical envelopes per stratum
- **Variant-level diff**: identifies gained, lost, and changed variants between runs
- **Drift classification**: distinguishes technical drift from biological variation
- **Root cause analysis**: evidence-based scoring with per-transition attribution
- **QC trending**: z-score monitoring of metrics across historical runs
- **Statistical Process Control (SPC)**: control chart computation with configurable sigma multipliers

### Multi-Caller Support

- **Ensemble analysis**: cross-caller concordance matrices and ensemble filtering (majority vote, intersection, union)
- **Audit trail**: timestamped event log for ingestion, baseline creation, locking, and verification

### Reporting

Reports are generated in three formats:

- **JSON**: machine-readable structured output with full per-stratum metrics, confidence intervals, MCC, acceptance results, BRAM posteriors, variant diffs, and QC trending data
- **HTML**: rendered from Jinja2 templates with color-coded verdicts, per-stratum tables (including CI bounds and MCC columns), BRAM risk tables, variant diff with drift classification, and QC trending
- **PDF**: print-ready rendering via weasyprint

## Architecture

Concorde is orchestrated by Snakemake with 11 rules executing in dependency order. All analytical logic is implemented as pure functions testable independently of the workflow engine.

```
Snakefile                    # Workflow orchestration (11 rules)
scripts/
  analysis/                  # 19 modules: stratification, baselines, acceptance,
                             #   BRAM, extended metrics, region annotation, RCA,
                             #   QC trending, SPC, audit, ensemble, verification
  ingestion/                 # Database ingestion with region annotation wiring
  matching/                  # Internal somatic variant matcher
  parsers/                   # 7 parsers (hap.py, RTG, som.py, VEP, nf-core QC)
  reporting/                 # JSON, HTML, PDF report generation
  validation/                # Pre-flight input validation
  validators/                # 31 validators across 8 modules
models/                      # 19 SQLAlchemy 2.0 ORM models
config/                      # 4 YAML configuration files
tests/                       # 705 tests across 8 directories
```

### Database Schema

19 SQLAlchemy ORM models stored in SQLite:

| Model | Purpose |
|-------|---------|
| Run | Pipeline run metadata (sample, caller, mode, ensemble) |
| Variant | Classified variants with genotype, VEP, somatic, and region annotations |
| Metric | Aggregate metrics with confidence intervals and MCC |
| StratifiedMetric | Per-stratum metrics with CIs, MCC, and genotype concordance |
| QCMetric | QC metrics (coverage, Ti/Tv, Het/Hom) |
| Baseline, BaselineEnvelope | Locked baselines with per-stratum statistical envelopes |
| FixtureRecord | Input fixture SHA-256 checksums |
| VerificationResult, VariantTransition | Run-vs-baseline diffs with drift classification |
| RootCauseEvidence | Evidence items for variant transitions |
| BayesianRiskAssessment, StratumPosterior | BRAM aggregate and per-stratum posteriors |
| AuditEvent | Timestamped audit trail |
| GeneSet, SoftwareVersion, LLMAnalysis | Ancillary metadata |

### Statistical Dependencies

All statistical computations use established numerical libraries:

| Computation | Library | Function |
|-------------|---------|----------|
| Wilson score CI quantiles | scipy.stats | `norm.ppf` |
| F1 bootstrap CI | scipy.stats | `bootstrap` (percentile, n=9999) |
| BRAM tail probability | scipy.stats | `beta.cdf` |
| BRAM prior std (Beta) | math | `sqrt` (exact Beta formula) |
| Empirical prior (sample variance) | numpy | `var(ddof=1)` |
| Risk aggregation | numpy | `prod`, `mean`, `max` |
| Precision / recall / F1 | Shared utility | `compute_classification_metrics` |

## Quick Start

```bash
# Install dependencies
pixi install

# Configure
vim config/config.yaml

# Validate inputs (dry run)
pixi run dry-run

# Execute pipeline
pixi run run
```

## Configuration

| File | Purpose |
|------|---------|
| `config/config.yaml` | Pipeline settings, input paths, region BED files, QC trending, SPC |
| `config/stratifications.yaml` | 11 stratification dimensions with bin definitions |
| `config/acceptance_criteria.yaml` | Three-tier thresholds and decision logic |
| `config/bram_config.yaml` | Bayesian risk assessment parameters |

## Testing

```bash
# Run all 705 tests
pixi run test

# Run specific test suite
.pixi/envs/default/bin/pytest tests/analysis/ -v
.pixi/envs/default/bin/pytest tests/reporting/ -v
```

## Documentation

| Document | Content |
|----------|---------|
| [docs/getting-started.md](docs/getting-started.md) | Installation, prerequisites, first run |
| [docs/configuration.md](docs/configuration.md) | Full configuration reference |
| [docs/pipeline-rules.md](docs/pipeline-rules.md) | Snakemake rule descriptions |
| [docs/stratification.md](docs/stratification.md) | Stratification engine and extended metrics |
| [docs/acceptance.md](docs/acceptance.md) | Three-tier acceptance criteria |
| [docs/bram.md](docs/bram.md) | Bayesian risk assessment methodology |
| [docs/reporting.md](docs/reporting.md) | Report formats and schema |
| [docs/models.md](docs/models.md) | Database schema reference |
| [docs/development.md](docs/development.md) | Developer guide and architecture |

## Dependencies

| Package | Purpose |
|---------|---------|
| Snakemake | Workflow orchestration |
| SQLAlchemy 2.0 | ORM and database management |
| pysam | VCF/BAM/FASTA parsing |
| pandas | Tabular data manipulation |
| scipy | Statistical inference (CIs, CDFs, bootstrap) |
| numpy | Numerical computation |
| Jinja2 | HTML report templating |
| weasyprint | PDF generation |
