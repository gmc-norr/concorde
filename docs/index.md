# Concorde Pipeline

> **Status: Under Active Development** — This project is not yet ready for production use. APIs, configuration formats, and database schemas may change without notice. Use at your own risk.

Concorde is a computational framework for the systematic validation and continuous verification of next-generation sequencing (NGS) variant calling pipelines. It provides statistically rigorous evaluation of variant calling accuracy through stratified performance metrics, frequentist confidence intervals, and Bayesian posterior inference.

## Core Capabilities

- **Truth-set concordance** using hap.py, RTG vcfeval, som.py, or an internal position+allele matcher
- **11-dimension stratification** with per-stratum precision, recall, F1, Wilson score CIs, bootstrap F1 CIs, MCC, and genotype concordance
- **Three-tier acceptance criteria** with risk-proportional thresholds (Critical, Standard, Informational)
- **Bayesian Risk Assessment (BRAM)** using Beta-Binomial conjugate posteriors for uncertainty-aware degradation detection
- **Continuous verification** with locked baselines, variant-level diffs, drift classification, and root cause analysis
- **QC trending and SPC** for longitudinal metric monitoring with z-score flagging and control charts
- **Multi-caller ensemble analysis** with cross-caller concordance and ensemble filtering
- **Audit trail** for regulatory traceability
- **Multi-format reporting**: JSON (machine-readable), HTML (Jinja2), PDF (weasyprint)

## Architecture

The pipeline is orchestrated by Snakemake with 10 rules. Analytical logic resides in pure functions (19 analysis modules) testable independently of the workflow engine. Data is persisted in an SQLite database via 16 SQLAlchemy 2.0 ORM models. All statistical computations delegate to scipy and numpy.

```
scripts/
  analysis/       # Stratification, baselines, acceptance, BRAM, extended metrics,
                  #   region annotation, RCA, QC trending, SPC, audit, ensemble
  ingestion/      # Database ingestion with region annotation wiring
  matching/       # Internal somatic variant matcher
  parsers/        # 7 parsers (hap.py, RTG, som.py, VEP, nf-core QC)
  reporting/      # JSON, HTML, PDF report generation
  validation/     # Pre-flight input validation
  validators/     # 31 validators across 8 modules
models/           # 16 SQLAlchemy 2.0 ORM models
config/           # 4 YAML configuration files
tests/            # 705 tests across 8 directories
```

## Documentation

- [Getting Started](getting-started.md) - Installation, prerequisites, first run
- [Configuration](configuration.md) - Full configuration reference (4 YAML files)
- [Pipeline Rules](pipeline-rules.md) - Snakemake rule descriptions
- [Stratification](stratification.md) - Stratification engine and extended metrics
- [Acceptance Criteria](acceptance.md) - Three-tier acceptance model
- [Bayesian Risk Assessment](bram.md) - BRAM methodology and integration
- [Reporting](reporting.md) - Report formats and schema
- [Database Models](models.md) - Schema reference (16 models)
- [Development](development.md) - Developer guide, testing, architecture
