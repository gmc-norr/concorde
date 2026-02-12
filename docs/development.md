# Developer Guide

## Getting Started

```bash
# Clone and install
git clone <repo-url>
cd concorde-pipeline
pixi install

# Verify setup
.pixi/envs/default/bin/pytest tests/ -v
```

### Key Commands

```bash
pixi run test        # Run all 705 tests
pixi run lint        # Ruff linting
pixi run format      # Ruff formatting
pixi run typecheck   # mypy type checking
pixi run security    # Bandit security scan
pixi run audit       # pip-audit dependency check
```

## Architecture

### Snakemake Script Pattern

Scripts use a `try/except NameError` pattern for testability:

```python
def my_function():
    """Core logic at module top level."""
    ...

try:
    from typing import TYPE_CHECKING
    if TYPE_CHECKING:
        from snakemake.script import Snakemake
        snakemake: Snakemake
    else:
        snakemake = snakemake  # type: ignore  # noqa: F821
    # Snakemake execution code here
except NameError:
    pass  # Not running via Snakemake
```

### Key Dependencies

| Package | Purpose |
|---------|---------|
| snakemake | Workflow orchestration |
| sqlalchemy | ORM for SQLite database (2.0 `Mapped[]` types) |
| pysam | VCF/BAM/FASTA file parsing |
| pandas | Data manipulation |
| scipy | Statistical inference (Wilson CIs, bootstrap CIs, BRAM CDFs) |
| numpy | Numerical computation (SE, stdev, aggregation) |
| jinja2 | HTML report templating |
| weasyprint | PDF generation |
| pyyaml | YAML configuration parsing |

### Statistical Computation Policy

All statistical computations delegate to established numerical libraries (scipy, numpy) rather than hand-rolled implementations:

| Computation | Library | Rationale |
|-------------|---------|-----------|
| Wilson CI z-quantile | `scipy.stats.norm.ppf` | Avoids manual z-table lookup |
| F1 confidence interval | `scipy.stats.bootstrap` | Calibrated percentile intervals for nonlinear statistics |
| BRAM tail probability | `scipy.stats.beta.cdf` | Beta CDF for conjugate posterior |
| BRAM prior std | `math.sqrt` | Exact Beta distribution formula |
| Sample variance (MoM prior) | `numpy.var(ddof=1)` | Bessel's correction, consistent across BRAM/SPC/trending |
| Risk aggregation | `numpy.prod`, `numpy.mean` | Numerical stability for product chains |
| Precision/recall/F1 | `compute_classification_metrics()` | Single shared utility with consistent None semantics |

### Module Organization

```
scripts/
  analysis/           # 19 modules
    acceptance.py         # Three-tier acceptance engine
    audit.py              # Audit trail recording
    baseline_manager.py   # Baseline locking, signing, envelopes
    bayesian_risk.py      # BRAM engine (Beta-Binomial conjugate)
    create_baseline.py    # Baseline creation from runs
    diff_classifier.py    # Drift vs biological classification
    ensemble.py           # Multi-caller concordance and filtering
    extended_metrics.py   # GT concordance, Ti/Tv, Het/Hom, CIs, MCC
    fixture_registry.py   # Input fixture checksum registry
    giab_stratifications.py  # GIAB BED file auto-discovery
    intersect_gene_sets.py   # Gene set BED intersection
    qc_trending.py        # Historical z-score trending
    region_annotator.py   # BED/FASTA region annotation
    root_cause.py         # Root cause evidence collection
    run_verification.py   # Verification workflow orchestrator
    spc.py                # Statistical Process Control charts
    stratification.py     # 11-dimension stratification engine
    variant_diff.py       # Variant-level diff computation
    verdicts.py           # Verdict constants (PASS, FAIL, BRAM_*, tiers)
  ingestion/
    ingest.py             # Main ingestion orchestrator
    ingest_helpers.py     # Row mapping, zygosity derivation, metric propagation
  matching/
    internal_matcher.py   # Position+allele somatic matcher
  parsers/
    parse_happy_vcf.py    # hap.py VCF parser (with indel_size, GT)
    parse_happy_metrics.py # hap.py summary metrics parser
    parse_rtg_vcf.py      # RTG vcfeval VCF parser (with indel_size, GT)
    parse_rtg_metrics.py  # RTG summary metrics parser
    parse_sompy_vcf.py    # som.py VCF parser
    parse_sompy_metrics.py # som.py summary metrics parser
    parse_vep_annotations.py  # VEP CSQ/ANN annotation parser
  reporting/
    generate_report.py    # Report orchestrator (queries DB, feeds templates)
    json_report.py        # JSON report generation
    html_report.py        # HTML report via Jinja2
    pdf_report.py         # PDF via weasyprint
    report_storage.py     # Report naming and checksums
    templates/
      report.html.j2      # HTML report template
  validation/
    validate_inputs.py    # Pre-flight validation orchestrator
  validators/             # 31 validators across 8 modules
  utils.py                # Shared utilities (safe_float/int, classification metrics)
  vcf_utils.py            # VCF helper functions
```

## Code Style

- **PEP 8** with ruff enforcement
- **Type hints** required for all function signatures
- **Google-style docstrings** for public functions
- **SQLAlchemy 2.0** `Mapped[]` type annotations for models
- Classes: `PascalCase`, functions: `snake_case`, constants: `UPPER_SNAKE_CASE`

### Import Order

```python
# 1. __future__ imports
from __future__ import annotations

# 2. Standard library
import logging
from pathlib import Path

# 3. Third-party
from sqlalchemy.orm import Mapped, mapped_column

# 4. Local
from models.variant import Variant
```

## Testing

### Test Organization (705 tests)

| Directory | Coverage |
|-----------|----------|
| `tests/analysis/` | Stratification, baselines, acceptance, BRAM, extended metrics, cross-validation, diff classifier, root cause, variant diff, region annotator, QC trending, SPC, audit, ensemble |
| `tests/ingestion/` | VEP merge, ingestion helpers, wiring integration tests |
| `tests/matching/` | Internal matcher (position+allele, VAF concordance) |
| `tests/models/` | All 19 ORM models including relationships and cascade deletes |
| `tests/parsers/` | som.py VCF/metrics, VEP annotation parsing |
| `tests/reporting/` | JSON schema, HTML rendering, PDF generation, BRAM reporting |
| `tests/validators/` | Config validation (germline, somatic, VEP), input files, data quality |

### Running Tests

```bash
# All tests
.pixi/envs/default/bin/pytest tests/ -v

# Single test with output
.pixi/envs/default/bin/pytest tests/analysis/test_bayesian_risk.py -v -s

# Keyword filter
.pixi/envs/default/bin/pytest tests/ -k "test_bram" -v

# Drop into debugger on failure
.pixi/envs/default/bin/pytest tests/ --pdb
```

### Expected Lint Errors

These are expected and acceptable:

- **E402** in test files: imports after `sys.path` modification
- **F821** in models: forward string references in SQLAlchemy relationships
- **F821** in Snakemake scripts: `snakemake` variable injected at runtime
- **F401** in `verify_validators.py`: intentional unused imports

## Adding New Features

### New Parser

1. Create `scripts/parsers/parse_<tool>.py` following existing patterns
2. Use pysam for VCF parsing with TYPE_CHECKING guard
3. Map tool-specific classifications to standard TP/FP/FN
4. Compute `indel_size` for INDEL variants
5. Extract `truth_gt` and `query_gt` where available
6. Create tests in `tests/parsers/test_parse_<tool>.py`
7. Add conditional rule in `Snakefile`

### New Stratification Dimension

1. Add dimension definition in `config/stratifications.yaml`
2. Add `_assign_<dimension>()` function in `scripts/analysis/stratification.py`
3. Add dimension block in `StratificationEngine.stratify()`
4. Ensure the Variant model has necessary fields (add with `nullable=True`)
5. Add tests in `tests/analysis/test_stratification.py`

### New Database Model

1. Create model in `models/<name>.py` using SQLAlchemy 2.0 `Mapped[]` types
2. Add relationships to related models
3. Import in `models/__init__.py` and add to `__all__`
4. Add tests in `tests/models/test_models.py`
5. Update ingestion scripts if needed

### New Validator

1. Add validation function in the appropriate `scripts/validators/` module
2. Register in `scripts/validation/validate_inputs.py`
3. Add tests in `tests/validators/`

## Debugging

### Pipeline Debugging

```bash
# Verbose Snakemake output
snakemake --cores all --configfile config/config.yaml --verbose --printshellcmds

# Keep going after errors
snakemake --cores all --configfile config/config.yaml --keep-going

# Debug DAG for specific target
snakemake --cores 1 --configfile config/config.yaml --debug-dag \
  results/NA12878_GATK_v1.0.0_decomposed/variants.tsv
```

### Test Debugging

```bash
# Run single test with verbose output
.pixi/envs/default/bin/pytest tests/analysis/test_bayesian_risk.py::TestBRAMEngine::test_assess_all_flag -v -s

# Run with print output visible
.pixi/envs/default/bin/pytest tests/ -k "test_name" -s
```

## Advanced Usage

### Baseline Management

Create a baseline from a validated run for continuous verification:

```python
from analysis.baseline_manager import BaselineManager

manager = BaselineManager(session)
baseline = manager.create_baseline(
    run_id=run.id,
    name="v4.5.0-baseline",
    description="Initial validated baseline"
)
# Auto-computes statistical envelopes from stratified metrics
```

### Continuous Verification

Compare a new run against an existing baseline to detect drift:

```python
from analysis.variant_diff import VariantDiffEngine
from analysis.diff_classifier import DiffClassifier

diff_engine = VariantDiffEngine()
transitions = diff_engine.compute_diff(baseline_variants, new_variants)

classifier = DiffClassifier()
classified = classifier.classify_transitions(transitions)
```

### Cluster Execution

```bash
snakemake --cores all --configfile config/config.yaml \
  --cluster "sbatch -p compute -c {threads} --mem={resources.mem_mb}" \
  --jobs 100
```

## References

### Development Tools

- [Snakemake](https://snakemake.readthedocs.io/)
- [SQLAlchemy 2.0](https://docs.sqlalchemy.org/)
- [pysam](https://pysam.readthedocs.io/)
- [Ruff](https://docs.astral.sh/ruff/)
- [pytest](https://docs.pytest.org/)
- [pixi](https://pixi.sh/)

### External Bioinformatics Tools

- [hap.py](https://github.com/Illumina/hap.py) — Germline variant benchmarking
- [RTG vcfeval](https://github.com/RealTimeGenomics/rtg-tools) — Variant comparison
- [som.py](https://github.com/Illumina/hap.py/blob/master/doc/sompy.md) — Somatic variant benchmarking
- [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) — Variant effect prediction
