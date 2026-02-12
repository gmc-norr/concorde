# Database Models

Concorde uses SQLAlchemy 2.0 ORM with 19 models stored in an SQLite database. All new columns use `nullable=True` for backward compatibility with existing databases.

## Model Reference

| Model | Table | Purpose |
|-------|-------|---------|
| `Run` | `runs` | Pipeline run metadata (sample, caller, mode, ensemble) |
| `Variant` | `variants` | Classified variants with genotype, VEP, somatic, and region annotations |
| `Metric` | `metrics` | Aggregate performance metrics with CIs and MCC |
| `QCMetric` | `qc_metrics` | QC metrics (coverage, Ti/Tv, Het/Hom) |
| `SoftwareVersion` | `software_versions` | Tool versions used |
| `LLMAnalysis` | `llm_analyses` | LLM-generated QC insights |
| `GeneSet` | `gene_sets` | Gene set membership annotations |
| `StratifiedMetric` | `stratified_metrics` | Per-stratum metrics with CIs, MCC, and genotype concordance |
| `Baseline` | `baselines` | Validation baselines |
| `BaselineEnvelope` | `baseline_envelopes` | Statistical envelopes per stratum |
| `FixtureRecord` | `fixture_records` | Input fixture checksums |
| `VerificationResult` | `verification_results` | Run-vs-baseline comparison results |
| `VariantTransition` | `variant_transitions` | Per-variant diff classifications |
| `RootCauseEvidence` | `root_cause_evidence` | Evidence items for drift explanation |
| `BayesianRiskAssessment` | `bayesian_risk_assessments` | BRAM aggregate scores and verdict |
| `StratumPosterior` | `stratum_posteriors` | Per-stratum posterior computation details |
| `AuditEvent` | `audit_events` | Timestamped audit trail for traceability |

## Core Models

### Run

Pipeline run metadata. Central model that most others reference via `run_id`.

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Auto-increment primary key |
| `sample` | String | Sample name |
| `pipeline_version` | String | Pipeline version |
| `caller` | String | Variant caller name |
| `mode` | String(20) | "germline" or "somatic" |
| `tumor_sample` | String | Tumor sample name (somatic only) |
| `normal_sample` | String | Normal sample name (somatic only) |
| `calling_mode` | String | "tumor_only" or "tumor_normal" (somatic only) |
| `callers_json` | Text | JSON list of callers (ensemble mode) |
| `ensemble_method` | String(50) | Ensemble method (majority_vote, intersection, union) |

### Variant

Classified variants with genotype, VEP, somatic, and region annotations.

**Core fields**: `chrom`, `pos`, `ref`, `alt`, `type`, `classification`, `dp`, `gq`, `af`, `qual`, `filter_status`, `gene`

**Genotype fields**: `truth_gt`, `query_gt`, `gt_concordant`

**VEP fields**: `consequence`, `gene_symbol`, `gene_id`, `transcript_id`, `hgvsc`, `hgvsp`, `exon`, `protein_position`, `sift_prediction`, `polyphen_prediction`, `impact`

**Somatic fields**: `tumor_dp`, `tumor_af`, `normal_dp`, `normal_af`, `somatic_quality`, `caller_filter`

**Stratification fields**: `zygosity`, `indel_size`

**Region annotation fields**: `gc_content`, `in_low_complexity`, `in_segdup`, `mappability_score`

### Metric

Aggregate performance metrics with confidence intervals.

| Field | Type | Description |
|-------|------|-------------|
| `precision` | Float | TP / (TP + FP) |
| `recall` | Float | TP / (TP + FN) |
| `f1` | Float | Harmonic mean of precision and recall |
| `precision_ci_lower` | Float | Wilson score CI lower bound |
| `precision_ci_upper` | Float | Wilson score CI upper bound |
| `recall_ci_lower` | Float | Wilson score CI lower bound |
| `recall_ci_upper` | Float | Wilson score CI upper bound |
| `f1_ci_lower` | Float | Bootstrap CI lower bound |
| `f1_ci_upper` | Float | Bootstrap CI upper bound |
| `mcc` | Float | Matthews Correlation Coefficient |
| `genotype_concordance` | Float | GT match rate among TPs |

## Analysis Models

### StratifiedMetric

Per-stratum metric computed by the stratification engine. Contains the same CI, MCC, and genotype concordance fields as Metric, plus:

| Field | Type | Description |
|-------|------|-------------|
| `dimension` | String | Stratification dimension name |
| `stratum` | String | Stratum within the dimension |
| `variant_type` | String | Variant type filter (default "ALL") |
| `total_variants` | Integer | Total variants in this stratum |
| `tp_count` | Integer | True positive count |
| `fp_count` | Integer | False positive count |
| `fn_count` | Integer | False negative count |
| `low_confidence` | Boolean | Below minimum variant count |

### Baseline

A validated reference point for continuous verification.

| Field | Type | Description |
|-------|------|-------------|
| `run_id` | Integer (FK) | Source validation run |
| `name` | String | Baseline identifier |
| `description` | String | Human-readable description |
| `is_locked` | Boolean | Prevent modifications |
| `signature` | String | Digital signature hash |

### BayesianRiskAssessment

BRAM aggregate results for a verification run.

| Field | Type | Description |
|-------|------|-------------|
| `run_id` | Integer (FK) | Verification run |
| `baseline_id` | Integer (FK) | Reference baseline |
| `aggregate_risk_score` | Float | Combined risk across strata |
| `mean_risk_score` | Float | Mean per-stratum risk |
| `flagged_stratum_count` | Integer | Number of flagged strata |
| `verdict` | String | BRAM_PASS, BRAM_FLAG, or BRAM_NOT_RUN |
| `aggregation_method` | String | max, weighted_mean, or product |

### StratumPosterior

Per-stratum posterior computation details from BRAM.

| Field | Type | Description |
|-------|------|-------------|
| `dimension` | String | Stratification dimension |
| `stratum` | String | Stratum name |
| `metric_name` | String | Metric being assessed |
| `delta_observed` | Float | Current - baseline metric value |
| `prior_mu` | Float | Prior mean |
| `prior_sigma` | Float | Prior standard deviation |
| `sigma_obs` | Float (nullable) | Unused in Beta model (set to None) |
| `posterior_mu` | Float | Beta posterior mean: alpha_post / (alpha_post + beta_post) |
| `posterior_sigma` | Float | Beta posterior standard deviation |
| `tail_probability` | Float | P(p_true < baseline - threshold \| data) |
| `risk_weight` | Float | Always 1.0 (retained for backward compatibility) |
| `weighted_risk` | Float | Equals tail_probability (since risk_weight=1.0) |
| `flagged` | Boolean | Whether this stratum was flagged |

### AuditEvent

Timestamped event log for regulatory traceability.

| Field | Type | Description |
|-------|------|-------------|
| `id` | Integer (PK) | Auto-increment primary key |
| `run_id` | Integer (FK) | Associated run (optional) |
| `baseline_id` | Integer (FK) | Associated baseline (optional) |
| `event_type` | String(50) | Event category (e.g., "ingestion", "baseline_created", "verification") |
| `actor` | String(255) | User or system identifier |
| `detail` | Text | Event details |
| `created_at` | DateTime | Timestamp (UTC) |

Indexed on `run_id`, `baseline_id`, `event_type`, and `created_at` for efficient querying.

## Relationships

```
Run ──< Variant
Run ──< Metric
Run ──< QCMetric
Run ──< SoftwareVersion
Run ──< LLMAnalysis
Run ──< GeneSet
Run ──< StratifiedMetric
Run ──< Baseline ──< BaselineEnvelope
                  ──< FixtureRecord
Run ──< VerificationResult ──< VariantTransition
Run ──< RootCauseEvidence
Run ──< BayesianRiskAssessment ──< StratumPosterior
Run ──< AuditEvent
```

All child models cascade delete when their parent Run (or intermediate parent) is deleted.
