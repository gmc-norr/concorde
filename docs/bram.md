# Bayesian Risk Assessment Module (BRAM)

BRAM provides uncertainty-aware degradation detection using a Beta-Binomial conjugate model. It augments (but never overrides) the deterministic acceptance decision, quantifying the posterior probability that a true performance degradation exceeds a clinically meaningful threshold.

## Why BRAM?

### The Limitations of Fixed Thresholds

Deterministic acceptance criteria answer a binary question: does the observed metric exceed a threshold? This works well when sample sizes are large and margins are clear, but it fails in exactly the scenarios where clinical validation is most consequential.

**Problem 1: Small strata produce noisy estimates.**

Consider a gene panel with 30 truth variants. A pipeline correctly calls 29 of 30 (96.7% sensitivity). Does this violate a 99% sensitivity requirement?

A deterministic check says yes -- 96.7% < 99%, therefore FAIL. But with N=30, the Wilson 95% confidence interval for true sensitivity extends from 83.3% to 99.4%. The observed miss may be sampling noise. BRAM computes the posterior probability that the true sensitivity has degraded meaningfully, accounting for the sample size. In this scenario, BRAM's posterior mean is 0.976, yielding a 45.8% tail probability for clinically meaningful degradation -- not flagged.

**Problem 2: Metrics above the threshold can still represent concerning regression.**

Suppose the acceptance threshold for INDEL sensitivity is 98.0%, and the current run shows 98.8% (baseline: 99.9%). Deterministic check: PASS. But BRAM detects that the pipeline has experienced a -1.1 percentage point degradation from baseline. With a large sample (N=5,000) and a Tier 1 threshold of 0.5%, BRAM's posterior is tightly concentrated around 0.988, and the probability that the true sensitivity is below 99.4% (= 99.9% - 0.5%) is effectively 100%. BRAM escalates to REVIEW_REQUIRED.

This is BRAM's primary value: it detects regression *from a validated baseline* that fixed thresholds miss because the metric is still nominally above the floor.

**Problem 3: Risk is not uniform across strata.**

A 0.5% sensitivity decline in ACMG SF v3.1 genes -- where missed variants directly affect clinical decisions -- carries different consequences than the same decline in intergenic MODIFIER variants. BRAM applies per-tier degradation thresholds (0.5% for Critical, 1% for Standard, 2% for Informational) so that critical panels are held to tighter standards.

**Problem 4: Prior runs provide information that fixed thresholds ignore.**

If eight previous verification runs show deltas of [-0.001, +0.002, -0.003, +0.001, -0.002, +0.003, -0.001, +0.002], the normal run-to-run variability is approximately +/-0.002. A new run with delta = -0.008 is unremarkable given this history. BRAM can incorporate this historical distribution as an empirical prior (via method-of-moments Beta fit), producing tighter posteriors and fewer false alarms.

### When BRAM Is Most Useful

| Scenario | Why BRAM Helps |
|----------|---------------|
| Small truth sets (N < 200) | Prior regularization prevents false alarms from noisy estimates |
| Continuous verification | Detects regression from baseline even when metrics remain above fixed thresholds |
| Tiered clinical panels | Per-tier degradation thresholds ensure critical panels receive proportionally more scrutiny |
| Mature pipelines (5+ runs) | Historical priors provide calibrated baselines for normal variability |
| Multi-stratum evaluation | Aggregation captures the global risk picture across dozens of strata |

### When BRAM Is Less Necessary

- **Initial validation** with no baseline: BRAM requires a baseline to compute deltas. On the first run, only deterministic thresholds apply.
- **Very large strata** (N > 10,000) with wide margin: When observed sensitivity is 99.95% against a 99% threshold, the posterior probability of meaningful degradation is effectively zero regardless of model.
- **Single-stratum evaluation**: If only global metrics are assessed (no stratification), the overhead of Bayesian computation adds little over a simple threshold check.

## Statistical Model

### Beta-Binomial Conjugate Model

BRAM uses a Beta-Binomial conjugate model. Sensitivity and precision are binomial proportions (counts of successes over trials), and the Beta distribution is the natural conjugate prior for binomial likelihoods. The conjugate posterior update is exact (no normal approximation is needed). The overall model does involve assumptions: the baseline is treated as a fixed reference (equivalent to assuming infinite baseline sample size — justified because baseline validations typically use large N, making their variance negligible relative to verification runs), the prior is subjective, and the historical Beta fit uses a method-of-moments approximation.

### Per-Stratum Computation

For each (stratum, metric) pair:

**1. Observation counts**

For sensitivity: successes = TP, trials = TP + FN. For precision: successes = TP, trials = TP + FP.

Raw counts are used when available. When only proportions are provided, counts are reconstructed as `successes = round(proportion * sample_size)`, which can introduce off-by-one rounding error.

**2. Beta prior centered at baseline**

The prior encodes "we expect the true proportion to be near the baseline, with confidence proportional to the concentration parameter kappa":

alpha_prior = baseline_value * kappa
beta_prior = (1 - baseline_value) * kappa

The prior mean equals the baseline value. Kappa controls how strongly the prior believes this -- kappa=20 means "20 pseudo-observations worth of confidence." The baseline is treated as a fixed reference point; its sampling uncertainty is not modeled, under the assumption that baselines are established from initial validation with sufficiently large sample sizes.

**3. Prior selection** (by priority):

1. Per-stratum override (custom `prior_concentration` from `bram_config.yaml`)
2. Empirical prior from historical data when `use_historical_priors: true` and at least `min_historical_runs` runs exist. Historical deltas are converted to proportions (`baseline + delta`), then fitted to a Beta distribution using the method of moments: `kappa = mu*(1-mu)/var - 1` (clamped to >= 2). Variance uses `ddof=1` (Bessel's correction), consistent with QC trending and SPC modules.
3. Default prior with `concentration: 20` when neither override nor sufficient history exists.

**4. Conjugate posterior update**

alpha_post = alpha_prior + successes
beta_post = beta_prior + (trials - successes)

This is the exact conjugate update for the Beta-Binomial model.

**5. Tail probability**

P(p_true < baseline - threshold | data) = Beta_CDF(baseline - threshold; alpha_post, beta_post)

Computed via `scipy.stats.beta.cdf`. This is the posterior probability that the true metric value has degraded by more than the threshold from the baseline.

**6. Per-tier degradation thresholds**

Each tier has its own degradation threshold:

| Tier | Threshold | Interpretation |
|------|-----------|----------------|
| Tier 1 (Critical) | 0.5% (0.005) | Critical panels: any clinically significant regression |
| Tier 2 (Standard) | 1.0% (0.01) | Standard coding variants: moderate regression |
| Tier 3 (Informational) | 2.0% (0.02) | Informational: only substantial degradation |

The tail probability is used directly as the risk score (no multiplicative weighting). A stratum is flagged when its tail probability exceeds the `alert_threshold` (default: 0.80).

**7. Aggregation**

Per-stratum tail probabilities are combined into an aggregate score:

| Method | Formula | Implementation |
|--------|---------|----------------|
| `max` | max(tail_probs) | `numpy.ndarray.max()` |
| `weighted_mean` | mean(tail_probs) | `numpy.mean()` |
| `product` | 1 - prod(1 - prob_i) | `numpy.prod()` for numerical stability |

The BRAM verdict is BRAM_FLAG if any stratum is flagged, BRAM_PASS otherwise.

## Worked Examples

### Example 1: BRAM Catches Regression That Deterministic Thresholds Miss

A lab validates an INDEL calling pipeline. The deterministic acceptance threshold for INDELs in the ACMG59 panel is sensitivity >= 98.0%.

| | Baseline | Verification |
|---|---------|-------------|
| Sensitivity | 99.9% | 98.8% |
| N (TP + FN) | 5,000 | 5,000 |

Deterministic check: 98.8% >= 98.0% -> **PASS**.

BRAM computation (Tier 1, threshold = 0.005):

| Step | Value |
|------|-------|
| delta_obs | -0.011 |
| successes | 4,940 / 5,000 |
| Prior | Beta(alpha=19.98, beta=0.02), mean=0.999, kappa=20 |
| Posterior | Beta(alpha=4959.98, beta=60.02) |
| Posterior mean | 0.988 |
| Posterior std | 0.00153 |
| Cutoff | 0.999 - 0.005 = 0.994 |
| P(p_true < 0.994) | **~100%** |
| Flagged? | **Yes** (>0.80) |

Final verdict: Deterministic PASS + BRAM_FLAG = **REVIEW_REQUIRED**.

Interpretation: Although sensitivity remains above the absolute threshold, the 1.1 percentage point regression from baseline is overwhelming. With N=5,000 the posterior is tightly concentrated around 0.988, far below the 0.994 cutoff. BRAM escalates so a clinical scientist can review whether the degradation is explained by a known pipeline change.

### Example 2: BRAM Absorbs Noise in Small Strata

A rare disease panel has only 30 truth variants in the HIGH-impact stratum. The pipeline misses 1 variant.

| | Baseline | Verification |
|---|---------|-------------|
| Sensitivity | 99.0% | 96.7% (29/30) |
| N | 30 | 30 |

Deterministic check: 96.7% < 99.0% -> **FAIL**.

BRAM computation (Tier 2, threshold = 0.01):

| Step | Value |
|------|-------|
| delta_obs | -0.023 |
| successes | 29 / 30 |
| Prior | Beta(alpha=19.80, beta=0.20), mean=0.990, kappa=20 |
| Posterior | Beta(alpha=48.80, beta=1.20) |
| Posterior mean | 0.976 |
| Posterior std | 0.02143 |
| Cutoff | 0.990 - 0.010 = 0.980 |
| P(p_true < 0.980) | 0.458 (45.8%) |
| Flagged? | **No** (0.458 < 0.80) |

Interpretation: With only 30 observations, the prior (kappa=20 pseudo-observations) contributes meaningfully to the posterior. The posterior is wide (std=0.021), reflecting genuine uncertainty. BRAM assigns a 45.8% tail probability that the true degradation exceeds 1% -- concerning but below the flag threshold. The deterministic FAIL still stands (BRAM never overrides FAIL), but the moderate tail probability provides context that the failure may be partly attributable to sampling noise.

### Example 3: Sample Size and the Direction of Evidence

BRAM's behavior with increasing sample size depends critically on whether the observed delta falls within or beyond the degradation threshold. This example illustrates both regimes using the default prior (kappa=20) and Tier 2 threshold (0.01).

**Panel A: Observed delta within threshold (delta = -0.008)**

The observed degradation is 0.8 percentage points -- less than the 1% Tier 2 threshold.

| N | Posterior Mean | Posterior Std | Tail Probability | Flagged? |
|---|---------------|---------------|-----------------|----------|
| 200 | 0.9908 | 0.00642 | 0.310 | No |
| 500 | 0.9923 | 0.00384 | 0.180 | No |
| 1,000 | 0.9912 | 0.00293 | 0.214 | No |
| 2,000 | 0.9911 | 0.00209 | 0.158 | No |
| 5,000 | 0.9910 | 0.00133 | 0.070 | No |

As N increases, the posterior mean converges toward the observed 0.991, which is *above* the cutoff of 0.989 (= 0.999 - 0.01). At large N the posterior tightens around 0.991, and the mass below 0.989 shrinks toward zero. More evidence *confirms* that the degradation is within acceptable limits.

**Panel B: Observed delta beyond threshold (delta = -0.020)**

The observed degradation is 2.0 percentage points -- well beyond the 1% threshold.

| N | Posterior Mean | Posterior Std | Tail Probability | Flagged? |
|---|---------------|---------------|-----------------|----------|
| 200 | 0.9817 | 0.00901 | 0.781 | No |
| 500 | 0.9807 | 0.00602 | 0.937 | **Yes** |
| 1,000 | 0.9794 | 0.00445 | 0.995 | **Yes** |
| 2,000 | 0.9792 | 0.00318 | ~1.000 | **Yes** |
| 5,000 | 0.9791 | 0.00202 | ~1.000 | **Yes** |

The posterior converges toward 0.979, which is well below the cutoff of 0.989. As N increases, the posterior tightens and the tail probability grows monotonically, reaching near-certainty at N=1,000. More evidence *strengthens* the flag.

The contrast between these two panels illustrates BRAM's key property: for observations within the threshold, evidence resolves toward "safe"; for observations beyond the threshold, evidence accumulates toward certainty of degradation.

### Example 4: Per-Tier Degradation Thresholds

Using a per-stratum override to set `degradation_threshold: 0.005` for the ACMG59 panel:

| | Baseline | Verification |
|---|---------|-------------|
| Sensitivity | 99.9% | 99.2% |
| N | 2,000 | 2,000 |

BRAM computation (delta = -0.007, posterior mean = 0.992, posterior std = 0.00197):

With default Tier 2 threshold (0.01): P(p_true < 0.989) = 7.1%. **Not flagged.**

With tighter Tier 1 threshold (0.005): P(p_true < 0.994) = 83.8%. **Flagged.**

The tighter threshold catches the 0.7 percentage point regression that the default threshold would dismiss.

```yaml
per_stratum_overrides:
  ACMG59:
    degradation_threshold: 0.005   # More sensitive detection for critical panel
```

### Example 5: Historical Priors Reduce False Alarms

A mature pipeline with 8 prior verification runs shows historical deltas of:
[-0.001, +0.002, -0.003, +0.001, -0.002, +0.003, -0.001, +0.002]

Method-of-moments Beta fit from these proportions yields kappa=537 (high confidence from stable history).

For a new run with delta = -0.008 at N=1,000 (p=0.991), Tier 2 threshold = 0.01:

| Prior Source | Posterior Mean | Posterior Std | Tail Prob |
|-------------|---------------|---------------|-----------|
| Default (kappa=20) | 0.9912 | 0.00293 | 0.214 |
| Empirical (kappa=537) | 0.9937 | 0.00202 | 0.022 |

The empirical prior, calibrated from the pipeline's demonstrated stability, reduces the tail probability from 21.4% to 2.2%. The high-concentration empirical prior (kappa=537, equivalent to 537 pseudo-observations) pulls the posterior mean toward the historical average, substantially reducing the risk score. This reflects the statistical principle that a pipeline with low historical variability should not be flagged for a delta that is within its established noise floor.

However, this behavior also means that historical priors can mask genuine regressions if the pipeline's stability changes. For this reason, BRAM's historical prior should be combined with periodic manual review, particularly after pipeline software updates. The `min_historical_runs` parameter (default: 5) prevents the empirical prior from activating before sufficient calibration data exists.

## Two-Stage Decision

BRAM never overrides a deterministic FAIL. It can only escalate a PASS or CONDITIONAL_PASS to REVIEW_REQUIRED.

| Deterministic | BRAM | Final Verdict |
|--------------|------|---------------|
| PASS | BRAM_PASS | PASS |
| PASS | BRAM_FLAG | REVIEW_REQUIRED |
| CONDITIONAL_PASS | BRAM_PASS | CONDITIONAL_PASS |
| CONDITIONAL_PASS | BRAM_FLAG | REVIEW_REQUIRED |
| REVIEW_REQUIRED | any | REVIEW_REQUIRED |
| FAIL | any | FAIL |

Key property: BRAM is asymmetric by design. It can raise concerns (escalate PASS to REVIEW_REQUIRED) but cannot dismiss them (a FAIL remains FAIL regardless of BRAM's assessment). This ensures that the deterministic safety floor is always preserved.

## Sample Size Behavior

BRAM's posterior adapts to the available evidence:

- **N < 100**: The prior dominates. BRAM acts as a regularizer, preventing false alarms from noisy small-sample estimates.
- **N = 100-500**: Prior and observation contribute roughly equally. Risk quantification becomes meaningful, but marginal degradations remain ambiguous.
- **N > 500**: The observation dominates. BRAM closely tracks the observed proportion with calibrated uncertainty bounds.
- **N > 5,000**: The posterior is tightly concentrated around the observed value. BRAM's assessment converges to near-certainty for degradations beyond the threshold.

For reliable flagging, aim for N >= 500 variants per stratum. Historical priors from 5+ verification runs provide the most stable risk estimates.

## Configuration

Key parameters in `config/bram_config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alert_threshold` | 0.80 | Flag strata with tail probability above this value |
| `aggregation_method` | `max` | Aggregation: `max`, `weighted_mean`, or `product` |
| `default_prior.concentration` | 20 | Prior effective sample size (kappa) |
| `use_historical_priors` | true | Use empirical priors from historical data |
| `min_historical_runs` | 5 | Minimum runs for empirical prior |
| `degradation_thresholds.tier_1` | 0.005 | Clinically meaningful degradation for critical panels |
| `degradation_thresholds.tier_2` | 0.01 | Clinically meaningful degradation for standard strata |
| `degradation_thresholds.tier_3` | 0.02 | Clinically meaningful degradation for informational strata |

### Tuning Guidance

| Parameter | Effect of Increasing | When to Adjust |
|-----------|---------------------|----------------|
| `alert_threshold` | Fewer flags, more permissive | Too many REVIEW_REQUIRED verdicts |
| `degradation_thresholds.*` | Fewer flags for small deltas | Only flag large regressions |
| `default_prior.concentration` | Stronger prior, slower convergence to data | Pipeline variability is well-characterized |
| Lower `degradation_thresholds.tier_1` | More sensitive flagging for critical strata | Critical panels need tighter monitoring |

### Per-Stratum Overrides

Override prior concentration or thresholds for specific strata:

```yaml
per_stratum_overrides:
  ACMG59:
    prior_concentration: 50       # Stronger prior for critical panels
    degradation_threshold: 0.003  # More sensitive detection (0.3%)
  INDEL:
    degradation_threshold: 0.02   # More tolerant for INDELs (2%)
```

### Aggregation Methods

| Method | Formula | Behavior | Recommended When |
|--------|---------|----------|-----------------|
| `max` | max(tail_probs) | Flags if *any* single stratum is concerning | Default. Conservative, catches single-stratum regressions. |
| `weighted_mean` | mean(tail_probs) | Averages risk across all strata | Many strata with uniform importance. Dilutes single-stratum signals. |
| `product` | 1 - prod(1 - prob_i) | Joint probability of at least one failure | Many independent strata. More sensitive than `max` when multiple strata show moderate risk simultaneously. |

## Report Integration

BRAM results appear in validation reports as:

- A **verdict banner** (BRAM_PASS / BRAM_FLAG) with aggregate and mean risk scores
- **Configuration summary**: alert threshold, degradation threshold, aggregation method
- **Per-stratum table**: dimension, stratum, metric name, delta observed, posterior mean/std, tail probability, degradation threshold, risk score, and flagged status

Flagged strata are highlighted in the HTML report for immediate visual identification.
