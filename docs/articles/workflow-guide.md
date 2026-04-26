# Bayesian meta-analysis workflow with bayesma

A comprehensive guide to conducting Bayesian meta-analyses using the
**bayesma** package, from data preparation through interpretation.

------------------------------------------------------------------------

## Overview

This workflow follows best practices for transparent, reproducible
meta-analysis. The stages are:

1.  **Data Preparation** — Clean and structure the data
2.  **Trustworthiness Assessment** (INSPECT-SR) — Evaluate study quality
    systematically
3.  **Preliminary Exploration** (PRIMED) — Understand data structure
    before modeling
4.  **Risk of Bias Assessment** (RoB) — Evaluate bias domains
5.  **Model Specification** — Define prior distributions and model
    structure
6.  **Model Fitting** — Fit Bayesian models (standard, bias-corrected,
    selection, robust)
7.  **Model Comparison** — Compare competing models via LOSO-CV or
    LOO-IC
8.  **Heterogeneity Assessment** — Evaluate between-study variation
9.  **Bias Assessment** — Detect and adjust for publication/reporting
    bias
10. **Subgroup & Moderation Analysis** — Explore effect modification
    (meta-regression)
11. **Sensitivity Analysis** — Assess robustness to modeling choices
12. **Interpretation & Reporting** — Summarize findings and communicate
    uncertainty

------------------------------------------------------------------------

## Stage 1: Data Preparation

### Input Structure

Prepare a data frame with **one row per study** (or per arm for
multi-arm studies). The exact columns depend on the likelihood:

#### For Binary Outcomes (Binomial)

``` r
data <- tibble::tibble(
  study = c("Author 2020", "Author 2021"),           # Study identifier
  event_int = c(10, 15),                              # Events in intervention
  event_ctrl = c(8, 12),                              # Events in control
  n_int = c(100, 120),                                # Sample size, intervention
  n_ctrl = c(100, 110)                                # Sample size, control
)
```

#### For Continuous Outcomes (Gaussian)

``` r
data <- tibble::tibble(
  study = c("Author 2020", "Author 2021"),
  mean_int = c(5.2, 6.1),                             # Intervention mean
  mean_ctrl = c(3.8, 4.5),                            # Control mean
  sd_int = c(1.5, 1.8),                               # Intervention SD
  sd_ctrl = c(1.4, 1.6),
  n_int = c(100, 120),
  n_ctrl = c(100, 110)
)
```

#### For Count Outcomes (Poisson)

``` r
data <- tibble::tibble(
  study = c("Author 2020", "Author 2021"),
  event_int = c(25, 30),                              # Count, intervention
  event_ctrl = c(15, 18),                             # Count, control
  n_int = c(100, 120),
  n_ctrl = c(100, 110)
)
```

#### For Multi-Arm Studies (One-Stage Only)

Add a `multi_arm` column to group arms within a study:

``` r
data <- tibble::tibble(
  study = rep(c("Author 2020", "Author 2020"), each = 2),  # Repeated study name
  arm = c("Arm1", "Arm2", "Arm1", "Arm2"),
  # ... outcome columns ...
  multi_arm = study                                   # Column grouping arms
)
```

### Data Quality Checks

Run basic validation:

``` r
# Check required columns
bayesma::check_statistics_consistency(data, likelihood = "binomial")
bayesma::check_n_consistency(data)

# Verify unique studies
stopifnot(!duplicated(data$study))
```

------------------------------------------------------------------------

## Stage 2: Trustworthiness Assessment (INSPECT-SR)

**When to use:** Assessing which studies are trustworthy enough to
include in analysis.

**What it does:** Systematically evaluates RCTs across 4 domains using
21 items, including automated statistical checks (Carlisle, GRIM, N
consistency, p-value verification).

### Step 1: Prepare Study Metadata

Create a data frame with study-level information and trustworthiness
judgements:

``` r
sr_data <- tibble::tibble(
  study = c("Smith 2020", "Jones 2021"),
  # D1: Post-publication notices (3 items)
  d1_1 = c("No concerns", "No concerns"),              # Retraction/EOC check
  d1_2 = c("No concerns", "Some concerns"),            # Withdrawal check
  d1_3 = c("No concerns", "No concerns"),              # Corrigendum check
  
  # D2: Conduct, governance, transparency (5 items)
  d2_1 = c("No concerns", "Some concerns"),            # Trial registration
  d2_2 = c("No concerns", "No concerns"),              # Ethical approval
  # ... d2_3, d2_4, d2_5 ...
  
  # D3: Text & publication (2 items)
  d3_1 = c("No concerns", "No concerns"),
  d3_2 = c("No concerns", "No concerns"),
  
  # D4: Results & statistics (11 items - some automated)
  d4_1 = c("No concerns", "No concerns"),              # Outcomes clearly stated
  # ... d4_2, d4_3 (automated), d4_4, ... d4_11 ...
  
  # List-columns for automated checks
  baseline = list(
    tibble::tibble(variable = "age", group = "int", mean = 45.2, sd = 12.1),
    tibble::tibble(variable = "age", group = "int", mean = 46.1, sd = 11.8)
  ),
  
  statistics = list(
    tibble::tibble(test_name = "t-test", test_statistic = 2.34, df = 198),
    tibble::tibble(test_name = "t-test", test_statistic = 2.10, df = 228)
  )
)
```

### Step 2: Run INSPECT-SR Assessment

``` r
inspect_results <- bayesma::inspect_sr(sr_data, studyvar = study)

# View summary
print(inspect_results)

# View per-check details and failures
bayesma::inspect_sr_table(inspect_results, only_failed = TRUE)
```

### Step 3: Visualize Trustworthiness

``` r
# Plot trustworthiness across studies
bayesma::inspect_plot(inspect_results)

# Summary across all studies
bayesma::inspect_summary_plot(inspect_results)
```

### Step 4: Filter for Analysis

``` r
# Option 1: Exclude serious concerns
meta_data <- bayesma::filter_trustworthy(
  data, 
  inspect_results,
  threshold = "some_concerns"  # Allow "No concerns" and "Some concerns"
)

# Option 2: Keep all, flag in sensitivity analysis
# ... proceed with all data, later check if results change without flagged studies
```

------------------------------------------------------------------------

## Stage 3: Preliminary Exploration (PRIMED)

**When to use:** Before fitting any model, explore your data structure
systematically.

**What it does:** Implements the PRIMED workflow (Pustejovsky, Zhang, &
Tipton) — describes dependence structure, explores moderators, inspects
auxiliary data, visualizes effect size distribution (last, to avoid
outcome-driven decisions).

### Full Workflow in One Call

``` r
primed_results <- bayesma::primed(
  data = meta_data,
  studyvar = study,
  moderators = c("intervention_type", "study_design"),  # Optional: variables to explore
  es_col = "yi",                                        # If pre-calculated effects
  se_col = "sei"                                        # If pre-calculated SEs
)

print(primed_results)
```

### Key Outputs

The [`primed()`](https://blmoran.github.io/bayesma/reference/primed.md)
function produces:

**1. Data Summary** — Number of studies, effect sizes, studies with
multiple effects, missing data

``` r
# Access via:
primed_results$summary
```

**2. Dependence Structure** — How effect sizes nest within studies

``` r
# Access via:
primed_results$dependence_table  # E.g., "52 effects from 18 studies"
```

**3. Study Characteristics** — Moderator frequencies and distributions

``` r
primed_results$moderator_tables
```

**4. Auxiliary Data Quality** — Standard errors, outcome ranges,
completeness

``` r
primed_results$auxiliary_inspection
```

**5. Effect Size Distribution** — Visualizations (histogram, Q-Q,
summary)

``` r
primed_results$es_distribution_plot
```

### Manual Exploration (if needed)

If you prefer granular control:

``` r
# Dependence structure
meta_data |>
  dplyr::group_by(study) |>
  dplyr::summarise(n_es = dplyr::n(), .groups = "drop") |>
  dplyr::count(n_es) |>
  dplyr::rename(n_effect_sizes_per_study = n_es, n_studies = n)

# Moderator frequencies
meta_data |>
  dplyr::count(intervention_type)

# SE ranges
meta_data |>
  dplyr::summarise(
    se_min = min(sei, na.rm = TRUE),
    se_mean = mean(sei, na.rm = TRUE),
    se_max = max(sei, na.rm = TRUE)
  )
```

------------------------------------------------------------------------

## Stage 4: Risk of Bias Assessment (RoB)

**When to use:** Evaluating methodological quality of included studies.

**What it does:** Visualizes risk of bias judgements across studies and
bias domains.

### Prepare RoB Data

``` r
rob_data <- tibble::tibble(
  study = meta_data$study,
  selection_bias = c("Low", "Unclear", "High"),         # Domain judgements
  performance_bias = c("Low", "Low", "Unclear"),
  detection_bias = c("Low", "Unclear", "High"),
  attrition_bias = c("Low", "Low", "Low"),
  reporting_bias = c("Low", "Low", "Unclear")
)
```

### Visualize RoB

``` r
# RoB traffic light plot
bayesma::rob_plot(rob_data, studyvar = study)

# RoB summary table
bayesma::rob_table(rob_data)
```

### Notes

- INSPECT-SR evaluates *trustworthiness* (is the data real?)
- RoB evaluates *bias risk* (could the design/conduct bias results?)
- Use both: combine trustworthiness filtering + RoB sensitivity analysis

------------------------------------------------------------------------

## Stage 5: Model Specification

**When to use:** Deciding what model to fit.

**What it does:** Specifies prior distributions, likelihood, random
effect structure, and bias-handling mechanisms.

### Standard Meta-Analysis

Fit a simple random-effects model:

``` r
spec <- bayesma::bayesma_spec(
  data = meta_data,
  studyvar = study,
  event_int = event_int,              # For binomial
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "random_effect",       # or "common_effect"
  stage = "two_stage",                # or "one_stage"
  re_dist = "normal",                 # or "t", "skew_normal", "mixture"
  
  # Priors
  mu_prior = bayesma::normal(0, 1),   # Pooled effect prior
  tau_prior = bayesma::half_normal(0, 0.5)  # Heterogeneity prior
)

# Inspect the specification
print(spec)
```

### With Bias Correction

``` r
spec <- bayesma::bayesma_spec(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "bias_corrected",      # Systematic bias adjustment
  stage = "two_stage",
  
  # Bias prior
  b_prior = bayesma::normal(0, 0.25)  # Prior on bias shift
)
```

### With Publication Bias (Selection Model)

``` r
spec <- bayesma::bayesma_spec(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "selection_weight",    # or "selection_copas"
  stage = "two_stage",
  p_cutoffs = c(0.025, 0.05),        # P-value thresholds for weighting
  
  # Selection priors
  selection_priors = bayesma::prior_weight_function()
)
```

### With Robust Outlier Mixture

``` r
spec <- bayesma::bayesma_spec(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "random_effect",
  re_dist = "mixture",               # Two-component robust mixture
  robust = TRUE,                     # Add outlier component
  n_components = 2,
  robust_weight = bayesma::beta(1, 9)  # Prior on outlier weight
)
```

### For Meta-Regression (Moderation)

``` r
spec <- bayesma::meta_reg_spec(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  moderators = ~ intervention_type + study_design,  # Formula of moderators
  
  # Priors
  mu_prior = bayesma::normal(0, 1),
  tau_prior = bayesma::half_normal(0, 0.5),
  # Beta priors for moderator coefficients (auto-generated)
  rescale_priors = 1  # Scale to outcome SD
)
```

### Robust Model Averaging (RoBMA)

Fit multiple bias-handling models simultaneously:

``` r
spec <- bayesma::robma_spec(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  method = "bridge",                  # or "ss" (spike-slab)
  bias_indicator = c("bias_corrected", "pet_peese", "selection_weight"),
  null_range = c(-0.1, 0.1)          # Define practically equivalent effect
)
```

### Prior Justification

Document your choices:

``` r
# Visualize priors vs. possible effect sizes
bayesma::bayesma_prior_density(
  spec,
  plot_range = c(-2, 2)
)

# Compare to published guidelines
bayesma::build_mu_prior_overlay(likelihood = "binomial")
```

------------------------------------------------------------------------

## Stage 6: Model Fitting

**When to use:** After specification, fit the model to data.

**What it does:** Compiles Stan code, samples from the posterior, and
extracts results.

### Standard Workflow

``` r
# One-call pipeline
fit <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "random_effect",
  stage = "two_stage",
  
  # MCMC settings
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = 1234,
  
  # Optional: inspect intermediate stages
  return_stage = "full"  # or "spec", "code", "data", "fit"
)

# Summary
print(fit)
```

### Modular Pipeline

For more control, use individual stages:

``` r
# Stage 1: Specification
spec <- bayesma::bayesma_spec(...)

# Stage 2: Generate Stan code
code <- bayesma::bayesma_stan_code(spec)
print(code)  # Inspect if needed

# Stage 3: Build Stan data
stan_data <- bayesma::bayesma_stan_data(spec)

# Stage 4: Compile and sample
fit <- bayesma::bayesma_fit(stan_data, code, chains = 4, ...)

# Stage 5: Extract results
effects <- bayesma::bayesma_extract(fit)

# Stage 6: Build output object
output <- bayesma::bayesma_output(fit, spec, effects)
```

### Meta-Regression Fitting

``` r
# One call
mreg_fit <- bayesma::meta_reg(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  moderators = ~ intervention_type + study_design,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Or modular
mreg_spec <- bayesma::meta_reg_spec(...)
mreg_code <- bayesma::meta_reg_stan_code(mreg_spec)
mreg_data <- bayesma::meta_reg_stan_data(mreg_spec)
mreg_fit <- bayesma::meta_reg_fit(mreg_data, mreg_code, chains = 4, ...)
mreg_effects <- bayesma::meta_reg_extract(mreg_fit)
mreg_output <- bayesma::meta_reg_output(mreg_fit, mreg_spec, mreg_effects)
```

### RoBMA Fitting

``` r
robma_fit <- bayesma::robma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  method = "bridge",
  bias_indicator = c("bias_corrected", "pet_peese", "selection_weight"),
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

print(robma_fit)
```

### Egger Regression for Publication Bias

``` r
egger_fit <- bayesma::egger(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "random_effect",
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Plot Egger regression
bayesma::egger_plot(egger_fit)
```

------------------------------------------------------------------------

## Stage 7: Model Comparison

**When to use:** Deciding between competing models.

**What it does:** Compares models using leave-one-study-out
cross-validation (LOSO-CV) for across-stage comparison, or LOO-IC for
within-stage comparison.

### Compare Multiple Models

``` r
fit_re <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  # ... specification for random-effects model
  model_type = "random_effect"
)

fit_ce <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  # ... specification for common-effect model
  model_type = "common_effect"
)

fit_robust <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  # ... specification for robust mixture model
  re_dist = "mixture",
  robust = TRUE
)

# Compare via LOSO-CV
comparison <- bayesma::compare_models(
  "Random-effects" = fit_re,
  "Common-effect" = fit_ce,
  "Robust" = fit_robust,
  data = meta_data,
  studyvar = study,
  criterion = "loso",  # or "loo"
  quiet = FALSE
)

print(comparison)
```

### Visualization

``` r
# LOSO-CRPS (continuous ranked probability score) comparison
bayesma::compare_plot(comparison, type = "loso_crps")

# Calibration plot (are prediction intervals accurate?)
bayesma::compare_plot(comparison, type = "calibration")

# Summary table
bayesma::compare_table(comparison)
```

### Key Metrics

- **LOSO-CRPS**: Lower is better. Averaged proper scoring rule.
- **LOSO-ELPD**: Higher is better. Expected log predictive density.
- **Calibration**: Empirical coverage should ≈ nominal (e.g., 95% CI
  should contain true value 95% of the time).
- **Miscalibration**: Mean \|empirical − nominal\| across coverage
  levels. Zero is perfect.

------------------------------------------------------------------------

## Stage 8: Heterogeneity Assessment

**When to use:** Understanding variation between studies.

**What it does:** Estimates and visualizes between-study variance (τ²)
and intra-class correlation.

### Extract Heterogeneity Estimates

``` r
# From standard meta-analysis
tau_summary <- fit |>
  bayesma::extract_summary() |>
  dplyr::filter(variable == "tau")

print(tau_summary)
# Shows: median, mad, quantiles (2.5%, 97.5%)

# Interpret
cat("Posterior median τ:", tau_summary$median, "\n")
cat("95% CrI for τ:", tau_summary$q2.5, "to", tau_summary$q97.5, "\n")
```

### Visualizations

``` r
# Prior-posterior plot for tau
bayesma::build_tau_prior_overlay(fit)

# Forest plot (includes tau visually)
bayesma::bayes_forest(fit, data = meta_data)
```

### Predictive Distribution

What to expect in a new study:

``` r
# Extract predictive draws for a new study
pred_draws <- fit$draws |>
  posterior::subset_draws(variable = "mu_new") |>
  as.numeric()

# Summarize prediction
cat("Posterior predictive median:", median(pred_draws), "\n")
cat("95% prediction interval:", 
    quantile(pred_draws, c(0.025, 0.975)), "\n")
```

### Heterogeneity Interpretation Guide

| τ estimate | Interpretation |
|----|----|
| ~0 | Minimal heterogeneity; effects nearly constant |
| Small (0-0.2) | Mild heterogeneity; common-effect reasonable as reference |
| Moderate (0.2-0.5) | Meaningful heterogeneity; random-effects preferred |
| Large (\>0.5) | High heterogeneity; consider moderator analysis |

------------------------------------------------------------------------

## Stage 9: Bias Assessment

**When to use:** Evaluating whether publication bias, small-study
effects, or systematic bias affect results.

**What it does:** Fits models with bias-handling mechanisms and compares
results to unbiased model.

### Publication Bias (Selection Models)

``` r
# Fit selection model
fit_selection <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "selection_weight",
  p_cutoffs = c(0.025, 0.05),
  chains = 4
)

# Compare to unbiased model
comparison_bias <- bayesma::compare_models(
  "No bias" = fit_re,
  "Selection bias" = fit_selection,
  data = meta_data,
  studyvar = study
)

print(comparison_bias)
```

### Small-Study Effects (PET-PEESE)

``` r
# PET-PEESE regression
fit_pet_peese <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "pet_peese",
  chains = 4
)

# Plot regression line
bayesma::overall_plot(fit_pet_peese, data = meta_data)
```

### Systematic Study Bias (Bias-Corrected Model)

``` r
# Bias-corrected model: assumes some studies have systematic bias
fit_bias_corrected <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "bias_corrected",
  b_prior = bayesma::normal(0, 0.25),
  chains = 4
)

# Compare
comparison_bias2 <- bayesma::compare_models(
  "No bias" = fit_re,
  "Bias-corrected" = fit_bias_corrected,
  data = meta_data,
  studyvar = study
)
```

### Robust Model Averaging (RoBMA)

Fit multiple bias models simultaneously and average:

``` r
robma_fit <- bayesma::robma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  method = "bridge",
  bias_indicator = c("none", "bias_corrected", "pet_peese", "selection_weight"),
  chains = 4
)

print(robma_fit)
# Shows posterior model probabilities + averaged effect estimate
```

### Egger Regression

Visual test for publication bias:

``` r
fit_egger <- bayesma::egger(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  chains = 4
)

# Plot with credible bands
bayesma::egger_plot(fit_egger, data = meta_data)

# Extract slope (publication bias indicator)
slope_summary <- fit_egger |>
  bayesma::extract_summary() |>
  dplyr::filter(variable == "beta_egger")

print(slope_summary)
```

### Funnel Plot

``` r
# Standard funnel plot
bayesma::funnel_plot(fit_re, data = meta_data)

# With publication bias correction
bayesma::funnel_plot(fit_selection, data = meta_data)
```

------------------------------------------------------------------------

## Stage 10: Subgroup & Moderation Analysis

**When to use:** Exploring effect modification (e.g., does efficacy vary
by age group, intervention type, or study quality?)

**What it does:** Fits meta-regression models with study-level
moderators.

### Model Specification

``` r
spec <- bayesma::meta_reg_spec(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  
  # Moderator formula (like lm/glm)
  moderators = ~ intervention_type + study_duration,
  
  # Priors
  mu_prior = bayesma::normal(0, 1),
  tau_prior = bayesma::half_normal(0, 0.5),
  # Beta priors auto-scaled to outcome variance
  rescale_priors = 1
)

print(spec)
```

### Fitting

``` r
mreg_fit <- bayesma::meta_reg(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  moderators = ~ intervention_type + study_duration,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

print(mreg_fit)
```

### Results

``` r
# Extract moderator coefficients
mreg_effects <- bayesma::extract_summary(mreg_fit) |>
  dplyr::filter(grepl("beta", variable))

print(mreg_effects)

# Visualize: predicted effects across moderator values
bayesma::metareg_mod_plot(
  mreg_fit,
  data = meta_data,
  moderator = "intervention_type"
)

# Or custom prediction plot
bayesma::sensitivity_plot(
  mreg_fit,
  type = "moderator",
  moderator = "study_duration"
)
```

### Continuous Moderators

For a continuous moderator (e.g., study year), plot the regression
surface:

``` r
# Data with new study years to predict
new_years <- tibble::tibble(study_year = seq(2010, 2025, by = 1))

# Predicted effects
pred <- bayesma::sensitivity_plot(
  mreg_fit,
  type = "moderator",
  moderator = "study_year",
  newdata = new_years
)

print(pred)
```

### Categorical Moderators

For categorical moderators (e.g., intervention type), compare subgroups:

``` r
# Summarize effects by subgroup
subgroup_summary <- mreg_fit$summary |>
  dplyr::filter(grepl("beta_intervention", variable)) |>
  dplyr::mutate(
    subgroup = gsub("beta_intervention_type\\[(.*)\\]", "\\1", variable)
  )

print(subgroup_summary)
```

------------------------------------------------------------------------

## Stage 11: Sensitivity Analysis

**When to use:** Assessing robustness to modeling assumptions and data
decisions.

**What it does:** Refits models under different assumptions and
summarizes how results change.

### Prior Sensitivity

Does the posterior depend on the prior?

``` r
# Fit with weak prior
spec_weak <- bayesma::bayesma_spec(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  mu_prior = bayesma::normal(0, 2),        # Wider prior
  tau_prior = bayesma::half_normal(0, 1)   # Wider prior
)
fit_weak <- bayesma::bayesma_fit(..., spec = spec_weak)

# Fit with strong prior
spec_strong <- bayesma::bayesma_spec(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  mu_prior = bayesma::normal(0, 0.5),      # Narrower prior
  tau_prior = bayesma::half_normal(0, 0.25)
)
fit_strong <- bayesma::bayesma_fit(..., spec = spec_strong)

# Compare posteriors
comparison_prior <- bayesma::compare_models(
  "Weak prior" = fit_weak,
  "Default prior" = fit_re,
  "Strong prior" = fit_strong,
  data = meta_data,
  studyvar = study
)

print(comparison_prior)
```

### Studies-Included Sensitivity

Does the conclusion hold when excluding studies?

``` r
# Remove outliers or studies with concerns
meta_subset <- meta_data |>
  dplyr::filter(!study %in% c("Outlier 2020", "Flagged 2021"))

fit_subset <- bayesma::bayesma(
  data = meta_subset,
  studyvar = study,
  # ... same specification as fit_re ...
)

comparison_subset <- bayesma::compare_models(
  "All studies" = fit_re,
  "Excluding outliers" = fit_subset,
  data = meta_data,
  studyvar = study
)
```

### ROBMA Sensitivity

Automated sensitivity across multiple bias assumptions:

``` r
# RoBMA includes implicit sensitivity via model averaging
robma_fit <- bayesma::robma(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  method = "bridge",
  bias_indicator = c("none", "bias_corrected", "pet_peese", "selection_weight")
)

print(robma_fit)
# Posterior model probabilities show sensitivity to bias model

# Visualize effect estimate under each model
bayesma::robma_table(robma_fit)
```

### Likelihood Sensitivity (One-Stage vs Two-Stage)

Does the stage of analysis matter?

``` r
fit_1s <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  stage = "one_stage",            # Full likelihood on raw data
  chains = 4
)

fit_2s <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  stage = "two_stage",            # Summary effect estimates
  chains = 4
)

comparison_stage <- bayesma::compare_models(
  "One-stage" = fit_1s,
  "Two-stage" = fit_2s,
  data = meta_data,
  studyvar = study
)

print(comparison_stage)
```

### Random Effects Distribution Sensitivity

``` r
fit_normal <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  re_dist = "normal"   # Symmetric
)

fit_t <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  re_dist = "t"        # Heavy-tailed, robust to outliers
)

fit_skew <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  re_dist = "skew_normal"  # Asymmetric
)

comparison_re <- bayesma::compare_models(
  "Normal" = fit_normal,
  "t-distribution" = fit_t,
  "Skew-normal" = fit_skew,
  data = meta_data,
  studyvar = study
)

print(comparison_re)
```

### Automated Sensitivity Plots

``` r
# Create sensitivity plot object
sens_plot <- bayesma::sensitivity_plot(
  fit_re,
  data = meta_data,
  # Include multiple sensitivity scenarios
  include_null_range = TRUE,
  include_prior_scales = c(0.5, 1, 2)
)

# Render with custom titles
sens_rendered <- bayesma::render_sensitivity_patchwork(
  sens_plot,
  title = "Sensitivity to Modeling Choices"
)

print(sens_rendered)
```

### For RoBMA: Robust Model Averaging Sensitivity

``` r
# RoBMA automatically averages over multiple bias assumptions
robma_fit <- bayesma::robma(
  data = meta_data,
  studyvar = study,
  # ... outcome columns ...
  method = "bridge",
  bias_indicator = c("none", "bias_corrected", "pet_peese", "selection_weight")
)

# Posterior model probabilities show which assumptions drive conclusions
robma_summary <- print(robma_fit)
# Inspect: "Model posteriors" section

# Forest plot showing weighted effect (averaged over models)
bayesma::bayes_forest(robma_fit, data = meta_data)
```

### ROBMA Sensitivity Plot

``` r
# Visualize robustness: how does effect estimate change with assumption?
robma_sensitivity <- bayesma::robma_sensitivity(
  robma_fit,
  data = meta_data,
  studyvar = study
)

# Default is to show effect estimate under each bias assumption
print(robma_sensitivity)

# Plot
bayesma::render_sensitivity_patchwork(robma_sensitivity)
```

------------------------------------------------------------------------

## Stage 12: Interpretation & Reporting

**When to use:** Summarizing findings for publication or stakeholders.

**What it does:** Extracts and formats results, confidence intervals,
and evidential summaries.

### Summary Statistics

``` r
# Get all summary statistics
summary_df <- fit_re |>
  bayesma::extract_summary()

# Pooled effect estimate
mu_summary <- summary_df |>
  dplyr::filter(variable == "mu")

cat("Posterior estimate of pooled effect:\n")
cat("  Median (95% CrI):", 
    round(mu_summary$median, 3), 
    "(", round(mu_summary$q2.5, 3), " to ", round(mu_summary$q97.5, 3), ")\n")

# Heterogeneity (τ)
tau_summary <- summary_df |>
  dplyr::filter(variable == "tau")

cat("\nBetween-study standard deviation:\n")
cat("  Median (95% CrI):", 
    round(tau_summary$median, 3), 
    "(", round(tau_summary$q2.5, 3), " to ", round(tau_summary$q97.5, 3), ")\n")
```

### Probability of Effect Direction & Magnitude

``` r
# P(θ > 0) — probability effect is positive
theta_draws <- fit_re$draws |>
  posterior::subset_draws(variable = "mu") |>
  as.numeric()

prob_positive <- mean(theta_draws > 0)
cat("P(effect > 0):", round(prob_positive, 3), "\n")

# P(effect in null range) — probability of negligible effect
null_range <- c(-0.1, 0.1)
prob_null <- mean(theta_draws >= null_range[1] & theta_draws <= null_range[2])
cat("P(effect in [−0.1, 0.1]):", round(prob_null, 3), "\n")
```

### Bayes Factor & Coefficient of Evidence

``` r
# Bayes factor: null (θ = 0) vs. alternative
bf_null <- bayesma::coefficient_evidence(fit_re)
print(bf_null)
# Shows BF10 and BF01

# Interpretation
if (bf_null$bf10 > 10) {
  cat("Strong evidence for non-zero effect\n")
} else if (bf_null$bf10 < 0.1) {
  cat("Strong evidence for zero effect\n")
} else {
  cat("Inconclusive evidence\n")
}
```

### Model Weights (RoBMA)

``` r
# Posterior model probabilities
robma_summary <- print(robma_fit)

# Shows: "No bias", "Bias-corrected", "PET-PEESE", "Selection model"
# E.g., 60% weight to "Selection model" means that model best describes data

# Weighted effect estimate across models
effect_weighted <- robma_fit |>
  bayesma::extract_summary() |>
  dplyr::filter(variable == "mu")

cat("Model-averaged effect estimate:\n")
cat("  Median (95% CrI):", 
    round(effect_weighted$median, 3), 
    "(", round(effect_weighted$q2.5, 3), " to ", round(effect_weighted$q97.5, 3), ")\n")
```

### Publication-Ready Tables

``` r
# Forest plot table
forest_table <- bayesma::bayes_forest(fit_re, data = meta_data)
print(forest_table)

# Model comparison table
comparison_table <- bayesma::compare_table(comparison)
print(comparison_table)

# RoBMA table (per-model results)
robma_table <- bayesma::robma_table(robma_fit)
print(robma_table)

# Sensitivity analysis table (if using ROBMA)
robma_sens <- bayesma::robma_sensitivity(robma_fit, data = meta_data, studyvar = study)
sensitivity_table <- bayesma::robma_table(robma_sens)
print(sensitivity_table)
```

### Diagnostic Plots

**MCMC Diagnostics:**

``` r
# Trace plots (visual inspection for convergence)
bayesplot::mcmc_trace(fit_re$fit$draws(), variables = c("mu", "tau"))

# Rhat (potential scale reduction factor; should be < 1.01)
fit_re$fit$summary(variables = c("mu", "tau")) |>
  dplyr::select(variable, rhat)

# Effective sample size
fit_re$fit$summary(variables = c("mu", "tau")) |>
  dplyr::select(variable, ess_bulk, ess_tail)
```

**Posterior Predictive Checks:**

``` r
# Do the posterior predictions match the observed data?
bayesma::pp_check(fit_re, type = "dens_overlay")
bayesma::pp_check(fit_re, type = "stat", stat = "median")
```

**Bias Diagnostics:**

``` r
# Pareto k values (influential observations)
bayesma::diagnostics(fit_re)

# If high k values, those studies are influential
```

### Publication Checklist

Ensure you report:

**Study Selection** — Search strategy, inclusion/exclusion, PRISMA
diagram

**Study Characteristics** — Table of included studies (N, design,
outcome, baseline characteristics)

**Risk of Bias** — RoB plot + summary, INSPECT-SR results if used

**Preliminary Exploration** — PRIMED tables + visualizations

**Main Result** — Pooled effect (point estimate + 95% CrI), τ
(heterogeneity)

**Model Justification** — Why random-effects? Why this likelihood? Prior
specification.

**Model Comparison** — LOSO-CV or LOO-IC results comparing main models

**Bias Assessment** — Selection models, PET-PEESE, egger plots,
robustness to model choice

**Subgroup Analysis** — Meta-regression results (if applicable)

**Sensitivity Analysis** — Prior sensitivity, studies-excluded
sensitivity, random-effects distribution sensitivity

**Interpretation** — Effect magnitude, uncertainty, implications

**Reproducibility** — R code, seed, data availability statement

------------------------------------------------------------------------

## Example Workflow End-to-End

Below is a condensed example for a binary outcome:

``` r
# Load data
data(dat.smith2021, package = "bayesma")
meta_data <- dat.smith2021

# Stage 1: Trustworthiness check (abbreviated)
sr_results <- bayesma::inspect_sr(meta_data, studyvar = study)
meta_data <- bayesma::filter_trustworthy(meta_data, sr_results)

# Stage 2: Preliminary exploration
primed_results <- bayesma::primed(meta_data, studyvar = study)
print(primed_results)

# Stage 3-4: RoB (external; use your own rob_data)
# bayesma::rob_plot(rob_data)

# Stage 5-6: Fit main model
fit_re <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "random_effect",
  stage = "two_stage",
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

print(fit_re)

# Fit bias-adjusted model
fit_selection <- bayesma::bayesma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  model_type = "selection_weight",
  stage = "two_stage",
  chains = 4
)

# Stage 7: Compare models
comparison <- bayesma::compare_models(
  "Random-effects" = fit_re,
  "Selection model" = fit_selection,
  data = meta_data,
  studyvar = study
)

print(comparison)

# Stage 8: Heterogeneity
tau <- bayesma::extract_summary(fit_re) |>
  dplyr::filter(variable == "tau")
print(tau)

# Stage 9: Bias assessment
bayesma::egger_plot(fit_selection, data = meta_data)
bayesma::funnel_plot(fit_re, data = meta_data)

# Stage 10: Meta-regression (if moderators available)
mreg_fit <- bayesma::meta_reg(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  moderators = ~ study_year + study_design,
  chains = 4
)

# Stage 11: Sensitivity via ROBMA
robma_fit <- bayesma::robma(
  data = meta_data,
  studyvar = study,
  event_int = event_int,
  event_ctrl = event_ctrl,
  n_int = n_int,
  n_ctrl = n_ctrl,
  likelihood = "binomial",
  method = "bridge",
  bias_indicator = c("none", "bias_corrected", "pet_peese", "selection_weight"),
  chains = 4
)

print(robma_fit)

# Stage 12: Results & reporting
summary_main <- bayesma::extract_summary(fit_re) |>
  dplyr::filter(variable %in% c("mu", "tau"))

print(summary_main)

# Forest plot
bayesma::bayes_forest(fit_re, data = meta_data)

# Model-averaged forest plot
bayesma::bayes_forest(robma_fit, data = meta_data)
```

------------------------------------------------------------------------

## Summary: Key Functions by Task

| Task | Function |
|----|----|
| **Data Integrity** | [`check_statistics_consistency()`](https://blmoran.github.io/bayesma/reference/check_statistics_consistency.md), [`check_n_consistency()`](https://blmoran.github.io/bayesma/reference/check_n_consistency.md) |
| **Trustworthiness** | [`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md), [`inspect_sr_table()`](https://blmoran.github.io/bayesma/reference/inspect_sr_table.md), [`inspect_plot()`](https://blmoran.github.io/bayesma/reference/inspect_plot.md), [`filter_trustworthy()`](https://blmoran.github.io/bayesma/reference/filter_trustworthy.md) |
| **Exploration** | [`primed()`](https://blmoran.github.io/bayesma/reference/primed.md) |
| **RoB** | [`rob_plot()`](https://blmoran.github.io/bayesma/reference/rob_plot.md), `rob_table()` |
| **Specification** | [`bayesma_spec()`](https://blmoran.github.io/bayesma/reference/bayesma_spec.md), [`meta_reg_spec()`](https://blmoran.github.io/bayesma/reference/meta_reg_spec.md), [`robma_spec()`](https://blmoran.github.io/bayesma/reference/robma_spec.md), [`egger_spec()`](https://blmoran.github.io/bayesma/reference/egger_spec.md) |
| **Fitting** | [`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md), [`meta_reg()`](https://blmoran.github.io/bayesma/reference/meta_reg.md), [`robma()`](https://blmoran.github.io/bayesma/reference/robma.md), [`egger()`](https://blmoran.github.io/bayesma/reference/egger.md) |
| **Extraction** | [`bayesma_extract()`](https://blmoran.github.io/bayesma/reference/bayesma_extract.md), [`meta_reg_extract()`](https://blmoran.github.io/bayesma/reference/meta_reg_extract.md), `extract_summary()` |
| **Comparison** | [`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md), [`compare_plot()`](https://blmoran.github.io/bayesma/reference/compare_plot.md), [`compare_table()`](https://blmoran.github.io/bayesma/reference/compare_table.md) |
| **Visualization** | `bayes_forest()`, [`funnel_plot()`](https://blmoran.github.io/bayesma/reference/funnel_plot.md), [`egger_plot()`](https://blmoran.github.io/bayesma/reference/egger_plot.md), [`metareg_mod_plot()`](https://blmoran.github.io/bayesma/reference/metareg_mod_plot.md), [`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md) |
| **Diagnostics** | [`diagnostics()`](https://blmoran.github.io/bayesma/reference/diagnostics.md), [`pp_check()`](https://blmoran.github.io/bayesma/reference/pp_check.md), [`bayesplot::mcmc_trace()`](https://mc-stan.org/bayesplot/reference/MCMC-traces.html) |
| **Reporting** | [`coefficient_evidence()`](https://blmoran.github.io/bayesma/reference/coefficient_evidence.md), [`robma_table()`](https://blmoran.github.io/bayesma/reference/robma_table.md), [`robma_sensitivity()`](https://blmoran.github.io/bayesma/reference/robma_sensitivity.md) |

------------------------------------------------------------------------

## References & Further Reading

- **PRIMED**: Pustejovsky & Tipton (2024). “Preliminary Investigation of
  Meta-analytic Databases”
- **INSPECT-SR**: Wilkinson et al. (2025). “INSPECT-SR: Trustworthiness
  of RCTs”
- **Bayesian Meta-Analysis**: Röver, Knapp, & Friede (2021).
  “Hartung-Knapp-Sidik-Jonkman approach”
- **RoBMA**: Bartoš & Maier (2023). “Robust Bayesian Model Averaging”
- **Stan for Meta-Analysis**: McElreath (2020). “Statistical Rethinking”

------------------------------------------------------------------------

**This workflow document accompanies the `bayesma` R package. For the
latest function documentation, see
[`?bayesma`](https://blmoran.github.io/bayesma/reference/bayesma.md),
[`?primed`](https://blmoran.github.io/bayesma/reference/primed.md),
[`?inspect_sr`](https://blmoran.github.io/bayesma/reference/inspect_sr.md),
and related help pages.**
