# Compare Multiple Bayesian Meta-Analysis Models

Compares fitted `bayesma` models using LOSO-CV (leave-one-study-out
cross-validation) as the primary ranking criterion, plus within-stage
LOO-CV for diagnostics. LOSO-CV evaluates all models on the same
question — "predict a held-out study" — on the effect-size scale, making
it valid for comparing across one-stage and two-stage models.

## Usage

``` r
compare_models(
  ...,
  data = NULL,
  studyvar = NULL,
  criterion = c("loso", "loo"),
  loso = TRUE,
  coverage_levels = c(0.5, 0.8, 0.9, 0.95),
  moment_match = FALSE,
  cores = getOption("mc.cores", 1),
  quiet = FALSE
)
```

## Arguments

- ...:

  Named or unnamed `bayesma` objects to compare. If unnamed, models are
  labeled "Model 1", "Model 2", etc.

- data:

  The original data frame used to fit all models. Required for LOSO-CV
  (the default).

- studyvar:

  Column identifying studies in `data` (unquoted). Required for LOSO-CV.

- criterion:

  Character. Primary criterion for ranking. `"loso"` (default):
  leave-one-study-out CRPS on the effect-size scale. Comparable across
  stages. `"loo"`: within-stage LOO-CV only.

- loso:

  Logical. Run LOSO-CV. Default `TRUE`. Set to `FALSE` for a fast
  within-stage-only comparison.

- coverage_levels:

  Numeric vector. Nominal coverage levels for calibration assessment in
  LOSO-CV. Default: `c(0.50, 0.80, 0.90, 0.95)`.

- moment_match:

  Logical. Attempt moment matching for LOO-CV. Default `FALSE`.

- cores:

  Integer. Number of cores for LOO computation. Default is
  `getOption("mc.cores", 1)`.

- quiet:

  Logical. Suppress progress messages. Default `FALSE`.

## Value

An object of class `"bayesma_comparison"` containing:

- comparison:

  Tibble with all comparison metrics

- loso_list:

  List of LOSO-CV results per model (if computed)

- loo_compare:

  List of
  [`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
  per stage group

- loo_list:

  List of LOO objects per model

- waic_list:

  List of WAIC objects per model

- diagnostics:

  LOO diagnostics (Pareto k)

- model_names:

  Character vector of model names

## Details

### Why LOSO-CV?

LOO-CV and WAIC operate on each model's native likelihood. One-stage
models condition on arm-level counts; two-stage models condition on
summary effect sizes. These likelihoods are on different scales, so
their ELPD values are **not comparable across stages**.

LOSO-CV sidesteps this by asking each model the same question: *"Given
all studies except study s, what is your predictive distribution for
study s's effect size?"*

Both model types answer on the same scale (the effect-size scale),
making the comparison valid.

### LOSO Metrics

- LOSO-CRPS:

  Continuous Ranked Probability Score averaged over held-out studies. A
  proper scoring rule that evaluates the full predictive distribution.
  Lower is better.

- LOSO-ELPD:

  Mean log predictive density at the held-out study. Higher is better.
  Analogous to LOO-ELPD but truly out-of-sample.

- LOSO-Coverage:

  At each nominal level, what proportion of held-out yi fall inside the
  corresponding prediction interval? Well-calibrated models have
  empirical ≈ nominal.

- LOSO-Miscalibration:

  Mean \|empirical − nominal\| across coverage levels. Zero is perfect.

### Within-Stage LOO

LOO-CV is still computed for within-stage diagnostics (Pareto k values
identify influential studies) and within-stage ranking when all models
share the same stage.

### Computational Cost

LOSO-CV refits each model S times (once per study). For M models with S
studies, this is M × S fits. Set `loso = FALSE` for a fast
within-stage-only comparison, or use `quiet = TRUE` to suppress per-fold
messages.

## See also

[`loo`](https://mc-stan.org/loo/reference/loo.html),
[`loo_compare`](https://mc-stan.org/loo/reference/loo_compare.html)

## Examples

``` r
if (FALSE) { # \dontrun{
comparison <- compare_models(
  "RE (two-stage)" = mod_2s_re,
  "RE (one-stage)"  = mod_1s_re,
  "Heavy-tailed"    = mod_2s_t,
  data = dat,
  studyvar = author
)
print(comparison)
compare_plot(comparison, type = "loso_crps")
compare_plot(comparison, type = "calibration")
compare_table(comparison)

# Fast within-stage only (no refitting)
compare_models(mod1, mod2, criterion = "loo", loso = FALSE)
} # }
```
