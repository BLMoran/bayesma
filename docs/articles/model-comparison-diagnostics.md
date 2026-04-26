# Model comparison and diagnostics

## Introduction

Fitting multiple models raises the question: which model best describes
the data? **bayesma** provides two complementary tools — formal
predictive model comparison via
[`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md)
and visual diagnostics via
[`compare_plot()`](https://blmoran.github.io/bayesma/reference/compare_plot.md)
and
[`ecdf_model_plot()`](https://blmoran.github.io/bayesma/reference/ecdf_model_plot.md).

## Leave-one-study-out cross-validation

The primary comparison criterion is leave-one-study-out cross-validation
(LOSO-CV). For each study $`i`$, the model fitted to the remaining
$`k - 1`$ studies is used to predict the withheld study’s effect. The
predictive quality is measured by the continuous ranked probability
score (CRPS), averaged over studies.

CRPS is a proper scoring rule that rewards sharpness (narrow predictive
distributions) and penalises miscalibration (predictive distributions
that miss the true value). For a predictive distribution $`F_i`$ and
observed value $`y_i`$:

``` math

\text{CRPS}(F_i, y_i) = \int_{-\infty}^{\infty} \left[F_i(x) - \mathbb{1}(x \geq y_i)\right]^2 dx
```

Lower CRPS is better. The total LOSO-CRPS is

``` math

\overline{\text{CRPS}} = \frac{1}{k} \sum_{i=1}^{k} \text{CRPS}(F_{-i}, y_i)
```

### Advantages of LOSO-CV for meta-analysis

- **Cross-stage comparability.** One-stage and two-stage models are
  fitted on different likelihoods but both produce predictive
  distributions on the effect-size scale. LOSO-CV on the effect scale
  makes them comparable.
- **Diagnostic value.** Per-study CRPS scores identify individual
  studies that a model predicts poorly — often the same studies that
  drive sensitivity analysis results.
- **No distributional assumption.** CRPS does not require the predictive
  distribution to be Gaussian.

## LOO-IC and WAIC

Within-stage model comparison (e.g., comparing RE distributions for the
same likelihood) can use LOO information criterion (LOO-IC) or the
widely applicable information criterion (WAIC). These approximate
leave-one-observation-out CV on the log-score scale.

``` math

\text{WAIC} = -2\left(\sum_i \log \overline{p}(y_i \mid \theta) - \sum_i \text{Var}_\theta[\log p(y_i \mid \theta)]\right)
```

LOO-IC is computed via Pareto-smoothed importance sampling and reports
Pareto $`\hat{k}`$ diagnostics for each observation. Studies with
$`\hat{k} > 0.7`$ are influential and their LOO contribution may be
unreliable; LOSO-CV is preferable for these cases.

## Running model comparison

``` r
fit_ce  <- bayesma(data, model_type = "common_effect")
fit_re  <- bayesma(data, model_type = "random_effect")
fit_t   <- bayesma(data, model_type = "random_effect", re_dist = "t")

comp <- compare_models(
  "Common-effect" = fit_ce,
  "RE Gaussian"   = fit_re,
  "RE Student-t"  = fit_t
)

compare_table(comp)
compare_plot(comp)
```

## Interpreting compare_table()

[`compare_table()`](https://blmoran.github.io/bayesma/reference/compare_table.md)
returns a ranked tibble with columns:

| Column        | Description                         |
|---------------|-------------------------------------|
| `model`       | Model label                         |
| `loso_crps`   | Mean LOSO-CRPS (lower = better)     |
| `delta_crps`  | Difference from best model          |
| `coverage_50` | Empirical 50% PI coverage           |
| `coverage_95` | Empirical 95% PI coverage           |
| `elpd_loo`    | LOO expected log-predictive density |

Differences in `loso_crps` less than 0.05 on the log-OR scale are
typically negligible. Prefer simpler models when differences are small.

## Interpreting compare_plot()

[`compare_plot()`](https://blmoran.github.io/bayesma/reference/compare_plot.md)
shows prediction interval coverage at each nominal level (50%, 80%, 90%,
95%). A perfectly calibrated model lies on the diagonal. Models above
the diagonal have intervals that are too wide (conservative); models
below the diagonal have intervals that are too narrow
(anti-conservative).

## Per-study diagnostics

The `$study_scores` element of a
[`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md)
result contains per-study CRPS and coverage statistics. Large per-study
CRPS values flag studies that are systematically mispredicted — a useful
diagnostic for identifying influential outliers or subgroups that
require separate modelling.

``` r
comp$study_scores |>
  dplyr::arrange(dplyr::desc(loso_crps)) |>
  head(5)
```

## When not to compare models

LOSO-CV is not an appropriate tool for selecting between a model with
and without a potential moderator — the moderator comparison is a
question about explanation, not prediction over a fixed set of studies.
Use posterior inclusion probabilities from
[`meta_reg()`](https://blmoran.github.io/bayesma/reference/meta_reg.md)
or [`robma()`](https://blmoran.github.io/bayesma/reference/robma.md) for
moderator selection.

Similarly, comparing a bias-corrected model to an uncorrected model with
CRPS conflates model fit with bias correction — a model that perfectly
reproduces biased data will outscore a correctly specified but
bias-adjusted model. See [Bias &
Heterogeneity](https://blmoran.github.io/bayesma/articles/bayesian-egger-test.md)
for appropriate diagnostics.
