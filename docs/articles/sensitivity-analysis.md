# Sensitivity analysis

## Introduction

A meta-analysis conclusion is only as credible as its robustness.
Sensitivity analysis systematically varies modelling choices and data
subsets to determine whether the pooled estimate is stable or dependent
on specific assumptions.

**bayesma** provides structured sensitivity analysis tools for prior
sensitivity, bias assumption sensitivity, and study-level influence.

## Prior sensitivity

The most common sensitivity analysis in Bayesian meta-analysis varies
the prior on between-study heterogeneity $`\tau`$.

``` r
fits <- list(
  narrow  = bayesma(data, tau_prior = half_normal(0, 0.1)),
  default = bayesma(data, tau_prior = half_cauchy(0, 0.5)),
  wide    = bayesma(data, tau_prior = half_cauchy(0, 1.0))
)

sensitivity_plot(fits, parameter = "mu")
```

[`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md)
shows how the posterior for $`\mu`$ (median and 95% CrI) changes across
prior specifications. A stable estimate across all priors indicates the
conclusion is not prior-driven.

### Prior specifications to test

| Prior                  | Interpretation                                 |
|------------------------|------------------------------------------------|
| `half_normal(0, 0.1)`  | Strong prior: most heterogeneity \< 0.2        |
| `half_normal(0, 0.25)` | Turner et al. (2015) pharmacological benchmark |
| `half_cauchy(0, 0.5)`  | Default weakly informative                     |
| `half_cauchy(0, 1.0)`  | Wide prior: heterogeneity up to 2+ plausible   |
| `uniform(0, 2)`        | Flat over plausible range                      |

## Bias assumption sensitivity

For publication bias models, sensitivity to the assumed bias magnitude
or selection mechanism is often more important than prior sensitivity.

### PET-PEESE vs selection models vs uncorrected

``` r
fits_bias <- list(
  uncorrected = bayesma(data, model_type = "random_effect"),
  pet_peese   = bayesma(data, model_type = "pet_peese"),
  selection   = bayesma(data, model_type = "selection_weight")
)

sensitivity_plot(fits_bias, parameter = "mu")
```

### Copas grid sensitivity

The Copas model’s selection parameters are not fully identified.
[`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md)
for Copas shows the pooled estimate across a grid of
$`(\gamma_0, \gamma_1)`$ values:

``` r
sensitivity_plot(fit_copas, type = "copas_grid")
```

A flat surface indicates the estimate is robust to the selection
mechanism assumption. A steep slope indicates high sensitivity.

### $`\pi_\text{bias}`$ sensitivity (mixture and BC-BNP models)

``` r
fits_pbias <- purrr::map(
  c(beta(1, 9), beta(1, 4), beta(1, 1)),
  \(p) bayesma(data, model_type = "mixture_model", p_bias_prior = p)
)

sensitivity_plot(fits_pbias, parameter = "mu", labels = c("skeptical", "default", "agnostic"))
```

## Study-level influence

**Leave-one-study-out (LOSO) sensitivity** reruns the model excluding
each study in turn, showing how the pooled estimate changes when each
study is removed:

``` r
loso_results <- purrr::map(
  seq_len(nrow(data)),
  \(i) bayesma(data[-i, ], model_type = "random_effect")
)

sensitivity_plot(loso_results, type = "loso", study_labels = data$study_label)
```

Studies with large influence (the estimate changes substantially when
they are removed) should be flagged and their characteristics examined.

## Likelihood sensitivity

For binary outcomes, sensitivity to the choice of summary measure can be
assessed by fitting both log-OR and log-RR models:

``` r
fit_lor <- bayesma(data, estimand = "OR")
fit_lrr <- bayesma(data, estimand = "RR")
compare_models(OR = fit_lor, RR = fit_lrr)
```

## Reporting sensitivity analyses

Report:

- Which priors were varied and the range tested.
- Whether conclusions changed qualitatively (from significant to
  non-significant, or vice versa).
- The specific study or assumption that most affects the conclusion.

Sensitivity analyses that change the direction of the estimate (positive
to null or null to positive) should be highlighted prominently; they
indicate that the conclusion is not robust.
