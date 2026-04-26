# bayesma: Bayesian Meta Analyses

## Overview

`bayesma` is an R package for conducting Bayesian meta-analyses using
Stan. It provides a comprehensive framework for model fitting,
publication bias adjustment, and visualisation, with support for various
likelihoods, heterogeneity structures, and bias-correction methods.

## Key Features

- **Flexible Model Specification**: One-stage (arm-level) or two-stage
  (effect size) models with Gaussian, binomial, or Poisson likelihoods
- **Heterogeneity Modelling**: Normal, Student-t, skew-normal, and
  finite mixture random-effects distributions
- **Publication Bias Adjustment**: Selection models (Copas,
  weight-function), regression-based methods (PET-PEESE), and robust
  Bayesian model averaging (RoBMA)
- **Extended Models**: Multi-arm trials, meta-regression, multivariate
  outcomes, and dose-response analyses
- **Comprehensive Diagnostics**: Posterior predictive checks, LOO-CV,
  WAIC, and convergence diagnostics (R̂, ESS)
- **Visualisation**: Forest plots, funnel plots, sensitivity analyses,
  and risk of bias summaries
- **Full Prior Control**: Sensible defaults with user-specified priors
  on all parameters

## Installation

You can install the development version of `bayesma` from GitHub:

``` r
remotes::install_github("BLMoran/bayesma")
```

`bayesma` requires `cmdstanr` and CmdStan:

``` r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Basic Usage

### One-Stage Common-Effect Meta-Analysis (Binary Outcome)

``` r
bayesma(
  data =  binary_outcome,
  study      = "Author",
  event_ctrl = "Event_Control",
  event_int  = "Event_Intervention",
  n_ctrl     = "N_Control",
  n_int      = "N_Intervention",
  likelihood = "binomial",
  model_type = "common_effect",
  stage = "one_stage",
  mu_prior = normal(0,1),
  tau_prior = half_cauchy(0, 0.5))
```

### One-Stage Random-Effect Meta-Analysis (Binary Outcome)

### Model Comparison

``` r
bayesma_model_comp <- compare_models(
  "One-Stage, Common Effect" = bayesma_one_stage_ce,
  "One-Stage, Random Effect" = bayesma_one_stage_re,
  "Two-Stage, Common Effect" = bayesma_two_stage_ce,
  "Two-Stage, Random Effect" = bayesma_two_stage_re,
  "One-Stage, Random Effect (adj Small Sample)" = bayesma_one_stage_small_samp
) 
```

## Features

### Analysis Framework

| Stage              | Options                                           |
|--------------------|---------------------------------------------------|
| **Likelihood**     | Binomial, Gaussian, Poisson                       |
| **Pooling**        | Common-effect, random-effects                     |
| **Data structure** | One-stage (arm-level) or two-stage (effect sizes) |

### Heterogeneity Modelling

Not all heterogeneity is the same. `bayesma` offers multiple approaches:

- **Distributional flexibility**: Normal, Student-t, skew-normal random
  effects
- **Mixture models**: Finite mixtures for multi-modal effect
  distributions
- **Robust estimation**: Observation-level outlier detection
- **Full prior control**: Specify priors on μ, τ, and all auxiliary
  parameters

### Publication Bias Adjustment

Three complementary strategies for small-study effects:

| Approach | Methods | Use when… |
|----|----|----|
| **Selection models** | Copas, weight-function | You suspect suppressed non-significant results |
| **Regression-based** | PET-PEESE | You want a simple bias-adjusted estimate |
| **Model averaging** | RoBMA | Mechanism is uncertain; want robust inference |

### Extended Models

- **Multi-arm trials**: Proper within-study correlation structure
- **Meta-regression**: Covariate-adjusted effects
- **Multivariate outcomes**: Correlated endpoints
- **Longitudinal**: Time-varying effects
- **Dose-response**: Flexible dose-response curves

### Visualisation

| Plot | Purpose |
|----|----|
| `forest_plot()` | Study effects with pooled estimate and prediction interval |
| [`overall_plot()`](https://blmoran.github.io/bayesma/reference/overall_plot.md) | Posterior densities for μ, τ, and predictive distribution |
| [`funnel_plot()`](https://blmoran.github.io/bayesma/reference/funnel_plot.md) | Funnel plot for visual bias assessment |
| [`sensitivity_plot()`](https://blmoran.github.io/bayesma/reference/sensitivity_plot.md) | Leave-one-out influence diagnostics |
| `ecdf_plot()` | Empirical CDF of posterior draws |
| [`rob_plot()`](https://blmoran.github.io/bayesma/reference/rob_plot.md) | Risk of bias summary (RoB2, ROBINS-I, ROBINS-E, QUADAS-2) |

### Diagnostics & Model Comparison

- Posterior predictive checks via `bayesplot`
- LOO-CV and WAIC for model comparison
- Bayes R² for explained heterogeneity
- Comprehensive convergence diagnostics (R̂, ESS, divergences)

## Workflow

A typical `bayesma` analysis follows this structure:

## Prior Specification

Sensible defaults with full customisation:

``` r
fit <- bayesma(
  ...,
  mu_prior = normal(0, 1),
  tau_prior = half_cauchy(0, 0.5),
  gamma_prior = normal(0, 10)
)
```

Available distributions:
[`normal()`](https://blmoran.github.io/bayesma/reference/priors.md),
[`half_normal()`](https://blmoran.github.io/bayesma/reference/priors.md),
[`half_cauchy()`](https://blmoran.github.io/bayesma/reference/priors.md),
[`half_student_t()`](https://blmoran.github.io/bayesma/reference/priors.md),
[`exponential()`](https://blmoran.github.io/bayesma/reference/priors.md),
[`uniform()`](https://blmoran.github.io/bayesma/reference/priors.md),
[`beta()`](https://blmoran.github.io/bayesma/reference/priors.md),
[`dirichlet()`](https://blmoran.github.io/bayesma/reference/priors.md).

## Requirements

- R ≥ 4.5.0
- [cmdstanr](https://mc-stan.org/cmdstanr/) and
  [CmdStan](https://mc-stan.org/users/interfaces/cmdstan)
- A C++ toolchain
  ([Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows,
  Xcode CLI on macOS)

## Related Packages

[metafor](https://www.metafor-project.org/) \| Comprehensive frequentist
meta-analysis \|  
[brms](https://paul-buerkner.github.io/brms/) \| General Bayesian
regression \|  
[RoBMA](https://fbartos.github.io/RoBMA/) \| Robust Bayesian model
averaging \|  
[bayesmeta](https://cran.r-project.org/package=bayesmeta) \| Bayesian
random-effects meta-analysis \|  
[meta](https://cran.r-project.org/package=meta) \| General frequentist
meta-analysis \|

## Citation

If you use bayesma in your research, please cite:

``` R
@Manual{,
  title = {bayesma: Bayesian Meta Analysis using Stan in R},
  author = {[Benjamin Moran and Thomas Payne]},
  year = {2026},
  note = {R package version 0.0.0.9000},
  url = {https://github.com/BLMoran/bayesma},
}
```

## Dependencies

bayesma depends on several R packages:

- [gt](https://gt.rstudio.com) for creating tables
- [gt](https://gt.rstudio.com) for creating tables
- [patchwork](https://patchwork.data-imaginist.com) for combining plots
- [ggplot2](https://ggplot2.tidyverse.org) for plotting
- [ggdist](https://mjskay.github.io/ggdist/) for density plotting
- [tidybayes](https://mjskay.github.io/tidybayes/index.html) for tidy
  workflow
- [dplyr](https://dplyr.tidyverse.org),
  [tidyr](https://tidyr.tidyverse.org),
  [purrr](https://purrr.tidyverse.org) for data manipulation
- [paletteer](https://emilhvitfeldt.github.io/paletteer/) for colour
  palettes
- [fontawesome](https://rstudio.github.io/fontawesome/) for risk of bias
  icons

## Feedback, Issues and Contributing

We welcome feedback, suggestions, issues and contributions. Please feel
free to contact either [Ben](mailto:ben.moran@newcastle.edu.au) or
[Tom](mailto:tompayne302@gmail.com) with any feedback. For any bugs,
please file it [here](https://github.com/BLMoran/bayesma/issues) with a
minimal code example to reproduce the issue. Pull requests can be made
[here](https://github.com/BLMoran/bayesma/pulls). Please note that the
bayesma project is released with a [Contributor Code of
Conduct](https://github.com/BLMoran/bayesma/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.

## License

This package is licensed under the GPL-3 License.

## Acknowledgments

bayesma builds upon the excellent work of the
[gt](https://gt.rstudio.com),
[patchwork](https://patchwork.data-imaginist.com),
[ggdist](https://mjskay.github.io/ggdist/),
[tidybayes](https://mjskay.github.io/tidybayes/index.html) and the
[tidyverse](https://www.tidyverse.org) suite of packages. Without the
work of the Stan and R community and their contributions, bayesma would
not be possible. \# bayesma
