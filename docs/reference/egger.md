# Egger's Regression Test for Small-Study Effects (Bayesian)

Thin orchestrator over the six-stage pipeline:
[`egger_spec()`](https://blmoran.github.io/bayesma/reference/egger_spec.md),
[`egger_stan_code()`](https://blmoran.github.io/bayesma/reference/egger_stan_code.md),
[`egger_stan_data()`](https://blmoran.github.io/bayesma/reference/egger_stan_data.md),
[`egger_fit()`](https://blmoran.github.io/bayesma/reference/egger_fit.md),
[`egger_extract()`](https://blmoran.github.io/bayesma/reference/egger_extract.md),
[`egger_output()`](https://blmoran.github.io/bayesma/reference/egger_output.md).

## Usage

``` r
egger(
  data,
  studyvar,
  n_ctrl,
  n_int,
  event_ctrl = NULL,
  event_int = NULL,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  likelihood = c("binomial", "gaussian", "poisson"),
  heterogeneity = c("multiplicative", "additive"),
  alpha_prior = NULL,
  beta_prior = NULL,
  kappa_prior = NULL,
  gamma_prior = NULL,
  d_prior = NULL,
  tau_prior = NULL,
  credible_level = 0.9,
  return_stage = c("full", "spec", "code", "data", "fit"),
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  adapt_delta = 0.95,
  seed = 1234,
  custom_model = NULL,
  custom_data = NULL,
  ...
)
```

## Arguments

- data:

  A data frame with one row per study.

- studyvar:

  Character. Column name of the study identifier.

- n_ctrl, n_int:

  Character. Column names of control and intervention sample sizes.

- event_ctrl, event_int:

  Character. Column names of event counts (binomial / Poisson
  likelihoods).

- mean_ctrl, mean_int, sd_ctrl, sd_int:

  Character. Column names of arm means and SDs (Gaussian likelihood).

- likelihood:

  Character. One of `"binomial"`, `"gaussian"`, `"poisson"`.

- heterogeneity:

  Character. `"multiplicative"` (default) or `"additive"`.

- alpha_prior:

  Prior on the intercept.

- beta_prior:

  Prior on the slope (the Egger coefficient).

- kappa_prior:

  Prior on the multiplicative heterogeneity coefficient.

- gamma_prior:

  Prior on the dispersion parameter.

- d_prior:

  Prior on the overdispersion parameter.

- tau_prior:

  Prior on the between-study SD for the additive heterogeneity
  parameterisation.

- credible_level:

  Numeric in `(0, 1)`. Credible-interval level for the summary. Default
  `0.90`.

- return_stage:

  Character. One of `"full"` (default), `"spec"`, `"code"`, `"data"`, or
  `"fit"`.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- custom_model:

  Optional character scalar of Stan code overriding the generated
  program.

- custom_data:

  Optional named list merged into the Stan data list.

- ...:

  Passed to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

An object of class `c("bayesma_egger", "bayesma")`.
