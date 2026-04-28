# Bayesian Meta-Regression

Bayesian Meta-Regression

## Usage

``` r
meta_reg(
  data,
  studyvar,
  yi = NULL,
  vi = NULL,
  mods,
  event_ctrl = NULL,
  event_int = NULL,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  n_ctrl = NULL,
  n_int = NULL,
  likelihood = c("binomial", "gaussian", "poisson"),
  model_type = c("random_effect", "common_effect"),
  stage = c("two_stage", "one_stage"),
  center = TRUE,
  scale = FALSE,
  small_sample = c("none", "t_approx", "hjsk"),
  mu_prior = NULL,
  tau_prior = NULL,
  gamma_prior = NULL,
  beta_prior = NULL,
  beta_priors = NULL,
  custom_model = NULL,
  custom_data = NULL,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = 1234,
  ...
)
```

## Arguments

- data:

  A data frame with one row per study.

- studyvar:

  Character. Column name of the study identifier.

- yi, vi:

  Character. Column names of pre-computed effect sizes and their
  sampling variances (two-stage only).

- mods:

  One-sided formula specifying moderators (e.g. `~ age + dose`).

- event_ctrl, event_int:

  Character. Column names of event counts for binomial / Poisson
  likelihoods.

- mean_ctrl, mean_int, sd_ctrl, sd_int:

  Character. Column names of arm means and SDs for the Gaussian
  likelihood.

- n_ctrl, n_int:

  Character. Column names of arm sample sizes.

- likelihood:

  Character. One of `"binomial"`, `"gaussian"`, `"poisson"`.

- model_type:

  Character. `"random_effect"` (default) or `"common_effect"`.

- stage:

  Character. `"two_stage"` (default) or `"one_stage"`.

- center:

  Logical. Mean-centre continuous moderators. Default `TRUE`.

- scale:

  Logical. Scale continuous moderators to unit SD. Default `FALSE`.

- small_sample:

  Character. Small-sample adjustment for two-stage models: `"none"`,
  `"t_approx"`, or `"hjsk"`.

- mu_prior:

  Prior on the intercept.

- tau_prior:

  Prior on the between-study SD (random-effects models).

- gamma_prior:

  Prior on the Gaussian arm-level intercept (one-stage only).

- beta_prior:

  Default prior for every regression coefficient.

- beta_priors:

  Named list of coefficient-specific priors, overriding `beta_prior` for
  those coefficients.

- custom_model:

  Optional character scalar of Stan code overriding the generated
  program.

- custom_data:

  Optional named list merged into the Stan data list.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- ...:

  Passed to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

An object of class `c("bayesma_metareg", "bayesma")`.
