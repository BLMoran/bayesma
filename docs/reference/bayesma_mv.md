# Run a Multivariate Bayesian Meta-Analysis in Stan

Thin orchestrator over the six-stage pipeline:
[`bayesma_mv_spec()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_spec.md),
[`bayesma_mv_stan_code()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_stan_code.md),
[`bayesma_mv_stan_data()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_stan_data.md),
[`bayesma_mv_fit()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_fit.md),
[`bayesma_mv_extract()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_extract.md),
[`bayesma_mv_output()`](https://blmoran.github.io/bayesma/reference/bayesma_mv_output.md).

## Usage

``` r
bayesma_mv(
  data,
  studyvar,
  mean_ctrl_1,
  mean_int_1,
  sd_ctrl_1,
  sd_int_1,
  n_ctrl_1,
  n_int_1,
  mean_ctrl_2,
  mean_int_2,
  sd_ctrl_2,
  sd_int_2,
  n_ctrl_2,
  n_int_2,
  outcome_labels = c("outcome_1", "outcome_2"),
  likelihood = c("gaussian"),
  stage = c("two_stage", "one_stage"),
  rho_within = 0.5,
  mu_prior = NULL,
  tau_prior = NULL,
  rho_between_prior = NULL,
  custom_model = NULL,
  custom_data = NULL,
  return_stage = c("full", "spec", "code", "data", "fit"),
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

  A data frame with one row per study containing arm-level data for both
  outcomes.

- studyvar:

  Column for study identifiers (unquoted).

- mean_ctrl_1, mean_int_1:

  Columns for control and intervention means for outcome 1 (unquoted).

- sd_ctrl_1, sd_int_1:

  Columns for control and intervention SDs for outcome 1 (unquoted).

- n_ctrl_1, n_int_1:

  Columns for control and intervention sample sizes for outcome 1
  (unquoted).

- mean_ctrl_2, mean_int_2:

  Columns for control and intervention means for outcome 2 (unquoted).

- sd_ctrl_2, sd_int_2:

  Columns for control and intervention SDs for outcome 2 (unquoted).

- n_ctrl_2, n_int_2:

  Columns for control and intervention sample sizes for outcome 2
  (unquoted).

- outcome_labels:

  Character vector of length 2 with labels for the two outcomes.
  Default: `c("outcome_1", "outcome_2")`.

- likelihood:

  Character. Currently only `"gaussian"`.

- stage:

  Character. `"two_stage"` (effect sizes computed then modelled) or
  `"one_stage"` (marginalised model).

- rho_within:

  Numeric scalar in `[-1, 1]`. Within-study correlation between the two
  outcomes, assumed known. Default: `0.5`.

- mu_prior:

  Prior on pooled effects. Either a single prior (applied to both
  outcomes) or a named list with elements matching `outcome_labels`.

- tau_prior:

  Prior on between-study SDs. Either a single prior or a named list.

- rho_between_prior:

  Prior on the between-study correlation. Default: `uniform(-1, 1)`.

- custom_model:

  Optional character scalar containing complete Stan code to override
  the generated program.

- custom_data:

  Optional named list merged into the Stan data list.

- return_stage:

  Character. One of `"full"` (default), `"spec"`, `"code"`, `"data"`, or
  `"fit"`.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- ...:

  Passed to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

An object of class `c("bayesma_mv", "bayesma")`.
