# Build a RoBMA specification object

Build a RoBMA specification object

## Usage

``` r
robma_spec(
  data,
  studyvar,
  event_ctrl = NULL,
  event_int = NULL,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  n_ctrl = NULL,
  n_int = NULL,
  likelihood = c("binomial", "gaussian", "poisson"),
  priors_effect = NULL,
  priors_effect_null = NULL,
  priors_heterogeneity = NULL,
  priors_heterogeneity_null = NULL,
  priors_bias = NULL,
  priors_bias_null = NULL,
  rescale_priors = 1,
  method = c("bridge", "ss"),
  bias_indicator = c("bias_corrected", "pet_peese", "selection_weight"),
  null_range = NULL,
  b_prior = NULL,
  p_bias_prior = NULL,
  p_cutoffs = c(0.025, 0.05),
  parallel = FALSE,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = 1234,
  quiet = FALSE,
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

- event_ctrl, event_int:

  Character. Column names of event counts (binomial / Poisson
  likelihoods).

- mean_ctrl, mean_int, sd_ctrl, sd_int:

  Character. Column names of arm means and SDs (Gaussian likelihood).

- n_ctrl, n_int:

  Character. Column names of arm sample sizes.

- likelihood:

  Character. One of `"binomial"`, `"gaussian"`, `"poisson"`.

- priors_effect, priors_effect_null, priors_heterogeneity,
  priors_heterogeneity_null, priors_bias, priors_bias_null:

  Lists of `bayesma_prior` objects for the effect, heterogeneity and
  publication-bias components (alternative and null). If `NULL`, the
  RoBMA defaults are used.

- rescale_priors:

  Numeric. Scale factor applied to default priors. Default `1`.

- method:

  Character. `"bridge"` (default) uses bridgesampling across the full
  model grid; `"ss"` uses a single spike-and-slab Stan model.

- bias_indicator:

  Character. Spike-and-slab bias mechanism: `"bias_corrected"`,
  `"pet_peese"`, or `"selection_weight"`.

- null_range:

  Numeric length-2 vector. Range treated as the point null for the
  effect component.

- b_prior:

  Prior on the `b` slope for spike-and-slab bias correction.

- p_bias_prior:

  Prior on the bias inclusion probability.

- p_cutoffs:

  Numeric vector of one-sided p-value cutoffs for selection-weight
  models. Default `c(0.025, 0.05)`.

- parallel:

  Logical. Fit the bridgesampling grid in parallel.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- quiet:

  Logical. Suppress per-step progress messages.

- custom_model:

  Optional character scalar of Stan code overriding the generated
  program (spike-and-slab only).

- custom_data:

  Optional named list merged into the Stan data list.

- ...:

  Additional arguments forwarded to `cmdstanr::CmdStanModel$sample()`.

## Value

An object of class `c("bayesma_robma_spec", "bayesma_spec")`.
