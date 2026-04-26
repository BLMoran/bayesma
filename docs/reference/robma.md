# Robust Bayesian Model Averaging for Meta-Analysis

`robma()` is a thin orchestrator over a six-stage pipeline. Each stage
is exported so users can pause for inspection or plug in their own Stan
program(s) via the `custom_model` argument.

## Usage

``` r
robma(
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
  return_stage = c("full", "spec", "code", "data", "fit"),
  format = TRUE,
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

  Numeric vector of length 2 giving the null range on the log scale
  (e.g., `c(-0.1, 0.1)` for log OR). Effects within this range are
  considered practically equivalent to zero. Defaults to `NULL` (point
  null at exactly zero). When specified, direction probabilities and the
  P(practically null) are computed relative to this range. For OR/RR,
  `c(-0.1, 0.1)` corresponds to OR/RR in `[0.905, 1.105]`.

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

  Optional Stan program(s) that override code generation. For
  `method = "bridge"`, a named list of character scalars keyed by model
  label. For `method = "ss"`, a single character scalar.

- custom_data:

  Optional Stan data overrides merged onto the auto-built data list(s).
  Same shape conventions as `custom_model`.

- return_stage:

  Character. One of `"full"` (default), `"spec"`, `"code"`, `"data"`, or
  `"fit"`. Returns the intermediate pipeline object instead of the final
  `bayesma_robma` object.

- format:

  Logical. If `TRUE` (default), auto-format generated Stan programs with
  `stanc --auto-format`.

- ...:

  Additional arguments forwarded to `cmdstanr::CmdStanModel$sample()`.

## Value

An object of class `c("bayesma_robma", "bayesma")`, or an intermediate
pipeline object if `return_stage` is not `"full"`.

## Details

The pipeline:

1.  [`robma_spec()`](https://blmoran.github.io/bayesma/reference/robma_spec.md)
    – validate arguments, resolve priors, build grid.

2.  [`robma_stan_code()`](https://blmoran.github.io/bayesma/reference/robma_stan_code.md)
    – generate Stan programs (one per component / indicator).

3.  [`robma_stan_data()`](https://blmoran.github.io/bayesma/reference/robma_stan_data.md)
    – build the corresponding Stan data lists.

4.  [`robma_fit()`](https://blmoran.github.io/bayesma/reference/robma_fit.md)
    – compile + sample + bridge sampling / PIP extraction.

5.  [`robma_extract()`](https://blmoran.github.io/bayesma/reference/robma_extract.md)
    – null-range probabilities + forest data frame.

6.  [`robma_output()`](https://blmoran.github.io/bayesma/reference/robma_output.md)
    – assemble the final `bayesma_robma` object.
