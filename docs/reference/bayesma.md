# Run a Bayesian Meta-Analysis in Stan

`bayesma()` is a thin orchestrator over a six-stage pipeline. Each stage
is exported so users can pause for inspection or plug in their own Stan
program via the `custom_model` argument.

## Usage

``` r
bayesma(
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
  model_type = c("random_effect", "common_effect", "bias_corrected", "bc_bnp",
    "selection_copas", "selection_weight", "pet_peese", "mixture_model"),
  stage = c("one_stage", "two_stage"),
  re_dist = c("normal", "t", "skew_normal", "mixture"),
  small_sample = c("none", "t_approx", "hjsk"),
  multi_arm = NULL,
  rho_prior = NULL,
  mu_prior = NULL,
  tau_prior = NULL,
  gamma_prior = NULL,
  nu_prior = NULL,
  alpha_prior = NULL,
  mixture_priors = NULL,
  b_prior = NULL,
  p_bias_prior = NULL,
  w_bias_prior = NULL,
  mu_beta_prior = NULL,
  tau_beta_prior = NULL,
  bnp_concentration_max = NULL,
  bnp_K_max = NULL,
  use_known_bias = FALSE,
  selection_priors = NULL,
  p_cutoffs = c(0.025, 0.05),
  n_components = 2L,
  robust = FALSE,
  robust_prior = NULL,
  robust_df = 4,
  robust_weight = NULL,
  custom_model = NULL,
  custom_data = NULL,
  estimand = NULL,
  cate_covariate = NULL,
  baseline_risk = NULL,
  re_min_k = NULL,
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

  A data frame with one row per study (or per arm for multi-arm
  studies).

- studyvar:

  Character. Column name for study identifiers.

- event_ctrl, event_int:

  Character (binomial/poisson). Event counts.

- mean_ctrl, mean_int:

  Character (gaussian). Arm means.

- sd_ctrl, sd_int:

  Character (gaussian). Arm SDs.

- n_ctrl, n_int:

  Character. Arm sample sizes.

- likelihood:

  Character. `"binomial"`, `"gaussian"`, or `"poisson"`.

- model_type:

  Character. `"random_effect"`, `"common_effect"`, `"bias_corrected"`,
  `"bc_bnp"` (Bias-Corrected Bayesian Non-Parametric; Verde 2025),
  `"selection_copas"` (Robust Bayesian Copas; Bai 2020),
  `"selection_weight"`, `"pet_peese"`, or `"mixture_model"` (Bayesian
  meta-analytic mixture; Maier 2024).

- stage:

  Character. `"one_stage"` or `"two_stage"`.

- re_dist:

  Character. `"normal"`, `"t"`, `"skew_normal"`, or `"mixture"`. For
  `model_type = "selection_copas"`, only `"normal"` and `"t"` are
  supported; pass `"t"` for the Bai (2020) Robust Bayesian Copas (RBC)
  formulation with heavy-tailed random effects (recommended; see
  `nu_prior`).

- small_sample:

  Character. `"none"`, `"t_approx"`, or `"hjsk"`.

- multi_arm:

  Character or NULL. Column name indicating multi-arm study grouping.
  Only used with `stage = "one_stage"`.

- rho_prior:

  Prior on the within-study correlation for multi-arm studies.

- mu_prior:

  Prior on the pooled effect mu.

- tau_prior:

  Prior on the RE standard deviation tau.

- gamma_prior:

  Prior on study baselines gamma. One-stage only.

- nu_prior:

  Prior on RE degrees of freedom nu. `re_dist = "t"` only.

- alpha_prior:

  Prior on RE skewness shape. `re_dist = "skew_normal"` only.

- mixture_priors:

  Named list of priors for mixture components.

- b_prior:

  Prior on bias shift B. `model_type = "bias_corrected"` only.

- p_bias_prior:

  Prior on proportion of biased studies. Used by `bias_corrected` and
  `bc_bnp`.

- w_bias_prior:

  Prior on per-study bias weight. `bias_corrected` only.

- mu_beta_prior:

  Prior on the bias-component location (DP base). Used by `bc_bnp` only.
  Default `uniform(-15, 15)`.

- tau_beta_prior:

  Prior on the bias-component scale (DP base). Used by `bc_bnp` only.
  Default `half_cauchy(0, 1)`.

- bnp_concentration_max:

  Numeric. Upper bound for the DP concentration
  `alpha ~ Uniform(0.5, bnp_concentration_max)`. `bc_bnp` only. Default
  `NULL` resolves from data via Verde (2025) eqn (20).

- bnp_K_max:

  Integer. Truncation level for the stick-breaking approximation.
  `bc_bnp` only. Default `NULL` resolves to
  `1 + ceiling(5 * bnp_concentration_max)`.

- use_known_bias:

  Logical. If `TRUE`, `data` must contain a `biased` column.

- selection_priors:

  Named list of priors for selection model parameters.

- p_cutoffs:

  Numeric vector. One-sided p-value cutpoints for the weight-function
  selection model.

- n_components:

  Integer. Number of mixture components.

- robust:

  Logical. Adds a two-component outlier mixture.

- robust_prior:

  Prior on the outlier component location.

- robust_df:

  Degrees of freedom for the outlier component.

- robust_weight:

  Prior on the mixing weight for the main component.

- custom_model:

  Optional Stan program (character scalar). When non-NULL, code
  generation is bypassed and the user's program is used verbatim.

- custom_data:

  Optional named list of Stan data entries that override or augment the
  automatically built data list.

- estimand:

  Character. The user-facing target effect: `"OR"`, `"RR"`, `"HR"`,
  `"IRR"`, `"MD"`, `"SMD"` (relative-effect / mean-difference estimands;
  behave as today), or `"RD"` / `"ARR"`, `"ATE"`, `"ATT"`, `"CATE"`
  (marginal estimands; computed from posterior draws via
  [`bayesma_marginal()`](https://blmoran.github.io/bayesma/reference/bayesma_marginal.md)
  and attached to the output as `$marginal`). Default `NULL` infers from
  `likelihood`.

- cate_covariate:

  Character. Column name of a study-level covariate to condition on for
  `estimand = "CATE"`. Required when CATE is requested.

- baseline_risk:

  Numeric, between 0 and 1, or `"study_mean"`. Reference baseline risk
  used to back-transform relative effects to absolute scale for
  two-stage binomial fits with `estimand` in
  `c("RD","ARR","ATE","ATT")`. Default `NULL` falls back to the observed
  control-arm event rate.

- re_min_k:

  Optional numeric. If set and `model_type = "random_effect"`, the model
  is downgraded to `"common_effect"` with a warning when the number of
  unique studies is below `re_min_k`. Default `NULL` (no enforcement).
  Useful for subgroup analyses where some strata have too few studies
  for stable random-effects estimation.

- return_stage:

  Character. One of `"full"` (default), `"spec"`, `"code"`, `"data"`, or
  `"fit"`. Returns the intermediate pipeline object instead of the final
  `bayesma` object.

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  MCMC settings.

- ...:

  Passed to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

A list of class `"bayesma"` (or the intermediate stage object when
`return_stage` is not `"full"`).

## Details

The pipeline:

1.  [`bayesma_spec()`](https://blmoran.github.io/bayesma/reference/bayesma_spec.md)
    – validate arguments, extract data, resolve priors.

2.  [`bayesma_stan_code()`](https://blmoran.github.io/bayesma/reference/bayesma_stan_code.md)
    – generate Stan code as named blocks.

3.  [`bayesma_stan_data()`](https://blmoran.github.io/bayesma/reference/bayesma_stan_data.md)
    – build the Stan data list.

4.  [`bayesma_fit()`](https://blmoran.github.io/bayesma/reference/bayesma_fit.md)
    – compile and sample.

5.  [`bayesma_extract()`](https://blmoran.github.io/bayesma/reference/bayesma_extract.md)
    – extract tidy effect components.

6.  [`bayesma_output()`](https://blmoran.github.io/bayesma/reference/bayesma_output.md)
    – assemble the final `bayesma` object.
