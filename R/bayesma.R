#' Run a Bayesian Meta-Analysis in Stan
#'
#' `bayesma()` is a thin orchestrator over a six-stage pipeline. Each stage is
#' exported so users can pause for inspection or plug in their own Stan program
#' via the `custom_model` argument.
#'
#' The pipeline:
#' 1. [bayesma_spec()] -- validate arguments, extract data, resolve priors.
#' 2. [bayesma_stan_code()] -- generate Stan code as named blocks.
#' 3. [bayesma_stan_data()] -- build the Stan data list.
#' 4. [bayesma_fit()] -- compile and sample.
#' 5. [bayesma_extract()] -- extract tidy effect components.
#' 6. [bayesma_output()] -- assemble the final `bayesma` object.
#'
#' @param data A data frame with one row per study (or per arm for multi-arm studies).
#' @param studyvar Character. Column name for study identifiers.
#' @param event_ctrl,event_int Character (binomial/poisson). Event counts.
#' @param mean_ctrl,mean_int Character (gaussian). Arm means.
#' @param sd_ctrl,sd_int Character (gaussian). Arm SDs.
#' @param n_ctrl,n_int Character. Arm sample sizes.
#' @param likelihood Character. `"binomial"`, `"gaussian"`, or `"poisson"`.
#' @param model_type Character. `"random_effect"`, `"common_effect"`,
#'   `"bias_corrected"`, `"bc_bnp"` (Bias-Corrected Bayesian Non-Parametric;
#'   Verde 2025), `"selection_copas"` (Robust Bayesian Copas; Bai 2020),
#'   `"selection_weight"`, `"pet_peese"`, or `"mixture_model"` (Bayesian
#'   meta-analytic mixture; Maier 2024).
#' @param stage Character. `"one_stage"` or `"two_stage"`.
#' @param re_dist Character. `"normal"`, `"t"`, `"skew_normal"`, or `"mixture"`.
#'   For `model_type = "selection_copas"`, only `"normal"` and `"t"` are
#'   supported; pass `"t"` for the Bai (2020) Robust Bayesian Copas (RBC)
#'   formulation with heavy-tailed random effects (recommended; see
#'   `nu_prior`).
#' @param small_sample Character. `"none"`, `"t_approx"`, or `"hjsk"`.
#' @param multi_arm Character or NULL. Column name indicating multi-arm study
#'   grouping. Only used with `stage = "one_stage"`.
#' @param rho_prior Prior on the within-study correlation for multi-arm studies.
#' @param mu_prior Prior on the pooled effect mu.
#' @param tau_prior Prior on the RE standard deviation tau.
#' @param gamma_prior Prior on study baselines gamma. One-stage only.
#' @param nu_prior Prior on RE degrees of freedom nu. `re_dist = "t"` only.
#' @param alpha_prior Prior on RE skewness shape. `re_dist = "skew_normal"` only.
#' @param mixture_priors Named list of priors for mixture components.
#' @param b_prior Prior on bias shift B. `model_type = "bias_corrected"` only.
#' @param p_bias_prior Prior on proportion of biased studies. Used by
#'   `bias_corrected` and `bc_bnp`.
#' @param w_bias_prior Prior on per-study bias weight. `bias_corrected` only.
#' @param mu_beta_prior Prior on the bias-component location (DP base). Used
#'   by `bc_bnp` only. Default `uniform(-15, 15)`.
#' @param tau_beta_prior Prior on the bias-component scale (DP base). Used
#'   by `bc_bnp` only. Default `half_cauchy(0, 1)`.
#' @param bnp_concentration_max Numeric. Upper bound for the DP concentration
#'   `alpha ~ Uniform(0.5, bnp_concentration_max)`. `bc_bnp` only. Default
#'   `NULL` resolves from data via Verde (2025) eqn (20).
#' @param bnp_K_max Integer. Truncation level for the stick-breaking
#'   approximation. `bc_bnp` only. Default `NULL` resolves to
#'   `1 + ceiling(5 * bnp_concentration_max)`.
#' @param use_known_bias Logical. If `TRUE`, `data` must contain a `biased` column.
#' @param selection_priors Named list of priors for selection model parameters.
#' @param p_cutoffs Numeric vector. One-sided p-value cutpoints for the
#'   weight-function selection model.
#' @param n_components Integer. Number of mixture components.
#' @param robust Logical. Adds a two-component outlier mixture.
#' @param robust_prior Prior on the outlier component location.
#' @param robust_df Degrees of freedom for the outlier component.
#' @param robust_weight Prior on the mixing weight for the main component.
#' @param custom_model Optional Stan program (character scalar). When non-NULL,
#'   code generation is bypassed and the user's program is used verbatim.
#' @param custom_data Optional named list of Stan data entries that override or
#'   augment the automatically built data list.
#' @param estimand Character. The user-facing target effect:
#'   `"OR"`, `"RR"`, `"HR"`, `"IRR"`, `"MD"`, `"SMD"` (relative-effect /
#'   mean-difference estimands; behave as today), or `"RD"` / `"ARR"`,
#'   `"ATE"`, `"ATT"`, `"CATE"` (marginal estimands; computed from posterior
#'   draws via [bayesma_marginal()] and attached to the output as `$marginal`).
#'   Default `NULL` infers from `likelihood`.
#' @param cate_covariate Character. Column name of a study-level covariate to
#'   condition on for `estimand = "CATE"`. Required when CATE is requested.
#' @param baseline_risk Numeric, between 0 and 1, or `"study_mean"`. Reference
#'   baseline risk used to back-transform relative effects to absolute scale
#'   for two-stage binomial fits with `estimand` in `c("RD","ARR","ATE","ATT")`.
#'   Default `NULL` falls back to the observed control-arm event rate.
#' @param re_min_k Optional numeric. If set and `model_type = "random_effect"`,
#'   the model is downgraded to `"common_effect"` with a warning when the
#'   number of unique studies is below `re_min_k`. Default `NULL` (no
#'   enforcement). Useful for subgroup analyses where some strata have too few
#'   studies for stable random-effects estimation.
#' @param return_stage Character. One of `"full"` (default), `"spec"`, `"code"`,
#'   `"data"`, or `"fit"`. Returns the intermediate pipeline object instead of
#'   the final `bayesma` object.
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed MCMC settings.
#' @param ... Passed to `cmdstanr::sample()`.
#'
#' @return A list of class `"bayesma"` (or the intermediate stage object when
#'   `return_stage` is not `"full"`).
#'
#' @export
bayesma <- function(
    data,
    studyvar,
    event_ctrl    = NULL,
    event_int     = NULL,
    mean_ctrl     = NULL,
    mean_int      = NULL,
    sd_ctrl       = NULL,
    sd_int        = NULL,
    n_ctrl        = NULL,
    n_int         = NULL,
    likelihood    = c("binomial", "gaussian", "poisson"),
    model_type    = c("random_effect", "common_effect", "bias_corrected",
                      "bc_bnp", "selection_copas", "selection_weight",
                      "pet_peese", "mixture_model"),
    stage         = c("one_stage", "two_stage"),
    re_dist       = c("normal", "t", "skew_normal", "mixture"),
    small_sample  = c("none", "t_approx", "hjsk"),
    multi_arm     = NULL,
    rho_prior     = NULL,
    mu_prior       = NULL,
    tau_prior      = NULL,
    gamma_prior    = NULL,
    nu_prior       = NULL,
    alpha_prior    = NULL,
    mixture_priors = NULL,
    b_prior        = NULL,
    p_bias_prior   = NULL,
    w_bias_prior   = NULL,
    mu_beta_prior  = NULL,
    tau_beta_prior = NULL,
    bnp_concentration_max = NULL,
    bnp_K_max      = NULL,
    use_known_bias = FALSE,
    selection_priors = NULL,
    p_cutoffs      = c(0.025, 0.05),
    n_components  = 2L,
    robust         = FALSE,
    robust_prior   = NULL,
    robust_df      = 4,
    robust_weight  = NULL,
    custom_model   = NULL,
    custom_data    = NULL,
    estimand       = NULL,
    cate_covariate = NULL,
    baseline_risk  = NULL,
    re_min_k       = NULL,
    return_stage   = c("full", "spec", "code", "data", "fit"),
    chains        = 4,
    iter_warmup   = 1000,
    iter_sampling = 1000,
    adapt_delta   = 0.95,
    seed          = 1234,
    ...
) {
  return_stage <- rlang::arg_match(return_stage)

  model_type <- enforce_re_min_k(
    model_type   = model_type,
    re_min_k     = re_min_k,
    data         = data,
    studyvar     = studyvar
  )

  spec <- bayesma_spec(
    data             = data,
    studyvar         = studyvar,
    event_ctrl       = event_ctrl,
    event_int        = event_int,
    mean_ctrl        = mean_ctrl,
    mean_int         = mean_int,
    sd_ctrl          = sd_ctrl,
    sd_int           = sd_int,
    n_ctrl           = n_ctrl,
    n_int            = n_int,
    likelihood       = likelihood,
    model_type       = model_type,
    stage            = stage,
    re_dist          = re_dist,
    small_sample     = small_sample,
    multi_arm        = multi_arm,
    rho_prior        = rho_prior,
    mu_prior         = mu_prior,
    tau_prior        = tau_prior,
    gamma_prior      = gamma_prior,
    nu_prior         = nu_prior,
    alpha_prior      = alpha_prior,
    mixture_priors   = mixture_priors,
    b_prior          = b_prior,
    p_bias_prior     = p_bias_prior,
    w_bias_prior     = w_bias_prior,
    mu_beta_prior    = mu_beta_prior,
    tau_beta_prior   = tau_beta_prior,
    bnp_concentration_max = bnp_concentration_max,
    bnp_K_max        = bnp_K_max,
    use_known_bias   = use_known_bias,
    selection_priors = selection_priors,
    p_cutoffs        = p_cutoffs,
    n_components     = n_components,
    robust           = robust,
    robust_prior     = robust_prior,
    robust_df        = robust_df,
    robust_weight    = robust_weight,
    custom_model     = custom_model,
    custom_data      = custom_data,
    estimand         = estimand,
    cate_covariate   = cate_covariate,
    baseline_risk    = baseline_risk
  )
  if (return_stage == "spec") return(spec)

  code <- bayesma_stan_code(spec)
  if (return_stage == "code") return(code)

  stan_data <- bayesma_stan_data(spec)
  if (return_stage == "data") return(stan_data)

  fit <- bayesma_fit(
    spec          = spec,
    code          = code,
    stan_data     = stan_data,
    chains        = chains,
    iter_warmup   = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta   = adapt_delta,
    seed          = seed,
    ...
  )
  if (return_stage == "fit") return(fit)

  effects <- bayesma_extract(fit, spec)
  out <- bayesma_output(spec, fit, effects)
  if (is_marginal_estimand(spec$estimand)) {
    out$marginal <- bayesma_marginal(fit, spec)
    out$pred_interval <- compute_marginal_pred_interval(fit, spec, out$marginal$draws)
  }
  out
}
