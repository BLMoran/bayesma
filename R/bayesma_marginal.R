#' Compute a marginal estimand from a bayesma fit
#'
#' Post-processes the posterior draws of a fitted [bayesma()] model to return
#' a marginal estimand on the natural scale: risk difference (RD/ARR), average
#' treatment effect (ATE), average treatment effect on the treated (ATT), or
#' a conditional average treatment effect (CATE).
#'
#' For relative-effect estimands (`"OR"`, `"RR"`, `"HR"`, `"IRR"`, `"MD"`,
#' `"SMD"`) this function is unnecessary — [bayesma_output()] already returns
#' the pooled effect on the appropriate scale.
#'
#' @section Methods by estimand and stage:
#'
#' \describe{
#'   \item{RD / ARR / ATE — binomial, one-stage}{Computed from posterior draws
#'     of the per-arm linear predictor. Mean over studies of
#'     `inv_logit(eta_int) - inv_logit(eta_ctrl)`.}
#'   \item{RD / ARR / ATE — binomial, two-stage}{Requires `baseline_risk`.
#'     Posterior draws of the pooled OR are back-transformed at the supplied
#'     baseline. Defaults to the unweighted mean of the observed control rates
#'     when `baseline_risk = NULL`.}
#'   \item{ATE — gaussian}{Equivalent to MD. Returns the pooled posterior on
#'     the absolute scale.}
#'   \item{ATT}{Computed as ATE weighted by treated-arm sample size
#'     (`n_int / sum(n_int)`). Without IPD this is an arm-size-weighted ATE,
#'     not a true causal ATT — interpret with caution.}
#'   \item{CATE}{Routes to [meta_reg()] with `moderators = cate_covariate`.
#'     Reports the meta-regression effect; the user is expected to evaluate
#'     it at a specific covariate value downstream.}
#' }
#'
#' @param fit A `bayesma_fit` object (or the cmdstanr fit inside one).
#' @param spec The matching `bayesma_spec`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`estimand`}{The estimand label.}
#'   \item{`draws`}{A numeric vector of posterior draws on the natural scale.}
#'   \item{`summary`}{A tibble with median, 95\% CrI, and posterior probability
#'     of being above zero.}
#' }
#'
#' @export
bayesma_marginal <- function(fit, spec) {
  estimand <- spec$estimand
  if (!is_marginal_estimand(estimand)) {
    cli::cli_abort(
      "{.arg estimand} {.val {estimand}} is not a marginal estimand. \\
       Use [bayesma_output()] for relative-effect estimands."
    )
  }

  cs_fit <- if (inherits(fit, "bayesma_fit")) fit$fit else fit

  draws <- switch(
    estimand,
    RD   = ,
    ARR  = ,
    ATE  = ,
    ATT  = compute_absolute_draws(cs_fit, spec, estimand),
    CATE = compute_cate_draws(cs_fit, spec)
  )

  list(
    estimand = estimand,
    draws    = draws,
    summary  = summarise_marginal(draws)
  )
}

#' @noRd
is_marginal_estimand <- function(x) {
  !is.null(x) && x %in% c("RD", "ARR", "ATE", "ATT", "CATE")
}

#' @noRd
resolve_estimand <- function(estimand, likelihood) {
  if (is.null(estimand)) {
    return(switch(
      likelihood,
      binomial = "OR",
      poisson  = "IRR",
      gaussian = "MD"
    ))
  }
  rlang::arg_match0(
    estimand,
    c("OR", "RR", "HR", "IRR", "MD", "SMD",
      "RD", "ARR", "ATE", "ATT", "CATE")
  )
}

#' @noRd
validate_estimand_args <- function(estimand, likelihood, cate_covariate,
                                   baseline_risk, data, stage) {
  binary_compatible <- c("OR", "RR", "HR", "RD", "ARR", "ATE", "ATT", "CATE")
  count_compatible  <- c("IRR", "ATE", "ATT", "CATE")
  cont_compatible   <- c("MD", "SMD", "ATE", "ATT", "CATE")

  ok <- switch(
    likelihood,
    binomial = estimand %in% binary_compatible,
    poisson  = estimand %in% count_compatible,
    gaussian = estimand %in% cont_compatible
  )
  if (!ok) {
    cli::cli_abort(
      "{.arg estimand} {.val {estimand}} is not compatible with \\
       {.arg likelihood} {.val {likelihood}}.",
      call = rlang::caller_env()
    )
  }

  if (estimand == "CATE" && is.null(cate_covariate)) {
    cli::cli_abort(
      "{.arg cate_covariate} is required when {.code estimand = \"CATE\"}.",
      call = rlang::caller_env()
    )
  }
  if (!is.null(cate_covariate) && !cate_covariate %in% names(data)) {
    cli::cli_abort(
      "{.arg cate_covariate} {.val {cate_covariate}} not found in {.arg data}.",
      call = rlang::caller_env()
    )
  }

  if (!is.null(baseline_risk)) {
    valid <- (is.numeric(baseline_risk) && length(baseline_risk) == 1 &&
              baseline_risk > 0 && baseline_risk < 1) ||
             identical(baseline_risk, "study_mean")
    if (!valid) {
      cli::cli_abort(
        "{.arg baseline_risk} must be a single number in (0, 1) or \\
         {.val study_mean}.",
        call = rlang::caller_env()
      )
    }
  }

  invisible(TRUE)
}

#' @noRd
compute_absolute_draws <- function(fit, spec, estimand) {
  if (spec$likelihood == "gaussian") {
    return(extract_pooled_draws(fit))
  }

  if (spec$likelihood != "binomial") {
    cli::cli_abort(
      "Estimand {.val {estimand}} not supported for likelihood \\
       {.val {spec$likelihood}}."
    )
  }

  if (spec$stage == "one_stage") {
    onestage_rd_draws(fit, spec, estimand)
  } else {
    twostage_rd_draws(fit, spec, estimand)
  }
}

#' @noRd
onestage_rd_draws <- function(fit, spec, estimand) {
  gamma_draws <- try_extract_arm_draws(fit, "gamma")
  mu_draws    <- extract_pooled_draws(fit)
  if (is.null(gamma_draws)) {
    cli::cli_warn(
      "One-stage fit does not expose per-study {.var gamma}; \\
       falling back to two-stage back-transform via {.arg baseline_risk}."
    )
    return(twostage_rd_draws(fit, spec, estimand))
  }

  S <- ncol(gamma_draws)
  n_draws <- nrow(gamma_draws)

  eps_draws <- if (spec$model_type == "random_effect") {
    try_extract_arm_draws(fit, "epsilon")
  } else NULL

  p_c <- stats::plogis(gamma_draws)
  mu_mat <- if (!is.null(eps_draws)) eps_draws + mu_draws
            else matrix(mu_draws, n_draws, S)
  p_i <- stats::plogis(gamma_draws + mu_mat)

  rd_per_study <- p_i - p_c
  weights <- arm_weights(spec, estimand)
  apply(rd_per_study, 1, stats::weighted.mean, w = weights)
}

#' @noRd
twostage_rd_draws <- function(fit, spec, estimand) {
  mu <- extract_pooled_draws(fit)
  baseline <- resolve_baseline(spec)
  or <- exp(mu)
  p1 <- baseline * or / (1 - baseline + baseline * or)
  p1 - baseline
}

#' @noRd
resolve_baseline <- function(spec) {
  br <- spec$baseline_risk
  if (is.null(br) || identical(br, "study_mean")) {
    rates <- spec$outcome_ctrl / spec$n_c
    return(mean(rates, na.rm = TRUE))
  }
  br
}

#' @noRd
arm_weights <- function(spec, estimand) {
  if (estimand == "ATT") return(spec$n_i / sum(spec$n_i))
  rep(1 / spec$S, spec$S)
}

#' @noRd
extract_pooled_draws <- function(fit) {
  d <- posterior::as_draws_df(fit$draws("mu"))
  as.numeric(d$mu)
}

#' @noRd
try_extract_arm_draws <- function(fit, var) {
  out <- tryCatch(
    posterior::as_draws_matrix(fit$draws(var)),
    error = function(e) NULL
  )
  if (is.null(out) || ncol(out) == 0) return(NULL)
  out
}

#' @noRd
compute_cate_draws <- function(fit, spec) {
  cli::cli_abort(c(
    "{.code estimand = \"CATE\"} requires a meta-regression fit.",
    i = "Use {.code meta_reg(moderators = cate_covariate)} and evaluate the \\
         posterior at the covariate value of interest."
  ))
}

#' @noRd
summarise_marginal <- function(draws) {
  tibble::tibble(
    median = stats::median(draws),
    lower  = stats::quantile(draws, 0.025, names = FALSE),
    upper  = stats::quantile(draws, 0.975, names = FALSE),
    p_gt_0 = mean(draws > 0)
  )
}
