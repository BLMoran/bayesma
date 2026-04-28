#' Compute a marginal estimand from a bayesma fit
#'
#' Post-processes the posterior draws of a fitted [bayesma()] model to return
#' a marginal estimand on the natural scale: risk difference (RD/ARR), average
#' treatment effect (ATE), average treatment effect on the treated (ATT), or
#' a conditional average treatment effect (CATE).
#'
#' For relative-effect estimands (`"OR"`, `"RR"`, `"HR"`, `"IRR"`, `"MD"`,
#' `"SMD"`) this function is unnecessary — [bayesma()] already returns
#' the pooled effect on the appropriate scale.
#'
#' @section Methods by estimand and stage:
#'
#' \describe{
#'   \item{RD / ARR / ATE — binomial, one-stage}{Computed via marginal
#'     standardisation (g-computation) over posterior draws of the per-study
#'     baseline logit (`gamma[s]`) and the pooled log-OR (`mu`), optionally
#'     shifted by study-level random effects (`epsilon[s]`). The per-study RD
#'     is `plogis(gamma[s] + mu + epsilon[s]) - plogis(gamma[s])`, averaged
#'     across studies weighted by the harmonic mean of arm sample sizes.
#'     This corresponds to a population-weighted ATE over the observed study
#'     mix and requires no external baseline assumption.}
#'   \item{RD / ARR / ATE — binomial, two-stage}{Back-transforms posterior
#'     draws of the pooled log-OR using a baseline risk drawn per-iteration
#'     from a Beta distribution fitted (method-of-moments) to the observed
#'     control-arm event rates. Propagates baseline uncertainty into the
#'     posterior RD. A fixed scalar `baseline_risk` bypasses this and uses
#'     the supplied value directly (old behaviour).}
#'   \item{ATE — gaussian}{Equivalent to MD. Returns the pooled posterior on
#'     the absolute scale.}
#'   \item{ATT}{One-stage: weighted by intervention-arm sample size
#'     (`n_int / sum(n_int)`). Two-stage: same back-transform as ATE but with
#'     intervention-size-weighted baseline. Without IPD this is an
#'     arm-size-weighted ATE on the treated, not a true causal ATT —
#'     interpret with caution.}
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


#' Marginal standardisation (g-computation) over the one-stage binomial model.
#'
#' Per-study RD is computed as:
#'   `plogis(gamma[s] + mu + epsilon[s]) - plogis(gamma[s])`
#'
#' where `gamma[s]` is the study baseline on the logit scale, `mu` is the pooled
#' log-OR, and `epsilon[s]` is the study-level random effect (zero for
#' common-effect models). This is mathematically equivalent to what
#' marginaleffects::avg_comparisons() does on a binomial(logit) model.
#'
#' Studies are weighted by the harmonic mean of arm sample sizes for ATE/RD/ARR
#' (precision-weighted, equivalent to inverse-variance weighting of the RD
#' estimator), and by intervention-arm size for ATT.
#'
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

  S       <- ncol(gamma_draws)
  n_draws <- nrow(gamma_draws)

  eps_draws <- if (spec$model_type == "random_effect") {
    try_extract_arm_draws(fit, "epsilon")
  } else {
    NULL
  }

  # Control arm probabilities: plogis(gamma[s])  [n_draws x S]
  p_c <- stats::plogis(gamma_draws)

  # Treatment arm linear predictor: gamma[s] + mu + epsilon[s]  [n_draws x S]
  mu_mat <- if (!is.null(eps_draws)) {
    eps_draws + matrix(mu_draws, nrow = n_draws, ncol = S)
  } else {
    matrix(mu_draws, nrow = n_draws, ncol = S)
  }
  p_i <- stats::plogis(gamma_draws + mu_mat)

  rd_per_study <- p_i - p_c  # [n_draws x S]

  # Weights: harmonic mean of arm sizes for ATE/RD/ARR (precision-weighted);
  # intervention-arm size for ATT.
  weights <- switch(
    estimand,
    ATE = ,
    RD  = ,
    ARR = {
      n_h <- 2 / (1 / spec$n_c + 1 / spec$n_i)
      n_h / sum(n_h)
    },
    ATT = spec$n_i / sum(spec$n_i),
    rep(1 / S, S)
  )

  apply(rd_per_study, 1, stats::weighted.mean, w = weights)
}


#' Back-transform two-stage pooled log-OR draws to the RD scale.
#'
#' When baseline_risk is NULL or "study_mean", the baseline is drawn
#' per-iteration from a Beta distribution fitted by method-of-moments to the
#' observed control-arm event rates. This propagates baseline uncertainty into
#' the posterior RD rather than conditioning on a fixed point estimate.
#'
#' When baseline_risk is a user-supplied scalar, that value is used directly
#' (fixed baseline, equivalent to the old behaviour).
#'
#' @noRd
twostage_rd_draws <- function(fit, spec, estimand) {
  mu      <- extract_pooled_draws(fit)
  n_draws <- length(mu)

  baseline <- resolve_baseline_draws(spec, n_draws)

  or <- exp(mu)
  p1 <- baseline * or / (1 - baseline + baseline * or)
  p1 - baseline
}


#' Resolve baseline risk as a vector of length n_draws.
#'
#' Fixed scalar  -> replicated to n_draws (no uncertainty propagated).
#' NULL / "study_mean" -> method-of-moments Beta fitted to observed control
#'   rates; one draw per posterior iteration so baseline uncertainty is
#'   propagated through to the RD posterior.
#'
#' @noRd
resolve_baseline_draws <- function(spec, n_draws) {
  br <- spec$baseline_risk

  # User-supplied fixed scalar: replicate as-is (old behaviour preserved)
  if (is.numeric(br) && length(br) == 1L) {
    return(rep(br, n_draws))
  }

  # Fit a Beta to the observed control-arm rates by method of moments,
  # then draw one baseline per posterior iteration.
  rates <- spec$outcome_ctrl / spec$n_c
  m     <- mean(rates, na.rm = TRUE)
  v     <- stats::var(rates, na.rm = TRUE)

  # Clamp variance so alpha and beta are strictly positive
  v_max <- m * (1 - m) * 0.99
  v     <- min(v, v_max)

  common <- m * (1 - m) / v - 1
  alpha  <- m * common
  beta_p <- (1 - m) * common

  stats::rbeta(n_draws, shape1 = alpha, shape2 = beta_p)
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

#' Compute a prediction interval for a new study on the marginal estimand scale.
#'
#' For binomial one-stage models: applies `mu_new` (the Stan-generated predicted
#' new study log-OR) at the per-draw average baseline logit across the existing
#' studies, then back-transforms to the marginal (ARR/RD) scale.
#'
#' For other likelihoods or when `mu_new` is unavailable: returns the log-OR
#' prediction interval as-is (same as the default `pred_interval`).
#'
#' @noRd
compute_marginal_pred_interval <- function(fit, spec, marginal_draws) {
  cs_fit <- if (inherits(fit, "bayesma_fit")) fit$fit else fit

  mu_new_draws <- tryCatch(
    as.numeric(posterior::as_draws_matrix(cs_fit$draws("mu_new"))),
    error = function(e) NULL
  )

  if (is.null(mu_new_draws)) return(NULL)

  if (spec$likelihood == "binomial") {
    gamma_rep <- if (spec$stage == "one_stage") {
      gamma_draws <- try_extract_arm_draws(cs_fit, "gamma")
      if (!is.null(gamma_draws)) as.numeric(rowMeans(gamma_draws)) else NULL
    } else {
      p_ctrl <- spec$outcome_ctrl / spec$n_c
      p_ctrl <- pmax(pmin(p_ctrl, 1 - 1e-6), 1e-6)
      rep(mean(stats::qlogis(p_ctrl)), length(mu_new_draws))
    }
    if (!is.null(gamma_rep)) {
      pred_arr <- stats::plogis(gamma_rep + mu_new_draws) - stats::plogis(gamma_rep)
      return(tibble::tibble(
        estimate = stats::median(pred_arr),
        lower    = stats::quantile(pred_arr, 0.025, names = FALSE),
        upper    = stats::quantile(pred_arr, 0.975, names = FALSE)
      ))
    }
  }

  tibble::tibble(
    estimate = stats::median(mu_new_draws),
    lower    = stats::quantile(mu_new_draws, 0.025, names = FALSE),
    upper    = stats::quantile(mu_new_draws, 0.975, names = FALSE)
  )
}
