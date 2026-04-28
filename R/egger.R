# egger modular pipeline
#
# Six stages mirroring bayesma_pipeline.R:
#   1. egger_spec()       -- validate + extract -> egger_spec
#   2. egger_stan_code()  -- spec -> named Stan blocks + full program
#   3. egger_stan_data()  -- spec -> cmdstanr data list
#   4. egger_fit()        -- compile + sample
#   5. egger_extract()    -- fit + spec -> tidy summaries
#   6. egger_output()     -- assemble final bayesma_egger object


# -----------------------------------------------------------------------------
# Stage 1: spec
# -----------------------------------------------------------------------------

#' @noRd
egger_spec <- function(
    data,
    studyvar,
    n_ctrl,
    n_int,
    event_ctrl    = NULL,
    event_int     = NULL,
    mean_ctrl     = NULL,
    mean_int      = NULL,
    sd_ctrl       = NULL,
    sd_int        = NULL,
    likelihood    = c("binomial", "gaussian", "poisson"),
    heterogeneity = c("multiplicative", "additive"),
    alpha_prior   = NULL,
    beta_prior    = NULL,
    kappa_prior   = NULL,
    gamma_prior   = NULL,
    d_prior       = NULL,
    tau_prior     = NULL,
    credible_level = 0.90,
    custom_model  = NULL,
    custom_data   = NULL
) {
  likelihood    <- rlang::arg_match(likelihood)
  heterogeneity <- rlang::arg_match(heterogeneity)

  if (!is.null(custom_model) && !is.character(custom_model)) {
    cli::cli_abort(
      "{.arg custom_model} must be a character scalar containing Stan code.",
      call = rlang::caller_env()
    )
  }
  if (!is.null(custom_data) && !is.list(custom_data)) {
    cli::cli_abort(
      "{.arg custom_data} must be a named list.",
      call = rlang::caller_env()
    )
  }

  extract_col <- function(d, var_name) {
    if (is.null(var_name)) return(NULL)
    val <- d[[var_name]]
    if (is.null(val)) {
      cli::cli_abort("Variable {.val {var_name}} not found in data.",
                     call = rlang::caller_env())
    }
    val
  }

  study_vec    <- extract_col(data, studyvar)
  S            <- length(study_vec)
  study_labels <- as.character(study_vec)
  n_c          <- extract_col(data, n_ctrl)
  n_i          <- extract_col(data, n_int)

  # ---- Compute effect sizes per likelihood ----
  r0 <- NULL; r1 <- NULL
  if (likelihood == "binomial") {
    if (is.null(event_ctrl) || is.null(event_int)) {
      cli::cli_abort(
        "{.val binomial} likelihood requires {.arg event_ctrl} and {.arg event_int}.",
        call = rlang::caller_env()
      )
    }
    r0 <- extract_col(data, event_ctrl)
    r1 <- extract_col(data, event_int)

    noev0 <- n_c - r0; noev1 <- n_i - r1
    zero_cells <- (r0 == 0) | (r1 == 0) | (noev0 == 0) | (noev1 == 0)
    r0s <- r0; r1s <- r1; noev0s <- noev0; noev1s <- noev1
    if (any(zero_cells)) {
      r0s[zero_cells]    <- r0[zero_cells] + 0.5
      r1s[zero_cells]    <- r1[zero_cells] + 0.5
      noev0s[zero_cells] <- noev0[zero_cells] + 0.5
      noev1s[zero_cells] <- noev1[zero_cells] + 0.5
    }
    yi      <- log(r1s / noev1s) - log(r0s / noev0s)
    sei     <- sqrt(1 / r0s + 1 / r1s + 1 / noev0s + 1 / noev1s)
    measure <- "log_or"

  } else if (likelihood == "gaussian") {
    if (is.null(mean_ctrl) || is.null(mean_int) ||
        is.null(sd_ctrl) || is.null(sd_int)) {
      cli::cli_abort(
        "{.val gaussian} likelihood requires {.arg mean_ctrl}, {.arg mean_int}, {.arg sd_ctrl}, {.arg sd_int}.",
        call = rlang::caller_env()
      )
    }
    m_c <- extract_col(data, mean_ctrl); m_i <- extract_col(data, mean_int)
    s_c <- extract_col(data, sd_ctrl);   s_i <- extract_col(data, sd_int)
    yi      <- m_i - m_c
    sei     <- sqrt(s_c^2 / n_c + s_i^2 / n_i)
    measure <- "mean_diff"
    cli::cli_warn(c(
      "For continuous outcomes, the Bayesian latent-SE method is not directly applicable.",
      "i" = "Using standard Bayesian Egger regression with observed standard errors."
    ))

  } else if (likelihood == "poisson") {
    if (is.null(event_ctrl) || is.null(event_int)) {
      cli::cli_abort(
        "{.val poisson} likelihood requires {.arg event_ctrl} and {.arg event_int}.",
        call = rlang::caller_env()
      )
    }
    r0      <- extract_col(data, event_ctrl)
    r1      <- extract_col(data, event_int)
    yi      <- log((r1 / n_i) / (r0 / n_c))
    sei     <- sqrt(1 / r0 + 1 / r1)
    measure <- "log_rr"
  }

  # ---- Resolve priors ----
  if (is.null(alpha_prior)) alpha_prior <- normal(0, 100)
  if (is.null(beta_prior))  beta_prior  <- normal(0, 100)
  if (is.null(kappa_prior)) kappa_prior <- uniform(0, 2)
  if (is.null(gamma_prior)) gamma_prior <- uniform(0, 2)
  if (is.null(d_prior))     d_prior     <- normal(0, 100)
  if (is.null(tau_prior))   tau_prior   <- uniform(0, 2)

  priors <- list(alpha = alpha_prior, beta = beta_prior,
                 kappa = kappa_prior, gamma = gamma_prior,
                 d = d_prior, tau = tau_prior)

  spec <- list(
    likelihood    = likelihood,
    heterogeneity = heterogeneity,
    S                = S,
    study_labels     = study_labels,
    n_c              = n_c,
    n_i              = n_i,
    r0               = r0,
    r1               = r1,
    yi               = yi,
    sei              = sei,
    measure          = measure,
    priors           = priors,
    credible_level   = credible_level,
    custom_model     = custom_model,
    custom_data      = custom_data
  )
  class(spec) <- c("egger_spec", "list")
  spec
}


#' @export
print.egger_spec <- function(x, ...) {
  cat("<egger_spec>\n",
      "  likelihood     : ", x$likelihood, "\n",
      "  heterogeneity  : ", x$heterogeneity, "\n",
      "  studies (S)    : ", x$S, "\n",
      "  credible_level : ", x$credible_level, "\n",
      sep = "")
  invisible(x)
}


# -----------------------------------------------------------------------------
# Stage 2: stan code
# -----------------------------------------------------------------------------

#' @noRd
egger_stan_code <- function(spec, format = TRUE) {
  if (!inherits(spec, "egger_spec")) {
    cli::cli_abort("{.arg spec} must be an {.cls egger_spec} object.")
  }

  raw <- if (!is.null(spec$custom_model)) {
    spec$custom_model
  } else if (spec$likelihood == "binomial") {
    generate_bayesian_egger_stan_binomial(spec$heterogeneity, spec$priors)
  } else {
    generate_bayesian_egger_stan_generic(spec$heterogeneity, spec$priors)
  }

  full   <- if (isTRUE(format)) format_stan_code(raw) else raw
  blocks <- parse_stan_blocks(full)
  out    <- c(blocks, list(full = full))
  class(out) <- c("egger_stan_code", "bayesma_stan_code", "list")
  out
}


# -----------------------------------------------------------------------------
# Stage 3: stan data
# -----------------------------------------------------------------------------

#' @noRd
egger_stan_data <- function(spec) {
  if (!inherits(spec, "egger_spec")) {
    cli::cli_abort("{.arg spec} must be an {.cls egger_spec} object.")
  }

  sd_list <- if (spec$likelihood == "binomial") {
    list(
      N         = spec$S,
      n0        = as.integer(spec$n_c),
      n1        = as.integer(spec$n_i),
      r0        = as.integer(spec$r0),
      r1        = as.integer(spec$r1),
      y         = spec$yi,
      upp_kappa = spec$priors$kappa$upper %||% 2,
      upp_gamma = spec$priors$gamma$upper %||% 2,
      upp_tau   = spec$priors$tau$upper %||% 2
    )
  } else {
    list(
      N         = spec$S,
      y         = spec$yi,
      se        = spec$sei,
      upp_kappa = spec$priors$kappa$upper %||% 2,
      upp_gamma = spec$priors$gamma$upper %||% 2
    )
  }

  if (!is.null(spec$custom_data)) {
    for (nm in names(spec$custom_data)) sd_list[[nm]] <- spec$custom_data[[nm]]
  }

  sd_list
}


# -----------------------------------------------------------------------------
# Stage 4: fit
# -----------------------------------------------------------------------------

#' @noRd
egger_fit <- function(spec,
                      code          = egger_stan_code(spec),
                      stan_data     = egger_stan_data(spec),
                      chains        = 4,
                      iter_warmup   = 2000,
                      iter_sampling = 4000,
                      adapt_delta   = 0.95,
                      seed          = 1234,
                      ...) {
  if (!inherits(spec, "egger_spec")) {
    cli::cli_abort("{.arg spec} must be an {.cls egger_spec} object.")
  }

  stan_program <- if (inherits(code, "egger_stan_code")) code$full
  else as.character(code)

  mod <- get_cmdstan_model_cached(stan_program)

  fit <- mod$sample(
    data = stan_data, chains = chains,
    parallel_chains = min(chains, parallel::detectCores()),
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
  )

  out <- list(fit = fit, stan_code = code, stan_data = stan_data)
  class(out) <- c("egger_fit", "list")
  out
}


# -----------------------------------------------------------------------------
# Stage 5: extract
# -----------------------------------------------------------------------------

#' @noRd
egger_extract <- function(fit, spec) {
  if (!inherits(fit, "egger_fit")) {
    cli::cli_abort("{.arg fit} must be an {.cls egger_fit} object.")
  }
  if (!inherits(spec, "egger_spec")) {
    cli::cli_abort("{.arg spec} must be an {.cls egger_spec} object.")
  }

  cmdstan_fit <- fit$fit
  alpha_lower <- (1 - spec$credible_level) / 2
  alpha_upper <- 1 - alpha_lower

  key_vars <- c("alpha", "beta")
  key_vars <- c(key_vars,
                if (spec$heterogeneity == "multiplicative") "kappa" else "gamma"
  )
  if (spec$likelihood == "binomial") key_vars <- c(key_vars, "d", "tau")

  summary_tbl <- cmdstan_fit$summary(variables = key_vars) |> tibble::as_tibble()
  draws       <- posterior::as_draws_df(cmdstan_fit$draws(variables = key_vars))

  # ---- Sigma (latent SE) summary for binomial ----
  sigma_summary <- NULL
  if (spec$likelihood == "binomial") {
    sigma_vars <- paste0("sigma_est[", seq_len(spec$S), "]")
    sigma_draws <- posterior::as_draws_df(cmdstan_fit$draws(variables = sigma_vars))
    sigma_summary <- purrr::map(seq_len(spec$S), function(i) {
      d <- as.numeric(sigma_draws[[paste0("sigma_est[", i, "]")]])
      tibble::tibble(
        study        = i,
        sigma_median = stats::median(d),
        sigma_lower  = stats::quantile(d, alpha_lower),
        sigma_upper  = stats::quantile(d, alpha_upper)
      )
    }) |> purrr::list_rbind()
  }

  # ---- Beta summary ----
  beta_draws <- as.numeric(draws$beta)
  beta_summary <- tibble::tibble(
    parameter     = "beta",
    estimate      = stats::median(beta_draws),
    sd            = stats::sd(beta_draws),
    mad           = stats::mad(beta_draws),
    lower         = stats::quantile(beta_draws, alpha_lower),
    upper         = stats::quantile(beta_draws, alpha_upper),
    prob_positive = mean(beta_draws > 0),
    prob_negative = mean(beta_draws < 0)
  )

  # ---- Conclusion ----
  conclusion <- build_egger_conclusion(beta_summary, spec$credible_level)

  eff <- list(
    summary       = summary_tbl,
    draws         = draws,
    beta_summary  = beta_summary,
    sigma_summary = sigma_summary,
    conclusion    = conclusion
  )
  class(eff) <- c("egger_effects", "list")
  eff
}


#' @noRd
build_egger_conclusion <- function(beta_summary, credible_level) {
  prob_positive <- beta_summary$prob_positive
  prob_negative <- beta_summary$prob_negative

  if (prob_negative > prob_positive) {
    direction       <- "negative"
    dominant_prob   <- prob_negative
    interpretation  <- "studies with larger SE tend to show more negative effects"
    bias_impl       <- "This pattern could suggest missing studies with positive/null results."
  } else {
    direction       <- "positive"
    dominant_prob   <- prob_positive
    interpretation  <- "studies with larger SE tend to show more positive effects"
    bias_impl       <- "This pattern could suggest missing studies with negative/null results (publication bias favoring positive findings)."
  }

  strength <- dplyr::case_when(
    dominant_prob >= 0.95 ~ "Strong evidence",
    dominant_prob >= 0.90 ~ "Substantial evidence",
    dominant_prob >= 0.80 ~ "Moderate evidence",
    dominant_prob >= 0.70 ~ "Weak evidence",
    .default = "Little evidence"
  )

  conclusion <- paste0(
    "Posterior probability of ", direction, " small-study effect: ",
    sprintf("%.1f%%", dominant_prob * 100), "\n\n",
    strength, " of ", direction, " small-study effects:\n",
    "  - ", interpretation, "\n",
    "  - ", bias_impl, "\n\n",
    sprintf("%.0f%% CrI for beta: [%.3f, %.3f]",
            credible_level * 100, beta_summary$lower, beta_summary$upper)
  )

  ci_excludes_zero <- (beta_summary$lower > 0) | (beta_summary$upper < 0)
  if (!ci_excludes_zero && dominant_prob >= 0.75) {
    conclusion <- paste0(
      conclusion, "\n\n",
      "Note: While the ", credible_level * 100, "% CrI includes zero, ",
      "there is still a ", sprintf("%.0f%%", dominant_prob * 100),
      " posterior probability\nof a ", direction, " effect. ",
      "Consider this probabilistic evidence rather than a binary decision."
    )
  }

  conclusion
}


# -----------------------------------------------------------------------------
# Stage 6: output
# -----------------------------------------------------------------------------

#' @noRd
egger_output <- function(spec, fit, effects) {
  if (!inherits(spec,    "egger_spec"))    cli::cli_abort("{.arg spec} must be {.cls egger_spec}.")
  if (!inherits(fit,     "egger_fit"))     cli::cli_abort("{.arg fit} must be {.cls egger_fit}.")
  if (!inherits(effects, "egger_effects")) cli::cli_abort("{.arg effects} must be {.cls egger_effects}.")

  code_full <- if (inherits(fit$stan_code, "egger_stan_code")) fit$stan_code$full
  else as.character(fit$stan_code)

  out <- list(
    fit           = fit$fit,
    summary       = effects$summary,
    draws         = effects$draws,
    beta_summary  = effects$beta_summary,
    sigma_summary = effects$sigma_summary,
    conclusion    = effects$conclusion,
    stan_code     = code_full,
    stan_data     = fit$stan_data,
    meta          = list(
      likelihood    = spec$likelihood,
      heterogeneity = spec$heterogeneity,
      measure          = spec$measure,
      study_labels     = spec$study_labels,
      priors           = spec$priors,
      credible_level   = spec$credible_level,
      yi               = spec$yi,
      sei              = spec$sei,
      n_studies        = spec$S
    )
  )
  class(out) <- "bayesma_egger"
  out
}


# -----------------------------------------------------------------------------
# Print / summary
# -----------------------------------------------------------------------------

#' @export
print.bayesma_egger <- function(x, digits = 3, ...) {
  print(summary(x), digits = digits, ...)
}

#' Summarise a fitted \code{bayesma_egger} model
#'
#' @param object A \code{bayesma_egger} object.
#' @param ... Currently unused.
#' @return An object of class \code{egger_summary}.
#' @keywords internal
#' @export
summary.bayesma_egger <- function(object, ...) {
  meta <- object$meta
  cred <- meta$credible_level
  al   <- (1 - cred) / 2
  au   <- 1 - al

  measure_label <- switch(
    meta$measure %||% "log_or",
    log_or    = "log-OR",
    arr       = "ARR",
    mean_diff = "MD",
    log_rr    = "log-IRR",
    meta$measure
  )

  md            <- tryCatch(object$fit$metadata(), error = function(e) NULL)
  chains        <- if (!is.null(md)) md$num_chains    else NA_integer_
  iter_warmup   <- if (!is.null(md)) md$iter_warmup   else NA_integer_
  iter_sampling <- if (!is.null(md)) md$iter_sampling else NA_integer_

  key_vars <- c("alpha", "beta",
    if (meta$heterogeneity == "multiplicative") "kappa" else "gamma")
  if (meta$likelihood == "binomial")
    key_vars <- c(key_vars, "d", "tau")

  drw <- posterior::as_draws_df(object$fit$draws(variables = key_vars))

  coef_mat <- do.call(rbind, purrr::map(key_vars, function(v) {
    d <- as.numeric(drw[[v]])
    matrix(
      c(stats::median(d), stats::mad(d),
        stats::quantile(d, al, names = FALSE),
        stats::quantile(d, au, names = FALSE),
        posterior::ess_bulk(drw[[v]]),
        posterior::ess_tail(drw[[v]]),
        posterior::rhat(drw[[v]])),
      nrow = 1,
      dimnames = list(v, c(
        "Estimate", "Est.Error",
        sprintf("Q%.1f", al * 100), sprintf("Q%.1f", au * 100),
        "Bulk_ESS", "Tail_ESS", "Rhat"
      ))
    )
  }))

  fmt_prior <- function(pr) if (is.null(pr)) NULL else format(pr)
  p <- meta$priors
  priors <- list(
    alpha = fmt_prior(p$alpha),
    beta  = fmt_prior(p$beta),
    het   = fmt_prior(p[[if (meta$heterogeneity == "multiplicative") "kappa" else "gamma"]]),
    het_name = if (meta$heterogeneity == "multiplicative") "kappa" else "gamma"
  )

  structure(
    list(
      likelihood    = meta$likelihood,
      measure       = measure_label,
      heterogeneity = meta$heterogeneity,
      n_studies     = meta$n_studies,
      chains        = chains,
      iter_warmup   = iter_warmup,
      iter_sampling = iter_sampling,
      total_draws   = chains * iter_sampling,
      credible_level = cred,
      priors        = priors,
      coef_mat      = coef_mat,
      beta_summary  = object$beta_summary,
      conclusion    = object$conclusion
    ),
    class = "egger_summary"
  )
}

#' @export
print.egger_summary <- function(x, digits = 3, ...) {
  # ---- Header ----
  cat(" Likelihood:", x$likelihood, "\n")
  cat("    Measure:", x$measure, "\n")
  cat("Heterogeneity:", x$heterogeneity, "\n")
  cat("    Studies:", x$n_studies, "\n")
  cat("      Draws:", paste0(
    x$chains, " chains, each with iter_sampling = ", x$iter_sampling,
    "; warmup = ", x$iter_warmup, ";\n",
    "             total post-warmup draws = ", x$total_draws
  ), "\n\n")

  # ---- Priors ----
  cat("Priors:\n")
  cat("  alpha:", x$priors$alpha, "\n")
  cat("   beta:", x$priors$beta,  "\n")
  cat(sprintf("  %s: %s\n", x$priors$het_name, x$priors$het))
  cat("\n")

  # ---- Coefficients ----
  cat("Regression Coefficients:\n")
  print_format(x$coef_mat, digits = digits)
  cat("\n")

  # ---- Egger conclusion ----
  cat(sprintf("Egger's Test (%g%% CrI):\n", x$credible_level * 100))
  cat(x$conclusion, "\n\n")

  cat(
    "Draws were sampled using CmdStan. For each parameter, Bulk_ESS\n",
    "and Tail_ESS are effective sample size measures, and Rhat is the\n",
    "potential scale reduction factor on split chains (at convergence, Rhat = 1).\n",
    sep = ""
  )

  invisible(x)
}


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

#' Egger's Regression Test for Small-Study Effects (Bayesian)
#'
#' @param data A data frame with one row per study.
#' @param studyvar Character. Column name of the study identifier.
#' @param n_ctrl,n_int Character. Column names of control and intervention
#'   sample sizes.
#' @param event_ctrl,event_int Character. Column names of event counts
#'   (binomial / Poisson likelihoods).
#' @param mean_ctrl,mean_int,sd_ctrl,sd_int Character. Column names of arm
#'   means and SDs (Gaussian likelihood).
#' @param likelihood Character. One of `"binomial"`, `"gaussian"`,
#'   `"poisson"`.
#' @param heterogeneity Character. `"multiplicative"` (default) or
#'   `"additive"`.
#' @param alpha_prior Prior on the intercept.
#' @param beta_prior Prior on the slope (the Egger coefficient).
#' @param kappa_prior Prior on the multiplicative heterogeneity coefficient.
#' @param gamma_prior Prior on the dispersion parameter.
#' @param d_prior Prior on the overdispersion parameter.
#' @param tau_prior Prior on the between-study SD for additive heterogeneity.
#' @param credible_level Numeric in `(0, 1)`. Credible-interval level.
#'   Default `0.90`.
#' @param custom_model Optional character scalar of Stan code overriding the
#'   generated program.
#' @param custom_data Optional named list merged into the Stan data list.
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed MCMC settings.
#' @param ... Passed to `cmdstanr::sample()`.
#' @return An object of class `c("bayesma_egger", "bayesma")`.
#' @export
egger <- function(
    data,
    studyvar,
    n_ctrl,
    n_int,
    event_ctrl     = NULL,
    event_int      = NULL,
    mean_ctrl      = NULL,
    mean_int       = NULL,
    sd_ctrl        = NULL,
    sd_int         = NULL,
    likelihood     = c("binomial", "gaussian", "poisson"),
    heterogeneity  = c("multiplicative", "additive"),
    alpha_prior    = NULL,
    beta_prior     = NULL,
    kappa_prior    = NULL,
    gamma_prior    = NULL,
    d_prior        = NULL,
    tau_prior       = NULL,
    credible_level = 0.90,
    chains         = 4,
    iter_warmup    = 2000,
    iter_sampling  = 4000,
    adapt_delta    = 0.95,
    seed           = 1234,
    custom_model   = NULL,
    custom_data    = NULL,
    ...
) {
  spec      <- egger_spec(
    data = data, studyvar = studyvar,
    n_ctrl = n_ctrl, n_int = n_int,
    event_ctrl = event_ctrl, event_int = event_int,
    mean_ctrl = mean_ctrl, mean_int = mean_int,
    sd_ctrl = sd_ctrl, sd_int = sd_int,
    likelihood = likelihood, heterogeneity = heterogeneity,
    alpha_prior = alpha_prior, beta_prior = beta_prior,
    kappa_prior = kappa_prior, gamma_prior = gamma_prior,
    d_prior = d_prior, tau_prior = tau_prior,
    credible_level = credible_level,
    custom_model = custom_model, custom_data = custom_data
  )
  code      <- egger_stan_code(spec)
  stan_data <- egger_stan_data(spec)
  fit       <- egger_fit(
    spec = spec, code = code, stan_data = stan_data,
    chains = chains, iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, adapt_delta = adapt_delta,
    seed = seed, ...
  )
  effects <- egger_extract(fit, spec)
  egger_output(spec, fit, effects)
}
