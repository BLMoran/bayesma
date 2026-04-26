# Internal helpers for meta_reg.
# Called by meta_reg_spec().

build_moderator_matrix <- function(data, mods, center, scale,
                                   call = rlang::caller_env()) {
  mf <- stats::model.frame(mods, data = data, na.action = stats::na.pass)

  if (anyNA(mf)) {
    n_missing <- sum(!stats::complete.cases(mf))
    cli::cli_abort(c(
      "Missing values found in moderator variables.",
      "x" = "{n_missing} studies have missing moderator values.",
      "i" = "Remove studies with missing moderators or impute values before fitting."
    ), call = call)
  }

  X <- stats::model.matrix(mods, data = mf)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  if (ncol(X) == 0) {
    cli::cli_abort(c(
      "No moderator variables found after processing formula.",
      "i" = "Check that your formula specifies at least one moderator."
    ), call = call)
  }

  continuous_cols <- apply(X, 2, function(col) {
    unique_vals <- unique(col)
    length(unique_vals) > 2 || any(!unique_vals %in% c(0, 1))
  })

  center_values        <- rep(0, ncol(X))
  scale_values         <- rep(1, ncol(X))
  names(center_values) <- colnames(X)
  names(scale_values)  <- colnames(X)

  if (center || scale) {
    for (j in seq_len(ncol(X))) {
      if (continuous_cols[j]) {
        if (center) {
          m <- mean(X[, j])
          X[, j] <- X[, j] - m
          center_values[j] <- m
        }
        if (scale) {
          s <- stats::sd(X[, j])
          if (s > 0) {
            X[, j] <- X[, j] / s
            scale_values[j] <- s
          }
        }
      }
    }
  }

  list(
    X               = X,
    formula         = mods,
    center          = center,
    scale           = scale,
    center_values   = center_values,
    scale_values    = scale_values,
    continuous_cols = continuous_cols
  )
}

validate_prior_args_mreg <- function(model_type, stage, likelihood,
                                     mu_prior, tau_prior, gamma_prior,
                                     beta_prior, beta_priors, coef_names,
                                     call = rlang::caller_env()) {
  is_re <- model_type == "random_effect"

  check_is_prior <- function(x, name) {
    if (!is.null(x) && !inherits(x, "bayesma_prior")) {
      cli::cli_abort(c(
        "{.arg {name}} must be a prior object (e.g., {.code normal(0, 1)}).",
        "i" = "Use {.fn normal}, {.fn half_normal}, {.fn half_cauchy}, etc."
      ), call = call)
    }
  }

  check_is_prior(mu_prior,    "mu_prior")
  check_is_prior(tau_prior,   "tau_prior")
  check_is_prior(gamma_prior, "gamma_prior")
  check_is_prior(beta_prior,  "beta_prior")

  if (!is.null(tau_prior) && !is_re) {
    cli::cli_warn("{.arg tau_prior} is only used in random-effects models. Ignoring.")
  }

  if (!is.null(beta_priors)) {
    if (!is.list(beta_priors)) {
      cli::cli_abort("{.arg beta_priors} must be a named list of priors.", call = call)
    }
    unknown <- setdiff(names(beta_priors), coef_names)
    if (length(unknown) > 0) {
      cli::cli_warn(c(
        "Unknown coefficient names in {.arg beta_priors}: {.val {unknown}}.",
        "i" = "Available coefficients: {.val {coef_names}}."
      ))
    }
    purrr::walk2(beta_priors, names(beta_priors), function(p, nm) {
      check_is_prior(p, paste0("beta_priors$", nm))
    })
  }

  if (!is.null(tau_prior) && tau_prior$family == "normal" && is_re) {
    cli::cli_warn(c(
      "{.arg tau_prior} is lower-bounded at 0 but prior is {.fn normal}.",
      "i" = "Consider {.fn half_normal}, {.fn half_cauchy}, or {.fn half_student_t}."
    ))
  }

  invisible(NULL)
}

resolve_priors_mreg <- function(model_type, stage, likelihood,
                                mu_prior, tau_prior, gamma_prior,
                                beta_prior, beta_priors, K, coef_names) {
  is_re <- model_type == "random_effect"

  if (is.null(mu_prior)) {
    mu_prior <- if (stage == "one_stage" && likelihood == "gaussian") {
      normal(0, 100)
    } else {
      normal(0, 10)
    }
  }

  if (is.null(tau_prior) && is_re) {
    tau_prior <- if (stage == "one_stage" && likelihood == "gaussian") {
      half_student_t(3, 0, 10)
    } else {
      half_student_t(3, 0, 2.5)
    }
  }

  if (is.null(gamma_prior) && stage == "one_stage") {
    gamma_prior <- if (likelihood == "gaussian") {
      normal(0, 100)
    } else {
      normal(0, 10)
    }
  }

  if (is.null(beta_prior)) beta_prior <- normal(0, 1)

  beta_prior_list <- purrr::map(coef_names, function(nm) {
    if (!is.null(beta_priors) && nm %in% names(beta_priors)) {
      beta_priors[[nm]]
    } else {
      beta_prior
    }
  })
  names(beta_prior_list) <- coef_names

  list(
    mu    = mu_prior,
    tau   = tau_prior,
    gamma = gamma_prior,
    beta  = beta_prior_list
  )
}

# Coefficient Evidence - Bayesian summaries for moderator coefficients

#' Bayesian Evidence for Meta-Regression Coefficients
#'
#' Computes Bayesian measures of evidence for moderator coefficients,
#' including probability of direction, null range (Region of Practical Equivalence),
#' and credible intervals.
#'
#' @param object A `bayesma_reg` object.
#' @param null_range Numeric vector of length 2 defining the null range, or NULL
#'   (default) to skip null range calculation. For standardized coefficients,
#'   `c(-0.1, 0.1)` is a common choice.
#' @param ci_level Numeric. Credible interval width (default: 0.95).
#'
#' @return A tibble with columns:
#'   - `term`: Coefficient name
#'   - `estimate`: Posterior median
#'   - `std_error`: Posterior SD
#'   - `ci_lower`, `ci_upper`: Credible interval bounds
#'   - `pd`: Probability of direction (max of P(β>0), P(β<0))
#'   - `direction`: Most probable direction ("positive" or "negative")
#'   - `p_null`: Proportion of posterior inside null range (if specified
#'   - `p_outside_null`: Proportion outside null range (if specified)
#'
#' @details
#' ## Probability of Direction (pd)
#'
#' The probability of direction is the proportion of the posterior on
#' the same side of zero as the median. It ranges from 0.5 (no evidence)
#' to 1 (strong evidence for direction). A pd > 0.95 is often considered
#' meaningful evidence, though this is not a formal threshold.
#'
#' ## null range (Region of Practical Equivalence)
#'
#' The null range defines a range of values considered "practically equivalent
#' to zero". The proportion of the posterior inside the null range indicates
#' whether the effect is practically meaningful:
#'
#' - `p_null` > 0.95: Effect is practically zero (accept null)
#' - `p_null` < 0.05: Effect is practically meaningful (reject null)
#' - Otherwise: Inconclusive
#'
#' Common null range choices:
#' - `c(-0.1, 0.1)` for standardized coefficients (Cohen's d scale)
#' - `c(-0.05, 0.05)` for log-OR (corresponds to OR 0.95-1.05)
#'
#' @examples
#' \dontrun{
#' fit <- meta_reg(data, studyvar = "author", yi = "yi", vi = "vi",
#'                 mods = ~ year + quality)
#'
#' # Basic evidence summary
#' coefficient_evidence(fit)
#'
#' # With null range
#' coefficient_evidence(fit, null_range = c(-0.1, 0.1))
#' }
#'
#' @export
coefficient_evidence <- function(object,
                                 null_range = NULL,
                                 ci_level = 0.95) {

  if (!inherits(object, "bayesma_reg")) {
    cli::cli_abort("{.arg object} must be a {.cls bayesma_reg} object.")
  }

  if (!is.null(null_range)) {
    if (!is.numeric(null_range) || length(null_range) != 2) {
      cli::cli_abort("{.arg null_range} must be a numeric vector of length 2.")
    }
    if (null_range[1] >= null_range[2]) {
      cli::cli_abort("{.arg null_range[1]} must be less than {.arg null_range[2]}.")
    }
  }

  coef_names <- object$meta$coef_names
  K <- length(coef_names)
  ci_probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)

  # Extract draws for each coefficient
  results <- purrr::map(seq_len(K), function(k) {
    vn <- paste0("beta[", k, "]")
    draws <- as.vector(
      posterior::subset_draws(object$fit$draws(vn), variable = vn)
    )

    # Basic summaries
    est <- stats::median(draws)
    se <- stats::sd(draws)
    ci <- stats::quantile(draws, probs = ci_probs)

    # Probability of direction
    p_positive <- mean(draws > 0)
    p_negative <- mean(draws < 0)
    pd <- max(p_positive, p_negative)
    direction <- if (p_positive >= p_negative) "positive" else "negative"

    # null range
    if (!is.null(null_range)) {
      p_null <- mean(draws >= null_range[1] & draws <= null_range[2])
      p_outside <- 1 - p_null
    } else {
      p_null <- NA_real_
      p_outside <- NA_real_
    }

    tibble::tibble(
      term = coef_names[k],
      estimate = est,
      std_error = se,
      ci_lower = ci[1],
      ci_upper = ci[2],
      pd = pd,
      direction = direction,
      p_null = p_null,
      p_outside_null = p_outside
    )
  }) |> purrr::list_rbind()

  # Add intercept
  mu_draws <- as.vector(
    posterior::subset_draws(object$fit$draws("mu"), variable = "mu")
  )

  mu_ci <- stats::quantile(mu_draws, probs = ci_probs)
  p_pos_mu <- mean(mu_draws > 0)
  p_neg_mu <- mean(mu_draws < 0)

  if (!is.null(null_range)) {
    p_null_mu <- mean(mu_draws >= null_range[1] & mu_draws <= null_range[2])
  } else {
    p_null_mu <- NA_real_
  }

  intercept_row <- tibble::tibble(
    term = "(Intercept)",
    estimate = stats::median(mu_draws),
    std_error = stats::sd(mu_draws),
    ci_lower = mu_ci[1],
    ci_upper = mu_ci[2],
    pd = max(p_pos_mu, p_neg_mu),
    direction = if (p_pos_mu >= p_neg_mu) "positive" else "negative",
    p_null = p_null_mu,
    p_outside_null = if (!is.null(null_range)) 1 - p_null_mu else NA_real_
  )

  results <- dplyr::bind_rows(intercept_row, results)

  # Store metadata
  attr(results, "null_range") <- null_range
  attr(results, "ci_level") <- ci_level

  class(results) <- c("bayesma_coef_evidence", class(results))
  results
}

#' @export
print.bayesma_coef_evidence <- function(x, digits = 3, ...) {
  cat("Bayesian Evidence for Meta-Regression Coefficients\n")
  cat(rep("-", 55), "\n", sep = "")

  null_range <- attr(x, "null_range")
  ci_level <- attr(x, "ci_level")

  cat(sprintf("Credible interval: %.0f%%\n", ci_level * 100))
  if (!is.null(null_range)) {
    cat(sprintf("null range: [%.3f, %.3f]\n", null_range[1], null_range[2]))
  }
  cat("\n")

  # Format for printing
  print_df <- x |>
    dplyr::mutate(
      estimate = round(.data$estimate, digits),
      std_error = round(.data$std_error, digits),
      ci = sprintf("[%.3f, %.3f]", .data$ci_lower, .data$ci_upper),
      pd = sprintf("%.1f%%", .data$pd * 100)
    ) |>
    dplyr::select("term", "estimate", "std_error", "ci", "pd", "direction")

  if (!is.null(null_range)) {
    print_df <- print_df |>
      dplyr::left_join(
        x |> dplyr::select("term", "p_null"),
        by = dplyr::join_by(term)
      ) |>
      dplyr::mutate(
        p_null = sprintf("%.1f%%", .data$p_null * 100)
      )
  }

  print(as.data.frame(print_df), row.names = FALSE)

  cat("\n")
  cat("pd = Probability of Direction (confidence in sign)\n")
  if (!is.null(null_range)) {
    cat("p_null = Probability of effect being practically zero\n")
  }

  invisible(x)
}
