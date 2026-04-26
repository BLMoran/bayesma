#' Extract log_lik matrix and chain_id from bayesma object
#'
#' @param model A bayesma object
#' @param aggregate_to_study Logical. If TRUE and model is one-stage,
#'   aggregate arm-level log_lik to study-level by summing.
#' @noRd
extract_log_lik <- function(model, aggregate_to_study = FALSE) {
  all_vars <- posterior::variables(model$fit$draws())
  log_lik_vars <- grep("^log_lik\\[", all_vars, value = TRUE)

  draws_df <- posterior::as_draws_df(
    model$fit$draws(variables = log_lik_vars)
  )

  chain_id <- draws_df$.chain

  log_lik_cols <- grep("^log_lik\\[", names(draws_df), value = TRUE)
  log_lik <- as.matrix(draws_df[, log_lik_cols, drop = FALSE])

  stage <- model$meta$stage %||% "two_stage"

  if (aggregate_to_study && stage == "one_stage") {
    S <- model$stan_data$S

    if (!is.null(S) && ncol(log_lik) == 2 * S) {
      log_lik_study <- purrr::map(seq_len(S), function(s) {
        log_lik[, s] + log_lik[, s + S]
      }) |>
        purrr::list_cbind()
      colnames(log_lik_study) <- paste0("log_lik[", seq_len(S), "]")
      log_lik <- log_lik_study
    } else if (!is.null(S)) {
      cli::cli_warn(
        "Could not aggregate one-stage log_lik to study level. Expected {2 * S} columns, got {ncol(log_lik)}."
      )
    }
  }

  list(log_lik = log_lik, chain_id = chain_id)
}


#' Extract observed y values from bayesma model
#'
#' @noRd
extract_observed_y <- function(model) {
  stage <- model$meta$stage %||% "two_stage"
  likelihood <- model$meta$likelihood %||% "gaussian"

  if (stage == "two_stage") {
    if (!is.null(model$meta$es)) return(model$meta$es$yi)
  } else {
    if (likelihood == "gaussian" && !is.null(model$stan_data$y)) {
      return(model$stan_data$y)
    }
    if (likelihood %in% c("binomial", "poisson") && !is.null(model$stan_data$events)) {
      return(as.numeric(model$stan_data$events))
    }
    if (!is.null(model$arm_data$outcome)) {
      return(as.numeric(model$arm_data$outcome))
    }
  }

  NULL
}


# ============================================================================
# Study-level prediction extraction
# ============================================================================

#' Extract study-level effect draws from any bayesma model
#'
#' Returns an S x n_draws matrix of theta_i draws regardless of stage/likelihood.
#'
#' @param model A bayesma object.
#' @return A list with theta_draws (S x n_draws matrix), yi, sei, S.
#' @noRd
extract_study_theta_draws <- function(model) {
  stage      <- model$meta$stage %||% "two_stage"
  model_type <- model$meta$model_type %||% "random_effect"
  is_re      <- model_type == "random_effect"
  es         <- model$meta$es

  if (is.null(es)) {
    cli::cli_abort(
      "Model does not contain stored effect sizes in {.code meta$es}.",
      call = rlang::caller_env()
    )
  }

  yi  <- es$yi
  sei <- es$sei
  S   <- length(yi)

  all_vars <- posterior::variables(model$fit$draws())

  mu_draws <- as.vector(
    posterior::subset_draws(model$fit$draws("mu"), variable = "mu")
  )
  n_draws <- length(mu_draws)

  extract_var_draws <- function(vn) {
    as.vector(posterior::subset_draws(model$fit$draws(vn), variable = vn))
  }

  theta_mat <- if (is_re && stage == "two_stage" && any(grepl("^theta\\[", all_vars))) {
    purrr::map(seq_len(S), ~ extract_var_draws(paste0("theta[", .x, "]"))) |>
      purrr::list_rbind()
  } else if (is_re && stage == "one_stage" && any(grepl("^epsilon\\[", all_vars))) {
    purrr::map(seq_len(S), function(i) {
      mu_draws + extract_var_draws(paste0("epsilon[", i, "]"))
    }) |>
      purrr::list_rbind()
  } else {
    matrix(rep(mu_draws, each = S), nrow = S, ncol = n_draws)
  }

  list(theta_draws = theta_mat, yi = yi, sei = sei, S = S)
}


# ============================================================================
# CRPS (Continuous Ranked Probability Score)
# ============================================================================

#' Compute CRPS via energy form
#'
#' @param y Scalar observed value.
#' @param draws Numeric vector of posterior predictive draws.
#' @return Scalar CRPS (lower is better).
#' @noRd
crps_sample <- function(y, draws) {
  mean(abs(draws - y)) - 0.5 * mean(abs(outer(draws, draws, "-")))
}


# ============================================================================
# LOSO-CV (Leave-One-Study-Out Cross-Validation)
# ============================================================================

#' Extract predictive draws for a held-out study from a refitted model
#'
#' Uses mu_new draws (RE models) or mu draws (CE models) as the prediction
#' for the held-out study's true effect.
#'
#' @param refit A bayesma object fitted on leave-one-out data.
#' @param yi_held Scalar. Observed effect size of the held-out study.
#' @param sei_held Scalar. Standard error of the held-out study.
#' @return A list with theta_draws, y_draws, yi, sei.
#' @noRd
extract_loso_predictive_draws <- function(refit, yi_held, sei_held) {
  all_vars <- posterior::variables(refit$fit$draws())

  mu_draws <- as.vector(
    posterior::subset_draws(refit$fit$draws("mu"), variable = "mu")
  )
  n_draws <- length(mu_draws)

  theta_draws <- if ("mu_new" %in% all_vars) {
    as.vector(
      posterior::subset_draws(refit$fit$draws("mu_new"), variable = "mu_new")
    )
  } else {
    mu_draws
  }

  y_draws <- stats::rnorm(n_draws, mean = theta_draws, sd = sei_held)

  list(
    theta_draws = theta_draws,
    y_draws     = y_draws,
    yi          = yi_held,
    sei         = sei_held
  )
}


#' Run LOSO-CV for a single bayesma model
#'
#' Refits the model S times, each time holding out one study, then evaluates
#' the predictive distribution against the held-out study on the effect-size
#' scale.
#'
#' @param model A bayesma object with stored call_args.
#' @param data The original data frame.
#' @param studyvar Character. Column name identifying studies.
#' @param coverage_levels Numeric vector. Nominal coverage levels.
#' @param max_draws Integer. Maximum draws for CRPS computation.
#' @param quiet Logical. Suppress per-study progress messages.
#' @return A list of class `"bayesma_loso"`.
#' @noRd
run_loso_cv <- function(model, data, studyvar,
                        coverage_levels = c(0.50, 0.80, 0.90, 0.95),
                        max_draws = 1000L,
                        quiet = FALSE) {

  if (is.null(model$meta$call_args)) {
    cli::cli_abort(c(
      "Cannot run LOSO-CV: model does not contain stored call arguments.",
      "i" = "Refit your model with the latest version of bayesma."
    ), call = rlang::caller_env())
  }

  es <- model$meta$es
  if (is.null(es)) {
    cli::cli_abort(
      "Model does not contain stored effect sizes in {.code meta$es}.",
      call = rlang::caller_env()
    )
  }

  study_col    <- data[[studyvar]]
  study_labels <- if (is.factor(study_col)) levels(study_col)
  else unique(as.character(study_col))
  S <- length(study_labels)

  inform <- if (quiet) \(...) invisible(NULL) else cli::cli_inform

  per_study <- purrr::map(seq_len(S), function(s) {
    inform("LOSO fold {s}/{S}: holding out {.val {study_labels[s]}}")

    data_minus_s <- data[study_col != study_labels[s], , drop = FALSE]

    refit_s <- tryCatch(
      refit_bayesma(model, newdata = data_minus_s),
      error = function(e) {
        cli::cli_warn(
          "Refit failed for study {.val {study_labels[s]}}: {conditionMessage(e)}"
        )
        NULL
      }
    )

    if (is.null(refit_s)) {
      return(tibble::tibble(
        study     = study_labels[s],
        study_idx = s,
        yi        = es$yi[s],
        sei       = es$sei[s],
        crps      = NA_real_,
        log_pd    = NA_real_
      ))
    }

    pred <- extract_loso_predictive_draws(refit_s, es$yi[s], es$sei[s])

    theta_sub <- if (length(pred$theta_draws) > max_draws) {
      sample(pred$theta_draws, max_draws)
    } else {
      pred$theta_draws
    }

    crps_s <- crps_sample(es$yi[s], theta_sub)

    log_pd_s <- log(mean(stats::dnorm(es$yi[s],
                                      mean = pred$theta_draws,
                                      sd   = es$sei[s])))

    coverage_cols <- purrr::map(coverage_levels, function(lvl) {
      alpha_lower <- (1 - lvl) / 2
      alpha_upper <- 1 - alpha_lower
      q <- stats::quantile(pred$theta_draws, probs = c(alpha_lower, alpha_upper))
      covered <- (es$yi[s] >= q[[1]]) & (es$yi[s] <= q[[2]])
      width   <- q[[2]] - q[[1]]
      nm_cov  <- paste0("covered_", gsub("\\.", "", as.character(lvl * 100)))
      nm_wid  <- paste0("width_", gsub("\\.", "", as.character(lvl * 100)))
      tibble::tibble(!!nm_cov := covered, !!nm_wid := width)
    }) |>
      purrr::list_cbind()

    dplyr::bind_cols(
      tibble::tibble(
        study     = study_labels[s],
        study_idx = s,
        yi        = es$yi[s],
        sei       = es$sei[s],
        crps      = crps_s,
        log_pd    = log_pd_s
      ),
      coverage_cols
    )
  }) |> purrr::list_rbind()

  cal_tbl <- purrr::map(coverage_levels, function(lvl) {
    cov_col <- paste0("covered_", gsub("\\.", "", as.character(lvl * 100)))
    wid_col <- paste0("width_", gsub("\\.", "", as.character(lvl * 100)))
    tibble::tibble(
      nominal             = lvl,
      empirical           = mean(per_study[[cov_col]], na.rm = TRUE),
      n_covered           = sum(per_study[[cov_col]], na.rm = TRUE),
      n_total             = sum(!is.na(per_study[[cov_col]])),
      mean_interval_width = mean(per_study[[wid_col]], na.rm = TRUE)
    )
  }) |> purrr::list_rbind()

  miscal <- mean(abs(cal_tbl$empirical - cal_tbl$nominal))

  out <- list(
    per_study      = per_study,
    mean_crps      = mean(per_study$crps, na.rm = TRUE),
    median_crps    = stats::median(per_study$crps, na.rm = TRUE),
    mean_log_pd    = mean(per_study$log_pd, na.rm = TRUE),
    calibration    = cal_tbl,
    miscalibration = miscal,
    n_failed       = sum(is.na(per_study$crps)),
    S              = S
  )
  class(out) <- "bayesma_loso"
  out
}


#' @export
print.bayesma_loso <- function(x, digits = 3, ...) {
  cli::cli_h2("LOSO Cross-Validation Results")
  cli::cli_bullets(c(
    "*" = "Studies: {x$S} ({x$n_failed} refit failures)",
    "*" = "Mean CRPS: {round(x$mean_crps, digits)}",
    "*" = "LOSO-ELPD: {round(x$mean_log_pd, digits)}",
    "*" = "Miscalibration: {round(x$miscalibration, digits)}"
  ))
  cli::cli_text("")
  cli::cli_h3("Calibration")
  print(
    x$calibration |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits))),
    n = Inf
  )
  invisible(x)
}
