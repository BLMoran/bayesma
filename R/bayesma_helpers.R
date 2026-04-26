# Helper Functions

#' Enforce a minimum number of studies for random-effects models
#'
#' If `re_min_k` is set and the number of unique studies is below it,
#' a random-effects model_type is downgraded to common_effect with a warning.
#' Other model types are passed through unchanged.
#'
#' @noRd
enforce_re_min_k <- function(model_type, re_min_k, data, studyvar) {
  if (is.null(re_min_k)) return(model_type)

  if (!is.numeric(re_min_k) || length(re_min_k) != 1 || re_min_k < 1) {
    cli::cli_abort(
      "{.arg re_min_k} must be a single positive number.",
      call = rlang::caller_env()
    )
  }

  resolved <- if (length(model_type) > 1) model_type[1] else model_type
  if (resolved != "random_effect") return(model_type)

  k <- length(unique(data[[studyvar]]))
  if (k >= re_min_k) return(model_type)

  cli::cli_warn(c(
    "!" = "Switching to {.val common_effect}: only {k} stud{?y/ies} \\
           (< {.arg re_min_k} = {re_min_k}).",
    "i" = "Random-effects estimation is unstable below {re_min_k} studies."
  ))
  "common_effect"
}

#' @noRd
fixef.bayesma <- function(x) {
  mu_row <- x$summary |>
    dplyr::filter(.data$variable == "mu")

  mat <- matrix(
    c(mu_row$median, mu_row$mad, mu_row$q5, mu_row$q95),
    nrow = 1,
    dimnames = list("Intercept", c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  )
  mat
}


#' Extract draws from a bayesma object in bayesfoRest-compatible format
#'
#' Reshapes the posterior draws from a bayesma fit into the long-format
#' tibble that bayesfoRest's forest_data_fn expects.
#'
#' Columns produced:
#'   - Author (character): study label or "Pooled Effect" / "Prediction"
#'   - b_Intercept (numeric): study-level effect draw (mu + random effect)
#'   - r_Author (numeric): random effect deviation from pooled
#'   - sd_Author__Intercept (numeric): tau draw
#'   - .chain, .iteration, .draw: MCMC identifiers
#'
#' @param x A bayesma object
#' @param estimand The estimand string. When a marginal binary estimand
#'   (RD/ARR/ATE/ATT) is requested and per-study gamma draws are available,
#'   `b_Intercept` is returned on the risk-difference scale.
#' @return A tibble of posterior draws in long format
#'
#' @noRd
extract_forest_draws <- function(x, estimand = NULL) {
  meta       <- x$meta
  S          <- length(meta$study_labels)
  stage      <- meta$stage
  model_type <- meta$model_type
  re_dist    <- meta$re_dist
  is_re      <- model_type == "random_effect"

  is_marginal_binary <- !is.null(estimand) &&
    is_marginal_estimand(estimand) &&
    identical(meta$likelihood, "binomial")

  raw_draws <- x$draws
  n_draws   <- nrow(raw_draws)
  mu_draws  <- as.numeric(raw_draws[["mu"]])

  # Per-study gamma (control-arm baseline logit) for marginal standardisation.
  # One-stage: modelled gamma[s] draws available from Stan.
  # Two-stage: no gamma parameter — derive from observed control-arm event rates
  #   stored in model$arm_data. These are fixed across draws (scalar per study).
  gamma_mat      <- NULL
  gamma_fixed    <- NULL  # length-S numeric vector for two-stage
  if (is_marginal_binary) {
    if (stage == "one_stage") {
      gamma_mat <- try_extract_arm_draws(x$fit, "gamma")  # n_draws x S
    } else if (stage == "two_stage" &&
               !is.null(meta$outcome_ctrl) && !is.null(meta$n_c)) {
      p_ctrl      <- meta$outcome_ctrl / meta$n_c
      p_ctrl      <- pmax(pmin(p_ctrl, 1 - 1e-6), 1e-6)
      gamma_fixed <- stats::qlogis(p_ctrl)
    }
  }

  # Pooled effect: use marginal$draws if available (already on RD scale)
  pooled_effect_draws <- if (is_marginal_binary && !is.null(x$marginal)) {
    x$marginal$draws
  } else {
    mu_draws
  }

  if (is_re && "tau" %in% names(raw_draws)) {
    tau_draws <- as.numeric(raw_draws[["tau"]])
  } else {
    tau_draws <- NA_real_
  }

  chain_col <- if (".chain" %in% names(raw_draws)) {
    as.integer(raw_draws[[".chain"]])
  } else {
    rep(1L, n_draws)
  }
  iteration_col <- if (".iteration" %in% names(raw_draws)) {
    as.integer(raw_draws[[".iteration"]])
  } else {
    seq_len(n_draws)
  }
  draw_col <- if (".draw" %in% names(raw_draws)) {
    as.integer(raw_draws[[".draw"]])
  } else {
    seq_len(n_draws)
  }

  # --- Study-level draws ---
  study_draws_list <- purrr::map(seq_len(S), function(i) {
    label <- meta$study_labels[i]

    if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
      if (stage == "two_stage") {
        theta_var <- paste0("theta[", i, "]")
        theta_i <- as.numeric(raw_draws[[theta_var]])
      } else {
        eps_var <- paste0("epsilon[", i, "]")
        eps_i <- as.numeric(raw_draws[[eps_var]])
        theta_i <- mu_draws + eps_i
      }
    } else {
      theta_i <- mu_draws
    }

    # Transform to marginal (RD) scale
    if (is_marginal_binary && !is.null(gamma_mat)) {
      gamma_i     <- as.numeric(gamma_mat[, i])
      b_intercept <- stats::plogis(gamma_i + theta_i) - stats::plogis(gamma_i)
    } else if (is_marginal_binary && !is.null(gamma_fixed)) {
      g           <- gamma_fixed[i]
      b_intercept <- stats::plogis(g + theta_i) - stats::plogis(g)
    } else {
      b_intercept <- theta_i
    }

    tibble::tibble(
      Author               = label,
      b_Intercept          = b_intercept,
      r_Author             = b_intercept - pooled_effect_draws,
      sd_Author__Intercept = tau_draws,
      .chain               = chain_col,
      .iteration           = iteration_col,
      .draw                = draw_col
    )
  })

  study_draws <- purrr::list_rbind(study_draws_list)

  # --- Pooled draws ---
  pooled_draws <- tibble::tibble(
    Author               = "Pooled Effect",
    b_Intercept          = pooled_effect_draws,
    r_Author             = 0,
    sd_Author__Intercept = tau_draws,
    .chain               = chain_col,
    .iteration           = iteration_col,
    .draw                = draw_col
  )

  # --- Prediction draws ---
  pred_draws <- NULL
  mu_new_draws <- if ("mu_new" %in% names(raw_draws)) {
    as.numeric(raw_draws[["mu_new"]])
  } else {
    tryCatch(
      as.numeric(posterior::subset_draws(x$fit$draws("mu_new"), variable = "mu_new")),
      error = function(e) NULL
    )
  }

  if (is_re && !is.null(mu_new_draws)) {
    if (is_marginal_binary) {
      # Predict ARR for a new study at the average baseline risk.
      # One-stage: per-draw mean of modelled gamma; two-stage: mean of observed logits.
      gamma_rep <- if (!is.null(gamma_mat)) {
        as.numeric(rowMeans(gamma_mat))
      } else if (!is.null(gamma_fixed)) {
        rep(mean(gamma_fixed), n_draws)
      } else {
        NULL
      }
      if (!is.null(gamma_rep)) {
        arr_pred   <- stats::plogis(gamma_rep + mu_new_draws) - stats::plogis(gamma_rep)
        pred_draws <- tibble::tibble(
          Author               = "Prediction",
          b_Intercept          = arr_pred,
          r_Author             = arr_pred - pooled_effect_draws,
          sd_Author__Intercept = tau_draws,
          .chain               = chain_col,
          .iteration           = iteration_col,
          .draw                = draw_col
        )
      }
    } else {
      pred_draws <- tibble::tibble(
        Author               = "Prediction",
        b_Intercept          = mu_new_draws,
        r_Author             = 0,
        sd_Author__Intercept = tau_draws,
        .chain               = chain_col,
        .iteration           = iteration_col,
        .draw                = draw_col
      )
    }
  }

  dplyr::bind_rows(study_draws, pooled_draws, pred_draws)
}


#' Refit a bayesma model on new (subset) data
#'
#' Uses the stored call_args from the original fit to re-run bayesma()
#' with different data.
#'
#' @param x A bayesma object
#' @param newdata A data frame (subset of original)
#' @return A new bayesma object
#'
#' @noRd
refit_bayesma <- function(x, newdata) {
  if (is.null(x$meta$call_args)) {
    cli::cli_abort(
      "Cannot refit: {.cls bayesma} object does not contain stored call arguments.",
      "i" = "Refit your model with the latest version of bayesma to store call args."
    )
  }

  args <- x$meta$call_args
  args$data <- newdata
  do.call(bayesma::bayesma, args)
}


# Internal: Emit Stan target += statement for a prior

#' @noRd
emit_prior_target <- function(prior, par_name) {
  if (is.null(prior)) return("")
  out <- switch(prior$family,
                normal = glue::glue(
                  "target += normal_lpdf({par_name} | {prior$mean}, {prior$sd});"),
                half_normal = glue::glue(
                  "target += normal_lpdf({par_name} | {prior$mean}, {prior$sd});"),
                half_cauchy = glue::glue(
                  "target += cauchy_lpdf({par_name} | {prior$location}, {prior$scale});"),
                half_student_t = glue::glue(
                  "target += student_t_lpdf({par_name} | {prior$df}, {prior$location}, {prior$scale});"),
                exponential = glue::glue(
                  "target += exponential_lpdf({par_name} | {prior$rate});"),
                beta = glue::glue(
                  "target += beta_lpdf({par_name} | {prior$alpha}, {prior$beta});"),
                scaled_inv_chi_sq = glue::glue(
                  "target += scaled_inv_chi_square_lpdf({par_name} | {prior$df}, {prior$scale});"),
                uniform = "",
                dirichlet = "",
                cli::cli_abort("Unknown prior family: {.val {prior$family}}")
  )
  as.character(out)
}

# Determine Stan parameter bounds string for a prior

#' @noRd
emit_prior_bounds <- function(prior, default_lower = NULL, default_upper = NULL) {
  if (is.null(prior)) {
    parts <- character()
    if (!is.null(default_lower)) parts <- c(parts, glue::glue("lower={default_lower}"))
    if (!is.null(default_upper)) parts <- c(parts, glue::glue("upper={default_upper}"))
    if (length(parts) == 0) return("")
    return(paste0("<", paste(parts, collapse = ", "), ">"))
  }
  if (prior$family == "uniform") {
    return(as.character(glue::glue("<lower={prior$lower}, upper={prior$upper}>")))
  }
  parts <- character()
  if (!is.null(default_lower)) parts <- c(parts, glue::glue("lower={default_lower}"))
  if (!is.null(default_upper)) parts <- c(parts, glue::glue("upper={default_upper}"))
  if (length(parts) == 0) return("")
  paste0("<", paste(parts, collapse = ", "), ">")
}



# Validate prior arguments

#' @noRd
# validate_prior_args() has been moved to validation.R



# Resolve priors (fill defaults for NULLs)

#' @noRd
resolve_priors <- function(stage, likelihood, model_type, re_dist,
                           mu_prior, tau_prior, gamma_prior,
                           nu_prior, alpha_prior, mixture_priors,
                           b_prior, p_bias_prior, w_bias_prior,
                           selection_priors,
                           mu_beta_prior = NULL,
                           tau_beta_prior = NULL) {

  is_re   <- model_type == "random_effect"
  is_bias <- model_type == "bias_corrected"
  is_bc_bnp <- model_type == "bc_bnp"
  is_selection <- model_type %in% c("selection_copas", "selection_weight", "pet_peese")

  if (is.null(mu_prior)) {
    mu_prior <- if (is_bias) {
      normal(0, 1)
    } else if (is_bc_bnp) {
      normal(0, 100)
    } else if (stage == "one_stage" && likelihood == "gaussian") {
      normal(0, 100)
    } else {
      normal(0, 10)
    }
  }

  if (is.null(tau_prior) && (is_re || is_bias || is_bc_bnp || is_selection)) {
    tau_prior <- if (is_bias) {
      scaled_inv_chi_sq(1, 0.5)
    } else if (is_bc_bnp || model_type == "selection_copas") {
      half_cauchy(0, 1)
    } else if (stage == "one_stage" && likelihood == "gaussian") {
      half_student_t(3, 0, 10)
    } else {
      half_student_t(3, 0, 2.5)
    }
  }

  if (is.null(gamma_prior) && stage == "one_stage") {
    gamma_prior <- if (likelihood == "gaussian") normal(0, 100)
    else normal(0, 10)
  }

  if (is.null(nu_prior) && (re_dist == "t" || model_type == "selection_copas")) {
    nu_prior <- exponential(0.1)
  }

  if (is.null(alpha_prior) && re_dist == "skew_normal") {
    alpha_prior <- normal(0, 5)
  }

  if (re_dist == "mixture" || model_type == "mixture_model") {
    defaults <- if (model_type == "mixture_model") {
      list(w = dirichlet(1), mu_k = normal(0, 10),
           tau_k = half_normal(0, 0.5))
    } else if (stage == "two_stage") {
      list(w = dirichlet(1), mu_k = normal(0, 10),
           tau_k = half_student_t(3, 0, 2.5))
    } else {
      list(w = dirichlet(1), delta_k = normal(0, 2.5),
           tau_k = half_student_t(3, 0, 2.5))
    }
    if (is.null(mixture_priors)) {
      mixture_priors <- defaults
    } else {
      for (nm in names(defaults)) {
        if (is.null(mixture_priors[[nm]])) mixture_priors[[nm]] <- defaults[[nm]]
      }
    }
  }

  # Bias-corrected defaults
  if (is_bias) {
    if (is.null(b_prior))      b_prior      <- uniform(0, 2)
    if (is.null(p_bias_prior)) p_bias_prior <- beta(1, 1)
    if (is.null(w_bias_prior)) w_bias_prior <- beta(0.5, 1)
  }

  if (is_bc_bnp) {
    if (is.null(p_bias_prior))  p_bias_prior  <- beta(0.5, 1)
    if (is.null(mu_beta_prior)) mu_beta_prior <- uniform(-15, 15)
    if (is.null(tau_beta_prior)) tau_beta_prior <- half_cauchy(0, 1)
  }

  if (is.null(selection_priors)) selection_priors <- list()

  if (model_type == "selection_copas") {
    if (is.null(selection_priors$gamma0)) selection_priors$gamma0 <- uniform(-2, 2)
    if (is.null(selection_priors$rho))    selection_priors$rho    <- uniform(-1, 1)
  }

  list(
    mu = mu_prior, tau = tau_prior, gamma = gamma_prior,
    nu = nu_prior, alpha = alpha_prior, mixture = mixture_priors,
    b = b_prior, p_bias = p_bias_prior, w_bias = w_bias_prior,
    mu_beta = mu_beta_prior, tau_beta = tau_beta_prior,
    selection = selection_priors
  )
}



# Generate robust mixture Stan code fragments

#' @noRd
emit_robust_parameters <- function(robust_config) {
  if (!robust_config$enabled) return(list(par = "", tp = "", data = ""))

  w_bnds <- emit_prior_bounds(robust_config$weight,
                              default_lower = 0, default_upper = 1)

  par <- glue::glue("  real{w_bnds} pi_main;  // mixing weight for main component")

  list(
    par = as.character(par),
    tp = "",
    data = ""
  )
}

#' @noRd
emit_robust_priors <- function(robust_config) {
  if (!robust_config$enabled) return("")

  w_tgt <- emit_prior_target(robust_config$weight, "pi_main")
  as.character(w_tgt)
}

#' @noRd
emit_robust_likelihood <- function(main_ll, y_expr, se_expr, robust_config) {
  if (!robust_config$enabled) return(main_ll)

  mu_out <- robust_config$prior$mean

  if (is.finite(robust_config$df)) {
    outlier_ll <- glue::glue(
      "student_t_lpdf({y_expr} | {robust_config$df}, {mu_out}, {se_expr} * 3)"
    )
  } else {
    outlier_ll <- glue::glue(
      "normal_lpdf({y_expr} | {mu_out}, {se_expr} * 3)"
    )
  }

  as.character(glue::glue(
    "log_sum_exp(log(pi_main) + {main_ll}, log1m(pi_main) + {outlier_ll})"
  ))
}



# Internal: Compute study-level effect sizes

#' @noRd
compute_effect_sizes <- function(outcome_ctrl, outcome_int, n_c, n_i,
                                 sd_c, sd_i, S, likelihood) {
  switch(likelihood,
         binomial = {
           a <- outcome_int;  b <- n_i - outcome_int
           c <- outcome_ctrl; d <- n_c - outcome_ctrl
           # Apply 0.5 continuity correction only to studies with zero cells
           cc <- dplyr::if_else(a == 0 | b == 0 | c == 0 | d == 0, 0.5, 0)
           yi  <- log(((a + cc) * (d + cc)) / ((b + cc) * (c + cc)))
           sei <- sqrt(1 / (a + cc) + 1 / (b + cc) +
                         1 / (c + cc) + 1 / (d + cc))
           list(yi = yi, sei = sei, measure = "log_or")
         },
         poisson = {
           yi  <- log((outcome_int / n_i) / (outcome_ctrl / n_c))
           sei <- sqrt(1 / outcome_int + 1 / outcome_ctrl)
           list(yi = yi, sei = sei, measure = "log_rr")
         },
         gaussian = {
           yi  <- outcome_int - outcome_ctrl
           sp  <- sqrt(((n_c - 1) * sd_c^2 + (n_i - 1) * sd_i^2) /
                         (n_c + n_i - 2))
           sei <- sp * sqrt(1 / n_c + 1 / n_i)
           list(yi = yi, sei = sei, measure = "mean_diff")
         }
  )
}



# Internal function to fit two-stage model

#' @noRd
fit_two_stage <- function(outcome_ctrl, outcome_int, n_c, n_i, sd_c, sd_i, S,
                          study_labels, likelihood, model_type,
                          re_dist, small_sample, priors, n_components,
                          robust_config = list(enabled = FALSE),
                          chains, iter_warmup, iter_sampling,
                          adapt_delta, seed, ...) {

  is_re <- model_type == "random_effect"
  use_t_likelihood <- small_sample %in% c("t_approx", "hjsk")

  es <- compute_effect_sizes(outcome_ctrl, outcome_int, n_c, n_i,
                             sd_c, sd_i, S, likelihood)

  stan_data <- list(S = S, y = es$yi, se = es$sei)

  if (use_t_likelihood) {
    n_total <- as.integer(n_c + n_i)
    if (any(n_total <= 2))
      cli::cli_abort("Student-t adjustment requires total sample size > 2.")
    stan_data$df <- as.numeric(n_total - 2)
  }

  if (re_dist == "mixture" && is_re) {
    stan_data$K <- as.integer(n_components)
    stan_data$prior_dirichlet_alpha <- rep(priors$mixture$w$alpha, n_components)
  }

  stan_code <- generate_stan_code_two_stage(
    model_type = model_type, re_dist = re_dist,
    use_t_likelihood = use_t_likelihood, priors = priors,
    robust_config = robust_config
  )

  mod <- get_cmdstan_model_cached(stan_code)

  fit <- mod$sample(
    data = stan_data, chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
  )

  build_output(
    fit = fit, likelihood = likelihood, model_type = model_type,
    re_dist = re_dist, study_labels = study_labels,
    stan_code = stan_code, stan_data = stan_data,
    arm_data = tibble::tibble(study = study_labels,
                              yi = es$yi, sei = es$sei),
    priors = priors, stage = "two_stage", es = es,
    n_components = n_components, robust_config = robust_config
  )
}



# Internal function to fit one-stage model

#' @noRd
fit_one_stage <- function(outcome_ctrl, outcome_int, n_c, n_i, sd_c, sd_i, S,
                          study_labels, likelihood, model_type,
                          re_dist, priors, n_components,
                          robust_config = list(enabled = FALSE),
                          multi_arm_config = list(enabled = FALSE),
                          chains, iter_warmup, iter_sampling,
                          adapt_delta, seed, ...) {

  is_re <- model_type == "random_effect"
  has_multi_arm <- isTRUE(multi_arm_config$enabled)

  arm_data <- tibble::tibble(
    study_id = rep(seq_len(S), times = 2),
    treat    = rep(c(0L, 1L), each = S),
    outcome  = c(outcome_ctrl, outcome_int),
    n        = c(n_c, n_i)
  )

  if (likelihood == "gaussian") {
    arm_data <- dplyr::mutate(arm_data,
                              sd = c(sd_c, sd_i), se = .data$sd / sqrt(.data$n)
    )
  }

  # Add multi-arm study identifiers to arm_data
  if (has_multi_arm) {
    arm_data$ma_study_id <- rep(multi_arm_config$ma_study_id, times = 2)
  }

  stan_data <- list(
    N = nrow(arm_data), S = S,
    treat = arm_data$treat, study = arm_data$study_id
  )

  if (likelihood == "binomial") {
    stan_data$events <- as.integer(arm_data$outcome)
    stan_data$n      <- as.integer(arm_data$n)
  } else if (likelihood == "gaussian") {
    stan_data$y  <- arm_data$outcome
    stan_data$se <- arm_data$se
  } else if (likelihood == "poisson") {
    stan_data$events   <- as.integer(arm_data$outcome)
    stan_data$exposure <- as.numeric(arm_data$n)
  }

  if (re_dist == "mixture" && is_re) {
    stan_data$K <- as.integer(n_components)
    stan_data$prior_dirichlet_alpha <- rep(priors$mixture$w$alpha, n_components)
  }

  # Add multi-arm data to stan_data
  if (has_multi_arm) {
    stan_data$n_ma_studies <- multi_arm_config$n_ma_studies
    stan_data$comp_to_ma <- multi_arm_config$ma_study_id
  }

  stan_code <- generate_stan_code_one_stage(
    likelihood = likelihood, model_type = model_type,
    re_dist = re_dist, priors = priors,
    robust_config = robust_config,
    multi_arm_config = multi_arm_config
  )

  mod <- get_cmdstan_model_cached(stan_code)

  fit <- mod$sample(
    data = stan_data, chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
  )

  es <- compute_effect_sizes(
    outcome_ctrl = arm_data$outcome[seq_len(S)],
    outcome_int  = arm_data$outcome[seq_len(S) + S],
    n_c = arm_data$n[seq_len(S)], n_i = arm_data$n[seq_len(S) + S],
    sd_c = sd_c, sd_i = sd_i, S = S, likelihood = likelihood
  )

  build_output(
    fit = fit, likelihood = likelihood, model_type = model_type,
    re_dist = re_dist, study_labels = study_labels,
    stan_code = stan_code, stan_data = stan_data,
    arm_data = arm_data, priors = priors,
    stage = "one_stage", es = es, n_components = n_components,
    robust_config = robust_config,
    multi_arm_config = multi_arm_config
  )
}

# Internal function to fit bias-corrected model (Jung)

#' @noRd
fit_bias_corrected <- function(yi, sei, S, study_labels, priors,
                               use_known_bias, known_bias,
                               chains, iter_warmup, iter_sampling,
                               adapt_delta, seed, ...) {

  p <- priors

  stan_data <- list(
    N              = S,
    y              = yi,
    se_y           = sei,
    use_known_bias = as.integer(use_known_bias),
    known_bias     = known_bias
  )

  stan_code <- generate_stan_code_bias_corrected(priors)

  mod <- get_cmdstan_model_cached(stan_code)

  fit <- mod$sample(
    data = stan_data, chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
  )

  # ---- Build output ----
  key_vars <- c("mu", "B", "tau", "inv_var", "p_bias")
  draw_vars <- c(key_vars, "mu_biased",
                 paste0("prob_biased[", seq_len(S), "]"))

  build_output_simple(
    fit = fit, key_vars = key_vars, draw_vars = draw_vars,
    yi = yi, sei = sei, study_labels = study_labels,
    pooled_label = "Pooled (unbiased)",
    stan_code = stan_code, stan_data = stan_data,
    extra_meta = list(use_known_bias = use_known_bias),
    extra_study_cols = list(
      prob_biased = extract_per_study_medians(fit, "prob_biased", S)
    )
  )
}


# Internal functin to fit weight function selection model (Vevea & Hedges 1995)

#' @noRd
fit_selection_weight <- function(yi, sei, S, study_labels, priors,
                                 p_cutoffs, chains, iter_warmup,
                                 iter_sampling, adapt_delta, seed, ...) {
  p <- priors

  K <- length(p_cutoffs) + 1L

  sp <- p$selection %||% list()
  if (is.null(sp$omega)) sp$omega <- dirichlet(rep(1, K))

  mu_tgt   <- emit_prior_target(p$mu, "mu")
  tau_tgt  <- emit_prior_target(p$tau, "tau")
  tau_bnds <- emit_prior_bounds(p$tau, default_lower = 0)

  stan_code <- paste0(
    "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=2> K;
  vector[K-1] p_cutoffs;
}
transformed data {
  // Pre-compute z-score boundaries for each p-value interval.
  // Clamp to avoid +/-Inf from inv_Phi at 0 or 1.
  vector[K+1] z_bounds;
  // z_bounds[1] corresponds to p_lower = 0 -> z_upper = inv_Phi(1) = +Inf
  // z_bounds[K+1] corresponds to p_upper = 1 -> z_lower = inv_Phi(0) = -Inf
  // We clamp to +/- 8 (beyond which Phi is effectively 0 or 1)
  z_bounds[1] = 8.0;      // inv_Phi(1 - 0) clamped
  z_bounds[K+1] = -8.0;   // inv_Phi(1 - 1) clamped
  for (k in 1:(K-1))
    z_bounds[k+1] = inv_Phi(1.0 - p_cutoffs[k]);
}
parameters {
  real mu;
  real", tau_bnds, " tau;
  vector<lower=0.01, upper=0.99>[K-1] omega_raw;
}
transformed parameters {
  vector[K] omega;
  vector[N] sigma;
  omega[1] = 1.0;
  for (k in 1:(K-1))
    omega[k+1] = omega_raw[k];
  for (i in 1:N)
    sigma[i] = sqrt(square(tau) + square(se[i]));
}
model {
  ", mu_tgt, "
  ", tau_tgt, "
  // Prior on weights: beta(1,1) = uniform on each
  for (k in 1:(K-1))
    target += beta_lpdf(omega_raw[k] | 1, 1);
  // Likelihood with weight function
  for (i in 1:N) {
    real z_i = y[i] / se[i];
    // One-sided p-value (assuming positive effects favoured)
    real p_val = 1.0 - Phi(z_i);
    // Determine weight interval
    real log_w;
    if (p_val < p_cutoffs[1]) {
      log_w = 0.0;  // log(omega[1]) = log(1) = 0
    } else if (K == 2) {
      log_w = log(omega[2]);
    } else {
      log_w = log(omega[K]);
      for (k in 1:(K-2)) {
        if (p_val >= p_cutoffs[k] && p_val < p_cutoffs[k+1])
          log_w = log(omega[k+1]);
      }
    }
    target += normal_lpdf(y[i] | mu, sigma[i]) + log_w;
    // Normalisation: sum of omega_k * P(y_i in interval k)
    {
      real norm_c = 0;
      for (k in 1:K) {
        // z_bounds[k] is the upper z boundary for interval k
        // z_bounds[k+1] is the lower z boundary for interval k
        real prob_k = Phi((z_bounds[k] * se[i] - mu) / sigma[i])
                    - Phi((z_bounds[k+1] * se[i] - mu) / sigma[i]);
        norm_c += omega[k] * fmax(prob_k, 1e-15);
      }
      target += -log(fmax(norm_c, 1e-15));
    }
  }
}
generated quantities {
  real pooled = mu;
  vector[K] weights = omega;
  real mu_new = normal_rng(mu, tau);
}")

  stan_data <- list(
    N = S, y = yi, se = sei,
    K = K, p_cutoffs = sort(p_cutoffs)
  )

  mod <- get_cmdstan_model_cached(stan_code)

  fit <- mod$sample(
    data = stan_data, chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
  )

  key_vars <- c("mu", "tau", paste0("omega[", seq_len(K), "]"))

  build_output_simple(
    fit = fit, key_vars = key_vars, draw_vars = key_vars,
    yi = yi, sei = sei, study_labels = study_labels,
    pooled_label = "Pooled (weight-adjusted)",
    stan_code = stan_code, stan_data = stan_data,
    extra_meta = list(p_cutoffs = p_cutoffs)
  )
}


# Internal: Fit PET-PEESE model

#' @noRd
fit_pet_peese <- function(yi, sei, S, study_labels, priors,
                          n_total,
                          chains, iter_warmup,
                          iter_sampling, adapt_delta, seed, ...) {
  p <- priors

  mu_tgt <- emit_prior_target(p$mu, "mu")

  sp <- p$selection %||% list()
  if (is.null(sp$beta_bias)) sp$beta_bias <- normal(0, 5)
  beta_tgt <- emit_prior_target(sp$beta_bias, "beta_bias")

  stan_data <- list(N = S, y = yi, se = sei,
                    n_total = as.array(n_total))

  build_stan <- function(predictor) {
    paste0(
      "data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] n_total;
}
transformed data {
  vector[N] inv_sqrt_n;
  vector[N] inv_n;
  for (i in 1:N) {
    inv_sqrt_n[i] = 1.0 / sqrt(n_total[i]);
    inv_n[i] = 1.0 / n_total[i];
  }
}
parameters {
  real mu;
  real beta_bias;
}
model {
  ", mu_tgt, "
  ", beta_tgt, "
  for (i in 1:N)
    target += normal_lpdf(y[i] | mu + beta_bias * ", predictor, "[i], se[i]);
}
generated quantities {
  real pooled = mu;
  real bias_slope = beta_bias;
}")
  }

# Compile and fit
run_model <- function(stan_code) {
  mod <- get_cmdstan_model_cached(stan_code)
  mod$sample(
    data = stan_data, chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
  )
}

# ---- Step 1: Fit PET ----
pet_code <- build_stan("inv_sqrt_n")
fit_pet  <- run_model(pet_code)

pet_draws <- as.vector(
  posterior::subset_draws(fit_pet$draws("mu"), variable = "mu")
)

# ---- Step 2: Posterior probability switching rule ----
# Determine direction from the data (sign of mean effect)
direction <- if (mean(yi) < 0) "negative" else "positive"
threshold <- 0.10

# P(effect in expected direction) per study-equivalent
prob_effect <- if (direction == "negative") {
  mean(pet_draws < 0)
} else {
  mean(pet_draws > 0)
}

# If strong evidence of effect from PET -> use PEESE for final estimate
use_peese <- prob_effect > (1 - threshold)

# ---- Step 3: Fit chosen model (or PEESE if switching) ----
if (use_peese) {
  peese_code <- build_stan("inv_n")
  fit_final  <- run_model(peese_code)
  recommended <- "peese"
  stan_code <- peese_code
} else {
  fit_final <- fit_pet
  recommended <- "pet"
  stan_code <- pet_code
}

  # ---- Build output ----
  key_vars <- c("mu", "beta_bias")

  build_output_simple(
    fit = fit_final, key_vars = key_vars, draw_vars = key_vars,
    yi = yi, sei = sei, study_labels = study_labels,
    pooled_label = paste0("Pooled (", toupper(recommended), "-adjusted)"),
    stan_code = stan_code, stan_data = stan_data,
    extra_meta = list(
      recommended = recommended,
      prob_effect = prob_effect,
      direction   = direction,
      threshold   = threshold,
      fit_pet     = fit_pet
    )
  )
}



# Helper: extract per-study posterior medians for an array variable

#' @noRd
extract_per_study_medians <- function(fit, var_prefix, S) {
  purrr::map_dbl(seq_len(S), function(i) {
    vn <- paste0(var_prefix, "[", i, "]")
    stats::median(as.vector(
      posterior::subset_draws(fit$draws(vn), variable = vn)
    ))
  })
}


# Build standardised output for bias/selection models

#' Shared output assembly for fit_bias_corrected, fit_selection_copas,
#' fit_selection_weight, and fit_pet_peese.
#'
#' @param fit CmdStanMCMC fit object
#' @param key_vars Character vector of parameter names to summarize
#' @param draw_vars Character vector of parameter names for draws extraction
#' @param yi Numeric vector of study-level effect sizes
#' @param sei Numeric vector of study-level standard errors
#' @param study_labels Character vector of study names
#' @param pooled_label Character label for the pooled row
#' @param stan_code Character. The Stan code used
#' @param stan_data List. The Stan data used
#' @param extra_meta Named list. Additional metadata for meta
#' @param extra_study_cols Named list of vectors to add as columns to study_rows
#' @return A bayesma object (list of class "bayesma")
#' @noRd
build_output_simple <- function(fit, key_vars, draw_vars,
                                yi, sei, study_labels,
                                pooled_label, stan_code, stan_data,
                                extra_meta = list(),
                                extra_study_cols = list()) {

  summary_tbl <- fit$summary(variables = key_vars) |> tibble::as_tibble()

  draws <- posterior::as_draws_df(fit$draws(variables = draw_vars))

  mu_draws <- as.vector(
    posterior::subset_draws(fit$draws("mu"), variable = "mu")
  )

  pooled_row <- tibble::tibble(
    study = pooled_label,
    estimate = stats::median(mu_draws),
    lower = stats::quantile(mu_draws, 0.025),
    upper = stats::quantile(mu_draws, 0.975),
    type = "pooled"
  )

  study_rows <- tibble::tibble(
    study = study_labels,
    estimate = yi,
    lower = yi - 1.96 * sei,
    upper = yi + 1.96 * sei,
    type = "study"
  )

  for (nm in names(extra_study_cols)) {
    study_rows[[nm]] <- extra_study_cols[[nm]]
  }

  forest_df <- dplyr::bind_rows(study_rows, pooled_row) |>
    dplyr::mutate(study = forcats::fct_inorder(.data$study))

  meta <- c(list(study_labels = study_labels), extra_meta)

  out <- list(
    fit = fit, summary = summary_tbl, forest_df = forest_df,
    draws = draws, stan_code = stan_code, stan_data = stan_data,
    meta = meta
  )
  class(out) <- "bayesma"
  out
}


# Function to build tidy output object

#' @noRd
build_output <- function(fit, likelihood, model_type, re_dist,
                         study_labels, stan_code, stan_data, arm_data,
                         priors, stage, es, n_components = 2L,
                         robust_config = list(enabled = FALSE),
                         multi_arm_config = list(enabled = FALSE)) {

  is_re <- model_type == "random_effect"
  use_robust <- isTRUE(robust_config$enabled)
  use_multi_arm <- isTRUE(multi_arm_config$enabled)
  S <- length(study_labels)

  effect_label <- switch(likelihood,
                         binomial = "log_or", gaussian = "mean_diff", poisson = "log_rr"
  )

  # Key parameters
  key_vars <- "mu"
  if (is_re) {
    key_vars <- switch(re_dist,
                       normal      = c(key_vars, "tau"),
                       t           = c(key_vars, "tau", "nu"),
                       skew_normal = c(key_vars, "tau", "alpha_skew"),
                       mixture     = {
                         if (stage == "two_stage") {
                           c(key_vars,
                             paste0("mu_k[", seq_len(n_components), "]"),
                             paste0("tau_k[", seq_len(n_components), "]"),
                             paste0("w[", seq_len(n_components), "]"))
                         } else {
                           c(key_vars, "tau",
                             paste0("delta_k[", seq_len(n_components), "]"),
                             paste0("tau_k[", seq_len(n_components), "]"),
                             paste0("w[", seq_len(n_components), "]"))
                         }
                       }
    )
    # Add rho for multi-arm models
    if (use_multi_arm && re_dist != "mixture") {
      key_vars <- c(key_vars, "rho")
    }
  }
  if (use_robust) key_vars <- c(key_vars, "pi_main")

  summary_tbl <- fit$summary(variables = key_vars) |> tibble::as_tibble()

  # Draws
  draw_vars <- key_vars
  if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
    study_pars <- if (stage == "two_stage") {
      paste0("theta[", seq_len(S), "]")
    } else {
      paste0("epsilon[", seq_len(S), "]")
    }
    draw_vars <- c(draw_vars, study_pars)
  } else if (is_re && re_dist == "mixture" && stage == "two_stage") {
    draw_vars <- c(draw_vars, paste0("cluster[", seq_len(S), "]"))
  }
  if (use_robust && is_re && re_dist %in% c("normal", "t", "skew_normal")) {
    draw_vars <- c(draw_vars,
                   paste0("prob_outlier[", seq_len(S), "]"))
  }
  # Include mu_new for prediction intervals (RE models only)
  if (is_re) {
    mu_new_available <- tryCatch({
      fit$draws("mu_new")
      TRUE
    }, error = function(e) FALSE)
    if (mu_new_available) {
      draw_vars <- c(draw_vars, "mu_new")
    }
  }
  draws <- posterior::as_draws_df(fit$draws(variables = draw_vars))

  # Pooled row
  pooled_draws <- as.vector(
    posterior::subset_draws(fit$draws("mu"), variable = "mu")
  )
  pooled_row <- tibble::tibble(
    study = "Pooled",
    estimate = stats::median(pooled_draws),
    lower = stats::quantile(pooled_draws, 0.025),
    upper = stats::quantile(pooled_draws, 0.975),
    type = "pooled"
  )

  # Study-level rows
  if (is_re && re_dist %in% c("normal", "t", "skew_normal")) {
    if (stage == "two_stage") {
      study_rows <- purrr::map(seq_len(S), function(i) {
        vn <- paste0("theta[", i, "]")
        d <- as.vector(
          posterior::subset_draws(fit$draws(vn), variable = vn)
        )
        tibble::tibble(study = study_labels[i],
                       estimate = stats::median(d),
                       lower = stats::quantile(d, 0.025),
                       upper = stats::quantile(d, 0.975),
                       type = "study")
      }) |> purrr::list_rbind()
    } else {
      study_rows <- purrr::map(seq_len(S), function(i) {
        vn <- paste0("epsilon[", i, "]")
        eps <- as.vector(
          posterior::subset_draws(fit$draws(vn), variable = vn)
        )
        eff <- pooled_draws + eps
        tibble::tibble(study = study_labels[i],
                       estimate = stats::median(eff),
                       lower = stats::quantile(eff, 0.025),
                       upper = stats::quantile(eff, 0.975),
                       type = "study")
      }) |> purrr::list_rbind()
    }
  } else {
    study_rows <- tibble::tibble(
      study = study_labels, estimate = es$yi,
      lower = es$yi - 1.96 * es$sei,
      upper = es$yi + 1.96 * es$sei,
      type = "study"
    )
  }

  # Add outlier probabilities if robust
  if (use_robust && is_re &&
      re_dist %in% c("normal", "t", "skew_normal")) {
    study_rows$prob_outlier <- extract_per_study_medians(fit, "prob_outlier", S)
  }

  forest_df <- dplyr::bind_rows(study_rows, pooled_row) |>
    dplyr::mutate(study = forcats::fct_inorder(.data$study),
                  effect_scale = effect_label)

  # Prediction interval
  pred_interval <- NULL
  if (is_re) {
    tryCatch({
      mn <- as.vector(
        posterior::subset_draws(fit$draws("mu_new"),
                                variable = "mu_new")
      )
      pred_interval <- tibble::tibble(
        estimate = stats::median(mn),
        lower = stats::quantile(mn, 0.025),
        upper = stats::quantile(mn, 0.975)
      )
    }, error = function(e) NULL)
  }

  out <- list(
    fit = fit, summary = summary_tbl, forest_df = forest_df,
    draws = draws, pred_interval = pred_interval,
    stan_code = stan_code, stan_data = stan_data, arm_data = arm_data,
    meta = list(
      likelihood = likelihood, model_type = model_type,
      re_dist = re_dist, stage = stage,
      study_labels = study_labels,
      priors = priors, effect_label = effect_label,
      robust = use_robust,
      multi_arm = use_multi_arm,
      es = es
    )
  )
  if (use_robust) out$meta$robust_config <- robust_config
  if (use_multi_arm) out$meta$multi_arm_config <- multi_arm_config
  class(out) <- "bayesma"
  out
}

#' Print Stan code from a fitted model
#'
#' @param fit A `CmdStanFit` object.
#' @returns The Stan code string (invisibly).
#' @export
stan_code <- function(fit) {
  base::cat(fit$stan_code)
  base::invisible(fit$stan_code)
}

#' Print a summary for a fitted model represented by a \code{bayesma} object
#'
#' @aliases print.bayesma_summary
#'
#' @param x An object of class \code{bayesma}.
#' @param digits The number of significant digits for printing out the summary;
#'   defaults to 2. Bulk_ESS and Tail_ESS are always rounded to integers.
#' @param ... Additional arguments passed to \code{\link{summary.bayesma}}.
#'
#' @seealso \code{\link{summary.bayesma}}
#'
#' @keywords internal
#' @export
print.bayesma <- function(x, digits = 2, ...) {
  print(summary(x, ...), digits = digits, ...)
}

#' @export
print.bayesma_summary <- function(x, digits = 2, ...) {
  digits <- as.integer(digits)

  # ---- Header ----
  cat(" Likelihood:", x$likelihood, "\n")
  cat("     Stage:", x$stage, "\n")
  cat("     Model:", x$model_type)
  if (!is.null(x$re_dist) && x$model_type == "random_effect") {
    cat(paste0(" (", x$re_dist, " random effects)"))
  }
  cat("\n")
  cat("  Estimand:", x$estimand, "\n")
  cat("    Studies:", x$n_studies, "\n")
  cat("      Draws:", paste0(
    x$chains, " chains, each with iter_sampling = ", x$iter_sampling,
    "; warmup = ", x$iter_warmup, ";\n",
    "             total post-warmup draws = ", x$total_draws
  ), "\n")
  cat("\n")

  # ---- Priors ----
  cat("Priors:\n")
  cat("  mu:", x$priors$mu, "\n")
  if (!is.null(x$priors$tau)) {
    cat("  tau:", x$priors$tau, "\n")
  }
  if (!is.null(x$priors$gamma)) {
    cat("  gamma (baseline):", x$priors$gamma, "\n")
  }
  if (!is.null(x$priors$nu)) {
    cat("  nu (df):", x$priors$nu, "\n")
  }
  if (!is.null(x$priors$alpha)) {
    cat("  alpha (skewness):", x$priors$alpha, "\n")
  }
  cat("\n")

  # ---- Pooled effect ----
  cat("Pooled Effect (", x$effect_label, "):\n", sep = "")
  print_format(x$pooled, digits = digits)
  cat("\n")

  # ---- Heterogeneity ----
  if (!is.null(x$heterogeneity)) {
    cat("Heterogeneity:\n")
    print_format(x$heterogeneity, digits = digits)
    cat("\n")
  }

  # ---- Marginal estimand ----
  if (!is.null(x$marginal)) {
    cat("Marginal Estimand (", x$marginal$estimand, "):\n", sep = "")
    print_format(x$marginal$table, digits = digits)
    cat("\n")
  }

  # ---- Convergence note ----
  cat(
    "Draws were sampled using CmdStan. For each parameter, Bulk_ESS\n",
    "and Tail_ESS are effective sample size measures, and Rhat is the\n",
    "potential scale reduction factor on split chains (at convergence, Rhat = 1).\n",
    sep = ""
  )

  invisible(x)
}

#' Summarise a fitted \code{bayesma} model
#'
#' Collects all relevant posterior summaries into a \code{bayesma_summary}
#' object that \code{print.bayesma_summary} renders in brms-style layout.
#'
#' @param object An object of class \code{bayesma}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{bayesma_summary} (a named list).
#'
#' @keywords internal
#' @export
summary.bayesma <- function(object, ...) {
  meta <- object$meta

  # ---- MCMC diagnostics ----
  fit      <- object$fit
  var_info <- fit$summary(variables = c("mu")) |> tibble::as_tibble()

  # Try to infer iter/warmup/chains from the fit metadata
  md <- tryCatch(fit$metadata(), error = function(e) NULL)
  chains        <- if (!is.null(md)) md$num_chains       else NA_integer_
  iter_warmup   <- if (!is.null(md)) md$iter_warmup      else NA_integer_
  iter_sampling <- if (!is.null(md)) md$iter_sampling    else NA_integer_
  total_draws   <- if (!is.null(md)) chains * iter_sampling else NA_integer_

  # ---- Estimand label ----
  estimand     <- meta$estimand %||% resolve_estimand(NULL, meta$likelihood)
  effect_label <- meta$effect_label %||% switch(
    meta$likelihood,
    binomial = "log-OR",
    poisson  = "log-IRR",
    gaussian = "MD"
  )

  # ---- Priors (format as strings) ----
  p         <- meta$priors
  fmt_prior <- function(pr) if (is.null(pr)) NULL else format(pr)
  priors <- list(
    mu    = fmt_prior(p$mu),
    tau   = fmt_prior(p$tau),
    gamma = fmt_prior(p$gamma),
    nu    = fmt_prior(p$nu),
    alpha = fmt_prior(p$alpha)
  )

  # ---- Pooled effect table ----
  mu_vars  <- fit$summary(variables = "mu") |> tibble::as_tibble()
  pooled   <- build_summary_table(mu_vars, row_label = "mu")

  # ---- Heterogeneity table ----
  is_re        <- meta$model_type == "random_effect"
  heterogeneity <- NULL

  if (is_re) {
    tau_vars <- tryCatch(
      fit$summary(variables = "tau") |> tibble::as_tibble(),
      error = function(e) NULL
    )
    if (!is.null(tau_vars)) {
      heterogeneity <- build_summary_table(tau_vars, row_label = "tau")

      # Add nu (t-dist df) if present
      nu_vars <- tryCatch(
        fit$summary(variables = "nu") |> tibble::as_tibble(),
        error = function(e) NULL
      )
      if (!is.null(nu_vars)) {
        heterogeneity <- rbind(
          heterogeneity,
          build_summary_table(nu_vars, row_label = "nu")
        )
      }

      # Add alpha_skew if present
      as_vars <- tryCatch(
        fit$summary(variables = "alpha_skew") |> tibble::as_tibble(),
        error = function(e) NULL
      )
      if (!is.null(as_vars)) {
        heterogeneity <- rbind(
          heterogeneity,
          build_summary_table(as_vars, row_label = "alpha_skew")
        )
      }
    }
  }

  # ---- Marginal estimand ----
  marginal_out <- NULL
  if (!is.null(object$marginal)) {
    mg <- object$marginal
    marginal_tbl <- matrix(
      c(mg$summary$median,
        mg$summary$upper - mg$summary$lower,   # approximate spread
        mg$summary$lower,
        mg$summary$upper,
        mg$summary$p_gt_0),
      nrow = 1,
      dimnames = list(
        mg$estimand,
        c("Median", "MAD", "Q2.5", "Q97.5", "P(>0)")
      )
    )
    marginal_out <- list(estimand = mg$estimand, table = marginal_tbl)
  }

  structure(
    list(
      likelihood    = meta$likelihood,
      stage         = meta$stage,
      model_type    = meta$model_type,
      re_dist       = meta$re_dist,
      estimand      = estimand,
      effect_label  = effect_label,
      n_studies     = length(meta$study_labels),
      chains        = chains,
      iter_warmup   = iter_warmup,
      iter_sampling = iter_sampling,
      total_draws   = total_draws,
      priors        = priors,
      pooled        = pooled,
      heterogeneity = heterogeneity,
      marginal      = marginal_out
    ),
    class = "bayesma_summary"
  )
}


# ---- Internal helpers ----

#' Build a formatted summary matrix from a cmdstanr summary tibble
#'
#' @param tbl Tibble from \code{fit$summary()}.
#' @param row_label Character label for the row name.
#' @return A 1-row matrix with columns Estimate, Est.Error, Q2.5, Q97.5,
#'   Bulk_ESS, Tail_ESS, Rhat.
#' @noRd
build_summary_table <- function(tbl, row_label) {
  mat <- matrix(
    c(
      tbl$median   %||% tbl$mean,
      tbl$mad      %||% tbl$sd,
      tbl$q2.5     %||% NA_real_,
      tbl$q97.5    %||% NA_real_,
      tbl$ess_bulk %||% NA_real_,
      tbl$ess_tail %||% NA_real_,
      tbl$rhat     %||% NA_real_
    ),
    nrow = 1,
    dimnames = list(
      row_label,
      c("Estimate", "Est.Error", "Q2.5", "Q97.5", "Bulk_ESS", "Tail_ESS", "Rhat")
    )
  )
  mat
}

#' Print a numeric matrix with controlled digits, matching brms style
#'
#' Bulk_ESS and Tail_ESS are always shown as integers.
#'
#' @param x Object coercible to matrix.
#' @param digits Number of decimal places for numeric columns.
#' @param no_digits Column names to show as integers.
#' @noRd
print_format <- function(x, digits = 2,
                         no_digits = c("Bulk_ESS", "Tail_ESS")) {
  x   <- as.matrix(x)
  fmt <- paste0("%.", digits, "f")
  out <- x

  for (i in seq_len(ncol(x))) {
    if (colnames(x)[i] %in% no_digits) {
      out[, i] <- sprintf("%.0f", as.numeric(x[, i]))
    } else {
      out[, i] <- sprintf(fmt, as.numeric(x[, i]))
    }
  }

  print(out, quote = FALSE, right = TRUE)
  invisible(x)
}
