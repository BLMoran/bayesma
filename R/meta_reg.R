# meta_reg modular pipeline
#
# Six user-inspectable stages mirroring bayesma_pipeline.R:
#   1. meta_reg_spec()       -- validate + extract -> meta_reg_spec
#   2. meta_reg_stan_code()  -- spec -> named list of Stan blocks + full program
#   3. meta_reg_stan_data()  -- spec -> cmdstanr data list
#   4. meta_reg_fit()        -- compile + sample
#   5. meta_reg_extract()    -- fit + spec -> tidy effect components
#   6. meta_reg_output()     -- assemble final bayesma_metareg object


# -----------------------------------------------------------------------------
# Stage 1: spec
# -----------------------------------------------------------------------------

#' Build a meta-regression specification object
#'
#' @param data A data frame with one row per study.
#' @param studyvar Character. Column name of the study identifier.
#' @param yi,vi Character. Column names of pre-computed effect sizes and
#'   their sampling variances (two-stage only).
#' @param mods One-sided formula specifying moderators
#'   (e.g. `~ age + dose`).
#' @param event_ctrl,event_int Character. Column names of event counts for
#'   binomial / Poisson likelihoods.
#' @param mean_ctrl,mean_int,sd_ctrl,sd_int Character. Column names of arm
#'   means and SDs for the Gaussian likelihood.
#' @param n_ctrl,n_int Character. Column names of arm sample sizes.
#' @param likelihood Character. One of `"binomial"`, `"gaussian"`,
#'   `"poisson"`.
#' @param model_type Character. `"random_effect"` (default) or
#'   `"common_effect"`.
#' @param stage Character. `"two_stage"` (default) or `"one_stage"`.
#' @param center Logical. Mean-centre continuous moderators. Default `TRUE`.
#' @param scale Logical. Scale continuous moderators to unit SD. Default
#'   `FALSE`.
#' @param small_sample Character. Small-sample adjustment for two-stage
#'   models: `"none"`, `"t_approx"`, or `"hjsk"`.
#' @param mu_prior Prior on the intercept (pooled effect at the reference
#'   moderator values).
#' @param tau_prior Prior on the between-study SD (random-effects models).
#' @param gamma_prior Prior on the Gaussian arm-level intercept
#'   (one-stage only).
#' @param beta_prior Default prior for every regression coefficient.
#' @param beta_priors Named list of coefficient-specific priors, overriding
#'   `beta_prior` for those coefficients.
#' @param custom_model Optional character scalar of Stan code overriding
#'   the generated program.
#' @param custom_data Optional named list merged into the Stan data list.
#' @return An object of class `"meta_reg_spec"`.
#' @export
meta_reg_spec <- function(
    data,
    studyvar,
    yi            = NULL,
    vi            = NULL,
    mods,
    event_ctrl    = NULL,
    event_int     = NULL,
    mean_ctrl     = NULL,
    mean_int      = NULL,
    sd_ctrl       = NULL,
    sd_int        = NULL,
    n_ctrl        = NULL,
    n_int         = NULL,
    likelihood    = c("binomial", "gaussian", "poisson"),
    model_type    = c("random_effect", "common_effect"),
    stage         = c("two_stage", "one_stage"),
    center        = TRUE,
    scale         = FALSE,
    small_sample  = c("none", "t_approx", "hjsk"),
    mu_prior      = NULL,
    tau_prior     = NULL,
    gamma_prior   = NULL,
    beta_prior    = NULL,
    beta_priors   = NULL,
    custom_model  = NULL,
    custom_data   = NULL
) {
  likelihood   <- rlang::arg_match(likelihood)
  model_type   <- rlang::arg_match(model_type)
  stage        <- rlang::arg_match(stage)
  small_sample <- rlang::arg_match(small_sample)

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

  if (missing(mods) || is.null(mods)) {
    cli::cli_abort(c(
      "{.arg mods} is required for meta-regression.",
      "i" = "Use {.code mods = ~ var1 + var2} to specify moderators.",
      "i" = "For meta-analysis without moderators, use {.fn bayesma}."
    ), call = rlang::caller_env())
  }
  if (!inherits(mods, "formula")) {
    cli::cli_abort(
      "{.arg mods} must be a one-sided formula (e.g., {.code ~ var1 + var2}).",
      call = rlang::caller_env()
    )
  }
  if (length(mods) == 3) {
    cli::cli_warn("Ignoring left-hand side of {.arg mods} formula.")
    mods <- mods[-2]
  }

  if (stage == "one_stage" && !is.null(yi)) {
    cli::cli_abort(c(
      "One-stage models require arm-level data, not pre-computed effect sizes.",
      "i" = "Use {.arg stage = 'two_stage'} with {.arg yi} and {.arg vi}."
    ), call = rlang::caller_env())
  }
  if (stage == "one_stage" && small_sample != "none") {
    cli::cli_warn("Small-sample adjustments are two-stage only. Ignoring.")
    small_sample <- "none"
  }
  if (stage == "two_stage" && !is.null(gamma_prior)) {
    cli::cli_warn("{.arg gamma_prior} is only used in one-stage models. Ignoring.")
  }

  extract_col <- function(d, var_name) {
    if (is.null(var_name)) return(NULL)
    val <- d[[var_name]]
    if (is.null(val)) {
      cli::cli_abort("Variable {.val {var_name}} not found in {.arg data}.",
                     call = rlang::caller_env())
    }
    val
  }

  study_vec    <- extract_col(data, studyvar)
  S            <- length(unique(study_vec))
  study_labels <- if (is.factor(study_vec)) levels(study_vec)
  else unique(as.character(study_vec))

  mod_result <- build_moderator_matrix(data = data, mods = mods,
                                       center = center, scale = scale)
  X          <- mod_result$X
  K          <- ncol(X)
  coef_names <- colnames(X)

  is_re <- model_type == "random_effect"
  validate_prior_args_mreg(model_type, stage, likelihood,
                           mu_prior, tau_prior, gamma_prior,
                           beta_prior, beta_priors, coef_names)

  priors <- resolve_priors_mreg(model_type, stage, likelihood,
                                mu_prior, tau_prior, gamma_prior,
                                beta_prior, beta_priors, K, coef_names)

  n_c <- extract_col(data, n_ctrl)
  n_i <- extract_col(data, n_int)

  if (stage == "two_stage" && !is.null(yi) && !is.null(vi)) {
    y_vec  <- extract_col(data, yi)
    v_vec  <- extract_col(data, vi)
    es     <- list(yi = y_vec, sei = sqrt(v_vec), vi = v_vec,
                   measure = "user_provided")
    outcome_ctrl <- NULL; outcome_int <- NULL
    sd_c <- NULL; sd_i <- NULL
  } else {
    if (likelihood %in% c("binomial", "poisson")) {
      outcome_ctrl <- extract_col(data, event_ctrl)
      outcome_int  <- extract_col(data, event_int)
      sd_c <- NULL; sd_i <- NULL
    } else {
      outcome_ctrl <- extract_col(data, mean_ctrl)
      outcome_int  <- extract_col(data, mean_int)
      sd_c <- extract_col(data, sd_ctrl)
      sd_i <- extract_col(data, sd_int)
    }
    es <- compute_effect_sizes(outcome_ctrl, outcome_int, n_c, n_i,
                               sd_c, sd_i, S, likelihood)
    es$vi <- es$sei^2
  }

  call_args <- list(
    studyvar = studyvar, yi = yi, vi = vi, mods = mods,
    event_ctrl = event_ctrl, event_int = event_int,
    mean_ctrl = mean_ctrl, mean_int = mean_int,
    sd_ctrl = sd_ctrl, sd_int = sd_int,
    n_ctrl = n_ctrl, n_int = n_int,
    likelihood = likelihood, model_type = model_type,
    stage = stage, center = center, scale = scale,
    small_sample = small_sample,
    mu_prior = mu_prior, tau_prior = tau_prior,
    gamma_prior = gamma_prior, beta_prior = beta_prior,
    beta_priors = beta_priors,
    custom_model = custom_model, custom_data = custom_data
  )

  spec <- list(
    likelihood   = likelihood,
    model_type   = model_type,
    stage        = stage,
    small_sample = small_sample,
    study_vec    = study_vec,
    study_labels = study_labels,
    S            = S,
    K            = K,
    X            = X,
    coef_names   = coef_names,
    mod_result   = mod_result,
    outcome_ctrl = outcome_ctrl,
    outcome_int  = outcome_int,
    n_c          = n_c,
    n_i          = n_i,
    sd_c         = sd_c,
    sd_i         = sd_i,
    es           = es,
    priors       = priors,
    custom_model = custom_model,
    custom_data  = custom_data,
    call_args    = call_args
  )
  class(spec) <- c("meta_reg_spec", "list")
  spec
}


#' @export
print.meta_reg_spec <- function(x, ...) {
  cat("<meta_reg_spec>\n",
      "  likelihood   : ", x$likelihood, "\n",
      "  model_type   : ", x$model_type, "\n",
      "  stage        : ", x$stage, "\n",
      "  small_sample : ", x$small_sample, "\n",
      "  studies (S)  : ", x$S, "\n",
      "  moderators (K): ", x$K, "\n",
      "  coefficients : ", paste(x$coef_names, collapse = ", "), "\n",
      sep = "")
  invisible(x)
}


# -----------------------------------------------------------------------------
# Stage 2: stan code
# -----------------------------------------------------------------------------

#' Generate Stan code for a meta-regression specification
#'
#' @param spec A `meta_reg_spec` object.
#' @param format Logical. If `TRUE` (default), run the generated program
#'   through Stan's `stanc --auto-format` for consistent indentation and
#'   spacing. Falls back to the raw program if the formatter is unavailable.
#' @return An object of class `"meta_reg_stan_code"`.
#' @export
meta_reg_stan_code <- function(spec, format = TRUE) {
  if (!inherits(spec, "meta_reg_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls meta_reg_spec} object.")
  }

  raw <- if (!is.null(spec$custom_model)) {
    spec$custom_model
  } else if (spec$stage == "two_stage") {
    use_t <- spec$small_sample %in% c("t_approx", "hjsk")
    generate_stan_code_mreg_two_stage(
      model_type       = spec$model_type,
      use_t_likelihood = use_t,
      priors           = spec$priors,
      K                = spec$K
    )
  } else {
    generate_stan_code_mreg_one_stage(
      likelihood = spec$likelihood,
      model_type = spec$model_type,
      priors     = spec$priors,
      K          = spec$K
    )
  }

  full   <- if (isTRUE(format)) format_stan_code(raw) else raw
  blocks <- parse_stan_blocks(full)
  out    <- c(blocks, list(full = full))
  class(out) <- c("meta_reg_stan_code", "bayesma_stan_code", "list")
  out
}


# -----------------------------------------------------------------------------
# Stage 3: stan data
# -----------------------------------------------------------------------------

#' Build the Stan data list for a meta-regression specification
#'
#' @param spec A `meta_reg_spec` object.
#' @return A named list suitable for `cmdstanr::CmdStanModel$sample(data = ...)`.
#' @export
meta_reg_stan_data <- function(spec) {
  if (!inherits(spec, "meta_reg_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls meta_reg_spec} object.")
  }

  sd_list <- if (spec$stage == "two_stage") {
    build_stan_data_mreg_two_stage(spec)
  } else {
    build_stan_data_mreg_one_stage(spec)
  }

  if (!is.null(spec$custom_data)) {
    for (nm in names(spec$custom_data)) sd_list[[nm]] <- spec$custom_data[[nm]]
  }

  sd_list
}


#' @noRd
build_stan_data_mreg_two_stage <- function(spec) {
  sd_list <- list(
    S  = spec$S,
    K  = spec$K,
    y  = spec$es$yi,
    se = spec$es$sei,
    X  = spec$X
  )

  if (spec$small_sample %in% c("t_approx", "hjsk")) {
    n_total <- as.integer(spec$n_c + spec$n_i)
    if (any(n_total <= 2)) {
      cli::cli_abort("Student-t adjustment requires total sample size > 2.")
    }
    sd_list$df <- as.numeric(n_total - 2)
  }

  sd_list
}


#' @noRd
build_stan_data_mreg_one_stage <- function(spec) {
  arm_data <- tibble::tibble(
    study_id = rep(seq_len(spec$S), times = 2),
    treat    = rep(c(0L, 1L), each = spec$S),
    outcome  = c(spec$outcome_ctrl, spec$outcome_int),
    n        = c(spec$n_c, spec$n_i)
  )
  if (spec$likelihood == "gaussian") {
    arm_data <- dplyr::mutate(arm_data,
                              sd = c(spec$sd_c, spec$sd_i), se = .data$sd / sqrt(.data$n)
    )
  }

  X_arm <- rbind(spec$X, spec$X)

  sd_list <- list(
    N     = nrow(arm_data),
    S     = spec$S,
    K     = spec$K,
    treat = arm_data$treat,
    study = arm_data$study_id,
    X     = X_arm
  )

  if (spec$likelihood == "binomial") {
    sd_list$events <- as.integer(arm_data$outcome)
    sd_list$n      <- as.integer(arm_data$n)
  } else if (spec$likelihood == "gaussian") {
    sd_list$y  <- arm_data$outcome
    sd_list$se <- arm_data$se
  } else if (spec$likelihood == "poisson") {
    sd_list$events   <- as.integer(arm_data$outcome)
    sd_list$exposure <- as.numeric(arm_data$n)
  }

  attr(sd_list, "arm_data") <- arm_data
  sd_list
}


# -----------------------------------------------------------------------------
# Stage 4: fit
# -----------------------------------------------------------------------------

#' Compile and sample a meta-regression model
#'
#' @param spec      A `meta_reg_spec` object.
#' @param code      A `meta_reg_stan_code` object. Defaults to
#'   `meta_reg_stan_code(spec)`.
#' @param stan_data A Stan data list. Defaults to `meta_reg_stan_data(spec)`.
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed MCMC settings.
#' @param ... Passed to `cmdstanr::CmdStanModel$sample()`.
#' @return An object of class `"meta_reg_fit"`.
#' @export
meta_reg_fit <- function(spec,
                         code          = meta_reg_stan_code(spec),
                         stan_data     = meta_reg_stan_data(spec),
                         chains        = 4,
                         iter_warmup   = 1000,
                         iter_sampling = 1000,
                         adapt_delta   = 0.95,
                         seed          = 1234,
                         ...) {
  if (!inherits(spec, "meta_reg_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls meta_reg_spec} object.")
  }

  stan_program <- if (inherits(code, "meta_reg_stan_code")) code$full
  else as.character(code)

  mod <- get_cmdstan_model_cached(stan_program)

  fit <- mod$sample(
    data = stan_data, chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
  )

  out <- list(fit = fit, stan_code = code, stan_data = stan_data)
  class(out) <- c("meta_reg_fit", "list")
  out
}


# -----------------------------------------------------------------------------
# Stage 5: extract
# -----------------------------------------------------------------------------

#' Extract tidy effect components from a meta-regression fit
#'
#' @param fit  A `meta_reg_fit` object.
#' @param spec A `meta_reg_spec` object.
#' @return An object of class `"meta_reg_effects"`.
#' @export
meta_reg_extract <- function(fit, spec) {
  if (!inherits(fit, "meta_reg_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls meta_reg_fit} object.")
  }
  if (!inherits(spec, "meta_reg_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls meta_reg_spec} object.")
  }

  cmdstan_fit <- fit$fit
  is_re <- spec$model_type == "random_effect"
  S     <- spec$S
  K     <- spec$K

  effect_label <- switch(spec$likelihood,
                         binomial = "log_or", gaussian = "mean_diff", poisson = "log_rr",
                         spec$es$measure
  )

  # ---- Key parameters ----
  key_vars <- c("mu", paste0("beta[", seq_len(K), "]"))
  if (is_re) key_vars <- c(key_vars, "tau")

  summary_tbl <- cmdstan_fit$summary(variables = key_vars) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      coef_name = c("(Intercept)", spec$coef_names,
                    if (is_re) "tau" else NULL)
    )

  # ---- Draws ----
  draw_vars <- key_vars
  if (is_re) {
    study_pars <- if (spec$stage == "two_stage") {
      paste0("theta[",   seq_len(S), "]")
    } else {
      paste0("epsilon[", seq_len(S), "]")
    }
    draw_vars <- c(draw_vars, study_pars)
    if (tryCatch({ cmdstan_fit$draws("mu_new"); TRUE },
                 error = function(e) FALSE)) {
      draw_vars <- c(draw_vars, "mu_new")
    }
  }

  draws <- posterior::as_draws_df(cmdstan_fit$draws(variables = draw_vars))

  for (k in seq_len(K)) {
    old_nm <- paste0("beta[", k, "]")
    new_nm <- paste0("beta_", spec$coef_names[k])
    if (old_nm %in% names(draws)) {
      names(draws)[names(draws) == old_nm] <- new_nm
    }
  }

  # ---- Coefficient table ----
  mu_draws <- as.vector(
    posterior::subset_draws(cmdstan_fit$draws("mu"), variable = "mu")
  )

  coef_tbl <- tibble::tibble(
    term      = c("(Intercept)", spec$coef_names),
    estimate  = NA_real_,
    std_error = NA_real_,
    q2.5      = NA_real_,
    q97.5     = NA_real_
  )
  coef_tbl$estimate[1]  <- stats::median(mu_draws)
  coef_tbl$std_error[1] <- stats::sd(mu_draws)
  coef_tbl$q2.5[1]      <- stats::quantile(mu_draws, 0.025)
  coef_tbl$q97.5[1]     <- stats::quantile(mu_draws, 0.975)

  for (k in seq_len(K)) {
    vn <- paste0("beta[", k, "]")
    b  <- as.vector(posterior::subset_draws(
      cmdstan_fit$draws(vn), variable = vn
    ))
    coef_tbl$estimate[k + 1]  <- stats::median(b)
    coef_tbl$std_error[k + 1] <- stats::sd(b)
    coef_tbl$q2.5[k + 1]      <- stats::quantile(b, 0.025)
    coef_tbl$q97.5[k + 1]     <- stats::quantile(b, 0.975)
  }

  # ---- Forest data ----
  pooled_row <- tibble::tibble(
    study    = "Pooled",
    estimate = coef_tbl$estimate[1],
    lower    = coef_tbl$q2.5[1],
    upper    = coef_tbl$q97.5[1],
    type     = "pooled"
  )

  if (is_re) {
    var_prefix <- if (spec$stage == "two_stage") "theta" else "epsilon"
    study_rows <- purrr::map(seq_len(S), function(i) {
      vn <- paste0(var_prefix, "[", i, "]")
      d  <- as.vector(posterior::subset_draws(
        cmdstan_fit$draws(vn), variable = vn
      ))
      eff <- if (var_prefix == "theta") d else mu_draws + d
      tibble::tibble(
        study    = spec$study_labels[i],
        estimate = stats::median(eff),
        lower    = stats::quantile(eff, 0.025),
        upper    = stats::quantile(eff, 0.975),
        type     = "study"
      )
    }) |> purrr::list_rbind()
  } else {
    study_rows <- tibble::tibble(
      study    = spec$study_labels,
      estimate = spec$es$yi,
      lower    = spec$es$yi - 1.96 * spec$es$sei,
      upper    = spec$es$yi + 1.96 * spec$es$sei,
      type     = "study"
    )
  }

  forest_df <- dplyr::bind_rows(study_rows, pooled_row) |>
    dplyr::mutate(study        = forcats::fct_inorder(.data$study),
                  effect_scale = effect_label)

  # ---- Tau summary ----
  tau_summary <- NULL
  if (is_re) {
    tau_d <- as.vector(posterior::subset_draws(
      cmdstan_fit$draws("tau"), variable = "tau"
    ))
    tau_summary <- tibble::tibble(
      estimate  = stats::median(tau_d),
      std_error = stats::sd(tau_d),
      q2.5      = stats::quantile(tau_d, 0.025),
      q97.5     = stats::quantile(tau_d, 0.975)
    )
  }

  # ---- Prediction interval ----
  pred_interval <- NULL
  if (is_re) {
    tryCatch({
      mn <- as.vector(posterior::subset_draws(
        cmdstan_fit$draws("mu_new"), variable = "mu_new"
      ))
      pred_interval <<- tibble::tibble(
        estimate = stats::median(mn),
        lower    = stats::quantile(mn, 0.025),
        upper    = stats::quantile(mn, 0.975)
      )
    }, error = function(e) NULL)
  }

  eff <- list(
    summary       = summary_tbl,
    draws         = draws,
    coefficients  = coef_tbl,
    tau           = tau_summary,
    forest_df     = forest_df,
    pred_interval = pred_interval,
    effect_label  = effect_label
  )
  class(eff) <- c("meta_reg_effects", "list")
  eff
}


# -----------------------------------------------------------------------------
# Stage 6: output
# -----------------------------------------------------------------------------

#' Assemble a `bayesma_metareg` object from pipeline outputs
#'
#' @param spec    A `meta_reg_spec` object.
#' @param fit     A `meta_reg_fit` object.
#' @param effects A `meta_reg_effects` object.
#' @return A list of class `c("bayesma_metareg", "bayesma")`.
#' @export
meta_reg_output <- function(spec, fit, effects) {
  if (!inherits(spec,    "meta_reg_spec"))    cli::cli_abort("{.arg spec} must be {.cls meta_reg_spec}.")
  if (!inherits(fit,     "meta_reg_fit"))     cli::cli_abort("{.arg fit} must be {.cls meta_reg_fit}.")
  if (!inherits(effects, "meta_reg_effects")) cli::cli_abort("{.arg effects} must be {.cls meta_reg_effects}.")

  code_full <- if (inherits(fit$stan_code, "meta_reg_stan_code")) fit$stan_code$full
  else as.character(fit$stan_code)

  arm_data <- attr(fit$stan_data, "arm_data")

  out <- list(
    fit           = fit$fit,
    summary       = effects$summary,
    coefficients  = effects$coefficients,
    tau           = effects$tau,
    forest_df     = effects$forest_df,
    draws         = effects$draws,
    pred_interval = effects$pred_interval,
    stan_code     = code_full,
    stan_data     = fit$stan_data,
    meta          = list(
      likelihood   = spec$likelihood,
      model_type   = spec$model_type,
      stage        = spec$stage,
      small_sample = spec$small_sample,
      study_labels = spec$study_labels,
      coef_names   = spec$coef_names,
      priors       = spec$priors,
      effect_label = effects$effect_label,
      es           = spec$es,
      mod_info     = spec$mod_result,
      call_args    = spec$call_args
    )
  )
  if (!is.null(arm_data)) out$arm_data <- arm_data

  class(out) <- c("bayesma_metareg", "bayesma")
  out
}


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

#' Bayesian Meta-Regression
#'
#' Thin orchestrator over the six-stage pipeline:
#' [meta_reg_spec()], [meta_reg_stan_code()], [meta_reg_stan_data()],
#' [meta_reg_fit()], [meta_reg_extract()], [meta_reg_output()].
#'
#' @inheritParams meta_reg_spec
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed MCMC settings.
#' @param return_stage Character. One of `"full"` (default), `"spec"`,
#'   `"code"`, `"data"`, or `"fit"`.
#' @param ... Passed to `cmdstanr::sample()`.
#' @return An object of class `c("bayesma_metareg", "bayesma")`.
#' @export
meta_reg <- function(
    data,
    studyvar,
    yi            = NULL,
    vi            = NULL,
    mods,
    event_ctrl    = NULL,
    event_int     = NULL,
    mean_ctrl     = NULL,
    mean_int      = NULL,
    sd_ctrl       = NULL,
    sd_int        = NULL,
    n_ctrl        = NULL,
    n_int         = NULL,
    likelihood    = c("binomial", "gaussian", "poisson"),
    model_type    = c("random_effect", "common_effect"),
    stage         = c("two_stage", "one_stage"),
    center        = TRUE,
    scale         = FALSE,
    small_sample  = c("none", "t_approx", "hjsk"),
    mu_prior      = NULL,
    tau_prior     = NULL,
    gamma_prior   = NULL,
    beta_prior    = NULL,
    beta_priors   = NULL,
    custom_model  = NULL,
    custom_data   = NULL,
    return_stage  = c("full", "spec", "code", "data", "fit"),
    chains        = 4,
    iter_warmup   = 1000,
    iter_sampling = 1000,
    adapt_delta   = 0.95,
    seed          = 1234,
    ...
) {
  return_stage <- rlang::arg_match(return_stage)

  spec <- meta_reg_spec(
    data = data, studyvar = studyvar, yi = yi, vi = vi, mods = mods,
    event_ctrl = event_ctrl, event_int = event_int,
    mean_ctrl = mean_ctrl, mean_int = mean_int,
    sd_ctrl = sd_ctrl, sd_int = sd_int,
    n_ctrl = n_ctrl, n_int = n_int,
    likelihood = likelihood, model_type = model_type,
    stage = stage, center = center, scale = scale,
    small_sample = small_sample,
    mu_prior = mu_prior, tau_prior = tau_prior,
    gamma_prior = gamma_prior, beta_prior = beta_prior,
    beta_priors = beta_priors,
    custom_model = custom_model, custom_data = custom_data
  )
  if (return_stage == "spec") return(spec)

  code <- meta_reg_stan_code(spec)
  if (return_stage == "code") return(code)

  stan_data <- meta_reg_stan_data(spec)
  if (return_stage == "data") return(stan_data)

  fit <- meta_reg_fit(
    spec = spec, code = code, stan_data = stan_data,
    chains = chains, iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, adapt_delta = adapt_delta,
    seed = seed, ...
  )
  if (return_stage == "fit") return(fit)

  effects <- meta_reg_extract(fit, spec)
  meta_reg_output(spec, fit, effects)
}
