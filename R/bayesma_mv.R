# bayesma_mv modular pipeline
#
# Six stages mirroring bayesma_pipeline.R:
#   1. bayesma_mv_spec()       -- validate + extract -> bayesma_mv_spec
#   2. bayesma_mv_stan_code()  -- spec -> named Stan blocks + full program
#   3. bayesma_mv_stan_data()  -- spec -> cmdstanr data list
#   4. bayesma_mv_fit()        -- compile + sample
#   5. bayesma_mv_extract()    -- fit + spec -> tidy effect components
#   6. bayesma_mv_output()     -- assemble final bayesma_mv object


# -----------------------------------------------------------------------------
# Stage 1: spec
# -----------------------------------------------------------------------------

#' Build a bivariate meta-analysis specification object
#'
#' @param data A data frame with one row per study containing arm-level data
#'   for both outcomes.
#' @param studyvar Column for study identifiers (unquoted).
#' @param mean_ctrl_1,mean_int_1 Columns for control and intervention means
#'   for outcome 1 (unquoted).
#' @param sd_ctrl_1,sd_int_1 Columns for control and intervention SDs for
#'   outcome 1 (unquoted).
#' @param n_ctrl_1,n_int_1 Columns for control and intervention sample sizes
#'   for outcome 1 (unquoted).
#' @param mean_ctrl_2,mean_int_2 Columns for control and intervention means
#'   for outcome 2 (unquoted).
#' @param sd_ctrl_2,sd_int_2 Columns for control and intervention SDs for
#'   outcome 2 (unquoted).
#' @param n_ctrl_2,n_int_2 Columns for control and intervention sample sizes
#'   for outcome 2 (unquoted).
#' @param outcome_labels Character vector of length 2 with labels for the
#'   two outcomes. Default: `c("outcome_1", "outcome_2")`.
#' @param likelihood Character. Currently only `"gaussian"`.
#' @param stage Character. `"two_stage"` (effect sizes computed then
#'   modelled) or `"one_stage"` (marginalised model).
#' @param rho_within Numeric scalar in `[-1, 1]`. Within-study correlation
#'   between the two outcomes, assumed known. Default: `0.5`.
#' @param mu_prior Prior on pooled effects. Either a single prior (applied
#'   to both outcomes) or a named list with elements matching
#'   `outcome_labels`.
#' @param tau_prior Prior on between-study SDs. Either a single prior or a
#'   named list.
#' @param rho_between_prior Prior on the between-study correlation.
#'   Default: `uniform(-1, 1)`.
#' @param custom_model Optional character scalar containing complete Stan
#'   code to override the generated program.
#' @param custom_data Optional named list merged into the Stan data list.
#' @return An object of class `"bayesma_mv_spec"`.
#' @export
bayesma_mv_spec <- function(
    data,
    studyvar,
    mean_ctrl_1, mean_int_1, sd_ctrl_1, sd_int_1, n_ctrl_1, n_int_1,
    mean_ctrl_2, mean_int_2, sd_ctrl_2, sd_int_2, n_ctrl_2, n_int_2,
    outcome_labels    = c("outcome_1", "outcome_2"),
    likelihood        = c("gaussian"),
    stage             = c("two_stage", "one_stage"),
    rho_within        = 0.5,
    mu_prior          = NULL,
    tau_prior         = NULL,
    rho_between_prior = NULL,
    custom_model      = NULL,
    custom_data       = NULL
) {
  likelihood <- rlang::arg_match(likelihood)
  stage      <- rlang::arg_match(stage)

  if (length(outcome_labels) != 2) {
    cli::cli_abort("{.arg outcome_labels} must be a character vector of length 2.",
                   call = rlang::caller_env())
  }
  if (!is.numeric(rho_within) || length(rho_within) != 1 ||
      rho_within < -1 || rho_within > 1) {
    cli::cli_abort("{.arg rho_within} must be a single number in [-1, 1].",
                   call = rlang::caller_env())
  }
  if (!is.null(custom_model) &&
      (!is.character(custom_model) || length(custom_model) != 1)) {
    cli::cli_abort(
      "{.arg custom_model} must be a character scalar containing Stan code.",
      call = rlang::caller_env()
    )
  }
  if (!is.null(custom_data) &&
      (!is.list(custom_data) || is.null(names(custom_data)))) {
    cli::cli_abort(
      "{.arg custom_data} must be a named list.",
      call = rlang::caller_env()
    )
  }

  extract_col <- function(d, quo, label) {
    tryCatch(
      rlang::eval_tidy(quo, d),
      error = function(e) {
        cli::cli_abort("Column {.val {rlang::as_label(quo)}} not found in data.",
                       call = rlang::caller_env())
      }
    )
  }

  studyvar_quo    <- rlang::enquo(studyvar)
  mc1_quo <- rlang::enquo(mean_ctrl_1); mi1_quo <- rlang::enquo(mean_int_1)
  sc1_quo <- rlang::enquo(sd_ctrl_1);   si1_quo <- rlang::enquo(sd_int_1)
  nc1_quo <- rlang::enquo(n_ctrl_1);    ni1_quo <- rlang::enquo(n_int_1)
  mc2_quo <- rlang::enquo(mean_ctrl_2); mi2_quo <- rlang::enquo(mean_int_2)
  sc2_quo <- rlang::enquo(sd_ctrl_2);   si2_quo <- rlang::enquo(sd_int_2)
  nc2_quo <- rlang::enquo(n_ctrl_2);    ni2_quo <- rlang::enquo(n_int_2)

  study_vec <- extract_col(data, studyvar_quo, "studyvar")
  if (!is.factor(study_vec)) study_vec <- factor(study_vec)
  study_labels <- levels(study_vec)
  S <- length(study_labels)

  es_1 <- compute_effect_sizes_mv(
    mean_ctrl = extract_col(data, mc1_quo), mean_int = extract_col(data, mi1_quo),
    sd_ctrl   = extract_col(data, sc1_quo), sd_int   = extract_col(data, si1_quo),
    n_ctrl    = extract_col(data, nc1_quo), n_int    = extract_col(data, ni1_quo),
    likelihood = likelihood
  )
  es_2 <- compute_effect_sizes_mv(
    mean_ctrl = extract_col(data, mc2_quo), mean_int = extract_col(data, mi2_quo),
    sd_ctrl   = extract_col(data, sc2_quo), sd_int   = extract_col(data, si2_quo),
    n_ctrl    = extract_col(data, nc2_quo), n_int    = extract_col(data, ni2_quo),
    likelihood = likelihood
  )

  priors <- resolve_priors_mv(outcome_labels, mu_prior, tau_prior,
                              rho_between_prior)

  col_args <- c("studyvar",
                "mean_ctrl_1", "mean_int_1", "sd_ctrl_1", "sd_int_1",
                "n_ctrl_1", "n_int_1",
                "mean_ctrl_2", "mean_int_2", "sd_ctrl_2", "sd_int_2",
                "n_ctrl_2", "n_int_2")
  quos <- list(studyvar_quo, mc1_quo, mi1_quo, sc1_quo, si1_quo,
               nc1_quo, ni1_quo, mc2_quo, mi2_quo, sc2_quo, si2_quo,
               nc2_quo, ni2_quo)
  call_args <- purrr::map(quos, rlang::as_label) |> stats::setNames(col_args)
  call_args <- c(call_args, list(
    outcome_labels = outcome_labels, likelihood = likelihood,
    stage = stage, rho_within = rho_within,
    mu_prior = mu_prior, tau_prior = tau_prior,
    rho_between_prior = rho_between_prior,
    custom_model = custom_model, custom_data = custom_data
  ))

  spec <- list(
    likelihood     = likelihood,
    stage          = stage,
    study_labels   = study_labels,
    S              = S,
    outcome_labels = outcome_labels,
    rho_within     = rho_within,
    es_1           = es_1,
    es_2           = es_2,
    priors         = priors,
    custom_model   = custom_model,
    custom_data    = custom_data,
    call_args      = call_args
  )
  class(spec) <- c("bayesma_mv_spec", "list")
  spec
}


#' @export
print.bayesma_mv_spec <- function(x, ...) {
  cat("<bayesma_mv_spec>\n",
      "  likelihood : ", x$likelihood, "\n",
      "  stage      : ", x$stage, "\n",
      "  studies (S): ", x$S, "\n",
      "  outcomes   : ", paste(x$outcome_labels, collapse = ", "), "\n",
      "  rho_within : ", x$rho_within, "\n",
      sep = "")
  invisible(x)
}


# -----------------------------------------------------------------------------
# Stage 2: stan code
# -----------------------------------------------------------------------------

#' Generate Stan code for a bivariate meta-analysis specification
#'
#' @param spec A `bayesma_mv_spec` object.
#' @param format Logical. If `TRUE` (default), run the generated program
#'   through Stan's `stanc --auto-format` for consistent indentation and
#'   spacing. Falls back to the raw program if the formatter is unavailable.
#' @return An object of class `"bayesma_mv_stan_code"`.
#' @export
bayesma_mv_stan_code <- function(spec, format = TRUE) {
  if (!inherits(spec, "bayesma_mv_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_mv_spec} object.")
  }

  raw <- if (!is.null(spec$custom_model)) {
    spec$custom_model
  } else if (spec$stage == "one_stage") {
    generate_stan_code_one_stage_mv(spec$outcome_labels, spec$priors)
  } else {
    generate_stan_code_two_stage_mv(spec$outcome_labels, spec$priors)
  }

  full   <- if (isTRUE(format)) format_stan_code(raw) else raw
  blocks <- parse_stan_blocks(full)
  out    <- c(blocks, list(full = full))
  class(out) <- c("bayesma_mv_stan_code", "bayesma_stan_code", "list")
  out
}


# -----------------------------------------------------------------------------
# Stage 3: stan data
# -----------------------------------------------------------------------------

#' Build the Stan data list for a bivariate meta-analysis specification
#'
#' @param spec A `bayesma_mv_spec` object.
#' @return A named list suitable for `cmdstanr::CmdStanModel$sample(data = ...)`.
#' @export
bayesma_mv_stan_data <- function(spec) {
  if (!inherits(spec, "bayesma_mv_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_mv_spec} object.")
  }

  sd_list <- list(
    S          = spec$S,
    y1         = spec$es_1$yi,
    y2         = spec$es_2$yi,
    se1        = spec$es_1$sei,
    se2        = spec$es_2$sei,
    rho_within = spec$rho_within
  )

  if (!is.null(spec$custom_data)) {
    for (nm in names(spec$custom_data)) sd_list[[nm]] <- spec$custom_data[[nm]]
  }

  sd_list
}


# -----------------------------------------------------------------------------
# Stage 4: fit
# -----------------------------------------------------------------------------

#' Compile and sample a bivariate meta-analysis model
#'
#' @param spec      A `bayesma_mv_spec` object.
#' @param code      A `bayesma_mv_stan_code` object. Defaults to
#'   `bayesma_mv_stan_code(spec)`.
#' @param stan_data A Stan data list. Defaults to `bayesma_mv_stan_data(spec)`.
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed MCMC settings.
#' @param ... Passed to `cmdstanr::CmdStanModel$sample()`.
#' @return An object of class `"bayesma_mv_fit"`.
#' @export
bayesma_mv_fit <- function(spec,
                           code          = bayesma_mv_stan_code(spec),
                           stan_data     = bayesma_mv_stan_data(spec),
                           chains        = 4,
                           iter_warmup   = 1000,
                           iter_sampling = 1000,
                           adapt_delta   = 0.95,
                           seed          = 1234,
                           ...) {
  if (!inherits(spec, "bayesma_mv_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_mv_spec} object.")
  }

  stan_program <- if (inherits(code, "bayesma_mv_stan_code")) code$full
  else as.character(code)

  mod <- get_cmdstan_model_cached(stan_program)

  fit <- mod$sample(
    data = stan_data, chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...
  )

  out <- list(fit = fit, stan_code = code, stan_data = stan_data)
  class(out) <- c("bayesma_mv_fit", "list")
  out
}


# -----------------------------------------------------------------------------
# Stage 5: extract
# -----------------------------------------------------------------------------

#' Extract tidy effect components from a bivariate meta-analysis fit
#'
#' @param fit  A `bayesma_mv_fit` object.
#' @param spec A `bayesma_mv_spec` object.
#' @return An object of class `"bayesma_mv_effects"`.
#' @export
bayesma_mv_extract <- function(fit, spec) {
  if (!inherits(fit, "bayesma_mv_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls bayesma_mv_fit} object.")
  }
  if (!inherits(spec, "bayesma_mv_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_mv_spec} object.")
  }

  cmdstan_fit <- fit$fit
  S <- spec$S

  key_vars <- c("mu1", "mu2", "tau1", "tau2", "rho_between")

  summary_tbl <- cmdstan_fit$summary(variables = key_vars) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      variable = dplyr::recode_values(
        .data$variable,
        "mu1"  ~ paste0("mu_",  spec$outcome_labels[1]),
        "mu2"  ~ paste0("mu_",  spec$outcome_labels[2]),
        "tau1" ~ paste0("tau_", spec$outcome_labels[1]),
        "tau2" ~ paste0("tau_", spec$outcome_labels[2])
      )
    )

  # ---- Draws ----
  draw_vars <- key_vars
  if (spec$stage == "two_stage") {
    for (k in 1:2) {
      for (s in seq_len(S)) {
        draw_vars <- c(draw_vars, paste0("theta[", k, ",", s, "]"))
      }
    }
  }
  draw_vars <- c(draw_vars, "mu_new[1]", "mu_new[2]")

  draws <- tryCatch(
    posterior::as_draws_df(cmdstan_fit$draws(variables = draw_vars)),
    error = function(e) {
      posterior::as_draws_df(cmdstan_fit$draws(variables = key_vars))
    }
  )

  # ---- Forest data ----
  forest_list <- purrr::map(1:2, function(k) {
    outcome_nm <- spec$outcome_labels[k]
    es_k <- if (k == 1) spec$es_1 else spec$es_2

    mu_k_draws <- as.vector(posterior::subset_draws(
      cmdstan_fit$draws(paste0("mu", k)), variable = paste0("mu", k)
    ))

    pooled_row <- tibble::tibble(
      study = "Pooled", outcome = outcome_nm,
      estimate = stats::median(mu_k_draws),
      lower = stats::quantile(mu_k_draws, 0.025),
      upper = stats::quantile(mu_k_draws, 0.975),
      type = "pooled"
    )

    if (spec$stage == "two_stage") {
      study_rows <- purrr::map(seq_len(S), function(s) {
        vn <- paste0("theta[", k, ",", s, "]")
        theta_ks <- tryCatch(
          as.vector(posterior::subset_draws(
            cmdstan_fit$draws(vn), variable = vn
          )),
          error = function(e) NULL
        )
        if (!is.null(theta_ks)) {
          tibble::tibble(
            study = spec$study_labels[s], outcome = outcome_nm,
            estimate = stats::median(theta_ks),
            lower = stats::quantile(theta_ks, 0.025),
            upper = stats::quantile(theta_ks, 0.975),
            type = "study"
          )
        } else {
          tibble::tibble(
            study = spec$study_labels[s], outcome = outcome_nm,
            estimate = es_k$yi[s],
            lower = es_k$yi[s] - 1.96 * es_k$sei[s],
            upper = es_k$yi[s] + 1.96 * es_k$sei[s],
            type = "study"
          )
        }
      }) |> purrr::list_rbind()
    } else {
      study_rows <- tibble::tibble(
        study = spec$study_labels, outcome = outcome_nm,
        estimate = es_k$yi,
        lower = es_k$yi - 1.96 * es_k$sei,
        upper = es_k$yi + 1.96 * es_k$sei,
        type = "study"
      )
    }

    dplyr::bind_rows(study_rows, pooled_row)
  })

  forest_df <- dplyr::bind_rows(forest_list) |>
    dplyr::mutate(
      study   = forcats::fct_inorder(.data$study),
      outcome = factor(.data$outcome, levels = spec$outcome_labels)
    )

  # ---- Prediction intervals ----
  pred_intervals <- purrr::map(1:2, function(k) {
    mn_var <- paste0("mu_new[", k, "]")
    mn_k <- tryCatch(
      as.vector(posterior::subset_draws(
        cmdstan_fit$draws(mn_var), variable = mn_var
      )),
      error = function(e) NULL
    )
    if (!is.null(mn_k)) {
      tibble::tibble(
        outcome  = spec$outcome_labels[k],
        estimate = stats::median(mn_k),
        lower    = stats::quantile(mn_k, 0.025),
        upper    = stats::quantile(mn_k, 0.975)
      )
    }
  }) |> purrr::list_rbind()

  eff <- list(
    summary        = summary_tbl,
    draws          = draws,
    forest_df      = forest_df,
    pred_intervals = pred_intervals
  )
  class(eff) <- c("bayesma_mv_effects", "list")
  eff
}


# -----------------------------------------------------------------------------
# Stage 6: output
# -----------------------------------------------------------------------------

#' Assemble a `bayesma_mv` object from pipeline outputs
#'
#' @param spec    A `bayesma_mv_spec` object.
#' @param fit     A `bayesma_mv_fit` object.
#' @param effects A `bayesma_mv_effects` object.
#' @return A list of class `c("bayesma_mv", "bayesma")`.
#' @export
bayesma_mv_output <- function(spec, fit, effects) {
  if (!inherits(spec,    "bayesma_mv_spec"))    cli::cli_abort("{.arg spec} must be {.cls bayesma_mv_spec}.")
  if (!inherits(fit,     "bayesma_mv_fit"))     cli::cli_abort("{.arg fit} must be {.cls bayesma_mv_fit}.")
  if (!inherits(effects, "bayesma_mv_effects")) cli::cli_abort("{.arg effects} must be {.cls bayesma_mv_effects}.")

  code_full <- if (inherits(fit$stan_code, "bayesma_mv_stan_code")) fit$stan_code$full
  else as.character(fit$stan_code)

  out <- list(
    fit            = fit$fit,
    summary        = effects$summary,
    forest_df      = effects$forest_df,
    draws          = effects$draws,
    pred_intervals = effects$pred_intervals,
    stan_code      = code_full,
    stan_data      = fit$stan_data,
    es             = list(outcome_1 = spec$es_1, outcome_2 = spec$es_2),
    meta           = list(
      S              = spec$S,
      outcome_labels = spec$outcome_labels,
      study_labels   = spec$study_labels,
      likelihood     = spec$likelihood,
      stage          = spec$stage,
      priors         = spec$priors,
      rho_within     = spec$rho_within,
      call_args      = spec$call_args
    )
  )
  class(out) <- c("bayesma_mv", "bayesma")
  out
}


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

#' Run a Multivariate Bayesian Meta-Analysis in Stan
#'
#' Thin orchestrator over the six-stage pipeline:
#' [bayesma_mv_spec()], [bayesma_mv_stan_code()], [bayesma_mv_stan_data()],
#' [bayesma_mv_fit()], [bayesma_mv_extract()], [bayesma_mv_output()].
#'
#' @inheritParams bayesma_mv_spec
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed MCMC settings.
#' @param return_stage Character. One of `"full"` (default), `"spec"`,
#'   `"code"`, `"data"`, or `"fit"`.
#' @param ... Passed to `cmdstanr::sample()`.
#' @return An object of class `c("bayesma_mv", "bayesma")`.
#' @export
bayesma_mv <- function(
    data,
    studyvar,
    mean_ctrl_1, mean_int_1, sd_ctrl_1, sd_int_1, n_ctrl_1, n_int_1,
    mean_ctrl_2, mean_int_2, sd_ctrl_2, sd_int_2, n_ctrl_2, n_int_2,
    outcome_labels    = c("outcome_1", "outcome_2"),
    likelihood        = c("gaussian"),
    stage             = c("two_stage", "one_stage"),
    rho_within        = 0.5,
    mu_prior          = NULL,
    tau_prior         = NULL,
    rho_between_prior = NULL,
    custom_model      = NULL,
    custom_data       = NULL,
    return_stage      = c("full", "spec", "code", "data", "fit"),
    chains            = 4,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    adapt_delta       = 0.95,
    seed              = 1234,
    ...
) {
  return_stage <- rlang::arg_match(return_stage)

  spec <- bayesma_mv_spec(
    data = data, studyvar = {{ studyvar }},
    mean_ctrl_1 = {{ mean_ctrl_1 }}, mean_int_1 = {{ mean_int_1 }},
    sd_ctrl_1 = {{ sd_ctrl_1 }}, sd_int_1 = {{ sd_int_1 }},
    n_ctrl_1 = {{ n_ctrl_1 }}, n_int_1 = {{ n_int_1 }},
    mean_ctrl_2 = {{ mean_ctrl_2 }}, mean_int_2 = {{ mean_int_2 }},
    sd_ctrl_2 = {{ sd_ctrl_2 }}, sd_int_2 = {{ sd_int_2 }},
    n_ctrl_2 = {{ n_ctrl_2 }}, n_int_2 = {{ n_int_2 }},
    outcome_labels = outcome_labels, likelihood = likelihood,
    stage = stage, rho_within = rho_within,
    mu_prior = mu_prior, tau_prior = tau_prior,
    rho_between_prior = rho_between_prior,
    custom_model = custom_model, custom_data = custom_data
  )
  if (return_stage == "spec") return(spec)

  code <- bayesma_mv_stan_code(spec)
  if (return_stage == "code") return(code)

  stan_data <- bayesma_mv_stan_data(spec)
  if (return_stage == "data") return(stan_data)

  fit <- bayesma_mv_fit(
    spec = spec, code = code, stan_data = stan_data,
    chains = chains, iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, adapt_delta = adapt_delta,
    seed = seed, ...
  )
  if (return_stage == "fit") return(fit)

  effects <- bayesma_mv_extract(fit, spec)
  bayesma_mv_output(spec, fit, effects)
}
