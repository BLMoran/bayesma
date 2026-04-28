#' Robust Bayesian Model Averaging for Meta-Analysis
#'
#' Fits a Robust Bayesian Meta-Analysis (RoBMA) model using Bayesian model
#' averaging across models with and without an effect, heterogeneity, and
#' publication bias. Use `stan_code(model)` to inspect the generated Stan
#' programs after fitting.
#'
#' @param data A data frame with one row per study.
#' @param studyvar Character. Column name of the study identifier.
#' @param event_ctrl,event_int Character. Column names of event counts
#'   (binomial / Poisson likelihoods).
#' @param mean_ctrl,mean_int,sd_ctrl,sd_int Character. Column names of arm
#'   means and SDs (Gaussian likelihood).
#' @param n_ctrl,n_int Character. Column names of arm sample sizes.
#' @param likelihood Character. One of `"binomial"`, `"gaussian"`, `"poisson"`.
#' @param priors_effect,priors_effect_null,priors_heterogeneity,priors_heterogeneity_null,priors_bias,priors_bias_null
#'   Lists of prior objects for the effect, heterogeneity, and publication-bias
#'   components (alternative and null). If `NULL`, RoBMA defaults are used.
#' @param rescale_priors Numeric. Scale factor applied to default priors.
#'   Default `1`.
#' @param method Character. `"bridge"` (default) uses bridge sampling across
#'   the full model grid; `"ss"` uses a single spike-and-slab Stan model.
#' @param bias_indicator Character. Spike-and-slab bias mechanism:
#'   `"bias_corrected"`, `"pet_peese"`, or `"selection_weight"`.
#' @param null_range Numeric vector of length 2 giving the null range on the
#'   log scale (e.g., `c(-0.1, 0.1)` for log OR). Effects within this range
#'   are considered practically equivalent to zero. Defaults to `NULL` (point
#'   null at exactly zero). For OR/RR, `c(-0.1, 0.1)` corresponds to OR/RR
#'   in `[0.905, 1.105]`.
#' @param b_prior Prior on the `b` slope for spike-and-slab bias correction.
#' @param p_bias_prior Prior on the bias inclusion probability.
#' @param p_cutoffs Numeric vector of one-sided p-value cutoffs for
#'   selection-weight models. Default `c(0.025, 0.05)`.
#' @param parallel Logical. Fit the bridgesampling grid in parallel.
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed MCMC settings.
#' @param quiet Logical. Suppress per-step progress messages.
#' @param custom_model Optional Stan program(s) that override code generation.
#'   For `method = "bridge"`, a named list of character scalars keyed by model
#'   label. For `method = "ss"`, a single character scalar.
#' @param custom_data Optional Stan data overrides merged onto the auto-built
#'   data list(s). Same shape conventions as `custom_model`.
#' @param format Logical. If `TRUE` (default), auto-format generated Stan
#'   programs with `stanc --auto-format`.
#' @param ... Additional arguments passed to `cmdstanr::CmdStanModel$sample()`.
#' @return An object of class `c("bayesma_robma", "bayesma")`.
#' @export
robma <- function(
    data, studyvar,
    event_ctrl = NULL, event_int = NULL,
    mean_ctrl = NULL, mean_int = NULL,
    sd_ctrl = NULL, sd_int = NULL,
    n_ctrl = NULL, n_int = NULL,
    likelihood = c("binomial", "gaussian", "poisson"),
    priors_effect = NULL,
    priors_effect_null = NULL,
    priors_heterogeneity = NULL,
    priors_heterogeneity_null = NULL,
    priors_bias = NULL,
    priors_bias_null = NULL,
    rescale_priors = 1,
    method = c("bridge", "ss"),
    bias_indicator = c("bias_corrected", "pet_peese", "selection_weight"),
    null_range = NULL,
    b_prior = NULL, p_bias_prior = NULL, p_cutoffs = c(0.025, 0.05),
    parallel = FALSE, chains = 4, iter_warmup = 1000,
    iter_sampling = 1000, adapt_delta = 0.95, seed = 1234,
    quiet = FALSE,
    custom_model = NULL,
    custom_data = NULL,
    format = TRUE,
    ...
) {
  spec <- robma_spec(
    data                      = data,
    studyvar                  = studyvar,
    event_ctrl                = event_ctrl, event_int = event_int,
    mean_ctrl                 = mean_ctrl,  mean_int  = mean_int,
    sd_ctrl                   = sd_ctrl,    sd_int    = sd_int,
    n_ctrl                    = n_ctrl,     n_int     = n_int,
    likelihood                = likelihood,
    priors_effect             = priors_effect,
    priors_effect_null        = priors_effect_null,
    priors_heterogeneity      = priors_heterogeneity,
    priors_heterogeneity_null = priors_heterogeneity_null,
    priors_bias               = priors_bias,
    priors_bias_null          = priors_bias_null,
    rescale_priors            = rescale_priors,
    method                    = method,
    bias_indicator            = bias_indicator,
    null_range                = null_range,
    b_prior                   = b_prior,
    p_bias_prior              = p_bias_prior,
    p_cutoffs                 = p_cutoffs,
    parallel                  = parallel,
    chains                    = chains,
    iter_warmup               = iter_warmup,
    iter_sampling             = iter_sampling,
    adapt_delta               = adapt_delta,
    seed                      = seed,
    quiet                     = quiet,
    custom_model              = custom_model,
    custom_data               = custom_data,
    ...
  )
  code      <- robma_stan_code(spec, format = format)
  stan_data <- robma_stan_data(spec)
  fit       <- robma_fit(
    spec          = spec,
    code          = code,
    stan_data     = stan_data,
    chains        = chains,
    iter_warmup   = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta   = adapt_delta,
    seed          = seed,
    parallel      = parallel,
    quiet         = quiet,
    ...
  )
  effects <- robma_extract(fit, spec)
  robma_output(spec, fit, effects)
}



# S3: print

#' @export
print.bayesma_robma <- function(x, ...) {
  method <- x$meta$method %||% "bridge"
  cli::cli_h1("Robust Bayesian Meta-Analysis (RoBMA)")
  if (method == "bridge") {
    cli::cli_alert_info(
      "{x$meta$n_models} component models fitted via bridge sampling")
  } else {
    bi <- x$meta$bias_indicator %||% "bias_corrected"
    cli::cli_alert_info(
      "Spike-and-slab joint model (bias indicator: {.val {bi}})")
  }
  cli::cli_text("")

  pooled    <- dplyr::filter(x$forest_df, .data$type == "pooled")
  eff_label <- x$meta$effect_label
  is_ratio  <- eff_label %in% c("log_or", "log_rr")

  cli::cli_h2("Model-Averaged Pooled Estimate")
  cli::cli_text("  {eff_label}: {round(pooled$estimate, 3)} ",
                "[{round(pooled$lower, 3)}, {round(pooled$upper, 3)}]")
  if (is_ratio) {
    ratio_name <- gsub("log_", "", eff_label)
    cli::cli_text("  {ratio_name}: {round(exp(pooled$estimate), 3)} ",
                  "[{round(exp(pooled$lower), 3)}, {round(exp(pooled$upper), 3)}]")
  }

  # ---- Posterior Probabilities (model-level) ----
  pp <- x$posterior_probs
  cli::cli_h2("Posterior Probabilities (model comparison)")
  cli::cli_text("  P(effect != 0):          {round(pp$effect, 3)}")
  cli::cli_text("  P(heterogeneity):        {round(pp$heterogeneity, 3)}")
  cli::cli_text("  P(publication bias):     {round(pp$bias, 3)}")

  # ---- Direction / null range probabilities ----
  nrp <- x$meta$null_range_probs
  nr  <- x$meta$null_range

  cli::cli_h2("Direction Probabilities (model-averaged posterior)")
  if (is.null(nr)) {
    cli::cli_text("  P(mu < 0):               {round(nrp$p_negative, 3)}")
    cli::cli_text("  P(mu > 0):               {round(nrp$p_positive, 3)}")
    if (nrp$p_null > 0.001) {
      cli::cli_text("  P(mu = 0):               {round(nrp$p_null, 3)}")
    }
  } else {
    nr_log <- nrp$null_range
    nr_nat <- nrp$null_range_natural
    if (is_ratio) {
      ratio_name <- gsub("log_", "", eff_label)
      cli::cli_text(paste0(
        "  Null range: ",
        ratio_name, " in [", round(nr_nat[1], 3), ", ", round(nr_nat[2], 3), "]",
        " (log scale: [", round(nr_log[1], 4), ", ", round(nr_log[2], 4), "])"))
    } else {
      cli::cli_text(paste0(
        "  Null range: [", round(nr_log[1], 3), ", ", round(nr_log[2], 3), "]"))
    }
    cli::cli_text("")
    if (is_ratio) {
      ratio_name <- gsub("log_", "", eff_label)
      cli::cli_text(paste0(
        "  P(", ratio_name, " < ", round(nr_nat[1], 3), "):          ",
        round(nrp$p_negative, 3)))
      cli::cli_text(paste0(
        "  P(practically null):     ", round(nrp$p_null, 3)))
      cli::cli_text(paste0(
        "  P(", ratio_name, " > ", round(nr_nat[2], 3), "):          ",
        round(nrp$p_positive, 3)))
    } else {
      cli::cli_text(paste0(
        "  P(mu < ", round(nr_log[1], 3), "):          ",
        round(nrp$p_negative, 3)))
      cli::cli_text(paste0(
        "  P(practically null):     ", round(nrp$p_null, 3)))
      cli::cli_text(paste0(
        "  P(mu > ", round(nr_log[2], 3), "):          ",
        round(nrp$p_positive, 3)))
    }
  }

  # ---- Inclusion BFs ----
  bf <- x$inclusion_bf
  format_bf <- function(val) {
    if (is.na(val)) return("NA")
    if (is.infinite(val)) return("Inf")
    round(val, 2)
  }
  cli::cli_h2("Inclusion Bayes Factors")
  cli::cli_text("  Effect (H1 vs H0):       BF = {format_bf(bf$effect)}")
  cli::cli_text("  Publication bias:        BF = {format_bf(bf$bias)}")
  cli::cli_text("  Heterogeneity:           BF = {format_bf(bf$heterogeneity)}")

  # ---- Component Models (bridge only) ----
  if (method == "bridge") {
    cli::cli_h2("Component Models")
    mt <- x$model_table
    purrr::walk(seq_len(nrow(mt)), function(i) {
      pp_str <- if (mt$post_prob[i] < 0.001 && mt$post_prob[i] > 0) {
        formatC(mt$post_prob[i], format = "e", digits = 1)
      } else {
        as.character(round(mt$post_prob[i], 3))
      }
      cli::cli_text("  {mt$model[i]}: {pp_str} ",
                    "(log ML = {round(mt$log_ml[i], 1)})")
    })
  }
  invisible(x)
}


# ============================================================================
# S3: summary
# ============================================================================

#' @export
summary.bayesma_robma <- function(object, ...) {
  method <- object$meta$method %||% "bridge"
  cli::cli_h1("RoBMA Summary")
  if (method == "bridge") {
    cli::cli_text("Method: Bridge sampling ({object$meta$n_models} models)")
  } else {
    bi <- object$meta$bias_indicator %||% "bias_corrected"
    cli::cli_text("Method: Spike-and-slab (bias indicator: {.val {bi}})")
  }
  cli::cli_text("")

  pooled    <- dplyr::filter(object$forest_df, .data$type == "pooled")
  eff_label <- object$meta$effect_label
  is_ratio  <- eff_label %in% c("log_or", "log_rr")

  cli::cli_h2("Model-Averaged Effect")
  cli::cli_text("  Estimate ({eff_label}): {round(pooled$estimate, 4)}")
  cli::cli_text("  95% CrI: [{round(pooled$lower, 4)}, {round(pooled$upper, 4)}]")
  if (is_ratio) {
    ratio_name <- toupper(gsub("log_", "", eff_label))
    cli::cli_text("  {ratio_name}: {round(exp(pooled$estimate), 4)} ",
                  "[{round(exp(pooled$lower), 4)}, {round(exp(pooled$upper), 4)}]")
  }

  # ---- Model comparison ----
  pp <- object$posterior_probs
  cli::cli_h2("Model Comparison")
  cli::cli_text("  P(effect != 0):          {round(pp$effect, 4)}")
  cli::cli_text("  P(heterogeneity):        {round(pp$heterogeneity, 4)}")
  cli::cli_text("  P(publication bias):     {round(pp$bias, 4)}")

  format_bf <- function(val) {
    if (is.na(val)) return("NA")
    if (is.infinite(val)) return("Inf")
    round(val, 2)
  }
  bf <- object$inclusion_bf
  cli::cli_text("")
  cli::cli_text("  Effect BF10:             {format_bf(bf$effect)}")
  cli::cli_text("  Bias BF:                 {format_bf(bf$bias)}")
  cli::cli_text("  Heterogeneity BF:        {format_bf(bf$heterogeneity)}")

  # ---- Direction / null range ----
  nrp <- object$meta$null_range_probs
  nr  <- object$meta$null_range

  cli::cli_h2("Direction Probabilities (model-averaged)")
  if (is.null(nr)) {
    cli::cli_text("  P(mu < 0) = {round(nrp$p_negative, 4)}")
    cli::cli_text("  P(mu > 0) = {round(nrp$p_positive, 4)}")
    if (nrp$p_null > 0.001) {
      cli::cli_text("  P(mu = 0) = {round(nrp$p_null, 4)}")
    }
  } else {
    if (is_ratio) {
      nr_nat <- nrp$null_range_natural
      ratio_name <- toupper(gsub("log_", "", eff_label))
      cli::cli_text("  Null range: [{round(nr[1], 3)}, {round(nr[2], 3)}] ",
                    "({ratio_name}: [{round(nr_nat[1], 3)}, {round(nr_nat[2], 3)}])")
    } else {
      cli::cli_text("  Null range: [{round(nr[1], 3)}, {round(nr[2], 3)}]")
    }
    cli::cli_text("")
    cli::cli_text("  P(harmful):              {round(nrp$p_negative, 4)}",
                  "  (mu < {round(nr[1], 3)})")
    cli::cli_text("  P(practically null):     {round(nrp$p_null, 4)}",
                  "  (mu in null range)")
    cli::cli_text("  P(beneficial):           {round(nrp$p_positive, 4)}",
                  "  (mu > {round(nr[2], 3)})")
  }

  cli::cli_text("")
  cli::cli_text("  Note: P(effect != 0) is a model comparison quantity")
  cli::cli_text("  (proportion of posterior weight on H1 models).")
  cli::cli_text("  Direction probabilities are from the model-averaged")
  cli::cli_text("  posterior draws and reflect parameter uncertainty.")

  # ---- Component models or SS parameters ----
  if (method == "bridge") {
    cli::cli_h2("Component Models")
    print(object$model_table)
  } else {
    cli::cli_h2("Spike-and-Slab Parameter Summaries")
    fit <- object$component_fits[[1]]$fit
    tryCatch({
      summ <- fit$summary()
      key_pars <- c("mu", "mu_raw", "tau", "tau_raw",
                    "pip_effect", "pip_hetero", "pip_bias")
      summ_filt <- dplyr::filter(summ, .data$variable %in% key_pars)
      print(summ_filt)
    }, error = function(e) {
      cli::cli_alert_warning("Could not retrieve parameter summaries.")
    })
  }
  invisible(object)
}


# ============================================================================
# Main bridge sampling function
# ============================================================================
# Phase 1: Build model specs (stan_code + stan_data + metadata)
#
# Each spec is a list with:
#   stan_code, stan_data, label, prior_weight,
#   is_effect_null, is_hetero_null, is_bias_null,
#   analytic (TRUE if no Stan model needed),
#   log_ml_override (numeric if analytic, NULL otherwise)
# ============================================================================

build_model_spec <- function(spec, es, S, n_c, n_i, likelihood) {

  is_effect_null <- isTRUE(spec$is_effect_null)
  is_hetero_null <- isTRUE(spec$is_hetero_null)
  is_bias_null   <- isTRUE(spec$is_bias_null)
  bias_type      <- if (is_bias_null) "none" else spec$bias_prior$type %||% "none"
  has_re         <- !is_hetero_null
  model_type     <- if (has_re) "random_effect" else "common_effect"

  base_meta <- list(
    label          = spec$label,
    prior_weight   = spec$prior_weight,
    is_effect_null = is_effect_null,
    is_hetero_null = is_hetero_null,
    is_bias_null   = is_bias_null
  )

  # =========================================================
  # CASE 1: Effect null + bias null
  # =========================================================
  if (is_effect_null && is_bias_null) {
    if (is_hetero_null) {
      # CE / H0 / no bias: analytic, no Stan model
      log_ml_val <- compute_log_ml_ce_null(es$yi, es$sei)
      return(c(base_meta, list(
        stan_code = NULL, stan_data = NULL,
        analytic = TRUE, log_ml_override = log_ml_val
      )))
    } else {
      # RE / H0 / no bias: needs Stan for bridge sampling
      pr <- resolve_priors("two_stage", likelihood, "random_effect", "normal",
                           spec$effect_prior, spec$hetero_prior,
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
      sc <- generate_stan_code_two_stage_null(
        model_type = "random_effect", re_dist = "normal", priors = pr)
      sd <- list(S = S, y = es$yi, se = es$sei)
      return(c(base_meta, list(
        stan_code = sc, stan_data = sd,
        analytic = FALSE, log_ml_override = NULL,
        # Store hetero_prior for quadrature fallback
        hetero_prior_for_fallback = spec$hetero_prior
      )))
    }
  }

  # =========================================================
  # CASE 2: Effect null + bias H1
  # =========================================================
  if (is_effect_null && !is_bias_null) {
    result <- switch(
      bias_type,
      pet = spec_component_pet_peese_h0(es, S, n_c, n_i, has_re),
      peese = spec_component_peese_h0(es, S, n_c, n_i, has_re),
      weight_function = {
        p_cuts <- sort(spec$bias_prior$parameters$steps %||% c(0.025, 0.05))
        spec_component_selection_weight_h0(es, S, has_re, p_cuts)
      },
      jung = spec_component_jung_h0(es, S, has_re),
      {
        # Unknown bias type: treat as no-bias analytic
        log_ml_val <- if (has_re) {
          compute_log_ml_re_null(es$yi, es$sei, spec$hetero_prior)
        } else {
          compute_log_ml_ce_null(es$yi, es$sei)
        }
        list(stan_code = NULL, stan_data = NULL,
             analytic = TRUE, log_ml_override = log_ml_val)
      }
    )
    return(c(base_meta, result))
  }

  # =========================================================
  # CASE 3: Effect H1 + bias null
  # =========================================================
  if (is_bias_null || bias_type == "none") {
    pr <- resolve_priors("two_stage", likelihood, model_type, "normal",
                         spec$effect_prior, spec$hetero_prior,
                         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
    sc <- generate_stan_code_two_stage(
      model_type = model_type, re_dist = "normal",
      use_t_likelihood = FALSE, priors = pr,
      robust_config = list(enabled = FALSE))
    sd <- list(S = S, y = es$yi, se = es$sei)
    return(c(base_meta, list(stan_code = sc, stan_data = sd)))
  }

  # =========================================================
  # CASE 4: Effect H1 + bias H1
  # =========================================================
  pr <- resolve_priors("two_stage", likelihood, model_type, "normal",
                       spec$effect_prior, spec$hetero_prior,
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       spec$bias_prior$parameters)

  result <- switch(
    bias_type,
    pet = spec_component_pet_peese(es, S, n_c, n_i, pr, has_re),
    peese = spec_component_pet_peese(es, S, n_c, n_i, pr, has_re),
    weight_function = {
      p_cuts <- sort(spec$bias_prior$parameters$steps %||% c(0.025, 0.05))
      if (has_re) {
        spec_component_selection_weight(es, S, pr, p_cuts)
      } else {
        spec_component_selection_weight_fe(es, S, pr, p_cuts)
      }
    },
    copas = spec_component_selection_copas(es, S, pr, has_re),
    jung = spec_component_bias_corrected(es, S, pr, has_re)
  )
  c(base_meta, result)
}


# ============================================================================
# Phase 2: Main robma_bridge — build, compile, sample, bridge
# ============================================================================

robma_bridge <- function(
    model_grid, es, S, n_c, n_i, study_labels, likelihood,
    parallel, chains, iter_warmup, iter_sampling, adapt_delta, seed,
    quiet = FALSE, ...
) {
  if (!requireNamespace("bridgesampling", quietly = TRUE))
    cli::cli_abort("{.pkg bridgesampling} is required for method = 'bridge'.")

  inform <- if (quiet) \(...) invisible(NULL) else cli::cli_alert_info
  inform_step <- if (quiet) \(...) invisible(NULL) else cli::cli_alert

  n_models <- length(model_grid)

  # ------------------------------------------------------------------
  # Step 1: Build all model specs (sequential — calls internal helpers)
  # ------------------------------------------------------------------
  inform("RoBMA (bridge): Building {n_models} model specifications...")
  model_specs <- purrr::imap(model_grid, function(spec, idx) {
    inform_step("  [{idx}/{n_models}] {spec$label}")
    tryCatch({
      out <- build_model_spec(spec, es, S, n_c, n_i, likelihood)
      # Ensure defaults for fields that may not be set by all code paths
      if (is.null(out$analytic))       out$analytic <- FALSE
      if (is.null(out$log_ml_override)) out$log_ml_override <- NULL
      out
    },
    error = function(e) {
      cli::cli_warn("Spec build failed for {.val {spec$label}}: {e$message}")
      NULL
    }
    )
  })

  # Drop failures
  spec_ok <- !purrr::map_lgl(model_specs, is.null)
  model_specs <- model_specs[spec_ok]
  if (length(model_specs) < 2)
    cli::cli_abort("RoBMA requires at least 2 valid model specifications.")

  # ------------------------------------------------------------------
  # Step 2: Compile unique Stan models (sequential, cached)
  # ------------------------------------------------------------------
  inform("Compiling unique Stan models...")

  # Build a code -> compiled model lookup
  stan_codes <- purrr::map_chr(model_specs, function(sp) {
    if (isTRUE(sp$analytic) || is.null(sp$stan_code)) return(NA_character_)
    sp$stan_code
  })

  unique_codes <- unique(stats::na.omit(stan_codes))
  inform("  {length(unique_codes)} unique Stan programs to compile")

  compiled_models <- purrr::map(
    seq_along(unique_codes),
    function(i) {
      inform_step("  Compiling model {i}/{length(unique_codes)}")
      stan_file <- cmdstanr::write_stan_file(unique_codes[[i]])
      mod <- cmdstanr::cmdstan_model(stan_file, force_recompile = TRUE,
                                     quiet = TRUE)
      mod
    }
  )
  names(compiled_models) <- unique_codes

  # Attach compiled model to each spec
  model_specs <- purrr::map(model_specs, function(sp) {
    if (isTRUE(sp$analytic) || is.null(sp$stan_code)) {
      sp$compiled_model <- NULL
    } else {
      sp$compiled_model <- compiled_models[[sp$stan_code]]
    }
    sp
  })

  # ------------------------------------------------------------------
  # Step 3: Sample all models (parallelisable — no internal helpers needed)
  # ------------------------------------------------------------------
  inform("Sampling {length(model_specs)} models...")

  sample_one_model <- function(sp, idx) {
    inform_step("  [{idx}/{length(model_specs)}] {sp$label}")

    # Analytic models: no sampling needed
    if (isTRUE(sp$analytic)) {
      return(list(
        fit = NULL,
        stan_code = sp$stan_code,
        stan_data = sp$stan_data,
        label = sp$label,
        prior_weight = sp$prior_weight,
        is_effect_null = sp$is_effect_null,
        is_hetero_null = sp$is_hetero_null,
        is_bias_null = sp$is_bias_null,
        log_ml_override = sp$log_ml_override,
        analytic = TRUE
      ))
    }

    # Sample from pre-compiled model
    tryCatch({
      fit <- sp$compiled_model$sample(
        data = sp$stan_data,
        chains = chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        adapt_delta = adapt_delta,
        seed = seed,
        refresh = 0,
        show_messages = FALSE,
        show_exceptions = FALSE,
        ...
      )

      list(
        fit = fit,
        stan_code = sp$stan_code,
        stan_data = sp$stan_data,
        label = sp$label,
        prior_weight = sp$prior_weight,
        is_effect_null = sp$is_effect_null,
        is_hetero_null = sp$is_hetero_null,
        is_bias_null = sp$is_bias_null,
        log_ml_override = sp$log_ml_override,
        analytic = FALSE,
        # Pass through for RE/H0/no-bias fallback
        hetero_prior_for_fallback = sp$hetero_prior_for_fallback
      )
    }, error = function(e) {
      cli::cli_warn("Sampling failed for {.val {sp$label}}: {e$message}")
      NULL
    })
  }

  if (parallel) {
    if (requireNamespace("mirai", quietly = TRUE)) {
      daemons_were_set <- tryCatch({
        s <- mirai::status()
        is.data.frame(s$daemons) && nrow(s$daemons) > 0L
      }, error = function(e) FALSE)

      if (!daemons_were_set) {
        n_cores <- max(1L, parallel::detectCores() - 1L)
        mirai::daemons(n_cores)
        on.exit(mirai::daemons(0), add = TRUE)
      }

      # Load bayesma + cmdstanr on daemons
      tryCatch(
        mirai::everywhere({
          loadNamespace("bayesma"); loadNamespace("cmdstanr")
        }),
        error = function(e) NULL
      )

      inform("  Parallel sampling via mirai...")
      results <- mirai::mirai_map(
        seq_along(model_specs),
        function(idx, model_specs, sample_one_model) {
          sample_one_model(model_specs[[idx]], idx)
        },
        .args = list(model_specs = model_specs,
                     sample_one_model = sample_one_model)
      )
      component_fits <- results[]
    } else if (.Platform$OS.type != "windows") {
      # mclapply fallback (macOS/Linux)
      n_cores <- max(1L, parallel::detectCores() - 1L)
      inform("  Using {n_cores} cores via mclapply...")
      component_fits <- parallel::mclapply(
        seq_along(model_specs),
        function(idx) sample_one_model(model_specs[[idx]], idx),
        mc.cores = n_cores,
        mc.set.seed = TRUE
      )
    } else {
      cli::cli_warn("Parallel not available (install {.pkg mirai}). Running sequentially.")
      component_fits <- purrr::imap(model_specs, sample_one_model)
    }
  } else {
    component_fits <- purrr::imap(model_specs, sample_one_model)
  }
  component_fits <- purrr::compact(component_fits)
  if (length(component_fits) < 2)
    cli::cli_abort("RoBMA requires at least 2 successfully fitted models.")

  # ------------------------------------------------------------------
  # Step 4: Bridge sampling for marginal likelihoods
  # ------------------------------------------------------------------
  inform("Computing marginal likelihoods...")
  bs_cache <- new.env(parent = emptyenv())
  log_mls <- purrr::map_dbl(component_fits, function(comp) {
    if (!is.null(comp$log_ml_override)) return(comp$log_ml_override)

    # RE/H0/no-bias: try bridge, fall back to quadrature
    fit <- comp$fit
    if (is.null(fit)) return(-Inf)

    fit_id <- rlang::obj_address(fit)
    if (exists(fit_id, envir = bs_cache)) return(get(fit_id, envir = bs_cache))

    log_ml_val <- tryCatch({
      invisible(utils::capture.output(
        bs <- bridgesampling::bridge_sampler(fit, silent = TRUE),
        type = "message"
      ))
      bs$logml
    }, error = function(e) {
      cli::cli_warn("Bridge sampling failed for {.val {comp$label}}")
      -Inf
    })

    # Quadrature fallback for RE/H0/no-bias
    if (!is.finite(log_ml_val) &&
        isTRUE(comp$is_effect_null) && isTRUE(comp$is_bias_null) &&
        !isTRUE(comp$is_hetero_null) &&
        !is.null(comp$hetero_prior_for_fallback)) {
      log_ml_val <- compute_log_ml_re_null(
        es$yi, es$sei, comp$hetero_prior_for_fallback)
      cli::cli_alert_info(
        "    {comp$label}: quadrature fallback log ML = {round(log_ml_val, 2)}")
    }

    assign(fit_id, log_ml_val, envir = bs_cache)
    log_ml_val
  })

  # ------------------------------------------------------------------
  # Step 5: Posterior probs, BFs, model averaging (unchanged logic)
  # ------------------------------------------------------------------
  prior_weights <- purrr::map_dbl(component_fits, ~ .x$prior_weight)
  post_probs    <- compute_posterior_probs(log_mls, prior_weights)
  names(post_probs) <- purrr::map_chr(component_fits, ~ .x$label)
  finite_mask <- is.finite(log_mls + log(prior_weights / sum(prior_weights)))

  is_h1      <- !purrr::map_lgl(component_fits, ~ .x$is_effect_null)
  has_bias   <- !purrr::map_lgl(component_fits, ~ .x$is_bias_null)
  has_hetero <- !purrr::map_lgl(component_fits, ~ .x$is_hetero_null)

  bf_effect <- compute_inclusion_bf(is_h1, post_probs, prior_weights, finite_mask)
  bf_bias   <- compute_inclusion_bf(has_bias, post_probs, prior_weights, finite_mask)
  bf_hetero <- compute_inclusion_bf(has_hetero, post_probs, prior_weights, finite_mask)

  # Model-averaged posterior
  cli::cli_alert_info("Computing model-averaged posterior...")
  n_total_draws <- iter_sampling * chains
  averaged_mu <- purrr::imap(component_fits, function(comp, idx) {
    n_from <- max(1, round(post_probs[idx] * n_total_draws))
    if (isTRUE(comp$is_effect_null))
      return(tibble::tibble(mu = rep(0, n_from), model = comp$label))
    tryCatch({
      mu_draws <- as.vector(posterior::subset_draws(
        comp$fit$draws("mu"), variable = "mu"))
      tibble::tibble(mu = sample(mu_draws, n_from, replace = TRUE),
                     model = comp$label)
    }, error = function(e) {
      tibble::tibble(mu = rep(0, n_from), model = comp$label)
    })
  }) |> purrr::list_rbind()

  model_table <- tibble::tibble(
    model = purrr::map_chr(component_fits, ~ .x$label),
    log_ml = log_mls, prior_weight = prior_weights,
    post_prob = as.numeric(post_probs),
    null_effect = !is_h1, has_bias = has_bias,
    has_heterogeneity = has_hetero
  ) |> dplyr::arrange(dplyr::desc(.data$post_prob))

  list(
    averaged_draws = averaged_mu, model_table = model_table,
    inclusion_bf = list(effect = bf_effect, bias = bf_bias,
                        heterogeneity = bf_hetero),
    posterior_probs = list(effect = sum(post_probs[is_h1]),
                           bias = sum(post_probs[has_bias]),
                           heterogeneity = sum(post_probs[has_hetero])),
    component_fits = component_fits, post_probs = post_probs,
    log_marginal_likelihoods = log_mls)
}


# ============================================================================
# Spec builders: H1 effect components (return stan_code + stan_data only)
# ============================================================================


# ---- H1 / PET-PEESE (FE or RE based on has_re) ----
spec_component_pet_peese <- function(es, S, n_c, n_i, priors, has_re) {
  mu_tgt   <- emit_prior_target(priors$mu, "mu")
  sp       <- priors$selection %||% list()
  if (is.null(sp$beta_bias)) sp$beta_bias <- normal(0, 5)
  beta_tgt <- emit_prior_target(sp$beta_bias, "beta_bias")

  if (has_re) {
    tau_tgt  <- emit_prior_target(priors$tau, "tau")
    tau_bnds <- emit_prior_bounds(priors$tau, default_lower = 0)

    stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] n_total;
}}
transformed data {{
  vector[N] inv_sqrt_n;
  for (i in 1:N)
    inv_sqrt_n[i] = 1.0 / sqrt(n_total[i]);
}}
parameters {{
  real mu;
  real{tau_bnds} tau;
  real beta_bias;
}}
model {{
  {mu_tgt}
  {tau_tgt}
  {beta_tgt}
  for (i in 1:N) {{
    real sigma_i = sqrt(square(tau) + square(se[i]));
    target += normal_lpdf(y[i] | mu + beta_bias * inv_sqrt_n[i], sigma_i);
  }}
}}
generated quantities {{
  real pooled = mu;
  real bias_slope = beta_bias;
  real mu_new = normal_rng(mu, tau);
}}")
  } else {
    stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] n_total;
}}
transformed data {{
  vector[N] inv_sqrt_n;
  for (i in 1:N)
    inv_sqrt_n[i] = 1.0 / sqrt(n_total[i]);
}}
parameters {{
  real mu;
  real beta_bias;
}}
model {{
  {mu_tgt}
  {beta_tgt}
  for (i in 1:N)
    target += normal_lpdf(y[i] | mu + beta_bias * inv_sqrt_n[i], se[i]);
}}
generated quantities {{
  real pooled = mu;
  real bias_slope = beta_bias;
}}")
  }

  stan_data <- list(N = S, y = es$yi, se = es$sei,
                    n_total = as.array(n_c + n_i))
  list(stan_code = as.character(stan_code), stan_data = stan_data)
}


# ---- H1 / Weight function (RE) ----
spec_component_selection_weight <- function(es, S, priors, p_cuts) {
  K        <- length(p_cuts) + 1L
  mu_tgt   <- emit_prior_target(priors$mu, "mu")
  tau_tgt  <- emit_prior_target(priors$tau, "tau")
  tau_bnds <- emit_prior_bounds(priors$tau, default_lower = 0)

  stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=2> K;
  vector[K-1] p_cutoffs;
}}
transformed data {{
  vector[K+1] z_bounds;
  z_bounds[1] = 8.0;
  z_bounds[K+1] = -8.0;
  for (k in 1:(K-1))
    z_bounds[k+1] = inv_Phi(1.0 - p_cutoffs[k]);
}}
parameters {{
  real mu;
  real{tau_bnds} tau;
  vector<lower=0.01, upper=0.99>[K-1] omega_raw;
}}
transformed parameters {{
  vector[K] omega;
  vector[N] sigma;
  omega[1] = 1.0;
  for (k in 1:(K-1))
    omega[k+1] = omega_raw[k];
  for (i in 1:N)
    sigma[i] = sqrt(square(tau) + square(se[i]));
}}
model {{
  {mu_tgt}
  {tau_tgt}
  for (k in 1:(K-1))
    target += beta_lpdf(omega_raw[k] | 1, 1);
  for (i in 1:N) {{
    real z_i = y[i] / se[i];
    real p_val = 1.0 - Phi(z_i);
    real log_w;
    if (p_val < p_cutoffs[1]) {{
      log_w = 0.0;
    }} else if (K == 2) {{
      log_w = log(omega[2]);
    }} else {{
      log_w = log(omega[K]);
      for (k in 1:(K-2)) {{
        if (p_val >= p_cutoffs[k] && p_val < p_cutoffs[k+1])
          log_w = log(omega[k+1]);
      }}
    }}
    target += normal_lpdf(y[i] | mu, sigma[i]) + log_w;
    {{
      real norm_c = 0;
      for (k in 1:K) {{
        real prob_k = Phi((z_bounds[k] * se[i] - mu) / sigma[i])
                    - Phi((z_bounds[k+1] * se[i] - mu) / sigma[i]);
        norm_c += omega[k] * fmax(prob_k, 1e-15);
      }}
      target += -log(fmax(norm_c, 1e-15));
    }}
  }}
}}
generated quantities {{
  real pooled = mu;
  vector[K] weights = omega;
  real mu_new = normal_rng(mu, tau);
}}")

  stan_data <- list(N = S, y = es$yi, se = es$sei, K = K, p_cutoffs = p_cuts)
  list(stan_code = as.character(stan_code), stan_data = stan_data)
}


# ---- H1 / Weight function (FE) ----
spec_component_selection_weight_fe <- function(es, S, priors, p_cuts) {
  K      <- length(p_cuts) + 1L
  mu_tgt <- emit_prior_target(priors$mu, "mu")

  stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=2> K;
  vector[K-1] p_cutoffs;
}}
transformed data {{
  vector[K+1] z_bounds;
  z_bounds[1] = 8.0;
  z_bounds[K+1] = -8.0;
  for (k in 1:(K-1))
    z_bounds[k+1] = inv_Phi(1.0 - p_cutoffs[k]);
}}
parameters {{
  real mu;
  vector<lower=0.01, upper=0.99>[K-1] omega_raw;
}}
transformed parameters {{
  vector[K] omega;
  omega[1] = 1.0;
  for (k in 1:(K-1))
    omega[k+1] = omega_raw[k];
}}
model {{
  {mu_tgt}
  for (k in 1:(K-1))
    target += beta_lpdf(omega_raw[k] | 1, 1);
  for (i in 1:N) {{
    real z_i = y[i] / se[i];
    real p_val = 1.0 - Phi(z_i);
    real log_w;
    if (p_val < p_cutoffs[1]) {{
      log_w = 0.0;
    }} else if (K == 2) {{
      log_w = log(omega[2]);
    }} else {{
      log_w = log(omega[K]);
      for (k in 1:(K-2)) {{
        if (p_val >= p_cutoffs[k] && p_val < p_cutoffs[k+1])
          log_w = log(omega[k+1]);
      }}
    }}
    target += normal_lpdf(y[i] | mu, se[i]) + log_w;
    {{
      real norm_c = 0;
      for (k in 1:K) {{
        real prob_k = Phi((z_bounds[k] * se[i] - mu) / se[i])
                    - Phi((z_bounds[k+1] * se[i] - mu) / se[i]);
        norm_c += omega[k] * fmax(prob_k, 1e-15);
      }}
      target += -log(fmax(norm_c, 1e-15));
    }}
  }}
}}
generated quantities {{
  real pooled = mu;
  vector[K] weights = omega;
}}")

  stan_data <- list(N = S, y = es$yi, se = es$sei, K = K, p_cutoffs = p_cuts)
  list(stan_code = as.character(stan_code), stan_data = stan_data)
}


# ---- H1 / Copas (FE or RE) ----
spec_component_selection_copas <- function(es, S, priors, has_re) {
  sp <- priors$selection %||% list()
  if (is.null(sp$gamma0)) sp$gamma0 <- normal(0, 2)
  if (is.null(sp$gamma1)) sp$gamma1 <- normal(0, 2)
  if (is.null(sp$rho))    sp$rho    <- uniform(-1, 1)
  mu_tgt   <- emit_prior_target(priors$mu, "mu")
  g0_tgt   <- emit_prior_target(sp$gamma0, "gamma0")
  g1_tgt   <- emit_prior_target(sp$gamma1, "gamma1")
  rho_tgt  <- emit_prior_target(sp$rho, "rho")
  rho_bnds <- emit_prior_bounds(sp$rho, default_lower = -1, default_upper = 1)
  g1_bnds  <- emit_prior_bounds(sp$gamma1, default_lower = 0)

  if (has_re) {
    tau_tgt  <- emit_prior_target(priors$tau, "tau")
    tau_bnds <- emit_prior_bounds(priors$tau, default_lower = 0)

    stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}}
parameters {{
  real mu;
  real{tau_bnds} tau;
  real gamma0;
  real{g1_bnds} gamma1;
  real{rho_bnds} rho;
}}
transformed parameters {{
  vector[N] sigma;
  vector[N] a;
  for (i in 1:N) {{
    sigma[i] = sqrt(square(tau) + square(se[i]));
    a[i] = gamma0 + gamma1 / se[i];
  }}
}}
model {{
  {mu_tgt}
  {tau_tgt}
  {g0_tgt}
  {g1_tgt}
  {rho_tgt}
  for (i in 1:N) {{
    real u_i = (y[i] - mu) / sigma[i];
    real sel_arg = (a[i] + rho * u_i) / sqrt(1 - square(rho));
    target += normal_lpdf(y[i] | mu, sigma[i]);
    target += log(Phi(sel_arg));
    target += -log(Phi(a[i] / sqrt(1 + square(rho) * square(tau)
              / (square(tau) + square(se[i])))));
  }}
}}
generated quantities {{
  real pooled = mu;
  real mu_new = normal_rng(mu, tau);
}}")
  } else {
    stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}}
parameters {{
  real mu;
  real gamma0;
  real{g1_bnds} gamma1;
  real{rho_bnds} rho;
}}
transformed parameters {{
  vector[N] a;
  for (i in 1:N)
    a[i] = gamma0 + gamma1 / se[i];
}}
model {{
  {mu_tgt}
  {g0_tgt}
  {g1_tgt}
  {rho_tgt}
  for (i in 1:N) {{
    real u_i = (y[i] - mu) / se[i];
    real sel_arg = (a[i] + rho * u_i) / sqrt(1 - square(rho));
    target += normal_lpdf(y[i] | mu, se[i]);
    target += log(Phi(sel_arg));
    target += -log(Phi(a[i]));
  }}
}}
generated quantities {{
  real pooled = mu;
}}")
  }

  stan_data <- list(N = S, y = es$yi, se = es$sei)
  list(stan_code = as.character(stan_code), stan_data = stan_data)
}


# ---- H1 / Jung bias-corrected (FE or RE) ----
spec_component_bias_corrected <- function(es, S, priors, has_re) {
  if (has_re) {
    stan_code <- generate_stan_code_bias_corrected(priors)
  } else {
    mu_tgt  <- emit_prior_target(priors$mu, "mu")
    b_prior <- priors$bias %||% priors$b %||% uniform(0, 2)
    b_tgt   <- emit_prior_target(b_prior, "B")
    b_bnds  <- emit_prior_bounds(b_prior, default_lower = 0)
    pb_prior <- priors$p_bias %||% beta(1, 1)
    pb_tgt   <- emit_prior_target(pb_prior, "p_bias")

    stan_code <- as.character(glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se_y;
  int<lower=0, upper=1> use_known_bias;
  array[N] int<lower=0, upper=1> known_bias;
}}
parameters {{
  real mu;
  real{b_bnds} B;
  real<lower=0, upper=1> p_bias;
}}
model {{
  {mu_tgt}
  {b_tgt}
  {pb_tgt}
  for (i in 1:N) {{
    if (use_known_bias == 1) {{
      if (known_bias[i] == 1) {{
        target += normal_lpdf(y[i] | mu + B, se_y[i]);
      }} else {{
        target += normal_lpdf(y[i] | mu, se_y[i]);
      }}
    }} else {{
      real lp_u = log1m(p_bias) + normal_lpdf(y[i] | mu, se_y[i]);
      real lp_b = log(p_bias + 1e-15)
                + normal_lpdf(y[i] | mu + B, se_y[i]);
      target += log_sum_exp(lp_u, lp_b);
    }}
  }}
}}
generated quantities {{
  real pooled = mu;
}}"))
  }

  stan_data <- list(N = S, y = es$yi, se_y = es$sei,
                    use_known_bias = 0L, known_bias = rep(0L, S))
  list(stan_code = stan_code, stan_data = stan_data)
}


# ============================================================================
# Spec builders: H0 effect + bias H1 components
# ============================================================================


# ---- H0 / PET (mu=0, beta_bias on 1/sqrt(n)) ----
spec_component_pet_peese_h0 <- function(es, S, n_c, n_i, has_re) {
  if (has_re) {
    stan_code <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] n_total;
}
transformed data {
  vector[N] inv_sqrt_n;
  for (i in 1:N)
    inv_sqrt_n[i] = 1.0 / sqrt(n_total[i]);
}
parameters {
  real<lower=0> tau;
  real beta_bias;
}
model {
  target += student_t_lpdf(tau | 3, 0, 2.5)
          - student_t_lccdf(0 | 3, 0, 2.5);
  target += normal_lpdf(beta_bias | 0, 5);
  for (i in 1:N) {
    real sigma_i = sqrt(square(tau) + square(se[i]));
    target += normal_lpdf(y[i] | beta_bias * inv_sqrt_n[i], sigma_i);
  }
}
generated quantities {
  real mu = 0.0;
  real bias_slope = beta_bias;
}"
  } else {
    stan_code <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] n_total;
}
transformed data {
  vector[N] inv_sqrt_n;
  for (i in 1:N)
    inv_sqrt_n[i] = 1.0 / sqrt(n_total[i]);
}
parameters {
  real beta_bias;
}
model {
  target += normal_lpdf(beta_bias | 0, 5);
  for (i in 1:N)
    target += normal_lpdf(y[i] | beta_bias * inv_sqrt_n[i], se[i]);
}
generated quantities {
  real mu = 0.0;
  real bias_slope = beta_bias;
}"
  }

stan_data <- list(N = S, y = es$yi, se = es$sei,
                  n_total = as.array(n_c + n_i))
list(stan_code = stan_code, stan_data = stan_data)
}


# ---- H0 / PEESE (mu=0, beta_bias on 1/n) ----
spec_component_peese_h0 <- function(es, S, n_c, n_i, has_re) {
  if (has_re) {
    stan_code <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] n_total;
}
transformed data {
  vector[N] inv_n;
  for (i in 1:N)
    inv_n[i] = 1.0 / n_total[i];
}
parameters {
  real<lower=0> tau;
  real beta_bias;
}
model {
  target += student_t_lpdf(tau | 3, 0, 2.5)
          - student_t_lccdf(0 | 3, 0, 2.5);
  target += normal_lpdf(beta_bias | 0, 5);
  for (i in 1:N) {
    real sigma_i = sqrt(square(tau) + square(se[i]));
    target += normal_lpdf(y[i] | beta_bias * inv_n[i], sigma_i);
  }
}
generated quantities {
  real mu = 0.0;
  real bias_slope = beta_bias;
}"
  } else {
    stan_code <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] n_total;
}
transformed data {
  vector[N] inv_n;
  for (i in 1:N)
    inv_n[i] = 1.0 / n_total[i];
}
parameters {
  real beta_bias;
}
model {
  target += normal_lpdf(beta_bias | 0, 5);
  for (i in 1:N)
    target += normal_lpdf(y[i] | beta_bias * inv_n[i], se[i]);
}
generated quantities {
  real mu = 0.0;
  real bias_slope = beta_bias;
}"
  }

stan_data <- list(N = S, y = es$yi, se = es$sei,
                  n_total = as.array(n_c + n_i))
list(stan_code = stan_code, stan_data = stan_data)
}


# ---- H0 / Weight function (mu=0, selection weights active) ----
spec_component_selection_weight_h0 <- function(es, S, has_re, p_cuts) {
  K <- length(p_cuts) + 1L

  if (has_re) {
    stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=2> K;
  vector[K-1] p_cutoffs;
}}
transformed data {{
  vector[K+1] z_bounds;
  z_bounds[1] = 8.0;
  z_bounds[K+1] = -8.0;
  for (k in 1:(K-1))
    z_bounds[k+1] = inv_Phi(1.0 - p_cutoffs[k]);
}}
parameters {{
  real<lower=0> tau;
  vector<lower=0.01, upper=0.99>[K-1] omega_raw;
}}
transformed parameters {{
  vector[K] omega;
  omega[1] = 1.0;
  for (k in 1:(K-1))
    omega[k+1] = omega_raw[k];
}}
model {{
  target += student_t_lpdf(tau | 3, 0, 2.5)
          - student_t_lccdf(0 | 3, 0, 2.5);
  for (k in 1:(K-1))
    target += beta_lpdf(omega_raw[k] | 1, 1);
  for (i in 1:N) {{
    real sigma_i = sqrt(square(tau) + square(se[i]));
    real z_i = y[i] / se[i];
    real p_val = 1.0 - Phi(z_i);
    real log_w;
    if (p_val < p_cutoffs[1]) {{
      log_w = 0.0;
    }} else if (K == 2) {{
      log_w = log(omega[2]);
    }} else {{
      log_w = log(omega[K]);
      for (k in 1:(K-2)) {{
        if (p_val >= p_cutoffs[k] && p_val < p_cutoffs[k+1])
          log_w = log(omega[k+1]);
      }}
    }}
    target += normal_lpdf(y[i] | 0.0, sigma_i) + log_w;
    {{
      real norm_c = 0;
      for (k in 1:K) {{
        real prob_k = Phi(z_bounds[k] * se[i] / sigma_i)
                    - Phi(z_bounds[k+1] * se[i] / sigma_i);
        norm_c += omega[k] * fmax(prob_k, 1e-15);
      }}
      target += -log(fmax(norm_c, 1e-15));
    }}
  }}
}}
generated quantities {{
  real mu = 0.0;
  vector[K] weights = omega;
}}")
  } else {
    stan_code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=2> K;
  vector[K-1] p_cutoffs;
}}
transformed data {{
  vector[K+1] z_bounds;
  z_bounds[1] = 8.0;
  z_bounds[K+1] = -8.0;
  for (k in 1:(K-1))
    z_bounds[k+1] = inv_Phi(1.0 - p_cutoffs[k]);
}}
parameters {{
  vector<lower=0.01, upper=0.99>[K-1] omega_raw;
}}
transformed parameters {{
  vector[K] omega;
  omega[1] = 1.0;
  for (k in 1:(K-1))
    omega[k+1] = omega_raw[k];
}}
model {{
  for (k in 1:(K-1))
    target += beta_lpdf(omega_raw[k] | 1, 1);
  for (i in 1:N) {{
    real z_i = y[i] / se[i];
    real p_val = 1.0 - Phi(z_i);
    real log_w;
    if (p_val < p_cutoffs[1]) {{
      log_w = 0.0;
    }} else if (K == 2) {{
      log_w = log(omega[2]);
    }} else {{
      log_w = log(omega[K]);
      for (k in 1:(K-2)) {{
        if (p_val >= p_cutoffs[k] && p_val < p_cutoffs[k+1])
          log_w = log(omega[k+1]);
      }}
    }}
    target += normal_lpdf(y[i] | 0.0, se[i]) + log_w;
    {{
      real norm_c = 0;
      for (k in 1:K) {{
        real prob_k = Phi(z_bounds[k]) - Phi(z_bounds[k+1]);
        norm_c += omega[k] * fmax(prob_k, 1e-15);
      }}
      target += -log(fmax(norm_c, 1e-15));
    }}
  }}
}}
generated quantities {{
  real mu = 0.0;
  vector[K] weights = omega;
}}")
  }

  stan_data <- list(N = S, y = es$yi, se = es$sei,
                    K = K, p_cutoffs = p_cuts)
  list(stan_code = as.character(stan_code), stan_data = stan_data)
}


# ---- H0 / Jung bias-corrected (mu=0, bias shift B active) ----
spec_component_jung_h0 <- function(es, S, has_re) {
  if (has_re) {
    stan_code <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}
parameters {
  real<lower=0> tau;
  real<lower=0> B;
  real<lower=0, upper=1> p_bias;
}
model {
  target += student_t_lpdf(tau | 3, 0, 2.5)
          - student_t_lccdf(0 | 3, 0, 2.5);
  target += uniform_lpdf(B | 0, 2);
  target += beta_lpdf(p_bias | 1, 1);
  for (i in 1:N) {
    real sigma_i = sqrt(square(tau) + square(se[i]));
    real lp_u = log1m(p_bias) + normal_lpdf(y[i] | 0.0, sigma_i);
    real lp_b = log(p_bias + 1e-15) + normal_lpdf(y[i] | B, sigma_i);
    target += log_sum_exp(lp_u, lp_b);
  }
}
generated quantities {
  real mu = 0.0;
}"
  } else {
    stan_code <- "
data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}
parameters {
  real<lower=0> B;
  real<lower=0, upper=1> p_bias;
}
model {
  target += uniform_lpdf(B | 0, 2);
  target += beta_lpdf(p_bias | 1, 1);
  for (i in 1:N) {
    real lp_u = log1m(p_bias) + normal_lpdf(y[i] | 0.0, se[i]);
    real lp_b = log(p_bias + 1e-15) + normal_lpdf(y[i] | B, se[i]);
    target += log_sum_exp(lp_u, lp_b);
  }
}
generated quantities {
  real mu = 0.0;
}"
  }

stan_data <- list(N = S, y = es$yi, se = es$sei)
list(stan_code = stan_code, stan_data = stan_data)
}
# ============================================================================
# Source together with robma_main.R and robma_bridge_part1/2/3.R
#
# All three SS generators use the full 8-combination log_mix:
#   Effect ON/OFF x Heterogeneity ON/OFF x Bias ON/OFF
#
# When effect is OFF and bias is ON, the bias mechanism operates
# with mu = 0. This captures "the apparent effect is entirely
# due to publication bias."
# ============================================================================


robma_spike_slab <- function(
    models, es, S, n_c, n_i, study_labels, likelihood,
    bias_indicator, mu_prior, tau_prior,
    b_prior, p_bias_prior, p_cutoffs,
    chains, iter_warmup, iter_sampling, adapt_delta, seed, ...
) {
  cli::cli_alert_info(paste0(
    "RoBMA (spike-and-slab via log_mix): bias indicator = ",
    bias_indicator))

  if (bias_indicator == "pet_peese") {
    cli::cli_alert_warning(paste0(
      "PET-PEESE has identifiability issues in the spike-and-slab ",
      "framework: the effect (mu) and bias slope (beta * 1/sqrt(n)) ",
      "are confounded. Consider bias_indicator = 'bias_corrected' or ",
      "method = 'bridge'."))
  }

  cli::cli_alert_info(paste0(
    "Note: The spike-and-slab method is a fast approximation. ",
    "It typically produces more conservative inclusion BFs than ",
    "bridge sampling. For publication-quality results, use ",
    "method = 'bridge'."))

  stan_result <- switch(bias_indicator,
                        bias_corrected   = ss_stan_jung(
                          es, S, mu_prior, tau_prior, b_prior, p_bias_prior),
                        pet_peese        = ss_stan_pet_peese(
                          es, S, n_c, n_i, mu_prior, tau_prior),
                        selection_weight = ss_stan_selection_weight(
                          es, S, mu_prior, tau_prior, p_cutoffs))

  stan_file <- cmdstanr::write_stan_file(as.character(stan_result$code))
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = stan_result$data, chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, ...)

  draws <- posterior::as_draws_df(fit$draws())
  mu_draws    <- as.vector(draws$mu)
  averaged_mu <- tibble::tibble(mu = mu_draws, model = "spike_slab_joint")

  pip_effect <- mean(as.vector(draws$pip_effect))
  pip_hetero <- mean(as.vector(draws$pip_hetero))
  pip_bias   <- mean(as.vector(draws$pip_bias))

  pip_to_bf <- function(pip) {
    pip <- max(min(pip, 1 - 1e-10), 1e-10)
    pip / (1 - pip)
  }

  n_div <- tryCatch({
    sum(fit$diagnostic_summary()$num_divergent)
  }, error = function(e) NA_integer_)

  if (!is.na(n_div) && n_div > 0) {
    n_total <- iter_sampling * chains
    pct <- round(100 * n_div / n_total, 1)
    cli::cli_alert_warning(paste0(
      n_div, " of ", n_total, " (", pct, "%) transitions diverged. ",
      "Results may be unreliable. Consider increasing adapt_delta ",
      "or iter_warmup."))
  }

  model_table <- tibble::tibble(
    model = paste0("Spike-and-slab (", bias_indicator, ")"),
    log_ml = NA_real_, prior_weight = 1, post_prob = 1,
    null_effect = NA, has_bias = TRUE, has_heterogeneity = TRUE,
    pip_effect = pip_effect, pip_hetero = pip_hetero,
    pip_bias = pip_bias)

  list(
    averaged_draws  = averaged_mu,
    model_table     = model_table,
    inclusion_bf    = list(
      effect = pip_to_bf(pip_effect),
      bias = pip_to_bf(pip_bias),
      heterogeneity = pip_to_bf(pip_hetero)),
    posterior_probs = list(
      effect = pip_effect,
      heterogeneity = pip_hetero,
      bias = pip_bias),
    component_fits  = list(list(
      fit = fit, label = "spike_slab_joint",
      null_effect = FALSE, has_bias = TRUE, has_hetero = TRUE)),
    post_probs               = c(spike_slab_joint = 1),
    log_marginal_likelihoods = NA_real_,
    pip = list(
      effect = pip_effect,
      heterogeneity = pip_hetero,
      bias = pip_bias))
}


# ============================================================================
# Jung bias-corrected
# ============================================================================

ss_stan_jung <- function(es, S, mu_prior, tau_prior,
                         b_prior, p_bias_prior) {
  mu_tgt   <- emit_prior_target(mu_prior, "mu_raw")
  tau_tgt  <- emit_prior_target(tau_prior, "tau_raw")
  tau_bnds <- emit_prior_bounds(tau_prior, default_lower = 0)
  b_bnds   <- emit_prior_bounds(b_prior, default_lower = 0)
  b_tgt    <- emit_prior_target(b_prior, "B_raw")
  pb_tgt   <- emit_prior_target(p_bias_prior, "p_bias")

  code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
}}
parameters {{
  real<lower=0, upper=1> pip_effect;
  real<lower=0, upper=1> pip_hetero;
  real<lower=0, upper=1> pip_bias;
  real mu_raw;
  real{tau_bnds} tau_raw;
  real{b_bnds} B_raw;
  real<lower=0, upper=1> p_bias;
}}
model {{
  target += beta_lpdf(pip_effect | 1, 1);
  target += beta_lpdf(pip_hetero | 1, 1);
  target += beta_lpdf(pip_bias | 1, 1);
  {mu_tgt}
  {tau_tgt}
  {b_tgt}
  {pb_tgt}

  for (i in 1:N) {{
    real s_h1 = sqrt(square(tau_raw) + square(se[i]));
    real s_h0 = se[i];

    // === Effect ON, Hetero ON, Bias ON ===
    real ll_e1_h1_b1;
    {{
      real lp_u = log1m(p_bias)
                + normal_lpdf(y[i] | mu_raw, s_h1);
      real lp_b = log(p_bias + 1e-15)
                + normal_lpdf(y[i] | mu_raw + B_raw, s_h1);
      ll_e1_h1_b1 = log_sum_exp(lp_u, lp_b);
    }}
    // === Effect ON, Hetero ON, Bias OFF ===
    real ll_e1_h1_b0 = normal_lpdf(y[i] | mu_raw, s_h1);

    // === Effect ON, Hetero OFF, Bias ON ===
    real ll_e1_h0_b1;
    {{
      real lp_u = log1m(p_bias)
                + normal_lpdf(y[i] | mu_raw, s_h0);
      real lp_b = log(p_bias + 1e-15)
                + normal_lpdf(y[i] | mu_raw + B_raw, s_h0);
      ll_e1_h0_b1 = log_sum_exp(lp_u, lp_b);
    }}
    // === Effect ON, Hetero OFF, Bias OFF ===
    real ll_e1_h0_b0 = normal_lpdf(y[i] | mu_raw, s_h0);

    // === Effect OFF, Hetero ON, Bias ON ===
    real ll_e0_h1_b1;
    {{
      real lp_u = log1m(p_bias)
                + normal_lpdf(y[i] | 0.0, s_h1);
      real lp_b = log(p_bias + 1e-15)
                + normal_lpdf(y[i] | B_raw, s_h1);
      ll_e0_h1_b1 = log_sum_exp(lp_u, lp_b);
    }}
    // === Effect OFF, Hetero ON, Bias OFF ===
    real ll_e0_h1_b0 = normal_lpdf(y[i] | 0.0, s_h1);

    // === Effect OFF, Hetero OFF, Bias ON ===
    real ll_e0_h0_b1;
    {{
      real lp_u = log1m(p_bias)
                + normal_lpdf(y[i] | 0.0, s_h0);
      real lp_b = log(p_bias + 1e-15)
                + normal_lpdf(y[i] | B_raw, s_h0);
      ll_e0_h0_b1 = log_sum_exp(lp_u, lp_b);
    }}
    // === Effect OFF, Hetero OFF, Bias OFF ===
    real ll_e0_h0_b0 = normal_lpdf(y[i] | 0.0, s_h0);

    // Marginalise: bias -> hetero -> effect
    real ll_e1_h1 = log_mix(pip_bias, ll_e1_h1_b1, ll_e1_h1_b0);
    real ll_e1_h0 = log_mix(pip_bias, ll_e1_h0_b1, ll_e1_h0_b0);
    real ll_e1    = log_mix(pip_hetero, ll_e1_h1, ll_e1_h0);

    real ll_e0_h1 = log_mix(pip_bias, ll_e0_h1_b1, ll_e0_h1_b0);
    real ll_e0_h0 = log_mix(pip_bias, ll_e0_h0_b1, ll_e0_h0_b0);
    real ll_e0    = log_mix(pip_hetero, ll_e0_h1, ll_e0_h0);

    target += log_mix(pip_effect, ll_e1, ll_e0);
  }}
}}
generated quantities {{
  real mu;
  if (bernoulli_rng(pip_effect)) {{
    mu = mu_raw;
  }} else {{
    mu = 0.0;
  }}
}}")

  list(code = code, data = list(N = S, y = es$yi, se = es$sei))
}


# ============================================================================
# PET-PEESE
# ============================================================================

ss_stan_pet_peese <- function(es, S, n_c, n_i,
                              mu_prior, tau_prior) {
  mu_tgt   <- emit_prior_target(mu_prior, "mu_raw")
  tau_tgt  <- emit_prior_target(tau_prior, "tau_raw")
  tau_bnds <- emit_prior_bounds(tau_prior, default_lower = 0)

  code <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector<lower=0>[N] inv_sqrt_n;
}}
parameters {{
  real<lower=0, upper=1> pip_effect;
  real<lower=0, upper=1> pip_hetero;
  real<lower=0, upper=1> pip_bias;
  real mu_raw;
  real{tau_bnds} tau_raw;
  real beta_bias_raw;
}}
model {{
  target += beta_lpdf(pip_effect | 1, 1);
  target += beta_lpdf(pip_hetero | 1, 1);
  target += beta_lpdf(pip_bias | 1, 1);
  {mu_tgt}
  {tau_tgt}
  target += normal_lpdf(beta_bias_raw | 0, 5);

  for (i in 1:N) {{
    real s_h1 = sqrt(square(tau_raw) + square(se[i]));
    real s_h0 = se[i];
    real bias_i = beta_bias_raw * inv_sqrt_n[i];

    // === Effect ON ===
    real ll_e1_h1_b1 = normal_lpdf(y[i] | mu_raw + bias_i, s_h1);
    real ll_e1_h1_b0 = normal_lpdf(y[i] | mu_raw, s_h1);
    real ll_e1_h0_b1 = normal_lpdf(y[i] | mu_raw + bias_i, s_h0);
    real ll_e1_h0_b0 = normal_lpdf(y[i] | mu_raw, s_h0);

    // === Effect OFF (mu=0, bias slope still active) ===
    real ll_e0_h1_b1 = normal_lpdf(y[i] | bias_i, s_h1);
    real ll_e0_h1_b0 = normal_lpdf(y[i] | 0.0, s_h1);
    real ll_e0_h0_b1 = normal_lpdf(y[i] | bias_i, s_h0);
    real ll_e0_h0_b0 = normal_lpdf(y[i] | 0.0, s_h0);

    // Marginalise: bias -> hetero -> effect
    real ll_e1_h1 = log_mix(pip_bias, ll_e1_h1_b1, ll_e1_h1_b0);
    real ll_e1_h0 = log_mix(pip_bias, ll_e1_h0_b1, ll_e1_h0_b0);
    real ll_e1    = log_mix(pip_hetero, ll_e1_h1, ll_e1_h0);

    real ll_e0_h1 = log_mix(pip_bias, ll_e0_h1_b1, ll_e0_h1_b0);
    real ll_e0_h0 = log_mix(pip_bias, ll_e0_h0_b1, ll_e0_h0_b0);
    real ll_e0    = log_mix(pip_hetero, ll_e0_h1, ll_e0_h0);

    target += log_mix(pip_effect, ll_e1, ll_e0);
  }}
}}
generated quantities {{
  real mu;
  if (bernoulli_rng(pip_effect)) {{
    mu = mu_raw;
  }} else {{
    mu = 0.0;
  }}
}}")

  list(
    code = code,
    data = list(
      N = S, y = es$yi, se = es$sei,
      inv_sqrt_n = 1 / sqrt(n_c + n_i)))
}


# ============================================================================
# Selection weight
# ============================================================================

ss_stan_selection_weight <- function(es, S, mu_prior,
                                     tau_prior, p_cutoffs) {
  K        <- length(p_cutoffs) + 1L
  mu_tgt   <- emit_prior_target(mu_prior, "mu_raw")
  tau_tgt  <- emit_prior_target(tau_prior, "tau_raw")
  tau_bnds <- emit_prior_bounds(tau_prior, default_lower = 0)

  # Split into two parts to stay within length limits
  code_part1 <- glue::glue("
data {{
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  int<lower=2> K;
  vector[K-1] p_cutoffs;
}}
transformed data {{
  vector[K+1] z_bounds;
  z_bounds[1] = 8.0;
  z_bounds[K+1] = -8.0;
  for (k in 1:(K-1))
    z_bounds[k+1] = inv_Phi(1.0 - p_cutoffs[k]);
}}
parameters {{
  real<lower=0, upper=1> pip_effect;
  real<lower=0, upper=1> pip_hetero;
  real<lower=0, upper=1> pip_bias;
  real mu_raw;
  real{tau_bnds} tau_raw;
  vector<lower=0.01, upper=0.99>[K-1] omega_raw;
}}
model {{
  target += beta_lpdf(pip_effect | 1, 1);
  target += beta_lpdf(pip_hetero | 1, 1);
  target += beta_lpdf(pip_bias | 1, 1);
  {mu_tgt}
  {tau_tgt}
  for (k in 1:(K-1))
    target += beta_lpdf(omega_raw[k] | 1, 1);

  for (i in 1:N) {{
    real s_h1 = sqrt(square(tau_raw) + square(se[i]));
    real s_h0 = se[i];

    // Build omega vector and compute selection weight
    real z_i = y[i] / se[i];
    real p_val = 1.0 - Phi(z_i);
    real log_w;
    if (p_val < p_cutoffs[1]) {{
      log_w = 0.0;
    }} else if (K == 2) {{
      log_w = log(omega_raw[1]);
    }} else {{
      log_w = log(omega_raw[K-1]);
      for (k in 1:(K-2)) {{
        if (p_val >= p_cutoffs[k] && p_val < p_cutoffs[k+1])
          log_w = log(omega_raw[k]);
      }}
    }}

    vector[K] omega_vec;
    omega_vec[1] = 1.0;
    for (k in 1:(K-1))
      omega_vec[k+1] = omega_raw[k];")

  code_part2 <- glue::glue("
    // === Effect ON, Bias ON: normalisation for mu_raw ===
    real nc_e1_h1 = 0;
    real nc_e1_h0 = 0;
    for (k in 1:K) {{
      real p_h1 = Phi((z_bounds[k] * se[i] - mu_raw) / s_h1)
                - Phi((z_bounds[k+1] * se[i] - mu_raw) / s_h1);
      real p_h0 = Phi((z_bounds[k] * se[i] - mu_raw) / s_h0)
                - Phi((z_bounds[k+1] * se[i] - mu_raw) / s_h0);
      nc_e1_h1 += omega_vec[k] * fmax(p_h1, 1e-15);
      nc_e1_h0 += omega_vec[k] * fmax(p_h0, 1e-15);
    }}
    real ll_e1_h1_b1 = normal_lpdf(y[i] | mu_raw, s_h1)
                     + log_w - log(fmax(nc_e1_h1, 1e-15));
    real ll_e1_h0_b1 = normal_lpdf(y[i] | mu_raw, s_h0)
                     + log_w - log(fmax(nc_e1_h0, 1e-15));

    // === Effect ON, Bias OFF ===
    real ll_e1_h1_b0 = normal_lpdf(y[i] | mu_raw, s_h1);
    real ll_e1_h0_b0 = normal_lpdf(y[i] | mu_raw, s_h0);

    // === Effect OFF, Bias ON: normalisation for mu=0 ===
    real nc_e0_h1 = 0;
    real nc_e0_h0 = 0;
    for (k in 1:K) {{
      real p_h1 = Phi(z_bounds[k] * se[i] / s_h1)
                - Phi(z_bounds[k+1] * se[i] / s_h1);
      real p_h0 = Phi(z_bounds[k])
                - Phi(z_bounds[k+1]);
      nc_e0_h1 += omega_vec[k] * fmax(p_h1, 1e-15);
      nc_e0_h0 += omega_vec[k] * fmax(p_h0, 1e-15);
    }}
    real ll_e0_h1_b1 = normal_lpdf(y[i] | 0.0, s_h1)
                     + log_w - log(fmax(nc_e0_h1, 1e-15));
    real ll_e0_h0_b1 = normal_lpdf(y[i] | 0.0, s_h0)
                     + log_w - log(fmax(nc_e0_h0, 1e-15));

    // === Effect OFF, Bias OFF ===
    real ll_e0_h1_b0 = normal_lpdf(y[i] | 0.0, s_h1);
    real ll_e0_h0_b0 = normal_lpdf(y[i] | 0.0, s_h0);

    // Marginalise: bias -> hetero -> effect
    real ll_e1_h1 = log_mix(pip_bias, ll_e1_h1_b1, ll_e1_h1_b0);
    real ll_e1_h0 = log_mix(pip_bias, ll_e1_h0_b1, ll_e1_h0_b0);
    real ll_e1    = log_mix(pip_hetero, ll_e1_h1, ll_e1_h0);

    real ll_e0_h1 = log_mix(pip_bias, ll_e0_h1_b1, ll_e0_h1_b0);
    real ll_e0_h0 = log_mix(pip_bias, ll_e0_h0_b1, ll_e0_h0_b0);
    real ll_e0    = log_mix(pip_hetero, ll_e0_h1, ll_e0_h0);

    target += log_mix(pip_effect, ll_e1, ll_e0);
  }}
}}
generated quantities {{
  real mu;
  if (bernoulli_rng(pip_effect)) {{
    mu = mu_raw;
  }} else {{
    mu = 0.0;
  }}
  vector[K] weights;
  weights[1] = 1.0;
  for (k in 1:(K-1))
    weights[k+1] = omega_raw[k];
}}")

  code <- paste0(as.character(code_part1), "\n",
                 as.character(code_part2))

  list(
    code = code,
    data = list(
      N = S, y = es$yi, se = es$sei,
      K = K, p_cutoffs = sort(p_cutoffs)))
}
