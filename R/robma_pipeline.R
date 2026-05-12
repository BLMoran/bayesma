# robma modular pipeline
#
# Six stages mirroring bayesma_pipeline.R:
#   1. robma_spec()      -- validate + extract -> bayesma_robma_spec
#   2. robma_stan_code() -- spec -> named Stan programs (per model label)
#   3. robma_stan_data() -- spec -> named cmdstanr data lists (per label)
#   4. robma_fit()       -- compile + sample + bridge / PIP extraction
#   5. robma_extract()   -- fit + spec -> null-range probs + forest df
#   6. robma_output()    -- assemble final bayesma_robma object


# -----------------------------------------------------------------------------
# Stage 1: robma_spec
# -----------------------------------------------------------------------------

#' @noRd
robma_spec <- function(
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
    b_prior = NULL,
    p_bias_prior = NULL,
    p_cutoffs = c(0.025, 0.05),
    horseshoe = FALSE,
    parallel = FALSE,
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    seed = 1234,
    quiet = FALSE,
    custom_model = NULL,
    custom_data = NULL,
    ...
) {
  likelihood     <- rlang::arg_match(likelihood)
  method         <- rlang::arg_match(method)
  bias_indicator <- rlang::arg_match(bias_indicator)

  priors_effect             <- priors_effect             %||% robma_default_priors("effect", rescale = rescale_priors)
  priors_effect_null        <- priors_effect_null        %||% robma_default_priors("effect", null = TRUE)
  priors_heterogeneity      <- priors_heterogeneity      %||% robma_default_priors("heterogeneity", rescale = rescale_priors)
  priors_heterogeneity_null <- priors_heterogeneity_null %||% robma_default_priors("heterogeneity", null = TRUE)
  priors_bias               <- priors_bias               %||% robma_default_priors("bias", rescale = rescale_priors)
  priors_bias_null          <- priors_bias_null          %||% robma_default_priors("bias", null = TRUE)
  b_prior      <- b_prior      %||% uniform(0, 2)
  p_bias_prior <- p_bias_prior %||% beta(1, 1)

  if (!is.null(null_range)) {
    if (!is.numeric(null_range) || length(null_range) != 2) {
      cli::cli_abort("{.arg null_range} must be a numeric vector of length 2.",
                     call = rlang::caller_env())
    }
    if (null_range[1] > null_range[2]) {
      cli::cli_abort("{.arg null_range[1]} must be <= {.arg null_range[2]}.",
                     call = rlang::caller_env())
    }
  }

  if (!is.logical(horseshoe) || length(horseshoe) != 1) {
    cli::cli_abort("{.arg horseshoe} must be a single logical value.",
                   call = rlang::caller_env())
  }
  if (isTRUE(horseshoe) && method != "ss") {
    cli::cli_abort("{.arg horseshoe} is only available for {.code method = 'ss'}.",
                   call = rlang::caller_env())
  }
  if (isTRUE(horseshoe) && bias_indicator != "bias_corrected") {
    cli::cli_abort(
      "{.arg horseshoe} is only supported with {.code bias_indicator = 'bias_corrected'}.",
      call = rlang::caller_env()
    )
  }

  validate_custom_model(custom_model, method)
  validate_custom_data(custom_data, method)

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
  S            <- length(study_vec)
  study_labels <- if (is.factor(study_vec)) levels(study_vec) else as.character(study_vec)
  n_c <- extract_col(data, n_ctrl)
  n_i <- extract_col(data, n_int)

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

  model_grid <- build_model_grid(
    priors_effect, priors_effect_null,
    priors_heterogeneity, priors_heterogeneity_null,
    priors_bias, priors_bias_null
  )

  if (length(model_grid) < 2) {
    cli::cli_abort("RoBMA requires at least 2 models in the grid.",
                   call = rlang::caller_env())
  }

  effect_label <- switch(likelihood,
                         binomial = "log_or",
                         gaussian = "mean_diff",
                         poisson  = "log_rr")

  dots <- list(...)
  call_args <- c(
    list(
      studyvar = studyvar,
      event_ctrl = event_ctrl, event_int = event_int,
      mean_ctrl = mean_ctrl,   mean_int = mean_int,
      sd_ctrl = sd_ctrl,       sd_int = sd_int,
      n_ctrl = n_ctrl,         n_int = n_int,
      likelihood = likelihood,
      method = method,
      bias_indicator = bias_indicator,
      priors_effect = priors_effect,
      priors_effect_null = priors_effect_null,
      priors_heterogeneity = priors_heterogeneity,
      priors_heterogeneity_null = priors_heterogeneity_null,
      priors_bias = priors_bias,
      priors_bias_null = priors_bias_null,
      rescale_priors = rescale_priors,
      null_range = null_range,
      b_prior = b_prior,
      p_bias_prior = p_bias_prior,
      p_cutoffs = p_cutoffs,
      horseshoe = horseshoe,
      parallel = parallel,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      seed = seed
    ),
    dots
  )

  spec <- list(
    likelihood       = likelihood,
    method           = method,
    bias_indicator   = bias_indicator,
    study_vec        = study_vec,
    study_labels     = study_labels,
    S                = S,
    outcome_ctrl     = outcome_ctrl,
    outcome_int      = outcome_int,
    n_c              = n_c,
    n_i              = n_i,
    sd_c             = sd_c,
    sd_i             = sd_i,
    es               = es,
    effect_label     = effect_label,
    model_grid       = model_grid,
    priors_effect             = priors_effect,
    priors_effect_null        = priors_effect_null,
    priors_heterogeneity      = priors_heterogeneity,
    priors_heterogeneity_null = priors_heterogeneity_null,
    priors_bias               = priors_bias,
    priors_bias_null          = priors_bias_null,
    rescale_priors   = rescale_priors,
    null_range       = null_range,
    b_prior          = b_prior,
    p_bias_prior     = p_bias_prior,
    p_cutoffs        = p_cutoffs,
    horseshoe        = horseshoe,
    parallel         = parallel,
    quiet            = quiet,
    custom_model     = custom_model,
    custom_data      = custom_data,
    call_args        = call_args
  )
  class(spec) <- c("bayesma_robma_spec", "bayesma_spec", "list")
  spec
}


#' @noRd
validate_custom_model <- function(custom_model, method,
                                  call = rlang::caller_env()) {
  if (is.null(custom_model)) return(invisible())

  if (is.character(custom_model) && length(custom_model) == 1) {
    return(invisible())
  }
  if (is.list(custom_model)) {
    if (is.null(names(custom_model)) || any(!nzchar(names(custom_model)))) {
      cli::cli_abort(
        "{.arg custom_model} list must be named by model label.",
        call = call
      )
    }
    ok <- purrr::map_lgl(custom_model,
                         \(x) is.character(x) && length(x) == 1)
    if (!all(ok)) {
      cli::cli_abort(
        "Each entry of {.arg custom_model} must be a character scalar.",
        call = call
      )
    }
    return(invisible())
  }
  cli::cli_abort(
    "{.arg custom_model} must be a character scalar or a named list of character scalars.",
    call = call
  )
}

#' @noRd
validate_custom_data <- function(custom_data, method,
                                 call = rlang::caller_env()) {
  if (is.null(custom_data)) return(invisible())
  if (!is.list(custom_data)) {
    cli::cli_abort("{.arg custom_data} must be a named list.", call = call)
  }
  if (is.null(names(custom_data)) || any(!nzchar(names(custom_data)))) {
    cli::cli_abort("{.arg custom_data} must be a named list.", call = call)
  }
  invisible()
}


#' @export
print.bayesma_robma_spec <- function(x, ...) {
  cat("<bayesma_robma_spec>\n",
      "  likelihood     : ", x$likelihood, "\n",
      "  method         : ", x$method, "\n",
      "  bias_indicator : ", x$bias_indicator, "\n",
      "  studies (S)    : ", x$S, "\n",
      "  n_models       : ", length(x$model_grid), "\n",
      "  custom_model   : ", !is.null(x$custom_model), "\n",
      sep = "")
  invisible(x)
}


# -----------------------------------------------------------------------------
# Stage 2: robma_stan_code
# -----------------------------------------------------------------------------

#' @noRd
robma_stan_code <- function(spec, format = TRUE) {
  if (!inherits(spec, "bayesma_robma_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_robma_spec} object.",
                   call = rlang::caller_env())
  }

  raw <- if (spec$method == "bridge") {
    robma_stan_code_bridge(spec)
  } else {
    robma_stan_code_ss(spec)
  }

  raw <- apply_custom_model(raw, spec$custom_model, spec$method)

  out <- purrr::map(raw, \(entry) finalise_code_entry(entry, format))
  class(out) <- c("bayesma_robma_stan_code", "list")
  out
}


#' @noRd
robma_stan_code_bridge <- function(spec) {
  purrr::map(spec$model_grid, function(g) {
    bs <- build_model_spec(g, spec$es, spec$S, spec$n_c, spec$n_i,
                           spec$likelihood)
    if (is.null(bs$analytic)) bs$analytic <- FALSE
    list(
      label           = bs$label,
      code            = bs$stan_code,
      analytic        = isTRUE(bs$analytic),
      log_ml_override = bs$log_ml_override,
      is_effect_null  = isTRUE(bs$is_effect_null),
      is_hetero_null  = isTRUE(bs$is_hetero_null),
      is_bias_null    = isTRUE(bs$is_bias_null),
      prior_weight    = bs$prior_weight,
      hetero_prior_for_fallback = bs$hetero_prior_for_fallback
    )
  }) |>
    stats::setNames(purrr::map_chr(spec$model_grid, \(g) g$label))
}


#' @noRd
robma_stan_code_ss <- function(spec) {
  mu_prior  <- spec$priors_effect[[1]]
  tau_prior <- spec$priors_heterogeneity[[1]]

  if (isTRUE(spec$horseshoe)) {
    res <- ss_stan_horseshoe(
      spec$es, spec$S, mu_prior, tau_prior,
      spec$b_prior, spec$p_bias_prior)
    label <- "horseshoe"
  } else {
    res <- switch(
      spec$bias_indicator,
      bias_corrected   = ss_stan_jung(
        spec$es, spec$S, mu_prior, tau_prior,
        spec$b_prior, spec$p_bias_prior),
      pet_peese        = ss_stan_pet_peese(
        spec$es, spec$S, spec$n_c, spec$n_i, mu_prior, tau_prior),
      selection_weight = ss_stan_selection_weight(
        spec$es, spec$S, mu_prior, tau_prior, spec$p_cutoffs)
    )
    label <- spec$bias_indicator
  }

  entry <- list(
    label    = label,
    code     = as.character(res$code),
    analytic = FALSE,
    ss_data  = res$data
  )
  stats::setNames(list(entry), label)
}


#' @noRd
apply_custom_model <- function(raw, custom_model, method) {
  if (is.null(custom_model)) return(raw)

  if (method == "ss") {
    if (is.character(custom_model) && length(custom_model) == 1) {
      key <- names(raw)[[1]]
      raw[[key]]$code     <- custom_model
      raw[[key]]$analytic <- FALSE
      return(raw)
    }
    if (is.list(custom_model)) {
      key <- names(raw)[[1]]
      if (!key %in% names(custom_model)) {
        cli::cli_abort(
          "{.arg custom_model} list has no entry for label {.val {key}}.",
          call = rlang::caller_env()
        )
      }
      raw[[key]]$code     <- custom_model[[key]]
      raw[[key]]$analytic <- FALSE
      return(raw)
    }
  }

  if (is.character(custom_model) && length(custom_model) == 1) {
    raw <- purrr::map(raw, function(entry) {
      entry$code     <- custom_model
      entry$analytic <- FALSE
      entry$log_ml_override <- NULL
      entry
    })
    return(raw)
  }

  if (is.list(custom_model)) {
    unknown <- setdiff(names(custom_model), names(raw))
    if (length(unknown) > 0) {
      cli::cli_abort(
        "Unknown label(s) in {.arg custom_model}: {.val {unknown}}.",
        call = rlang::caller_env()
      )
    }
    for (nm in names(custom_model)) {
      raw[[nm]]$code     <- custom_model[[nm]]
      raw[[nm]]$analytic <- FALSE
      raw[[nm]]$log_ml_override <- NULL
    }
    return(raw)
  }
  raw
}


#' @noRd
finalise_code_entry <- function(entry, format) {
  if (isTRUE(entry$analytic) || is.null(entry$code)) {
    entry$full   <- NA_character_
    entry$blocks <- NULL
    return(entry)
  }

  full <- if (isTRUE(format)) format_stan_code(entry$code) else entry$code
  entry$full   <- full
  entry$blocks <- parse_stan_blocks(full)
  entry
}


#' @export
print.bayesma_robma_stan_code <- function(x, ...) {
  for (nm in names(x)) {
    cat("\n=== ", nm, " ===\n", sep = "")
    entry <- x[[nm]]
    if (isTRUE(entry$analytic)) {
      cat("(analytic -- no Stan program)\n")
    } else {
      cat(entry$full, "\n", sep = "")
    }
  }
  invisible(x)
}


# -----------------------------------------------------------------------------
# Stage 3: robma_stan_data
# -----------------------------------------------------------------------------

#' @noRd
robma_stan_data <- function(spec) {
  if (!inherits(spec, "bayesma_robma_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_robma_spec} object.",
                   call = rlang::caller_env())
  }

  out <- if (spec$method == "bridge") {
    purrr::map(spec$model_grid, function(g) {
      bs <- build_model_spec(g, spec$es, spec$S, spec$n_c, spec$n_i,
                             spec$likelihood)
      bs$stan_data
    }) |>
      stats::setNames(purrr::map_chr(spec$model_grid, \(g) g$label))
  } else {
    mu_prior  <- spec$priors_effect[[1]]
    tau_prior <- spec$priors_heterogeneity[[1]]
    if (isTRUE(spec$horseshoe)) {
      res   <- ss_stan_horseshoe(
        spec$es, spec$S, mu_prior, tau_prior,
        spec$b_prior, spec$p_bias_prior)
      label <- "horseshoe"
    } else {
      res <- switch(
        spec$bias_indicator,
        bias_corrected   = ss_stan_jung(
          spec$es, spec$S, mu_prior, tau_prior,
          spec$b_prior, spec$p_bias_prior),
        pet_peese        = ss_stan_pet_peese(
          spec$es, spec$S, spec$n_c, spec$n_i, mu_prior, tau_prior),
        selection_weight = ss_stan_selection_weight(
          spec$es, spec$S, mu_prior, tau_prior, spec$p_cutoffs)
      )
      label <- spec$bias_indicator
    }
    stats::setNames(list(res$data), label)
  }

  out <- overlay_custom_data(out, spec$custom_data, spec$method)
  class(out) <- c("bayesma_robma_stan_data", "list")
  out
}


#' @noRd
overlay_custom_data <- function(stan_data_list, custom_data, method) {
  if (is.null(custom_data)) return(stan_data_list)

  is_per_label <- all(purrr::map_lgl(custom_data, is.list)) &&
    !is.null(names(custom_data)) &&
    all(names(custom_data) %in% names(stan_data_list))

  if (is_per_label) {
    for (nm in names(custom_data)) {
      for (k in names(custom_data[[nm]])) {
        stan_data_list[[nm]][[k]] <- custom_data[[nm]][[k]]
      }
    }
    return(stan_data_list)
  }

  if (method == "ss") {
    key <- names(stan_data_list)[[1]]
    for (k in names(custom_data)) {
      stan_data_list[[key]][[k]] <- custom_data[[k]]
    }
    return(stan_data_list)
  }

  for (nm in names(stan_data_list)) {
    if (is.null(stan_data_list[[nm]])) next
    for (k in names(custom_data)) {
      stan_data_list[[nm]][[k]] <- custom_data[[k]]
    }
  }
  stan_data_list
}


#' @export
print.bayesma_robma_stan_data <- function(x, ...) {
  cat("<bayesma_robma_stan_data> (", length(x), " entr",
      if (length(x) == 1) "y" else "ies", ")\n", sep = "")
  for (nm in names(x)) {
    cat("  ", nm, ": ", sep = "")
    if (is.null(x[[nm]])) {
      cat("(analytic -- no data)\n")
    } else {
      cat(length(x[[nm]]), " element(s) [",
          paste(names(x[[nm]]), collapse = ", "), "]\n", sep = "")
    }
  }
  invisible(x)
}


# -----------------------------------------------------------------------------
# Stage 4: robma_fit
# -----------------------------------------------------------------------------

#' @noRd
robma_fit <- function(spec, code, stan_data,
                      chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                      adapt_delta = 0.95, seed = 1234,
                      parallel = FALSE, quiet = FALSE, ...) {
  if (!inherits(spec, "bayesma_robma_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_robma_spec} object.",
                   call = rlang::caller_env())
  }

  fit <- if (spec$method == "bridge") {
    robma_fit_bridge(spec, code, stan_data,
                     chains, iter_warmup, iter_sampling,
                     adapt_delta, seed, parallel, quiet, ...)
  } else {
    robma_fit_ss(spec, code, stan_data,
                 chains, iter_warmup, iter_sampling,
                 adapt_delta, seed, quiet, ...)
  }
  class(fit) <- c("bayesma_robma_fit", "list")
  fit
}


#' @noRd
robma_fit_bridge <- function(spec, code, stan_data,
                             chains, iter_warmup, iter_sampling,
                             adapt_delta, seed, parallel, quiet, ...) {
  if (!requireNamespace("bridgesampling", quietly = TRUE)) {
    cli::cli_abort(
      "{.pkg bridgesampling} is required for {.code method = 'bridge'}.",
      call = rlang::caller_env()
    )
  }

  inform      <- if (quiet) \(...) invisible(NULL) else cli::cli_alert_info
  inform_step <- if (quiet) \(...) invisible(NULL) else cli::cli_alert

  labels <- names(code)
  n_models <- length(labels)

  model_specs <- purrr::map(labels, function(nm) {
    entry <- code[[nm]]
    list(
      label           = entry$label %||% nm,
      stan_code       = entry$full,
      stan_data       = stan_data[[nm]],
      analytic        = isTRUE(entry$analytic),
      log_ml_override = entry$log_ml_override,
      is_effect_null  = isTRUE(entry$is_effect_null),
      is_hetero_null  = isTRUE(entry$is_hetero_null),
      is_bias_null    = isTRUE(entry$is_bias_null),
      prior_weight    = entry$prior_weight %||% 1,
      hetero_prior_for_fallback = entry$hetero_prior_for_fallback
    )
  }) |>
    stats::setNames(labels)

  inform("Compiling unique Stan models...")
  stan_codes <- purrr::map_chr(model_specs, function(sp) {
    if (isTRUE(sp$analytic) || is.null(sp$stan_code) ||
        is.na(sp$stan_code)) return(NA_character_)
    sp$stan_code
  })
  unique_codes <- unique(stats::na.omit(stan_codes))
  inform("  {length(unique_codes)} unique Stan programs to compile")

  compiled_models <- purrr::map(
    seq_along(unique_codes),
    function(i) {
      inform_step("  Compiling model {i}/{length(unique_codes)}")
      stan_file <- cmdstanr::write_stan_file(unique_codes[[i]])
      cmdstanr::cmdstan_model(stan_file, force_recompile = TRUE, quiet = TRUE)
    }
  )
  names(compiled_models) <- unique_codes

  model_specs <- purrr::map(model_specs, function(sp) {
    sp$compiled_model <- if (isTRUE(sp$analytic) || is.null(sp$stan_code) ||
                             is.na(sp$stan_code)) {
      NULL
    } else {
      compiled_models[[sp$stan_code]]
    }
    sp
  })

  inform("Sampling {length(model_specs)} models...")
  sample_one_model <- function(sp, idx) {
    inform_step("  [{idx}/{length(model_specs)}] {sp$label}")
    if (isTRUE(sp$analytic)) {
      return(list(
        fit = NULL, stan_code = sp$stan_code, stan_data = sp$stan_data,
        label = sp$label, prior_weight = sp$prior_weight,
        is_effect_null = sp$is_effect_null,
        is_hetero_null = sp$is_hetero_null,
        is_bias_null   = sp$is_bias_null,
        log_ml_override = sp$log_ml_override, analytic = TRUE
      ))
    }
    tryCatch({
      fit <- sp$compiled_model$sample(
        data          = sp$stan_data,
        chains        = chains,
        iter_warmup   = iter_warmup,
        iter_sampling = iter_sampling,
        adapt_delta   = adapt_delta,
        seed          = seed,
        refresh       = 0,
        show_messages = FALSE,
        show_exceptions = FALSE,
        ...
      )
      list(
        fit = fit, stan_code = sp$stan_code, stan_data = sp$stan_data,
        label = sp$label, prior_weight = sp$prior_weight,
        is_effect_null = sp$is_effect_null,
        is_hetero_null = sp$is_hetero_null,
        is_bias_null   = sp$is_bias_null,
        log_ml_override = sp$log_ml_override, analytic = FALSE,
        hetero_prior_for_fallback = sp$hetero_prior_for_fallback
      )
    }, error = function(e) {
      cli::cli_warn("Sampling failed for {.val {sp$label}}: {e$message}")
      NULL
    })
  }

  component_fits <- if (parallel) {
    run_sampling_parallel(model_specs, sample_one_model, inform)
  } else {
    purrr::imap(model_specs, sample_one_model)
  }
  component_fits <- purrr::compact(component_fits)
  if (length(component_fits) < 2) {
    cli::cli_abort("RoBMA requires at least 2 successfully fitted models.",
                   call = rlang::caller_env())
  }

  inform("Computing marginal likelihoods...")
  bs_cache <- new.env(parent = emptyenv())
  log_mls <- purrr::map_dbl(component_fits, function(comp) {
    if (!is.null(comp$log_ml_override)) return(comp$log_ml_override)
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

    if (!is.finite(log_ml_val) &&
        isTRUE(comp$is_effect_null) && isTRUE(comp$is_bias_null) &&
        !isTRUE(comp$is_hetero_null) &&
        !is.null(comp$hetero_prior_for_fallback)) {
      log_ml_val <- compute_log_ml_re_null(
        spec$es$yi, spec$es$sei, comp$hetero_prior_for_fallback)
      cli::cli_alert_info(
        "    {comp$label}: quadrature fallback log ML = {round(log_ml_val, 2)}")
    }

    assign(fit_id, log_ml_val, envir = bs_cache)
    log_ml_val
  })

  prior_weights <- purrr::map_dbl(component_fits, ~ .x$prior_weight)
  post_probs    <- compute_posterior_probs(log_mls, prior_weights)
  names(post_probs) <- purrr::map_chr(component_fits, ~ .x$label)
  finite_mask <- is.finite(log_mls + log(prior_weights / sum(prior_weights)))

  is_h1      <- !purrr::map_lgl(component_fits, ~ .x$is_effect_null)
  has_bias   <- !purrr::map_lgl(component_fits, ~ .x$is_bias_null)
  has_hetero <- !purrr::map_lgl(component_fits, ~ .x$is_hetero_null)

  bf_effect     <- compute_inclusion_bf(is_h1, post_probs, prior_weights, finite_mask)
  bf_bias       <- compute_inclusion_bf(has_bias, post_probs, prior_weights, finite_mask)
  bf_hetero     <- compute_inclusion_bf(has_hetero, post_probs, prior_weights, finite_mask)
  bf_per_mech   <- compute_per_mechanism_bfs(component_fits, post_probs,
                                             prior_weights, finite_mask)

  cli::cli_alert_info("Computing model-averaged posterior...")
  n_total_draws <- iter_sampling * chains
  averaged_mu <- purrr::imap(component_fits, function(comp, idx) {
    n_from <- max(1, round(post_probs[idx] * n_total_draws))
    if (isTRUE(comp$is_effect_null)) {
      return(tibble::tibble(mu = rep(0, n_from), model = comp$label))
    }
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
    averaged_draws = averaged_mu,
    model_table    = model_table,
    inclusion_bf   = list(effect = bf_effect, bias = bf_bias,
                          heterogeneity = bf_hetero,
                          by_mechanism = bf_per_mech),
    posterior_probs = list(effect = sum(post_probs[is_h1]),
                           bias = sum(post_probs[has_bias]),
                           heterogeneity = sum(post_probs[has_hetero])),
    component_fits = component_fits,
    post_probs     = post_probs,
    log_marginal_likelihoods = log_mls
  )
}


#' @noRd
run_sampling_parallel <- function(model_specs, sample_one_model, inform) {
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
    return(results[])
  }

  if (.Platform$OS.type != "windows") {
    n_cores <- max(1L, parallel::detectCores() - 1L)
    inform("  Using {n_cores} cores via mclapply...")
    return(parallel::mclapply(
      seq_along(model_specs),
      function(idx) sample_one_model(model_specs[[idx]], idx),
      mc.cores = n_cores,
      mc.set.seed = TRUE
    ))
  }

  cli::cli_warn(
    "Parallel not available (install {.pkg mirai}). Running sequentially."
  )
  purrr::imap(model_specs, sample_one_model)
}


#' @noRd
robma_fit_ss <- function(spec, code, stan_data,
                         chains, iter_warmup, iter_sampling,
                         adapt_delta, seed, quiet, ...) {
  inform <- if (quiet) \(...) invisible(NULL) else cli::cli_alert_info
  warn   <- if (quiet) \(...) invisible(NULL) else cli::cli_alert_warning

  if (isTRUE(spec$horseshoe)) {
    inform("RoBMA (regularised horseshoe effect prior): bias indicator = bias_corrected")
  } else {
    inform("RoBMA (spike-and-slab via log_mix): bias indicator = {spec$bias_indicator}")
  }

  if (spec$bias_indicator == "pet_peese") {
    warn(paste0(
      "PET-PEESE has identifiability issues in the spike-and-slab ",
      "framework: the effect (mu) and bias slope (beta * 1/sqrt(n)) ",
      "are confounded. Consider bias_indicator = 'bias_corrected' or ",
      "method = 'bridge'."))
  }
  inform(paste0(
    "Note: The spike-and-slab method is a fast approximation. ",
    "It typically produces more conservative inclusion BFs than ",
    "bridge sampling. For publication-quality results, use method = 'bridge'."
  ))

  key <- names(code)[[1]]
  stan_program <- code[[key]]$full
  data_list    <- stan_data[[key]]

  stan_file <- cmdstanr::write_stan_file(stan_program)
  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data          = data_list,
    chains        = chains,
    iter_warmup   = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta   = adapt_delta,
    seed          = seed,
    refresh       = 0,
    show_messages = FALSE,
    show_exceptions = FALSE,
    ...
  )

  draws <- posterior::as_draws_df(fit$draws())
  mu_draws    <- as.vector(draws$mu)
  averaged_mu <- tibble::tibble(mu = mu_draws, model = "spike_slab_joint")

  use_horseshoe <- isTRUE(spec$horseshoe)
  pip_effect <- if (use_horseshoe) {
    mean(as.vector(draws$pip_effect_approx))
  } else {
    mean(as.vector(draws$pip_effect))
  }
  pip_hetero <- mean(as.vector(draws$pip_hetero))
  pip_bias   <- mean(as.vector(draws$pip_bias))

  pip_to_bf <- function(pip) {
    pip <- max(min(pip, 1 - 1e-10), 1e-10)
    pip / (1 - pip)
  }

  n_div <- tryCatch(
    sum(fit$diagnostic_summary()$num_divergent),
    error = function(e) NA_integer_
  )
  if (!is.na(n_div) && n_div > 0) {
    n_total <- iter_sampling * chains
    pct <- round(100 * n_div / n_total, 1)
    cli::cli_alert_warning(paste0(
      n_div, " of ", n_total, " (", pct, "%) transitions diverged. ",
      "Results may be unreliable. Consider increasing adapt_delta ",
      "or iter_warmup."
    ))
  }

  model_label <- if (use_horseshoe) {
    "Horseshoe (bias_corrected)"
  } else {
    paste0("Spike-and-slab (", spec$bias_indicator, ")")
  }

  model_table <- tibble::tibble(
    model = model_label,
    log_ml = NA_real_, prior_weight = 1, post_prob = 1,
    null_effect = NA, has_bias = TRUE, has_heterogeneity = TRUE,
    pip_effect = pip_effect, pip_hetero = pip_hetero, pip_bias = pip_bias
  )

  list(
    averaged_draws  = averaged_mu,
    model_table     = model_table,
    inclusion_bf    = list(
      effect        = pip_to_bf(pip_effect),
      bias          = pip_to_bf(pip_bias),
      heterogeneity = pip_to_bf(pip_hetero)
    ),
    posterior_probs = list(
      effect        = pip_effect,
      heterogeneity = pip_hetero,
      bias          = pip_bias
    ),
    component_fits  = list(list(
      fit = fit, label = "spike_slab_joint",
      null_effect = FALSE, has_bias = TRUE, has_hetero = TRUE
    )),
    post_probs               = c(spike_slab_joint = 1),
    log_marginal_likelihoods = NA_real_,
    pip = list(effect = pip_effect, heterogeneity = pip_hetero, bias = pip_bias)
  )
}


# -----------------------------------------------------------------------------
# Stage 5: robma_extract
# -----------------------------------------------------------------------------

#' @noRd
robma_extract <- function(fit, spec) {
  if (!inherits(spec, "bayesma_robma_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls bayesma_robma_spec} object.",
                   call = rlang::caller_env())
  }
  if (!inherits(fit, "bayesma_robma_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls bayesma_robma_fit} object.",
                   call = rlang::caller_env())
  }

  mu_avg <- fit$averaged_draws$mu

  null_range_probs <- compute_null_range_probs(
    mu_avg, spec$null_range, spec$effect_label
  )

  study_rows <- tibble::tibble(
    study    = spec$study_labels,
    estimate = spec$es$yi,
    lower    = spec$es$yi - 1.96 * spec$es$sei,
    upper    = spec$es$yi + 1.96 * spec$es$sei,
    type     = "study"
  )
  pooled_row <- tibble::tibble(
    study    = "Pooled (RoBMA)",
    estimate = stats::median(mu_avg),
    lower    = stats::quantile(mu_avg, 0.025),
    upper    = stats::quantile(mu_avg, 0.975),
    type     = "pooled"
  )
  forest_df <- dplyr::bind_rows(study_rows, pooled_row) |>
    dplyr::mutate(study = forcats::fct_inorder(.data$study),
                  effect_scale = spec$effect_label)

  out <- list(
    null_range_probs = null_range_probs,
    forest_df        = forest_df,
    ma_posterior     = as.numeric(mu_avg)
  )
  class(out) <- c("bayesma_robma_effects", "list")
  out
}


# -----------------------------------------------------------------------------
# Stage 6: robma_output
# -----------------------------------------------------------------------------

#' @noRd
robma_output <- function(spec, fit, effects) {
  result <- fit
  result$meta <- list(
    likelihood     = spec$likelihood,
    model_type     = "robma",
    stage          = "two_stage",
    method         = spec$method,
    bias_indicator = if (spec$method == "ss") spec$bias_indicator else NA_character_,
    horseshoe      = isTRUE(spec$horseshoe),
    study_labels   = spec$study_labels,
    effect_label   = spec$effect_label,
    es             = spec$es,
    n_models       = length(fit$component_fits),
    robust         = FALSE,
    null_range     = spec$null_range,
    null_range_probs = effects$null_range_probs,
    priors = list(
      effect             = spec$priors_effect,
      effect_null        = spec$priors_effect_null,
      heterogeneity      = spec$priors_heterogeneity,
      heterogeneity_null = spec$priors_heterogeneity_null,
      bias               = spec$priors_bias,
      bias_null          = spec$priors_bias_null
    ),
    call_args = spec$call_args
  )
  result$forest_df    <- effects$forest_df
  result$ma_posterior <- effects$ma_posterior

  class(result) <- c("bayesma_robma", "bayesma")
  result
}
