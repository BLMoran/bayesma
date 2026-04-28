#' Run RoBMA Models Across Multiple Prior Specifications
#'
#' @description
#' Pre-computes RoBMA models for each prior specification provided. These fits
#' are stored and can later be attached to a bayesma object for use in
#' `sensitivity_plot()`. This separation allows the computationally expensive
#' RoBMA fitting to be done once and reused.
#'
#' @param data A data frame containing the study data.
#' @param priors A named list of prior specifications. Each element must be a list
#'   with at least `mu_prior` and (optionally) `tau_prior`, and may include
#'   `name` used for display.
#'   Example:
#'   \code{list(
#'     prior_1 = list(name = "Vague", mu_prior = normal(0, 10), tau_prior = half_cauchy(0, 1)),
#'     prior_2 = list(name = "Regularising", mu_prior = normal(0, 1), tau_prior = half_cauchy(0, 0.5))
#'   )}
#' @param robma_template Either:
#'   \itemize{
#'     \item A fitted RoBMA object (class `bayesma_robma`) with stored `$meta$call_args`
#'       to use as a template for refitting across priors, OR
#'     \item `NULL` (default), in which case `robma_args` must be provided
#'   }
#' @param robma_args A list of arguments to pass to `robma()` when `robma_template`
#'   is `NULL`. Must include all required arguments except `priors_effect` and
#'   `priors_heterogeneity`, which will be set from the `priors` argument.
#' @param parallel Logical. If TRUE, uses parallel processing for refits.
#' @param workers Optional integer. Number of parallel workers.
#' @param seed Logical. If TRUE (default), uses parallel-safe seeding.
#' @param .progress Logical. If TRUE, displays a progress bar.
#'
#' @return A `bayesma_robma_sensitivity` object: a named list of RoBMA fits keyed
#'   by prior IDs, with additional metadata in `$meta`.
#'
#' @details
#' The returned object can be attached to a bayesma model using
#' `attach_robma_sensitivity()` or passed directly to `sensitivity_plot()`.
#'
#' @examples
#' \dontrun{
#' # Define priors for sensitivity analysis
#' priors <- list(
#'   vague = list(
#'     name = "Vague",
#'     mu_prior = normal(0, 10),
#'     tau_prior = half_cauchy(0, 1)
#'   ),
#'   informative = list(
#'     name = "Informative",
#'     mu_prior = normal(0, 0.5),
#'     tau_prior = half_cauchy(0, 0.3)
#'   )
#' )
#'
#' # Option 1: Using a template RoBMA fit
#' robma_sens <- robma_sensitivity(
#'   data = my_data,
#'   priors = priors,
#'   robma_template = my_robma_fit,
#'   parallel = TRUE
#' )
#'
#' # Option 2: Specifying robma arguments directly
#' robma_sens <- robma_sensitivity(
#'   data = my_data,
#'   priors = priors,
#'   robma_args = list(
#'     studyvar = "study_id",
#'     event_ctrl = "events_c",
#'     event_int = "events_i",
#'     n_ctrl = "n_c",
#'     n_int = "n_i",
#'     likelihood = "binomial"
#'   ),
#'   parallel = TRUE
#' )
#'
#' # Attach to bayesma model
#' model <- attach_robma_sensitivity(model, robma_sens)
#'
#' # Now sensitivity_plot can use the pre-computed RoBMA fits
#' sensitivity_plot(model, data, priors, estimand = "OR", incl_robma = TRUE)
#' }
#'
#' @export
robma_sensitivity <- function(
    data,
    priors,
    robma_template = NULL,
    robma_args = NULL,
    parallel = FALSE,
    workers = NULL,
    seed = TRUE,
    .progress = TRUE
) {


  # 1) Validate inputs

  if (!is.list(priors) || length(priors) == 0) {
    cli::cli_abort("{.arg priors} must be a non-empty list.")
  }

  prior_ids <- names(priors)
  if (is.null(prior_ids) || any(prior_ids == "")) {
    cli::cli_abort("{.arg priors} must be a named list (each prior must have an ID).")
  }

  # Validate each prior has mu_prior
  purrr::walk(prior_ids, function(pid) {
    if (is.null(priors[[pid]]$mu_prior)) {
      cli::cli_abort(c(
        "Each prior must include {.val mu_prior}.",
        "i" = "Missing mu_prior for prior id: {.val {pid}}"
      ))
    }
  })

  # Validate we have either template or args

  if (is.null(robma_template) && is.null(robma_args)) {
    cli::cli_abort(c(
      "Either {.arg robma_template} or {.arg robma_args} must be provided.",
      "i" = "Use {.arg robma_template} with an existing RoBMA fit, or",
      "i" = "Use {.arg robma_args} with a list of arguments for {.fn robma}."
    ))
  }

  if (!is.null(robma_template)) {
    if (!inherits(robma_template, "bayesma_robma")) {
      cli::cli_abort("{.arg robma_template} must be a {.cls bayesma_robma} object.")
    }
    if (is.null(robma_template$meta$call_args)) {
      cli::cli_abort(c(
        "RoBMA template does not contain stored call arguments.",
        "i" = "Ensure the RoBMA fit was created with a recent version that stores {.code meta$call_args}."
      ))
    }
  }

  # 2) Build prior label map & check prior equality

  prior_name_map <- purrr::map_chr(priors, ~ .x$name %||% "")
  prior_label_map <- purrr::map_chr(prior_ids, function(pid) {
    nm <- prior_name_map[[pid]]
    if (!is.null(nm) && nm != "") nm else pid
  })
  names(prior_label_map) <- prior_ids

  priors_match <- function(p1, p2) {
    if (is.null(p1) && is.null(p2)) return(TRUE)
    if (is.null(p1) || is.null(p2)) return(FALSE)
    if (!inherits(p1, "bayesma_prior") || !inherits(p2, "bayesma_prior")) return(FALSE)
    identical(unclass(p1), unclass(p2))
  }


  # 3) Extract original priors for comparison

  orig_mu <- NULL
  orig_tau <- NULL

  if (!is.null(robma_template)) {
    # RoBMA uses lists of priors; compare against the first element
    if (length(robma_template$meta$call_args$priors_effect) > 0) {
      orig_mu <- robma_template$meta$call_args$priors_effect[[1]]
    }
    if (length(robma_template$meta$call_args$priors_heterogeneity) > 0) {
      orig_tau <- robma_template$meta$call_args$priors_heterogeneity[[1]]
    }
  }


  # 4) Define task runner

  run_single_robma <- function(
    prior_id,
    .priors,
    .data,
    .robma_template,
    .robma_args,
    .orig_mu,
    .orig_tau,
    .priors_match_fn
  ) {
    ps <- .priors[[prior_id]]
    mu_prior <- ps$mu_prior
    tau_prior <- ps$tau_prior

    # Check if we can reuse the template (priors match)
    if (!is.null(.robma_template)) {
      mu_matches <- !is.null(mu_prior) && .priors_match_fn(mu_prior, .orig_mu)
      tau_matches <- is.null(tau_prior) || .priors_match_fn(tau_prior, .orig_tau)

      if (mu_matches && tau_matches) {
        cli::cli_alert_success("Prior {.val {prior_id}}: reusing template (priors match)")
        return(.robma_template)
      }
    }

    # Need to fit a new model
    cli::cli_alert_info("Prior {.val {prior_id}}: fitting RoBMA model...")

    if (!is.null(.robma_template)) {
      # Use template's call_args
      ca <- .robma_template$meta$call_args
      ca$data <- .data

      # Replace effect priors: apply mu_prior to all H1 effect components
      if (!is.null(mu_prior)) {
        ca$priors_effect <- purrr::map(ca$priors_effect, ~ mu_prior)
      }
      if (!is.null(tau_prior)) {
        ca$priors_heterogeneity <- purrr::map(ca$priors_heterogeneity, ~ tau_prior)
      }
    } else {
      # Use robma_args - need to be careful here

      ca <- .robma_args

      # Don't double-set data if already in robma_args
      if (!"data" %in% names(ca) || is.null(ca$data)) {
        ca$data <- .data
      }

      # Set effect priors - robma() expects a list of priors for model averaging
      # We set all H1 priors to the user's mu_prior
      ca$priors_effect <- list(mu_prior)
      if (!is.null(tau_prior)) {
        ca$priors_heterogeneity <- list(tau_prior)
      }

      # Debug: show what we're calling
      cli::cli_alert_info("  Effect prior: {format(mu_prior)}")
      if (!is.null(tau_prior)) {
        cli::cli_alert_info("  Heterogeneity prior: {format(tau_prior)}")
      }
    }

    fit <- tryCatch(
      {
        # Use robma directly (we're inside the bayesma package)
        do.call(robma, ca)
      },
      error = function(e) {
        cli::cli_alert_danger("Failed to fit RoBMA with prior {.val {prior_id}}")
        cli::cli_alert_warning("Error: {e$message}")
        NULL
      }
    )

    fit
  }


  # 5) Run all RoBMA fits

  cli::cli_h2("Running RoBMA sensitivity analysis")
  cli::cli_alert_info("Fitting {length(prior_ids)} RoBMA model{?s} across prior specifications")

  # Always run sequentially first to ensure errors are visible

  # Parallel execution can swallow errors
  if (isTRUE(parallel)) {
    cli::cli_alert_info("Running in parallel mode with {workers %||% 'auto'} workers")

    robma_fits <- bayesma_future_map(
      .x = prior_ids,
      .f = run_single_robma,
      parallel = TRUE,
      workers = workers,
      seed = seed,
      .priors = priors,
      .data = data,
      .robma_template = robma_template,
      .robma_args = robma_args,
      .orig_mu = orig_mu,
      .orig_tau = orig_tau,
      .priors_match_fn = priors_match
    )
    names(robma_fits) <- prior_ids
  } else {
    # Sequential execution - errors will be visible
    robma_fits <- purrr::map(
      prior_ids,
      ~ run_single_robma(
        prior_id = .x,
        .priors = priors,
        .data = data,
        .robma_template = robma_template,
        .robma_args = robma_args,
        .orig_mu = orig_mu,
        .orig_tau = orig_tau,
        .priors_match_fn = priors_match
      ),
      .progress = .progress
    )
    names(robma_fits) <- prior_ids
  }

  # Check for failures
  failed <- purrr::map_lgl(robma_fits, is.null)
  if (any(failed)) {
    cli::cli_warn(c(
      "Some RoBMA fits failed:",
      "x" = paste(prior_ids[failed], collapse = ", ")
    ))
  }

  n_success <- sum(!failed)
  cli::cli_alert_success("Completed: {n_success}/{length(prior_ids)} RoBMA fits")


  # 6) Build output object

  result <- structure(
    robma_fits,
    class = c("bayesma_robma_sensitivity", "list"),
    meta = list(
      prior_ids = prior_ids,
      prior_label_map = prior_label_map,
      priors = priors,
      n_fits = n_success,
      n_failed = sum(failed),
      created_at = Sys.time()
    )
  )

  result
}


#' Attach RoBMA Sensitivity Fits to a bayesma Object
#'
#' @description
#' Attaches pre-computed RoBMA sensitivity fits to a bayesma object.
#' These fits will be used by `sensitivity_plot()` when `incl_robma = TRUE`.
#'
#' @param model A `bayesma` object.
#' @param robma_sensitivity A `bayesma_robma_sensitivity` object created by
#'   `robma_sensitivity()`.
#'
#' @return The modified `bayesma` object with `$robma_sensitivity` attached.
#'
#' @export
attach_robma_sensitivity <- function(model, robma_sensitivity) {
  if (!inherits(model, "bayesma")) {
    cli::cli_abort("{.arg model} must be a {.cls bayesma} object.")
  }

  if (!inherits(robma_sensitivity, "bayesma_robma_sensitivity")) {
    cli::cli_abort(
      "{.arg robma_sensitivity} must be a {.cls bayesma_robma_sensitivity} object ",
      "created by {.fn robma_sensitivity}."
    )
  }

  model$robma_sensitivity <- robma_sensitivity
  cli::cli_alert_success(
    "Attached RoBMA sensitivity fits ({attr(robma_sensitivity, 'meta')$n_fits} models)"
  )

  model
}


#' Check if RoBMA Sensitivity Fits are Available
#'
#' @description
#' Checks whether a bayesma object has pre-computed RoBMA sensitivity fits attached.
#'
#' @param model A `bayesma` object.
#'
#' @return Logical. TRUE if robma_sensitivity is attached, FALSE otherwise.
#'
#' @export
has_robma_sensitivity <- function(model) {
  !is.null(model$robma_sensitivity) &&
    inherits(model$robma_sensitivity, "bayesma_robma_sensitivity")
}


#' Print Method for bayesma_robma_sensitivity
#'
#' @param x A `bayesma_robma_sensitivity` object.
#' @param ... Additional arguments (unused).
#'
#' @export
#' @noRd
print.bayesma_robma_sensitivity <- function(x, ...) {
  meta <- attr(x, "meta")

  cli::cli_h1("RoBMA Sensitivity Analysis Fits")
  cli::cli_alert_info("{meta$n_fits} RoBMA models fitted across prior specifications")

  if (meta$n_failed > 0) {
    cli::cli_alert_warning("{meta$n_failed} fit(s) failed")
  }

  cli::cli_text("")
  cli::cli_h2("Prior Specifications")

  purrr::iwalk(meta$priors, function(p, id) {
    label <- meta$prior_label_map[[id]] %||% id
    mu_str <- if (!is.null(p$mu_prior)) format(p$mu_prior) else "default"
    tau_str <- if (!is.null(p$tau_prior)) format(p$tau_prior) else "default"

    status <- if (is.null(x[[id]])) {
      cli::col_red("\u2717")
    } else {
      cli::col_green("\u2713")
    }

    cli::cli_text("  {status} {.val {label}} ({id}): \u03BC ~ {mu_str}, \u03C4 ~ {tau_str}")
  })

  cli::cli_text("")
  cli::cli_text("Created: {meta$created_at}")

  invisible(x)
}


# Internal: Parallel mapping helper

#' @noRd
bayesma_future_map <- function(.x, .f, parallel = FALSE, workers = NULL, seed = TRUE, ...) {
  if (!isTRUE(parallel)) {
    return(purrr::map(.x, .f, ...))
  }

  dots <- list(...)

  # mirai (preferred)
  if (requireNamespace("mirai", quietly = TRUE)) {
    teardown <- bayesma_ensure_daemons(workers)
    if (teardown) on.exit(mirai::daemons(0), add = TRUE)

    results <- mirai::mirai_map(
      .x,
      function(.x_i, .f, .dots) do.call(.f, c(list(.x_i), .dots)),
      .args = list(.f = .f, .dots = dots)
    )
    return(as.list(results[]))
  }

  # mclapply fallback (macOS/Linux only)
  if (.Platform$OS.type != "windows") {
    n_cores <- workers %||% max(1L, parallel::detectCores() - 1L)
    results <- parallel::mclapply(.x, function(.x_i) {
      do.call(.f, c(list(.x_i), dots))
    }, mc.cores = n_cores, mc.set.seed = TRUE)
    return(results)
  }

  cli::cli_warn("Parallel not available (install {.pkg mirai}). Running sequentially.")
  purrr::map(.x, .f, ...)
}
