# ============================================================================
# NOTE: The main sensitivity_plot() function is defined in sensitivity_plot.R.
# This file contains helper functions used by sensitivity_plot():
#   - bayesma_ensure_daemons()
#   - bayesma_future_map_dfr()
#   - bayesma_future_pmap_dfr()
#   - robma_to_sensitivity_draws()
#   - summarise_sensitivity_posteriors()
# ============================================================================


# ============================================================================
# Internal: parallel mapping helpers (mirai > mclapply > sequential)
# ============================================================================

#' Set up mirai daemons if not already configured
#' @noRd
bayesma_ensure_daemons <- function(workers) {
  daemons_were_set <- tryCatch({
    s <- mirai::status()
    is.data.frame(s$daemons) && nrow(s$daemons) > 0L
  }, error = function(e) FALSE)

  if (!daemons_were_set) {
    n_cores <- workers %||% max(1L, parallel::detectCores() - 1L)
    mirai::daemons(n_cores)
    tryCatch(mirai::everywhere(library(bayesma)), error = function(e) NULL)
    return(TRUE)
  }
  FALSE
}

#' @noRd
bayesma_future_map_dfr <- function(.x, .f, parallel = FALSE, workers = NULL, seed = TRUE, ...) {
  if (!isTRUE(parallel)) {
    return(purrr::map(.x, .f, ...) |> purrr::list_rbind())
  }

  if (requireNamespace("mirai", quietly = TRUE)) {
    teardown <- bayesma_ensure_daemons(workers)
    if (teardown) on.exit(mirai::daemons(0), add = TRUE)

    results <- mirai::mirai_map(
      .x,
      function(.x_i, .f, ...) .f(.x_i, ...),
      .args = list(.f = .f, ...)
    )
    return(dplyr::bind_rows(results[]))
  }

  # Fallback to mclapply on macOS/Linux
  if (.Platform$OS.type != "windows") {
    n_cores <- workers %||% max(1L, parallel::detectCores() - 1L)
    results <- parallel::mclapply(.x, .f, ..., mc.cores = n_cores, mc.set.seed = TRUE)
    return(dplyr::bind_rows(results))
  }

  cli::cli_warn("Parallel not available (install {.pkg mirai}). Running sequentially.")
  purrr::map(.x, .f, ...) |> purrr::list_rbind()
}

#' @noRd
bayesma_future_pmap_dfr <- function(.l, .f, parallel = FALSE, workers = NULL, seed = TRUE, ...) {
  if (!isTRUE(parallel)) {
    return(purrr::pmap(.l, .f, ...) |> purrr::list_rbind())
  }

  if (requireNamespace("mirai", quietly = TRUE)) {
    teardown <- bayesma_ensure_daemons(workers)
    if (teardown) on.exit(mirai::daemons(0), add = TRUE)

    n_tasks <- length(.l[[1]])
    task_args <- purrr::map(seq_len(n_tasks), function(i) {
      purrr::map(.l, ~ .x[[i]])
    })

    results <- mirai::mirai_map(
      task_args,
      function(args, .f) do.call(.f, args),
      .args = list(.f = .f)
    )
    return(dplyr::bind_rows(results[]))
  }

  # Fallback to mclapply on macOS/Linux
  if (.Platform$OS.type != "windows") {
    n_cores <- workers %||% max(1L, parallel::detectCores() - 1L)
    n_tasks <- length(.l[[1]])
    results <- parallel::mclapply(seq_len(n_tasks), function(i) {
      args <- purrr::map(.l, ~ .x[[i]])
      do.call(.f, args)
    }, mc.cores = n_cores, mc.set.seed = TRUE)
    return(dplyr::bind_rows(results))
  }

  cli::cli_warn("Parallel not available (install {.pkg mirai}). Running sequentially.")
  purrr::pmap(.l, .f, ...) |> purrr::list_rbind()
}


# ============================================================================
# Internal: RoBMA draw extraction
# ============================================================================
#' @noRd
robma_to_sensitivity_draws <- function(robma_fit,
                                       measure,
                                       prior,
                                       prior_label,
                                       section_label = "RoBMA") {

  x <- NULL

  # Preferred alias if present
  if (!is.null(robma_fit$ma_posterior)) {
    x <- robma_fit$ma_posterior
  }

  # Your bayesma_robma output: averaged_draws$mu
  if (is.null(x) && !is.null(robma_fit$averaged_draws)) {
    if (is.data.frame(robma_fit$averaged_draws) && "mu" %in% names(robma_fit$averaged_draws)) {
      x <- robma_fit$averaged_draws$mu
    } else if (is.list(robma_fit$averaged_draws) && !is.null(robma_fit$averaged_draws$mu)) {
      x <- robma_fit$averaged_draws$mu
    }
  }

  # Fallback: bayesma-style draws
  if (is.null(x) && !is.null(robma_fit$draws$mu)) {
    x <- robma_fit$draws$mu
  }

  if (is.null(x)) {
    cli::cli_abort(c(
      "Could not find model-averaged posterior draws in {.arg robma_fit}.",
      "i" = "Expected $ma_posterior, $averaged_draws$mu, or $draws$mu."
    ))
  }

  x <- as.numeric(x)

  # RoBMA returns mu on log scale for ratio measures in your implementation
  if (measure %in% c("OR", "RR", "HR", "IRR")) {
    x <- exp(x)
  }

  tibble::tibble(
    x = x,
    prior = prior,
    prior_label = prior_label,
    section_label = section_label
  )
}


# ============================================================================
# Internal: Posterior summaries for the sensitivity table
# ============================================================================
#' @noRd
summarise_sensitivity_posteriors <- function(draws, null_value, null_range) {

  section_levels <- unique(draws$section_label)
  prior_levels   <- unique(draws$prior)

  out <- draws |>
    dplyr::mutate(
      section_label = factor(section_label, levels = section_levels),
      prior         = factor(prior, levels = prior_levels)
    ) |>
    dplyr::summarise(
      median = stats::median(x),
      l95    = stats::quantile(x, 0.025),
      u95    = stats::quantile(x, 0.975),

      # Use names expected by sensitivity_table_right()
      pr_benefit    = round(mean(x < null_value) * 100, 1),
      pr_no_benefit = round(mean(x > null_value) * 100, 1),

      .by = c(section_label, prior, prior_label)
    ) |>
    dplyr::arrange(section_label, prior)

  if (!is.null(null_range) && length(null_range) == 2) {
    out <- out |>
      dplyr::mutate(
        pr_benefit_null_range    = round(mean(x < null_range[1]) * 100, 1),
        pr_no_benefit_null_range = round(mean(x > null_range[2]) * 100, 1)
      )
  } else {
    out <- out |>
      dplyr::mutate(
        pr_benefit_null_range    = NA_real_,
        pr_no_benefit_null_range = NA_real_
      )
  }

  out
}
