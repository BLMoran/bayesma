#' Posterior and Prior Predictive Checks for bayesma Objects
#'
#' @description
#' Compares observed data to replicated datasets drawn from either the
#' posterior predictive distribution (standard use) or the prior predictive
#' distribution (when the model was fitted with `sample_prior = TRUE`).
#'
#' **Posterior predictive check** (`sample_prior = FALSE`, the default): draws
#' `y_rep` from `p(y_rep | y)`. Agreement between `y` and `y_rep` indicates the
#' fitted model can reproduce the observed data.
#'
#' **Prior predictive check** (`sample_prior = TRUE`): draws `y_rep` from
#' `p(y_rep)`, integrating over the prior without conditioning on the data.
#' This is a tool for prior *elicitation and sanity-checking* — verifying that
#' the prior does not place substantial mass on data values that are impossible
#' or implausible given domain knowledge, before the data have been seen.
#'
#' @section Prior predictive checks and double-dipping:
#' It is legitimate to use a prior predictive check to verify that a prior is
#' coherent on the observable scale (e.g. that it does not imply impossible
#' event counts). It is **not** legitimate to iteratively adjust priors until
#' `y_rep` matches the observed `y`: doing so smuggles the data into the prior,
#' inflating posterior confidence. Priors should be specified from external
#' knowledge, expert elicitation, or independent reference data — not from the
#' analysis dataset itself.
#'
#' @param object A `bayesma` object. If fitted with `sample_prior = TRUE`,
#'   a prior predictive check is produced; otherwise a posterior predictive
#'   check.
#' @param type Character. Plot type passed to `bayesplot::ppc_*`. One of
#'   `"dens_overlay"` (default), `"hist"`, `"stat"`, `"scatter"`,
#'   `"bars"` (discrete only), `"ecdf_overlay"`, `"ribbon"`.
#' @param ndraws Integer. Number of draws to display. Default `100`.
#' @param ... Additional arguments passed to the underlying
#'   `bayesplot::ppc_*` function.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' # Posterior predictive check
#' fit <- bayesma(data, likelihood = "binomial", ...)
#' pp_check(fit)
#' pp_check(fit, type = "stat", stat = "mean")
#'
#' # Prior predictive check — use to verify priors are sensible,
#' # not to calibrate them against the analysis data
#' prior_fit <- bayesma(data, likelihood = "binomial", ..., sample_prior = TRUE)
#' pp_check(prior_fit)
#' }
#'
#' @export
pp_check <- function(object, type = "dens_overlay", ndraws = 100, ...) {

  # ---- Determine observed data and likelihood ----
  likelihood <- object$meta$likelihood %||% "gaussian"
  stage      <- object$meta$stage %||% "two_stage"

  y <- extract_observed(object, stage, likelihood)
  if (is.null(y)) {
    stop(
      "Could not extract observed data for posterior predictive check.\n",
      "Ensure object$stan_data contains the observed outcome.",
      call. = FALSE
    )
  }

  # ---- Verify Stan output files still exist ----
  check_fit_accessible(object$fit)

  # ---- Extract y_rep ----
  all_vars   <- posterior::variables(object$fit$draws())
  y_rep_vars <- grep("^y_rep\\[", all_vars, value = TRUE)

  if (length(y_rep_vars) == 0) {
    stop(
      "No 'y_rep' found in posterior draws.\n",
      "Add a generated quantities block with y_rep to your Stan model.\n",
      "See ?bayesma for models that support pp_check().",
      call. = FALSE
    )
  }

  y_rep <- posterior::as_draws_matrix(
    object$fit$draws(variables = y_rep_vars)
  )
  y_rep <- as.matrix(y_rep)  # strip draws_matrix class for bayesplot

  # ---- Subsample for readability ----
  if (nrow(y_rep) > ndraws) {
    idx   <- sample.int(nrow(y_rep), ndraws)
    y_rep <- y_rep[idx, , drop = FALSE]
  }

  # ---- Check dimension agreement ----
  if (ncol(y_rep) != length(y)) {
    stop(
      sprintf(
        "Dimension mismatch: y has %d observations but y_rep has %d columns.\n",
        length(y), ncol(y_rep)
      ),
      "Check that y_rep in generated quantities matches the observed data.",
      call. = FALSE
    )
  }

  # ---- Determine if data is actually discrete ----
  is_discrete <- is.integer(y) && all(y == round(y), na.rm = TRUE)

  # ---- Validate type ----
  continuous_types <- c("dens_overlay", "hist", "ecdf_overlay",
                        "stat", "scatter", "ribbon")
  discrete_types   <- c("bars", "hist", "stat", "rootogram")

  if (is_discrete && type == "dens_overlay") {
    message("Note: y_rep is on a continuous scale. Using 'dens_overlay'.")
  }

  if (!is_discrete && type %in% c("bars", "rootogram")) {
    message(sprintf("Type '%s' requires discrete data but y is continuous. ", type),
            "Switching to 'dens_overlay'.")
    type <- "dens_overlay"
  }

  # ---- Dispatch to bayesplot ----
  ppc_fn <- switch(type,
                   dens_overlay = bayesplot::ppc_dens_overlay,
                   hist         = bayesplot::ppc_hist,
                   ecdf_overlay = bayesplot::ppc_ecdf_overlay,
                   stat         = bayesplot::ppc_stat,
                   scatter      = bayesplot::ppc_scatter_avg,
                   bars         = bayesplot::ppc_bars,
                   ribbon       = bayesplot::ppc_ribbon,
                   rootogram    = bayesplot::ppc_rootogram,
                   {
                     stop(
                       sprintf("Unknown pp_check type: '%s'.\n", type),
                       sprintf("Available types: %s",
                               paste(c(continuous_types, discrete_types), collapse = ", ")),
                       call. = FALSE
                     )
                   }
  )

  check_label <- if (isTRUE(object$meta$prior_predictive)) {
    "Prior Predictive Check"
  } else {
    "Posterior Predictive Check"
  }

  p <- ppc_fn(y, y_rep, ...) +
    ggplot2::ggtitle(check_label) +
    ggplot2::theme_minimal(base_family = "") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )

  p
}


# Internal function to extract observed data vector

#' @noRd
extract_observed <- function(object, stage, likelihood) {

  # Two-stage: y is always in stan_data
  if (stage == "two_stage") {
    return(object$stan_data$y)
  }

  # One-stage gaussian: y in stan_data
  if (likelihood == "gaussian") {
    return(object$stan_data$y)
  }

  # One-stage binomial/poisson: events in stan_data
  if (likelihood %in% c("binomial", "poisson")) {
    return(object$stan_data$events)
  }

  NULL
}
