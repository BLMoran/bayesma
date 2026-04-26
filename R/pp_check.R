#' Posterior Predictive Check for bayesma Objects
#'
#' @description
#' Produces a posterior predictive check plot comparing observed data to
#' replicated datasets drawn from the posterior predictive distribution.
#' Provides a posterior predictive check interface analogous to standard Bayesian workflow.
#'
#' @param object A \code{bayesma} object.
#' @param type Character. Plot type passed to \code{bayesplot::ppc_*}.
#'   One of \code{"dens_overlay"} (default), \code{"hist"}, \code{"stat"},
#'   \code{"scatter"}, \code{"bars"} (discrete only), \code{"ecdf_overlay"},
#'   \code{"ribbon"}.
#' @param ndraws Number of posterior draws to use. Default 100.
#' @param ... Additional arguments passed to the underlying \code{bayesplot::ppc_*} function.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' fit <- bayesma(data, likelihood = "binomial", ...)
#' pp_check(fit)
#' pp_check(fit, type = "stat", stat = "mean")
#' pp_check(fit, type = "scatter")
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

  p <- ppc_fn(y, y_rep, ...) +
    ggplot2::ggtitle("Posterior Predictive Check") +
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
