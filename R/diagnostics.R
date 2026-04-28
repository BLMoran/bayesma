#' Single-Page Model Diagnostics for bayesma Objects
#'
#' @description
#' Produces a single-page diagnostic summary for a Bayesian meta-analysis model
#' fitted via \code{cmdstanr}. Combines visual diagnostics (trace plots,
#' posterior predictive check, Rhat, ESS) with a tabular summary of MCMC health
#' indicators. All output is arranged on a single page using \code{patchwork}.
#'
#' @param object A \code{bayesma} object containing a \code{CmdStanMCMC} fit
#'   in \code{object$fit}.
#' @param pars Character vector of parameter names for trace/ACF plots.
#'   If \code{NULL}, sensible defaults are chosen based on the model stage.
#' @param ndraws Number of posterior draws for the posterior predictive check.
#'   Default is 100.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list of diagnostic values. Prints a composite
#'   \code{patchwork} plot as a side effect.
#'
#' @details
#' The diagnostic page contains six panels:
#' \describe{
#'   \item{Top-left}{Trace plots for key parameters.}
#'   \item{Top-right}{Posterior predictive check.}
#'   \item{Middle-left}{Rhat values for all parameters.}
#'   \item{Middle-right}{Effective sample size ratios for all parameters.}
#'   \item{Bottom-left}{Autocorrelation function for key parameters.}
#'   \item{Bottom-right}{MCMC diagnostics summary table.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- bayesma(data, likelihood = "binomial", ...)
#' diagnostics(fit)
#' }
#'
#' @export
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}


#' @rdname diagnostics
#' @export
diagnostics.bayesma <- function(object, pars = NULL, ndraws = 100, ...) {

  # ---- validate object ----
  if (is.null(object$fit)) {
    stop("object$fit is NULL. Expected a CmdStanMCMC object.", call. = FALSE)
  }

  # ---- resolve stage and measure ----
  stage   <- object$meta$stage   %||% "unknown"
  measure <- object$meta$effect_label %||% "unknown"

  # ---- convert draws to array for bayesplot ----
  draws_array <- posterior::as_draws_array(object$fit$draws())

  # ---- resolve default parameters ----
  if (is.null(pars)) {
    pars <- default_pars(object, draws_array)
  }

  # ---- extract diagnostics ----
  diagnostics <- extract_diagnostics(object)

  # ---- build panels ----
  p_trace   <- panel_trace(draws_array, pars)
  p_ppc     <- panel_ppc(object, ndraws)
  p_rhat    <- panel_rhat(object$summary)
  p_neff    <- panel_neff(object$summary)
  p_acf     <- panel_acf(draws_array, pars)
  p_table   <- panel_diagnostics_table(diagnostics)

  # ---- compose single page ----
  composed <- patchwork::wrap_plots(
    p_trace,   p_ppc,
    p_rhat,    p_neff,
    p_acf,     p_table,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title    = "bayesma: Model Diagnostics",
      subtitle = paste0("Stage: ", stage, " | Measure: ", measure),
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(family = "", face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(family = "", size = 10, colour = "grey40")
      )
    )

  print(composed)

  invisible(object)
}
