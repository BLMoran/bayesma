# ============================================================================
# interpret(): Comprehensive narrative interpretation of bayesma analyses
# ============================================================================

#' Interpret a Bayesian meta-analysis workflow
#'
#' Generates a comprehensive narrative interpretation across one or more
#' bayesma fits -- overall effects, heterogeneity, publication bias, model
#' averaging, sensitivity, model comparison, and convergence diagnostics.
#' The function auto-detects what is provided and assembles only the
#' relevant sections.
#'
#' @param ... One or more fitted bayesma objects (named or unnamed). Accepted
#'   classes: `bayesma`, `bayesma_mv`, `bayesma_robma`, `bayesma_egger`,
#'   `bayesma_metareg`, `bayesma_robma_sensitivity`, `bayesma_comparison`.
#' @param null_range Optional length-2 numeric vector for direction/ROPE
#'   probabilities. Defaults to `c(-0.1, 0.1)` on the natural scale.
#' @param effect_label Optional character override for the effect label used
#'   in narrative text (e.g. "log_or"). Inherited from fits when `NULL`.
#' @param credible_level Credible interval width used in summaries.
#'   Default `0.95`.
#' @param quiet If `TRUE`, suppresses progress messages during assembly.
#'
#' @return An object of class `bayesma_interpretation` -- a list with one
#'   element per detected section plus a `meta` slot. The `print()` method
#'   renders the full narrative report.
#'
#' @examples
#' \dontrun{
#' fit  <- bayesma::bayesma(data = dat, yi = "yi", sei = "sei")
#' rob  <- bayesma::robma(data = dat, yi = "yi", sei = "sei")
#' egg  <- bayesma::egger(data = dat, yi = "yi", sei = "sei")
#' sens <- bayesma::robma_sensitivity(data = dat, priors = my_priors)
#' interpret(fit, rob, egg, sens)
#' }
#' @export
interpret <- function(...,
                      null_range = c(-0.1, 0.1),
                      effect_label = NULL,
                      credible_level = 0.95,
                      quiet = FALSE) {

  fits <- rlang::list2(...)
  if (length(fits) == 0) {
    cli::cli_abort("Provide at least one fitted bayesma object via {.arg ...}.")
  }

  inform <- if (quiet) \(...) invisible() else cli::cli_inform

  classified <- classify_fits(fits)

  inform("Step 1: Interpreting overall effects")
  effects_section <- interpret_effects_section(classified, credible_level)

  inform("Step 2: Interpreting heterogeneity")
  hetero_section <- interpret_heterogeneity_section(classified, credible_level)

  inform("Step 3: Interpreting publication bias")
  bias_section <- interpret_bias_section(classified, credible_level)

  inform("Step 4: Interpreting model averaging (RoBMA)")
  robma_section <- interpret_robma_section(classified)

  inform("Step 5: Interpreting model comparison")
  comparison_section <- interpret_comparison_section(classified)

  inform("Step 6: Interpreting sensitivity analysis")
  sensitivity_section <- interpret_sensitivity_section(classified, null_range)

  inform("Step 7: Interpreting convergence diagnostics")
  diagnostics_section <- interpret_diagnostics_section(classified)

  inform("Step 8: Synthesising overall conclusions")
  conclusions <- interpret_overall_conclusion(
    effects_section, hetero_section, bias_section,
    robma_section, sensitivity_section, diagnostics_section
  )

  out <- list(
    effects       = effects_section,
    heterogeneity = hetero_section,
    bias          = bias_section,
    robma         = robma_section,
    comparison    = comparison_section,
    sensitivity   = sensitivity_section,
    diagnostics   = diagnostics_section,
    conclusions   = conclusions,
    meta = list(
      classified     = classified,
      null_range     = null_range,
      effect_label   = effect_label %||% classified$effect_label,
      credible_level = credible_level,
      n_fits         = length(fits)
    )
  )
  class(out) <- c("bayesma_interpretation", "list")
  out
}


# ============================================================================
# Classification: identify which fits are which
# ============================================================================

#' @noRd
classify_fits <- function(fits) {
  classes <- purrr::map(fits, class)

  pick_first <- function(predicate) {
    idx <- which(purrr::map_lgl(fits, predicate))
    if (length(idx) == 0) NULL else fits[[idx[1]]]
  }

  pick_all <- function(predicate) {
    idx <- which(purrr::map_lgl(fits, predicate))
    if (length(idx) == 0) list() else fits[idx]
  }

  primary <- pick_first(\(x) inherits(x, "bayesma") &&
                          !inherits(x, c("bayesma_robma", "bayesma_metareg",
                                         "bayesma_mv")))
  mv      <- pick_first(\(x) inherits(x, "bayesma_mv"))
  robma   <- pick_first(\(x) inherits(x, "bayesma_robma"))
  egger   <- pick_first(\(x) inherits(x, "bayesma_egger"))
  meta_r  <- pick_first(\(x) inherits(x, "bayesma_metareg"))
  sens    <- pick_first(\(x) inherits(x, "bayesma_robma_sensitivity"))
  cmp     <- pick_first(\(x) inherits(x, "bayesma_comparison"))

  bayesian_fits <- pick_all(\(x) inherits(x, "bayesma") ||
                              inherits(x, "bayesma_egger"))

  effect_label <- robma$meta$effect_label %||%
    primary$meta$effect_label %||%
    mv$meta$effect_label %||%
    meta_r$meta$effect_label

  list(
    primary       = primary,
    mv            = mv,
    robma         = robma,
    egger         = egger,
    meta_reg      = meta_r,
    sensitivity   = sens,
    comparison    = cmp,
    bayesian_fits = bayesian_fits,
    effect_label  = effect_label
  )
}


# ============================================================================
# Section: Overall effects
# ============================================================================

#' @noRd
interpret_effects_section <- function(c, credible_level) {

  source_fit <- c$robma %||% c$primary %||% c$mv %||% c$meta_reg
  if (is.null(source_fit)) return(NULL)

  pooled <- pooled_estimate(source_fit)
  if (is.null(pooled)) return(NULL)

  is_ratio <- isTRUE(c$effect_label %in% c("log_or", "log_rr"))
  ratio_name <- if (is_ratio) gsub("log_", "", c$effect_label) else NULL

  ci_excludes_zero <- (pooled$lower > 0) | (pooled$upper < 0)
  direction <- if (pooled$estimate > 0) "positive" else "negative"

  narrative <- paste0(
    "The pooled estimate is ", round(pooled$estimate, 3),
    " (", credible_level * 100, "% CrI: [",
    round(pooled$lower, 3), ", ", round(pooled$upper, 3), "])"
  )
  if (is_ratio) {
    narrative <- paste0(
      narrative, ". On the ratio scale: ", ratio_name, " = ",
      round(exp(pooled$estimate), 3),
      " [", round(exp(pooled$lower), 3), ", ",
      round(exp(pooled$upper), 3), "]"
    )
  }
  narrative <- paste0(
    narrative,
    if (ci_excludes_zero) {
      paste0(". The CrI excludes zero, indicating credible evidence of a ",
             direction, " effect.")
    } else {
      ". The CrI includes zero -- the effect is not unambiguously different from null."
    }
  )

  list(
    source        = class_label(source_fit),
    pooled        = pooled,
    is_ratio      = is_ratio,
    ratio_name    = ratio_name,
    direction     = direction,
    excludes_zero = ci_excludes_zero,
    narrative     = narrative
  )
}


# ============================================================================
# Section: Heterogeneity
# ============================================================================

#' @noRd
interpret_heterogeneity_section <- function(c, credible_level) {

  source_fit <- c$primary %||% c$meta_reg %||% c$mv
  if (is.null(source_fit)) return(NULL)

  tau_summary <- extract_tau_summary(source_fit)
  if (is.null(tau_summary)) return(NULL)

  pred_int <- source_fit$pred_interval %||% source_fit$pred_intervals

  level <- dplyr::case_when(
    tau_summary$median <  0.10 ~ "low",
    tau_summary$median <  0.30 ~ "moderate",
    tau_summary$median <  0.60 ~ "substantial",
    .default = "high"
  )

  narrative <- paste0(
    "Between-study heterogeneity (tau) has a posterior median of ",
    round(tau_summary$median, 3),
    " (", credible_level * 100, "% CrI: [",
    round(tau_summary$lower, 3), ", ",
    round(tau_summary$upper, 3),
    "]) -- ", level, " heterogeneity."
  )

  if (!is.null(pred_int) && nrow(pred_int) > 0) {
    pi_lower <- pred_int$lower[1]
    pi_upper <- pred_int$upper[1]
    narrative <- paste0(
      narrative,
      " The ", credible_level * 100, "% prediction interval for a new study is [",
      round(pi_lower, 3), ", ", round(pi_upper, 3), "]."
    )
  }

  list(
    tau            = tau_summary,
    level          = level,
    pred_interval  = pred_int,
    narrative      = narrative
  )
}


# ============================================================================
# Section: Publication bias / small-study effects
# ============================================================================

#' @noRd
interpret_bias_section <- function(c, credible_level) {

  parts <- list()

  if (!is.null(c$egger)) {
    egg <- c$egger
    beta_summary <- egg$beta_summary
    if (!is.null(beta_summary)) {
      dom_prob <- max(beta_summary$prob_positive, beta_summary$prob_negative)
      direction <- if (beta_summary$prob_negative > beta_summary$prob_positive) {
        "negative"
      } else {
        "positive"
      }
      strength <- evidence_label(dom_prob)

      parts$egger <- list(
        beta_estimate = beta_summary$estimate,
        beta_lower    = beta_summary$lower,
        beta_upper    = beta_summary$upper,
        prob_positive = beta_summary$prob_positive,
        prob_negative = beta_summary$prob_negative,
        direction     = direction,
        strength      = strength,
        narrative     = paste0(
          "Egger's regression: beta = ", round(beta_summary$estimate, 3),
          " [", round(beta_summary$lower, 3), ", ",
          round(beta_summary$upper, 3), "]. ",
          strength, " of ", direction, " small-study effects ",
          "(P = ", sprintf("%.1f%%", dom_prob * 100), ")."
        ),
        conclusion = egg$conclusion
      )
    }
  }

  if (!is.null(c$robma)) {
    rob <- c$robma
    pp_bias <- rob$posterior_probs$bias
    bf_bias <- rob$inclusion_bf$bias
    parts$robma_bias <- list(
      post_prob_bias  = pp_bias,
      inclusion_bf    = bf_bias,
      bf_label        = bf_label(bf_bias),
      narrative       = paste0(
        "RoBMA inclusion analysis: P(publication bias) = ",
        round(pp_bias, 3), " (BF = ", format_bf(bf_bias), "; ",
        bf_label(bf_bias), ")."
      )
    )
  }

  if (length(parts) == 0) return(NULL)

  combined_narrative <- paste(
    purrr::map_chr(parts, "narrative"),
    collapse = " "
  )

  list(
    egger     = parts$egger,
    robma     = parts$robma_bias,
    narrative = combined_narrative
  )
}


# ============================================================================
# Section: RoBMA model averaging
# ============================================================================

#' @noRd
interpret_robma_section <- function(c) {
  if (is.null(c$robma)) return(NULL)

  rob <- c$robma
  pp  <- rob$posterior_probs
  bf  <- rob$inclusion_bf
  nrp <- rob$meta$null_range_probs

  effect_evidence <- bf_label(bf$effect)
  hetero_evidence <- bf_label(bf$heterogeneity)
  bias_evidence   <- bf_label(bf$bias)

  narrative_lines <- c(
    paste0("RoBMA averaged across ",
           rob$meta$n_models %||% length(rob$component_fits), " models."),
    paste0("Effect inclusion: P(H1) = ", round(pp$effect, 3),
           ", BF = ", format_bf(bf$effect), " -- ", effect_evidence, "."),
    paste0("Heterogeneity inclusion: P = ", round(pp$heterogeneity, 3),
           ", BF = ", format_bf(bf$heterogeneity), " -- ", hetero_evidence, "."),
    paste0("Bias inclusion: P = ", round(pp$bias, 3),
           ", BF = ", format_bf(bf$bias), " -- ", bias_evidence, ".")
  )

  if (!is.null(nrp)) {
    narrative_lines <- c(
      narrative_lines,
      paste0("Direction probabilities (model-averaged posterior): ",
             "P(positive) = ", round(nrp$p_positive, 3),
             ", P(negative) = ", round(nrp$p_negative, 3),
             if (!is.null(nrp$p_null)) {
               paste0(", P(practically null) = ", round(nrp$p_null, 3))
             } else "", ".")
    )
  }

  list(
    posterior_probs = pp,
    inclusion_bf    = bf,
    null_range_probs = nrp,
    method          = rob$meta$method %||% "bridge",
    n_models        = rob$meta$n_models %||% length(rob$component_fits),
    effect_evidence = effect_evidence,
    hetero_evidence = hetero_evidence,
    bias_evidence   = bias_evidence,
    narrative       = paste(narrative_lines, collapse = " ")
  )
}


# ============================================================================
# Section: Model comparison (LOO / LOSO)
# ============================================================================

#' @noRd
interpret_comparison_section <- function(c) {
  if (is.null(c$comparison)) return(NULL)

  cmp <- c$comparison
  comparison_tbl <- cmp$comparison
  if (is.null(comparison_tbl) || nrow(comparison_tbl) == 0) return(NULL)

  best <- comparison_tbl |>
    dplyr::arrange(.data$rank) |>
    dplyr::slice(1)

  criterion <- cmp$criterion %||% "loo"
  narrative <- paste0(
    "Model comparison via ", toupper(criterion),
    ": best-supported model is ", best$model,
    if (!is.null(best$elpd_loo)) {
      paste0(" (ELPD = ", round(best$elpd_loo, 2), ")")
    } else "",
    "."
  )

  list(
    comparison_tbl = comparison_tbl,
    best_model     = best$model,
    criterion      = criterion,
    narrative      = narrative
  )
}


# ============================================================================
# Section: Sensitivity analysis (across priors)
# ============================================================================

#' @noRd
interpret_sensitivity_section <- function(c, null_range) {
  if (is.null(c$sensitivity)) return(NULL)

  sens <- c$sensitivity
  prior_names <- names(sens)
  prior_names <- prior_names[!prior_names %in% "meta"]
  if (length(prior_names) == 0) return(NULL)

  per_prior <- purrr::map(prior_names, \(name) {
    fit <- sens[[name]]
    pooled <- pooled_estimate(fit)
    if (is.null(pooled)) return(NULL)
    tibble::tibble(
      prior     = name,
      estimate  = pooled$estimate,
      lower     = pooled$lower,
      upper     = pooled$upper,
      pp_effect = fit$posterior_probs$effect %||% NA_real_
    )
  }) |>
    purrr::compact() |>
    purrr::list_rbind()

  if (nrow(per_prior) == 0) return(NULL)

  est_range <- range(per_prior$estimate)
  est_diff  <- diff(est_range)
  signs <- sign(per_prior$estimate)
  signs_consistent <- length(unique(signs[signs != 0])) <= 1

  robust <- est_diff < 0.10 && signs_consistent

  narrative <- paste0(
    "Sensitivity analysis across ", nrow(per_prior),
    " prior specifications: estimates range from ",
    round(est_range[1], 3), " to ", round(est_range[2], 3),
    " (span = ", round(est_diff, 3), "). ",
    if (robust) {
      "Conclusions are robust to prior choice."
    } else if (signs_consistent) {
      "Sign is consistent across priors, but magnitude varies -- interpret point estimates with caution."
    } else {
      "Sign of effect is sensitive to prior choice -- conclusions are NOT robust."
    }
  )

  list(
    per_prior        = per_prior,
    estimate_range   = est_range,
    signs_consistent = signs_consistent,
    robust           = robust,
    narrative        = narrative
  )
}


# ============================================================================
# Section: Convergence diagnostics
# ============================================================================

#' @noRd
interpret_diagnostics_section <- function(c) {

  fits <- c$bayesian_fits
  if (length(fits) == 0) return(NULL)

  summaries <- purrr::map(fits, \(f) {
    s <- f$summary
    if (is.null(s) || !"rhat" %in% names(s)) return(NULL)
    tibble::tibble(
      model        = class_label(f),
      max_rhat     = max(s$rhat, na.rm = TRUE),
      min_ess_bulk = if ("ess_bulk" %in% names(s)) {
        min(s$ess_bulk, na.rm = TRUE)
      } else NA_real_,
      min_ess_tail = if ("ess_tail" %in% names(s)) {
        min(s$ess_tail, na.rm = TRUE)
      } else NA_real_,
      n_params     = nrow(s)
    )
  }) |>
    purrr::compact() |>
    purrr::list_rbind()

  if (nrow(summaries) == 0) return(NULL)

  worst_rhat <- max(summaries$max_rhat, na.rm = TRUE)
  worst_ess  <- min(summaries$min_ess_bulk, na.rm = TRUE)

  rhat_status <- dplyr::case_when(
    worst_rhat < 1.01 ~ "excellent",
    worst_rhat < 1.05 ~ "acceptable",
    worst_rhat < 1.10 ~ "borderline",
    .default = "poor"
  )
  ess_status <- dplyr::case_when(
    is.na(worst_ess)  ~ "not assessed",
    worst_ess >= 1000 ~ "excellent",
    worst_ess >= 400  ~ "acceptable",
    worst_ess >= 100  ~ "borderline",
    .default = "poor"
  )

  narrative <- paste0(
    "Convergence diagnostics across ", nrow(summaries),
    " fit(s): worst Rhat = ", round(worst_rhat, 3),
    " (", rhat_status, "), worst bulk ESS = ",
    if (is.na(worst_ess)) "NA" else round(worst_ess, 0),
    " (", ess_status, "). ",
    if (rhat_status %in% c("excellent", "acceptable") &&
        ess_status %in% c("excellent", "acceptable")) {
      "Chains have mixed adequately."
    } else {
      "Some parameters show convergence issues -- inspect traces, increase iterations, or revise priors before relying on inference."
    }
  )

  list(
    per_model    = summaries,
    worst_rhat   = worst_rhat,
    worst_ess    = worst_ess,
    rhat_status  = rhat_status,
    ess_status   = ess_status,
    narrative    = narrative
  )
}


# ============================================================================
# Synthesised conclusions
# ============================================================================

#' @noRd
interpret_overall_conclusion <- function(eff, het, bias, robma, sens, diag) {
  bullets <- character()

  if (!is.null(eff)) {
    bullets <- c(bullets, paste0("Effect: ", eff$narrative))
  }
  if (!is.null(het)) {
    bullets <- c(bullets, paste0("Heterogeneity: ", het$narrative))
  }
  if (!is.null(bias)) {
    bullets <- c(bullets, paste0("Bias: ", bias$narrative))
  }
  if (!is.null(robma)) {
    bullets <- c(bullets, paste0("Model averaging: ", robma$narrative))
  }
  if (!is.null(sens)) {
    bullets <- c(bullets, paste0("Sensitivity: ", sens$narrative))
  }
  if (!is.null(diag)) {
    bullets <- c(bullets, paste0("Diagnostics: ", diag$narrative))
  }

  bullets
}


# ============================================================================
# Print method
# ============================================================================

#' @export
print.bayesma_interpretation <- function(x, ...) {
  cli::cli_h1("Bayesian Meta-Analysis Interpretation")
  cli::cli_text("Synthesised across {x$meta$n_fits} fit{?s}.")
  cli::cli_text("")

  if (!is.null(x$effects)) {
    cli::cli_h2("1. Overall Effect")
    cli::cli_text(x$effects$narrative)
    cli::cli_text("")
  }

  if (!is.null(x$heterogeneity)) {
    cli::cli_h2("2. Heterogeneity")
    cli::cli_text(x$heterogeneity$narrative)
    cli::cli_text("")
  }

  if (!is.null(x$bias)) {
    cli::cli_h2("3. Publication Bias / Small-study Effects")
    cli::cli_text(x$bias$narrative)
    if (!is.null(x$bias$egger$conclusion)) {
      cli::cli_text("")
      cli::cli_text(x$bias$egger$conclusion)
    }
    cli::cli_text("")
  }

  if (!is.null(x$robma)) {
    cli::cli_h2("4. RoBMA Model Averaging")
    cli::cli_text(x$robma$narrative)
    cli::cli_text("")
  }

  if (!is.null(x$comparison)) {
    cli::cli_h2("5. Model Comparison")
    cli::cli_text(x$comparison$narrative)
    cli::cli_text("")
  }

  if (!is.null(x$sensitivity)) {
    cli::cli_h2("6. Sensitivity Analysis")
    cli::cli_text(x$sensitivity$narrative)
    cli::cli_text("")
    print(x$sensitivity$per_prior)
    cli::cli_text("")
  }

  if (!is.null(x$diagnostics)) {
    cli::cli_h2("7. Convergence Diagnostics")
    cli::cli_text(x$diagnostics$narrative)
    cli::cli_text("")
  }

  if (length(x$conclusions) > 0) {
    cli::cli_h2("Summary")
    purrr::walk(x$conclusions, \(b) cli::cli_li(b))
  }

  invisible(x)
}


# ============================================================================
# Helpers
# ============================================================================

#' @noRd
pooled_estimate <- function(fit) {
  fd <- fit$forest_df
  if (is.null(fd) || !"type" %in% names(fd)) return(NULL)
  pooled <- dplyr::filter(fd, .data$type == "pooled")
  if (nrow(pooled) == 0) return(NULL)
  list(
    estimate = pooled$estimate[1],
    lower    = pooled$lower[1],
    upper    = pooled$upper[1]
  )
}

#' @noRd
extract_tau_summary <- function(fit) {
  s <- fit$summary
  if (is.null(s)) return(NULL)
  tau_row <- dplyr::filter(s, .data$variable == "tau")
  if (nrow(tau_row) == 0) return(NULL)
  list(
    median = tau_row$median[1],
    lower  = tau_row$q5[1] %||% tau_row$lower[1] %||% NA_real_,
    upper  = tau_row$q95[1] %||% tau_row$upper[1] %||% NA_real_
  )
}

#' @noRd
class_label <- function(fit) {
  cls <- class(fit)
  cls[!cls %in% c("list", "S7_object")][1]
}

#' @noRd
evidence_label <- function(prob) {
  dplyr::case_when(
    prob >= 0.99 ~ "Very strong evidence",
    prob >= 0.95 ~ "Strong evidence",
    prob >= 0.90 ~ "Substantial evidence",
    prob >= 0.80 ~ "Moderate evidence",
    prob >= 0.70 ~ "Weak evidence",
    .default = "Little evidence"
  )
}

#' @noRd
bf_label <- function(bf) {
  if (is.na(bf)) return("BF unavailable")
  if (is.infinite(bf)) return("decisive evidence for inclusion")
  if (bf <= 0) return("decisive evidence against inclusion")
  log_bf <- log(bf)
  dplyr::case_when(
    log_bf >  log(100) ~ "extreme evidence for inclusion",
    log_bf >  log(30)  ~ "very strong evidence for inclusion",
    log_bf >  log(10)  ~ "strong evidence for inclusion",
    log_bf >  log(3)   ~ "moderate evidence for inclusion",
    log_bf >  log(1)   ~ "anecdotal evidence for inclusion",
    log_bf >  log(1/3) ~ "anecdotal evidence against inclusion",
    log_bf >  log(1/10) ~ "moderate evidence against inclusion",
    log_bf >  log(1/30) ~ "strong evidence against inclusion",
    .default = "very strong evidence against inclusion"
  )
}

#' @noRd
format_bf <- function(val) {
  if (is.na(val))         return("NA")
  if (is.infinite(val))   return("Inf")
  if (val < 0.01)         return(formatC(val, format = "e", digits = 1))
  if (val > 1000)         return(formatC(val, format = "e", digits = 1))
  as.character(round(val, 2))
}
