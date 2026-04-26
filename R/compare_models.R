#' Compare Multiple Bayesian Meta-Analysis Models
#'
#' @description
#' Compares fitted \code{bayesma} models using LOSO-CV (leave-one-study-out
#' cross-validation) as the primary ranking criterion, plus within-stage
#' LOO-CV for diagnostics. LOSO-CV evaluates all models on the same
#' question — "predict a held-out study" — on the effect-size scale,
#' making it valid for comparing across one-stage and two-stage models.
#'
#' @param ... Named or unnamed \code{bayesma} objects to compare. If unnamed,
#'   models are labeled "Model 1", "Model 2", etc.
#' @param data The original data frame used to fit all models. Required for
#'   LOSO-CV (the default).
#' @param studyvar Column identifying studies in `data` (unquoted). Required
#'   for LOSO-CV.
#' @param criterion Character. Primary criterion for ranking.
#'   \code{"loso"} (default): leave-one-study-out CRPS on the effect-size
#'   scale. Comparable across stages.
#'   \code{"loo"}: within-stage LOO-CV only.
#' @param loso Logical. Run LOSO-CV. Default \code{TRUE}. Set to \code{FALSE}
#'   for a fast within-stage-only comparison.
#' @param coverage_levels Numeric vector. Nominal coverage levels for
#'   calibration assessment in LOSO-CV. Default: \code{c(0.50, 0.80, 0.90, 0.95)}.
#' @param moment_match Logical. Attempt moment matching for LOO-CV.
#'   Default \code{FALSE}.
#' @param cores Integer. Number of cores for LOO computation.
#'   Default is \code{getOption("mc.cores", 1)}.
#' @param quiet Logical. Suppress progress messages. Default \code{FALSE}.
#'
#' @return An object of class \code{"bayesma_comparison"} containing:
#' \describe{
#'   \item{comparison}{Tibble with all comparison metrics}
#'   \item{loso_list}{List of LOSO-CV results per model (if computed)}
#'   \item{loo_compare}{List of \code{loo::loo_compare()} per stage group}
#'   \item{loo_list}{List of LOO objects per model}
#'   \item{waic_list}{List of WAIC objects per model}
#'   \item{diagnostics}{LOO diagnostics (Pareto k)}
#'   \item{model_names}{Character vector of model names}
#' }
#'
#' @details
#' ## Why LOSO-CV?
#'
#' LOO-CV and WAIC operate on each model's native likelihood. One-stage
#' models condition on arm-level counts; two-stage models condition on
#' summary effect sizes. These likelihoods are on different scales, so
#' their ELPD values are **not comparable across stages**.
#'
#' LOSO-CV sidesteps this by asking each model the same question:
#' *"Given all studies except study s, what is your predictive distribution
#' for study s's effect size?"*
#'
#' Both model types answer on the same scale (the effect-size scale),
#' making the comparison valid.
#'
#' ## LOSO Metrics
#'
#' \describe{
#'   \item{LOSO-CRPS}{Continuous Ranked Probability Score averaged over
#'     held-out studies. A proper scoring rule that evaluates the full
#'     predictive distribution. Lower is better.}
#'   \item{LOSO-ELPD}{Mean log predictive density at the held-out study.
#'     Higher is better. Analogous to LOO-ELPD but truly out-of-sample.}
#'   \item{LOSO-Coverage}{At each nominal level, what proportion of
#'     held-out yi fall inside the corresponding prediction interval?
#'     Well-calibrated models have empirical ≈ nominal.}
#'   \item{LOSO-Miscalibration}{Mean |empirical − nominal| across
#'     coverage levels. Zero is perfect.}
#' }
#'
#' ## Within-Stage LOO
#'
#' LOO-CV is still computed for within-stage diagnostics (Pareto k values
#' identify influential studies) and within-stage ranking when all models
#' share the same stage.
#'
#' ## Computational Cost
#'
#' LOSO-CV refits each model S times (once per study). For M models with
#' S studies, this is M × S fits. Set \code{loso = FALSE} for a fast
#' within-stage-only comparison, or use \code{quiet = TRUE} to suppress
#' per-fold messages.
#'
#' @examples
#' \dontrun{
#' comparison <- compare_models(
#'   "RE (two-stage)" = mod_2s_re,
#'   "RE (one-stage)"  = mod_1s_re,
#'   "Heavy-tailed"    = mod_2s_t,
#'   data = dat,
#'   studyvar = author
#' )
#' print(comparison)
#' compare_plot(comparison, type = "loso_crps")
#' compare_plot(comparison, type = "calibration")
#' compare_table(comparison)
#'
#' # Fast within-stage only (no refitting)
#' compare_models(mod1, mod2, criterion = "loo", loso = FALSE)
#' }
#'
#' @seealso
#' \code{\link[loo]{loo}}, \code{\link[loo]{loo_compare}}
#'
#' @export
compare_models <- function(
    ...,
    data              = NULL,
    studyvar          = NULL,
    criterion         = c("loso", "loo"),
    loso              = TRUE,
    coverage_levels   = c(0.50, 0.80, 0.90, 0.95),
    moment_match      = FALSE,
    cores             = getOption("mc.cores", 1),
    quiet             = FALSE
) {
  inform <- if (quiet) \(...) invisible(NULL) else cli::cli_inform

  dots      <- rlang::list2(...)
  criterion <- rlang::arg_match(criterion)

  # ---- Separate models from any stray named args ----
  is_bayesma <- purrr::map_lgl(dots, ~ inherits(.x, "bayesma"))
  models     <- dots[is_bayesma]

  if (length(models) < 2) {
    cli::cli_abort(c(
      "At least 2 {.cls bayesma} objects are required for comparison.",
      "i" = "Pass models as the first positional arguments."
    ), call = rlang::caller_env())
  }

  if (criterion == "loso" && !loso) {
    cli::cli_abort(c(
      '{.arg criterion} = "loso" requires {.arg loso} = TRUE.',
      "i" = 'Use {.code criterion = "loo"} with {.code loso = FALSE} for within-stage comparison.'
    ), call = rlang::caller_env())
  }

  # ---- Tidyeval: resolve studyvar ----
  studyvar_quo <- rlang::enquo(studyvar)

  if (loso) {
    if (is.null(data)) {
      cli::cli_abort(c(
        "LOSO-CV requires {.arg data}.",
        "i" = "Provide the original data frame used to fit the models.",
        "i" = "Or set {.code loso = FALSE, criterion = \"loo\"} for within-stage only."
      ), call = rlang::caller_env())
    }
    if (rlang::quo_is_null(studyvar_quo)) {
      cli::cli_abort(c(
        "LOSO-CV requires {.arg studyvar}.",
        "i" = "Provide the bare column name identifying studies (e.g. {.code studyvar = author})."
      ), call = rlang::caller_env())
    }

    studyvar_nm <- tidyselect_column_name(data, {{ studyvar }})
  } else {
    studyvar_nm <- NULL
  }

  # ---- Model names ----
  model_names <- names(models)
  if (is.null(model_names) || any(model_names == "")) {
    model_names <- purrr::map_chr(seq_along(models), function(i) {
      nm <- names(models)[i]
      if (is.null(nm) || nm == "") paste0("Model ", i) else nm
    })
  }
  names(models) <- model_names

  # ---- Stage detection ----
  stages <- purrr::map_chr(models, ~ .x$meta$stage %||% "unknown")
  has_mixed_stages <- length(unique(stages)) > 1

  if (has_mixed_stages && !loso) {
    cli::cli_warn(c(
      "Models span both one-stage and two-stage likelihoods.",
      "!" = "LOO/WAIC are {.strong not comparable} across stages.",
      "i" = "Set {.code loso = TRUE} for valid cross-stage comparison."
    ))
  }

  # ---- Check log_lik availability ----
  has_log_lik <- purrr::map_lgl(models, function(m) {
    all_vars <- posterior::variables(m$fit$draws())
    any(grepl("^log_lik", all_vars))
  })

  if (!all(has_log_lik)) {
    cli::cli_warn(c(
      "Models without {.code log_lik}: {.val {model_names[!has_log_lik]}}.",
      "i" = "LOO/WAIC will be skipped for these models."
    ))
  }

  # ---- LOO-CV (within-stage diagnostics) ----
  inform("Computing LOO-CV...")

  loo_waic <- purrr::map(model_names[has_log_lik], function(nm) {
    m       <- models[[nm]]
    ll_data <- extract_log_lik(m, aggregate_to_study = FALSE)
    r_eff   <- loo::relative_eff(exp(ll_data$log_lik),
                                 chain_id = ll_data$chain_id, cores = cores)

    loo_obj <- if (moment_match) {
      tryCatch(
        loo::loo(ll_data$log_lik, r_eff = r_eff, cores = cores,
                 moment_match = TRUE),
        error = function(e) {
          cli::cli_warn("Moment matching failed for {.val {nm}}.")
          loo::loo(ll_data$log_lik, r_eff = r_eff, cores = cores)
        }
      )
    } else {
      loo::loo(ll_data$log_lik, r_eff = r_eff, cores = cores)
    }

    waic_obj <- loo::waic(ll_data$log_lik)
    list(loo = loo_obj, waic = waic_obj)
  }) |>
    stats::setNames(model_names[has_log_lik])

  null_entries <- purrr::map(model_names[!has_log_lik], ~ NULL) |>
    stats::setNames(model_names[!has_log_lik])

  loo_list  <- c(purrr::map(loo_waic, "loo"),  null_entries)
  waic_list <- c(purrr::map(loo_waic, "waic"), null_entries)

  # ---- LOO compare within stage groups ----
  loo_compare_list <- unique(stages) |>
    purrr::map(function(stg) {
      stg_nms <- model_names[stages == stg & has_log_lik]
      if (length(stg_nms) >= 2) loo::loo_compare(loo_list[stg_nms]) else NULL
    }) |>
    stats::setNames(unique(stages)) |>
    purrr::compact()

  # ---- LOO diagnostics ----
  diagnostics_df <- purrr::map(loo_list[has_log_lik], function(loo_obj) {
    if (is.null(loo_obj)) return(NULL)
    pk <- loo::pareto_k_values(loo_obj)
    tibble::tibble(
      n_high_k      = sum(pk > 0.7),
      n_very_high_k = sum(pk > 1.0),
      max_k         = max(pk),
      prop_ok       = mean(pk <= 0.7)
    )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows(.id = "model")

  # ---- LOSO-CV ----
  loso_list <- NULL
  if (loso) {
    study_vec <- dplyr::pull(data, {{ studyvar }})
    n_studies <- length(unique(study_vec))
    n_total   <- length(models) * n_studies
    inform("Running LOSO-CV: {length(models)} models \u00d7 {n_studies} studies = {n_total} refits")

    loso_list <- purrr::map2(models, model_names, function(m, nm) {
      inform("LOSO-CV for {.val {nm}}...")
      run_loso_cv(
        model           = m,
        data            = data,
        studyvar        = studyvar_nm,
        coverage_levels = coverage_levels,
        quiet           = quiet
      )
    }) |>
      stats::setNames(model_names)
  }

  # ---- Build comparison table ----
  comparison_df <- purrr::map(model_names, function(nm) {
    row <- tibble::tibble(model = nm, stage = stages[[nm]])

    if (!is.null(loo_list[[nm]])) {
      loo_est  <- loo_list[[nm]]$estimates
      waic_est <- waic_list[[nm]]$estimates
      row$elpd_loo     <- loo_est["elpd_loo", "Estimate"]
      row$se_elpd_loo  <- loo_est["elpd_loo", "SE"]
      row$p_loo        <- loo_est["p_loo", "Estimate"]
    } else {
      row$elpd_loo <- row$se_elpd_loo <- row$p_loo <- NA_real_
    }

    if (!is.null(loso_list)) {
      ls <- loso_list[[nm]]
      row$loso_crps   <- ls$mean_crps
      row$loso_elpd   <- ls$mean_log_pd
      row$loso_miscal <- ls$miscalibration

      row_95 <- ls$calibration |>
        dplyr::filter(.data$nominal == 0.95)
      row$loso_cover_95 <- if (nrow(row_95) > 0) row_95$empirical else NA_real_
      row$loso_width_95 <- if (nrow(row_95) > 0) row_95$mean_interval_width else NA_real_
    }

    row
  }) |> purrr::list_rbind()

  # ---- Within-stage ELPD differences ----
  comparison_df$elpd_diff <- NA_real_
  comparison_df$se_diff   <- NA_real_

  elpd_diff_updates <- purrr::map(names(loo_compare_list), function(stg) {
    lc_df       <- as.data.frame(loo_compare_list[[stg]])
    lc_df$model <- rownames(lc_df)
    tibble::tibble(
      model     = lc_df$model,
      elpd_diff = lc_df$elpd_diff,
      se_diff   = lc_df$se_diff
    )
  }) |> purrr::list_rbind()

  if (nrow(elpd_diff_updates) > 0) {
    comparison_df <- comparison_df |>
      dplyr::select(-"elpd_diff", -"se_diff") |>
      dplyr::left_join(elpd_diff_updates, by = "model")
  }

  # ---- Merge diagnostics ----
  if (nrow(diagnostics_df) > 0) {
    comparison_df <- comparison_df |>
      dplyr::left_join(diagnostics_df, by = "model")
  }

  # ---- Rank ----
  comparison_df <- switch(criterion,
                          loso = comparison_df |>
                            dplyr::arrange(.data$loso_crps) |>
                            dplyr::mutate(rank = dplyr::row_number()),
                          loo = comparison_df |>
                            dplyr::mutate(
                              rank = dplyr::row_number(dplyr::desc(.data$elpd_loo)),
                              .by = "stage"
                            )
  )

  out <- list(
    comparison   = comparison_df,
    loso_list    = loso_list,
    loo_compare  = loo_compare_list,
    loo_list     = loo_list,
    waic_list    = waic_list,
    diagnostics  = diagnostics_df,
    model_names  = model_names,
    stages       = stages,
    mixed_stages = has_mixed_stages,
    criterion    = criterion,
    n_models     = length(models)
  )
  class(out) <- "bayesma_comparison"
  out
}


# ============================================================================
# Internal: resolve a tidyselect column to a string name
# ============================================================================

#' @noRd
tidyselect_column_name <- function(data, var) {
  nm <- rlang::as_name(rlang::enquo(var))
  if (!nm %in% names(data)) {
    cli::cli_abort(
      "Column {.val {nm}} not found in {.arg data}.",
      call = rlang::caller_env()
    )
  }
  nm
}


# ============================================================================
# Print
# ============================================================================

#' @export
print.bayesma_comparison <- function(x, digits = 3, ...) {

  cli::cli_h1("Bayesian Meta-Analysis Model Comparison")
  cli::cli_text("{.strong {x$n_models}} models | Ranking: {.val {x$criterion}}")

  if (!is.null(x$loso_list)) {
    cli::cli_text("")
    cli::cli_h2("LOSO Cross-Validation (cross-stage, effect-size scale)")

    loso_df <- x$comparison |>
      dplyr::arrange(.data$rank) |>
      dplyr::select(dplyr::any_of(c(
        "rank", "model", "stage",
        "loso_crps", "loso_elpd", "loso_cover_95",
        "loso_width_95", "loso_miscal"
      ))) |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits)))

    print(loso_df, n = Inf)

    cli::cli_text("")
    cli::cli_alert_info("LOSO-CRPS: lower is better (primary ranking criterion).")
    cli::cli_alert_info("LOSO-ELPD: higher is better.")
    cli::cli_alert_info("Coverage near 0.95 and miscalibration near 0 indicate well-calibrated predictions.")
  }

  cli::cli_text("")
  purrr::walk(unique(x$stages), function(stg) {
    stg_df <- x$comparison |>
      dplyr::filter(.data$stage == stg)

    if (nrow(stg_df) < 2) return(invisible(NULL))

    cli::cli_h2("{stg} models (within-stage LOO)")

    display <- stg_df |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits))) |>
      dplyr::select(dplyr::any_of(c(
        "model", "elpd_loo", "se_elpd_loo",
        "elpd_diff", "se_diff", "p_loo", "n_high_k"
      )))
    print(display, n = Inf)
  })

  best <- x$comparison |>
    dplyr::filter(.data$rank == min(.data$rank, na.rm = TRUE)) |>
    dplyr::pull("model")
  cli::cli_text("")
  cli::cli_alert_success("Best model: {.val {best}}")

  if (!is.null(x$diagnostics) && nrow(x$diagnostics) > 0 &&
      any(x$diagnostics$n_high_k > 0, na.rm = TRUE)) {
    cli::cli_alert_warning("Some models have observations with high Pareto k (> 0.7).")
  }

  invisible(x)
}


# ============================================================================
# Plots
# ============================================================================

#' Plot Model Comparison Results
#'
#' @param x A \code{bayesma_comparison} object.
#' @param type Character. Plot type:
#'   \code{"loso_crps"} (default): per-study CRPS from LOSO-CV.
#'   \code{"calibration"}: LOSO calibration curves.
#'   \code{"elpd"}: within-stage ELPD comparison.
#'   \code{"pareto_k"}: LOO Pareto k diagnostics.
#' @param ... Additional arguments (unused).
#' @return A \code{ggplot2} object.
#' @export
compare_plot <- function(x, type = c("loso_crps", "calibration",
                                     "elpd", "pareto_k"), ...) {
  type <- rlang::arg_match(type)

  switch(type,
         loso_crps   = plot_loso_crps(x),
         calibration = plot_loso_calibration(x),
         elpd        = plot_elpd_comparison(x),
         pareto_k    = plot_pareto_k(x)
  )
}


#' @noRd
plot_loso_crps <- function(x) {
  if (is.null(x$loso_list)) {
    cli::cli_abort("No LOSO data. Use {.code loso = TRUE}.")
  }

  crps_df <- purrr::map2(x$loso_list, names(x$loso_list), function(ls, nm) {
    ls$per_study |>
      dplyr::select("study", "crps") |>
      dplyr::mutate(model = nm)
  }) |> purrr::list_rbind()

  summary_df <- crps_df |>
    dplyr::summarise(mean_crps = mean(.data$crps, na.rm = TRUE), .by = "model") |>
    dplyr::mutate(model = forcats::fct_reorder(.data$model, .data$mean_crps))

  crps_df <- crps_df |>
    dplyr::mutate(model = factor(.data$model, levels = levels(summary_df$model)))

  ggplot2::ggplot(crps_df, ggplot2::aes(x = .data$model, y = .data$crps)) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 2) +
    ggplot2::geom_crossbar(
      data = summary_df,
      ggplot2::aes(x = .data$model, y = .data$mean_crps,
                   ymin = .data$mean_crps, ymax = .data$mean_crps),
      width = 0.5, colour = "firebrick", linewidth = 0.8
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "LOSO-CV: Per-Study CRPS",
      subtitle = "Lower is better. Red bar = mean. Truly out-of-sample.",
      x = NULL, y = "CRPS"
    ) +
    ggplot2::theme_minimal()
}


#' @noRd
plot_loso_calibration <- function(x) {
  if (is.null(x$loso_list)) {
    cli::cli_abort("No LOSO data. Use {.code loso = TRUE}.")
  }

  cal_df <- purrr::map2(x$loso_list, names(x$loso_list), function(ls, nm) {
    ls$calibration |> dplyr::mutate(model = nm)
  }) |> purrr::list_rbind()

  ggplot2::ggplot(cal_df, ggplot2::aes(x = .data$nominal,
                                       y = .data$empirical,
                                       colour = .data$model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(),
                                limits = c(0, 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                limits = c(0, 1)) +
    ggplot2::labs(
      title = "LOSO-CV: Prediction Calibration",
      subtitle = "Points on the diagonal = well-calibrated. Truly out-of-sample.",
      x = "Nominal coverage",
      y = "Empirical coverage",
      colour = "Model"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}


#' @noRd
plot_elpd_comparison <- function(x) {
  plot_df <- x$comparison |>
    dplyr::filter(!is.na(.data$elpd_loo)) |>
    dplyr::mutate(model = forcats::fct_reorder(.data$model, .data$elpd_loo))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$elpd_loo, y = .data$model)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = .data$elpd_loo - 2 * .data$se_elpd_loo,
        xmax = .data$elpd_loo + 2 * .data$se_elpd_loo
      ),
      height = 0.2
    ) +
    ggplot2::facet_wrap(~ .data$stage, scales = "free_x") +
    ggplot2::labs(
      title = "ELPD Comparison (within stage only)",
      subtitle = "Higher is better. Do not compare across stage panels.",
      x = "ELPD (LOO-CV)", y = NULL
    ) +
    ggplot2::theme_minimal()
}


#' @noRd
plot_pareto_k <- function(x) {
  pareto_df <- purrr::map2(x$loo_list, names(x$loo_list), function(loo_obj, nm) {
    if (is.null(loo_obj)) return(NULL)
    k_vals <- loo::pareto_k_values(loo_obj)
    tibble::tibble(
      model       = nm,
      observation = seq_along(k_vals),
      pareto_k    = k_vals
    )
  }) |>
    purrr::compact() |>
    purrr::list_rbind()

  ggplot2::ggplot(pareto_df, ggplot2::aes(x = .data$observation,
                                          y = .data$pareto_k)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed", colour = "orange") +
    ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed", colour = "red") +
    ggplot2::facet_wrap(~ .data$model, scales = "free_x") +
    ggplot2::labs(
      title = "Pareto k Diagnostics",
      subtitle = "> 0.7 (orange): influential; > 1.0 (red): very problematic",
      x = "Observation", y = "Pareto k"
    ) +
    ggplot2::theme_minimal()
}


# ============================================================================
# Table
# ============================================================================

#' Create a Table for Model Comparison Results
#'
#' Produces a publication-ready \code{gt} table with LOSO-CV metrics
#' (when available) and within-stage LOO diagnostics.
#'
#' @param x A \code{bayesma_comparison} object.
#' @param digits Integer. Decimal places. Default 2.
#' @param include_loo Logical. Include within-stage LOO columns. Default \code{TRUE}.
#' @param include_diagnostics Logical. Include Pareto k column. Default \code{TRUE}.
#'
#' @return A \code{gt} table object.
#'
#' @examples
#' \dontrun{
#' compare_table(comparison)
#' compare_table(comparison) |> gt::gtsave("comparison.html")
#' }
#'
#' @export
compare_table <- function(
    x,
    digits = 2,
    include_loo = TRUE,
    include_diagnostics = TRUE
) {
  if (!inherits(x, "bayesma_comparison")) {
    cli::cli_abort("{.arg x} must be a {.cls bayesma_comparison} object.")
  }

  fmt_1 <- paste0("%.", digits, "f")
  fmt_2 <- paste0("%.", digits, "f (%.", digits, "f)")

  tbl_df <- x$comparison |> dplyr::arrange(.data$rank)

  select_cols <- c("rank", "model", "stage")

  # ---- LOSO columns ----
  if ("loso_crps" %in% names(tbl_df)) {
    tbl_df <- tbl_df |>
      dplyr::mutate(
        loso_crps_fmt   = sprintf(fmt_1, .data$loso_crps),
        loso_elpd_fmt   = sprintf(fmt_1, .data$loso_elpd),
        loso_cover_fmt  = sprintf("%.0f%%", .data$loso_cover_95 * 100),
        loso_width_fmt  = sprintf(fmt_1, .data$loso_width_95),
        loso_miscal_fmt = sprintf(fmt_1, .data$loso_miscal)
      )
    select_cols <- c(select_cols,
                     "loso_crps_fmt", "loso_elpd_fmt",
                     "loso_cover_fmt", "loso_width_fmt", "loso_miscal_fmt")
  }

  # ---- LOO columns ----
  if (include_loo) {
    tbl_df <- tbl_df |>
      dplyr::mutate(
        elpd_loo_fmt = dplyr::if_else(
          is.na(.data$elpd_loo), "\u2014",
          sprintf(fmt_2, .data$elpd_loo, .data$se_elpd_loo)
        ),
        elpd_diff_fmt = dplyr::case_when(
          is.na(.data$elpd_diff) ~ "\u2014",
          .data$elpd_diff == 0   ~ "\u2014",
          .default = sprintf(fmt_2, .data$elpd_diff, .data$se_diff)
        )
      )
    select_cols <- c(select_cols, "elpd_loo_fmt", "elpd_diff_fmt")
  }

  # ---- Diagnostics column ----
  if (include_diagnostics && "n_high_k" %in% names(tbl_df)) {
    tbl_df <- tbl_df |>
      dplyr::mutate(
        pareto_k_fmt = dplyr::if_else(
          is.na(.data$n_high_k), "\u2014",
          dplyr::if_else(.data$n_high_k == 0, "\u2714",
                         sprintf("%d high", .data$n_high_k))
        )
      )
    select_cols <- c(select_cols, "pareto_k_fmt")
  }

  tbl_df <- tbl_df |> dplyr::select(dplyr::all_of(select_cols))

  # ---- Column labels ----
  col_labels <- list(
    rank  = "Rank",
    model = "Model",
    stage = "Stage"
  )
  if ("loso_crps_fmt" %in% select_cols) {
    col_labels$loso_crps_fmt   <- "CRPS"
    col_labels$loso_elpd_fmt   <- "ELPD"
    col_labels$loso_cover_fmt  <- "95% Cov."
    col_labels$loso_width_fmt  <- "95% Width"
    col_labels$loso_miscal_fmt <- "Miscal."
  }
  if ("elpd_loo_fmt" %in% select_cols) {
    col_labels$elpd_loo_fmt  <- "ELPD (SE)"
    col_labels$elpd_diff_fmt <- "\u0394ELPD (SE)"
  }
  if ("pareto_k_fmt" %in% select_cols) {
    col_labels$pareto_k_fmt <- "Pareto k"
  }

  # ---- Build gt ----
  gt_tbl <- gt::gt(tbl_df) |>
    gt::cols_label(.list = col_labels) |>
    gt::tab_header(
      title = "Bayesian Meta-Analysis Model Comparison",
      subtitle = sprintf("Ranked by %s", toupper(x$criterion))
    )

  loso_cols <- intersect(
    c("loso_crps_fmt", "loso_elpd_fmt", "loso_cover_fmt",
      "loso_width_fmt", "loso_miscal_fmt"),
    select_cols
  )
  if (length(loso_cols) > 0) {
    gt_tbl <- gt_tbl |>
      gt::tab_spanner(label = "LOSO-CV (out-of-sample)",
                      columns = dplyr::all_of(loso_cols))
  }

  loo_cols <- intersect(c("elpd_loo_fmt", "elpd_diff_fmt"), select_cols)
  if (length(loo_cols) > 0) {
    gt_tbl <- gt_tbl |>
      gt::tab_spanner(label = "Within-Stage LOO",
                      columns = dplyr::all_of(loo_cols))
  }

  gt_tbl <- gt_tbl |>
    gt::tab_source_note(
      source_note = gt::md(paste0(
        "*LOSO*: leave-one-study-out CV on the effect-size scale (comparable across stages). ",
        "*CRPS*: lower is better. *ELPD*: higher is better. ",
        "*Miscal.*: mean |empirical \u2212 nominal| (0 = perfect). ",
        "*Within-stage LOO*: only comparable within the same stage."
      ))
    )

  if (x$mixed_stages && length(loo_cols) > 0) {
    gt_tbl <- gt_tbl |>
      gt::tab_footnote(
        footnote = "ELPD and \u0394ELPD are only comparable within the same stage.",
        locations = gt::cells_column_spanners(spanners = "Within-Stage LOO")
      )
  }

  # ---- Styling ----
  best_row <- which(tbl_df$rank == min(tbl_df$rank))[1]
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = list(gt::cell_text(weight = "bold"),
                   gt::cell_fill(color = "#E8F4E8")),
      locations = gt::cells_body(rows = best_row)
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    ) |>
    gt::cols_align(align = "center", columns = -"model") |>
    gt::cols_align(align = "left", columns = "model") |>
    gt::tab_options(
      table.font.size                    = gt::px(12),
      heading.title.font.size            = gt::px(14),
      heading.title.font.weight          = "bold",
      heading.subtitle.font.size         = gt::px(11),
      column_labels.font.size            = gt::px(11),
      column_labels.border.bottom.color  = "grey40",
      column_labels.border.bottom.width  = gt::px(2),
      table_body.border.bottom.color     = "grey40",
      table.border.top.color             = "grey40",
      table.border.top.width             = gt::px(2),
      table.border.bottom.width          = gt::px(2),
      data_row.padding                   = gt::px(6)
    )

  if ("pareto_k_fmt" %in% select_cols && "n_high_k" %in% names(x$comparison)) {
    ranked_comp <- x$comparison |> dplyr::arrange(.data$rank)
    high_k_rows <- which(ranked_comp$n_high_k > 0)
    if (length(high_k_rows) > 0) {
      gt_tbl <- gt_tbl |>
        gt::tab_style(
          style = gt::cell_text(color = "#CC0000"),
          locations = gt::cells_body(columns = "pareto_k_fmt",
                                     rows = high_k_rows)
        )
    }
  }

  gt_tbl
}
