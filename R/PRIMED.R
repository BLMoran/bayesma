#' PRIMED: Preliminary Investigation of Meta-analytic Databases
#'
#' Runs the four-step exploratory/preliminary data analysis workflow for
#' meta-analysis of dependent effect sizes as described in Pustejovsky,
#' Zhang, & Tipton (2026).
#'
#' @param data A data frame containing the meta-analytic database.
#' @param es_col Character. Name of the column containing effect size estimates.
#' @param se_col Character. Name of the column containing standard errors.
#' @param study_col Character. Name of the column identifying studies.
#' @param sample_col Character. Name of the column identifying samples within
#'   studies. If NULL (default), assumes one sample per study.
#' @param n_col Character. Name of the column containing total sample sizes.
#'   If NULL, sample size plots are skipped.
#' @param moderators Character vector of column names to examine as potential
#'   moderators. If NULL (default), no moderator analysis is performed.
#' @param es_type Character. Type of effect size: "SMD" (standardized mean
#'   difference) or "correlation". Affects how scaled SEs and weights are
#'   computed. Default is "SMD".
#' @param df_col Character. Name of the column containing degrees of freedom
#'   (used for scaled SE calculation when es_type = "SMD"). If NULL, scaled
#'   SEs are not computed.
#' @param rho_values Numeric vector. Assumed within-sample correlations for
#'   ISC weight calculations. Default is c(0.1, 0.3, 0.5, 0.7, 0.9).
#' @param fence_multiplier Numeric. Multiplier of the IQR for outlier fences
#'   in effect size density plots. Default is 3 (following Tukey's conventions
#'   as described in the paper).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{summary}{A list of summary statistics about the database.}
#'     \item{plots}{A named list of ggplot objects for each workflow step.}
#'     \item{tables}{A named list of summary tables (tibbles).}
#'   }
#'
#' @details
#' The PRIMED workflow proceeds in four steps:
#' \enumerate{
#'   \item \strong{Data structure}: Counts observations at each level and
#'     describes the dependence structure (effects per sample, samples per
#'     study, sample size distributions).
#'   \item \strong{Moderators}: Examines marginal distributions, missingness,
#'     and hierarchical (within- vs between-sample) structure of covariates.
#'   \item \strong{Standard errors & weights}: Inspects SE distributions
#'     within samples, computes scaled SEs (for SMDs), and calculates ISC
#'     weights under varying assumed correlations.
#'   \item \strong{Effect size distribution}: Visualises marginal and
#'     sample-level densities with outlier fences, and produces a hierarchical
#'     forest plot of dependent effect sizes.
#' }
#'
#' @examples
#' \dontrun{
#' results <- primed(
#'   data = my_meta_data,
#'   es_col = "g",
#'   se_col = "se",
#'   study_col = "study",
#'   sample_col = "sample_id",
#'   n_col = "n_total",
#'   moderators = c("intervention_type", "mean_age", "pct_female"),
#'   es_type = "SMD",
#'   df_col = "df"
#' )
#'
#' # View all step-1 plots
#' results$plots$step1_es_per_sample
#' results$plots$step1_samples_per_study
#'
#' # Access summary statistics
#' results$summary
#' }
#'
#' @import dplyr tidyr ggplot2 purrr
#' @importFrom tibble tibble
#' @importFrom stats median quantile sd
#' @export
primed <- function(data,
                   es_col,
                   se_col,
                   study_col,
                   sample_col = NULL,
                   n_col = NULL,
                   moderators = NULL,
                   es_type = c("SMD", "correlation"),
                   df_col = NULL,
                   rho_values = c(0.1, 0.3, 0.5, 0.7, 0.9),
                   fence_multiplier = 3) {

  es_type <- rlang::arg_match(es_type)

  # --- Input validation ---
  required_cols <- c(es_col, se_col, study_col)
  optional_cols <- c(sample_col, n_col, df_col)
  all_cols <- c(required_cols, optional_cols[!is.null(optional_cols)])
  missing_cols <- setdiff(all_cols, names(data))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Column(s) not found in {.arg data}: {.val {missing_cols}}.")
  }

  # Create internal working copy with standardised names
  df <- data |>
    dplyr::mutate(
      .study = .data[[study_col]],
      .es = .data[[es_col]],
      .se = .data[[se_col]]
    )

  if (!is.null(sample_col)) {
    df <- dplyr::mutate(df, .sample = .data[[sample_col]])
  } else {
    df <- dplyr::mutate(df, .sample = .data[[study_col]])
  }

  if (!is.null(n_col)) {
    df <- dplyr::mutate(df, .n = .data[[n_col]])
  }

  if (!is.null(df_col)) {
    df <- dplyr::mutate(df, .df = .data[[df_col]])
  }

  # Unique identifier for sample nested in study
  df <- df |>
    dplyr::mutate(.sample_id = interaction(.study, .sample, drop = TRUE))

  # Initialise output containers
  plots <- list()
  tables <- list()
  summary_stats <- list()

  # =========================================================================

  # STEP 1: Describe the Amount of Data and Dependence Structure
  # =========================================================================

  step1 <- .primed_step1(df, sample_col, n_col)
  plots <- c(plots, step1$plots)
  tables <- c(tables, step1$tables)
  summary_stats <- c(summary_stats, step1$summary)

  # =========================================================================
  # STEP 2: Explore Study Characteristics and Potential Moderators
  # =========================================================================

  if (!is.null(moderators)) {
    step2 <- .primed_step2(df, moderators, study_col)
    plots <- c(plots, step2$plots)
    tables <- c(tables, step2$tables)
  }

  # =========================================================================
  # STEP 3: Inspect Standard Errors and Other Auxiliary Data
  # =========================================================================

  step3 <- .primed_step3(df, es_type, df_col, rho_values)
  plots <- c(plots, step3$plots)
  tables <- c(tables, step3$tables)

  # =========================================================================
  # STEP 4: Visualize the Distribution of Effect Size Estimates
  # =========================================================================

  step4 <- .primed_step4(df, fence_multiplier)
  plots <- c(plots, step4$plots)
  tables <- c(tables, step4$tables)

  structure(
    list(
      summary = summary_stats,
      plots = plots,
      tables = tables
    ),
    class = "primed"
  )
}


# ---------------------------------------------------------------------------
# Print method
# ---------------------------------------------------------------------------

#' @export
print.primed <- function(x, ...) {
  cat("PRIMED Workflow Results\n")
  cat("======================\n\n")

  cat("Summary:\n")
  cat("  Studies (J):", x$summary$n_studies, "\n")
  cat("  Samples (M):", x$summary$n_samples, "\n")
  cat("  Effect sizes (K):", x$summary$n_effects, "\n")
  cat("  Effects per sample: median =", x$summary$median_es_per_sample,
      ", range =", x$summary$range_es_per_sample[1], "-",
      x$summary$range_es_per_sample[2], "\n")
  if (!is.null(x$summary$median_samples_per_study)) {
    cat("  Samples per study: median =", x$summary$median_samples_per_study,
        ", range =", x$summary$range_samples_per_study[1], "-",
        x$summary$range_samples_per_study[2], "\n")
  }
  cat("\nAvailable plots (", length(x$plots), "):\n")
  purrr::walk(names(x$plots), ~ cat("  -", .x, "\n"))
  cat("\nAvailable tables (", length(x$tables), "):\n")
  purrr::walk(names(x$tables), ~ cat("  -", .x, "\n"))
  invisible(x)
}


# ===========================================================================
# STEP 1 — Data structure & dependence
# ===========================================================================

.primed_step1 <- function(df, sample_col, n_col) {

  plots <- list()
  tables <- list()
  summary_out <- list()

  # --- Effects per sample ---
  es_per_sample <- df |>
    dplyr::count(.sample_id, name = "k") |>
    dplyr::pull(k)

  # --- Samples per study ---
  samples_per_study <- df |>
    dplyr::distinct(.study, .sample_id) |>
    dplyr::count(.study, name = "m") |>
    dplyr::pull(m)

  has_multi_sample <- !is.null(sample_col)

  summary_out$n_studies <- dplyr::n_distinct(df$.study)
  summary_out$n_samples <- dplyr::n_distinct(df$.sample_id)
  summary_out$n_effects <- nrow(df)
  summary_out$median_es_per_sample <- stats::median(es_per_sample)
  summary_out$range_es_per_sample <- range(es_per_sample)
  if (has_multi_sample) {
    summary_out$median_samples_per_study <- stats::median(samples_per_study)
    summary_out$range_samples_per_study <- range(samples_per_study)
  }

  # Effects-per-sample histogram
  es_per_sample_df <- tibble::tibble(k = es_per_sample)
  plots$step1_es_per_sample <- ggplot2::ggplot(es_per_sample_df, ggplot2::aes(x = k)) +
    ggplot2::geom_histogram(
      binwidth = 1, fill = "#2c5f7c", colour = "white", boundary = 0.5
    ) +
    ggplot2::labs(
      title = paste0(
        "Effect Size Estimates per Sample (M = ",
        summary_out$n_samples, " samples)"
      ),
      x = "Effect Size Estimates per Sample",
      y = "count"
    ) +
    ggplot2::theme_minimal()

  # Samples-per-study histogram (only if sample != study)
  if (has_multi_sample) {
    sps_df <- tibble::tibble(m = samples_per_study)
    plots$step1_samples_per_study <- ggplot2::ggplot(sps_df, ggplot2::aes(x = m)) +
      ggplot2::geom_histogram(
        binwidth = 1, fill = "#c4973b", colour = "white", boundary = 0.5
      ) +
      ggplot2::labs(
        title = paste0(
          "Samples per Study (J = ",
          summary_out$n_studies, " studies)"
        ),
        x = "Samples per Study",
        y = "count"
      ) +
      ggplot2::theme_minimal()
  }

  # Data structure summary table
  tables$step1_structure <- tibble::tibble(
    Level = c("Studies (J)", "Samples (M)", "Effect sizes (K)"),
    Count = c(summary_out$n_studies, summary_out$n_samples, summary_out$n_effects)
  )

  # --- Sample size distribution ---
  if (!is.null(n_col)) {
    sample_sizes <- df |>
      dplyr::distinct(.sample_id, .keep_all = TRUE)

    plots$step1_sample_sizes <- ggplot2::ggplot(
      sample_sizes,
      ggplot2::aes(x = .n)
    ) +
      ggplot2::geom_density(fill = "#2c5f7c", alpha = 0.4) +
      ggplot2::geom_rug(sides = "b") +
      ggplot2::labs(
        title = "Primary Study Sample Sizes",
        x = "Total Sample Size",
        y = NULL
      ) +
      ggplot2::theme_minimal()

    tables$step1_sample_sizes <- sample_sizes |>
      dplyr::summarise(
        Mean = mean(.n, na.rm = TRUE),
        SD = stats::sd(.n, na.rm = TRUE),
        Min = min(.n, na.rm = TRUE),
        Q1 = stats::quantile(.n, 0.25, na.rm = TRUE),
        Median = stats::median(.n, na.rm = TRUE),
        Q3 = stats::quantile(.n, 0.75, na.rm = TRUE),
        Max = max(.n, na.rm = TRUE)
      )
  }

  list(plots = plots, tables = tables, summary = summary_out)
}


# ===========================================================================
# STEP 2 — Moderators
# ===========================================================================

.primed_step2 <- function(df, moderators, study_col) {

  plots <- list()
  tables <- list()

  present_mods <- intersect(moderators, names(df))
  if (length(present_mods) == 0) {
    message("No moderator columns found in data. Skipping Step 2.")
    return(list(plots = plots, tables = tables))
  }

  # Classify moderators
  mod_classes <- purrr::map_chr(present_mods, ~ {
    vals <- df[[.x]]
    if (is.factor(vals) || is.character(vals) || is.logical(vals)) {
      "categorical"
    } else {
      "continuous"
    }
  }) |>
    rlang::set_names(present_mods)

  cat_mods <- names(mod_classes[mod_classes == "categorical"])
  cont_mods <- names(mod_classes[mod_classes == "continuous"])

  # --- 2.1 Marginal distributions ---

  # Categorical moderators: paired bar chart (studies & effects)
  plots$step2_categorical <- purrr::map(cat_mods, ~ {
    mod_name <- .x
    # Count studies and effects per level
    counts <- df |>
      dplyr::mutate(.mod_val = as.character(.data[[mod_name]])) |>
      dplyr::mutate(.mod_val = dplyr::if_else(
        is.na(.mod_val), "Not Reported", .mod_val
      )) |>
      dplyr::summarise(
        studies = dplyr::n_distinct(.study),
        effects = dplyr::n(),
        .by = .mod_val
      ) |>
      tidyr::pivot_longer(
        cols = c(studies, effects),
        names_to = "type",
        values_to = "count"
      ) |>
      dplyr::mutate(
        .mod_val = forcats::fct_reorder(.mod_val, count, .fun = max)
      )

    ggplot2::ggplot(counts, ggplot2::aes(x = count, y = .mod_val, fill = type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_text(
        ggplot2::aes(label = count),
        position = ggplot2::position_dodge(width = 0.9),
        hjust = -0.2, size = 3
      ) +
      ggplot2::scale_fill_manual(values = c(studies = "#2c5f7c", effects = "#c4973b")) +
      ggplot2::labs(title = mod_name, x = NULL, y = NULL, fill = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
  }) |>
    rlang::set_names(paste0("step2_cat_", cat_mods))

  # Un-nest the named list so plots sit at top level
  plots <- c(plots[!names(plots) %in% "step2_categorical"], plots$step2_categorical)

  # Continuous moderators: density + rug

  plots_cont <- purrr::map(cont_mods, ~ {
    mod_name <- .x
    vals <- df[[mod_name]]
    pct_miss <- round(100 * mean(is.na(vals)), 1)
    n_obs <- sum(!is.na(vals))

    ggplot2::ggplot(df, ggplot2::aes(x = .data[[mod_name]])) +
      ggplot2::geom_density(fill = "#2c5f7c", alpha = 0.4) +
      ggplot2::geom_rug(sides = "b") +
      ggplot2::labs(
        title = paste0(mod_name, " (n = ", n_obs, ", ", pct_miss, "% missing)"),
        x = mod_name, y = NULL
      ) +
      ggplot2::theme_minimal()
  }) |>
    rlang::set_names(paste0("step2_cont_", cont_mods))

  plots <- c(plots, plots_cont)

  # --- Missingness summary ---
  miss_tbl <- tibble::tibble(
    moderator = present_mods,
    n_missing = purrr::map_int(present_mods, ~ sum(is.na(df[[.x]]))),
    pct_missing = round(100 * n_missing / nrow(df), 1)
  )
  tables$step2_missingness <- miss_tbl

  # --- 2.2 Hierarchical structure for effect-level categorical mods ---
  # Co-occurrence tables for categorical moderators
  if (length(cat_mods) > 0) {
    tables$step2_cooccurrence <- purrr::map(cat_mods, ~ {
      .primed_cooccurrence(df, .x)
    }) |>
      rlang::set_names(paste0("cooccurrence_", cat_mods))
  }

  # Group-mean centering for continuous moderators
  if (length(cont_mods) > 0) {
    tables$step2_decomposition <- purrr::map(cont_mods, ~ {
      .primed_decompose_continuous(df, .x)
    }) |>
      rlang::set_names(paste0("decomposition_", cont_mods))
  }

  list(plots = plots, tables = tables)
}


# ---------------------------------------------------------------------------
# Co-occurrence table for a categorical moderator
# ---------------------------------------------------------------------------
.primed_cooccurrence <- function(df, mod_col) {

  df_mod <- df |>
    dplyr::mutate(.mod = as.character(.data[[mod_col]])) |>
    dplyr::filter(!is.na(.mod))

  levels_used <- sort(unique(df_mod$.mod))

  # For each study, which levels are present?
  study_levels <- df_mod |>
    dplyr::summarise(levels = list(unique(.mod)), .by = .study)

  # Build co-occurrence matrix
  co_mat <- matrix(
    "", nrow = length(levels_used), ncol = length(levels_used),
    dimnames = list(levels_used, levels_used)
  )

  purrr::walk(seq_len(nrow(study_levels)), function(i) {
    lvls <- study_levels$levels[[i]]
    combos <- expand.grid(a = lvls, b = lvls, stringsAsFactors = FALSE) |>
      dplyr::filter(a <= b)  # upper triangle + diagonal
    # We just mark presence; actual counting below
  })

  # Simpler approach: count studies and effects per pair
  cross_list <- purrr::map(levels_used, function(lev_a) {
    purrr::map_chr(levels_used, function(lev_b) {
      # Studies that have both levels
      studies_a <- df_mod |>
        dplyr::filter(.mod == lev_a) |>
        dplyr::pull(.study) |>
        unique()
      studies_b <- df_mod |>
        dplyr::filter(.mod == lev_b) |>
        dplyr::pull(.study) |>
        unique()
      shared <- intersect(studies_a, studies_b)
      n_studies <- length(shared)
      if (n_studies == 0) return("")

      if (lev_a == lev_b) {
        # Diagonal: count total effects for this level in shared studies
        n_effects <- df_mod |>
          dplyr::filter(.study %in% shared, .mod == lev_a) |>
          nrow()
      } else {
        # Off-diagonal: count effects for lev_a in shared studies
        n_effects <- df_mod |>
          dplyr::filter(.study %in% shared, .mod == lev_a) |>
          nrow()
      }
      paste0(n_studies, " (", n_effects, ")")
    })
  })

  co_tbl <- purrr::list_rbind(cross_list) |>
    as.data.frame(stringsAsFactors = FALSE)
  names(co_tbl) <- levels_used
  co_tbl <- tibble::tibble(level = levels_used) |>
    dplyr::bind_cols(tibble::as_tibble(co_tbl))

  co_tbl
}


# ---------------------------------------------------------------------------
# Group-mean decomposition for a continuous moderator
# ---------------------------------------------------------------------------
.primed_decompose_continuous <- function(df, mod_col) {

  df_mod <- df |>
    dplyr::filter(!is.na(.data[[mod_col]])) |>
    dplyr::mutate(
      .sample_mean = mean(.data[[mod_col]], na.rm = TRUE),
      .centered = .data[[mod_col]] - .sample_mean,
      .by = .sample_id
    )

  .five_num_summary <- function(x, label) {
    tibble::tibble(
      Component = label,
      `% Missing` = round(100 * mean(is.na(df[[mod_col]])), 2),
      Mean = mean(x, na.rm = TRUE),
      SD = stats::sd(x, na.rm = TRUE),
      Min = min(x, na.rm = TRUE),
      Q1 = stats::quantile(x, 0.25, na.rm = TRUE),
      Median = stats::median(x, na.rm = TRUE),
      Q3 = stats::quantile(x, 0.75, na.rm = TRUE),
      Max = max(x, na.rm = TRUE)
    )
  }

  dplyr::bind_rows(
    .five_num_summary(df_mod[[mod_col]], "Effect-Level"),
    .five_num_summary(df_mod$.sample_mean, "Sample-Level Average"),
    .five_num_summary(df_mod$.centered, "Sample-Mean-Centered")
  )
}


# ===========================================================================
# STEP 3 — Standard errors & weights
# ===========================================================================

.primed_step3 <- function(df, es_type, df_col, rho_values) {

  plots <- list()
  tables <- list()

  # --- SE distribution by sample ---
  # Sort samples by median SE for readability
  sample_order <- df |>
    dplyr::summarise(med_se = stats::median(.se, na.rm = TRUE), .by = .sample_id) |>
    dplyr::arrange(med_se) |>
    dplyr::pull(.sample_id)

  df_se <- df |>
    dplyr::mutate(
      .sample_id = factor(.sample_id, levels = sample_order)
    )

  # Only plot if reasonable number of samples (< 80)
  if (dplyr::n_distinct(df$.sample_id) <= 80) {
    plots$step3_se_by_sample <- ggplot2::ggplot(
      df_se,
      ggplot2::aes(x = .se, y = .sample_id, colour = .sample_id)
    ) +
      ggplot2::geom_point(alpha = 0.7, show.legend = FALSE) +
      ggplot2::labs(
        title = "Standard Errors by Sample",
        x = "Standard Error",
        y = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6)
      )
  }

  # --- Scaled SE (for SMDs) ---
  if (es_type == "SMD" && !is.null(df_col)) {
    df <- df |>
      dplyr::mutate(
        .scaled_se = sqrt(pmax(.se^2 - (.es^2) / (2 * .df), 0))
      )

    scaled_order <- df |>
      dplyr::summarise(med = stats::median(.scaled_se, na.rm = TRUE), .by = .sample_id) |>
      dplyr::arrange(med) |>
      dplyr::pull(.sample_id)

    df_scaled <- df |>
      dplyr::mutate(.sample_id = factor(.sample_id, levels = scaled_order))

    if (dplyr::n_distinct(df$.sample_id) <= 80) {
      plots$step3_scaled_se <- ggplot2::ggplot(
        df_scaled,
        ggplot2::aes(x = .scaled_se, y = .sample_id, colour = .sample_id)
      ) +
        ggplot2::geom_point(alpha = 0.7, show.legend = FALSE) +
        ggplot2::labs(
          title = "Scaled Standard Errors by Sample",
          x = "Scaled SE",
          y = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 6)
        )
    }
  }

  # --- ISC Weights ---
  sample_stats <- df |>
    dplyr::summarise(
      k = dplyr::n(),
      sigma = stats::median(.se, na.rm = TRUE),
      .by = c(.study, .sample_id)
    )

  isc_weights <- purrr::map(rho_values, function(rho) {
    sample_stats |>
      dplyr::mutate(
        w = k / (sigma^2 * ((k - 1) * rho + 1)),
        w_norm = w / sum(w),
        rho = rho
      )
  }) |> purrr::list_rbind()

  # Sort by average weight
  weight_order <- isc_weights |>
    dplyr::summarise(mean_w = mean(w_norm), .by = .sample_id) |>
    dplyr::arrange(mean_w) |>
    dplyr::pull(.sample_id)

  isc_weights <- isc_weights |>
    dplyr::mutate(
      .sample_id = factor(.sample_id, levels = weight_order),
      rho = factor(rho)
    )

  if (dplyr::n_distinct(df$.sample_id) <= 80) {
    plots$step3_isc_weights <- ggplot2::ggplot(
      isc_weights,
      ggplot2::aes(x = w_norm, y = .sample_id, colour = rho, group = rho)
    ) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_line(alpha = 0.4) +
      ggplot2::labs(
        title = "Inverse Sampling Covariance (ISC) Weights",
        x = "Normalized Weight",
        y = NULL,
        colour = "Assumed\nCorrelation"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6)
      )
  }

  tables$step3_isc_weights <- isc_weights |>
    dplyr::select(.study, .sample_id, k, sigma, rho, w_norm) |>
    tidyr::pivot_wider(names_from = rho, values_from = w_norm, names_prefix = "rho_")

  list(plots = plots, tables = tables)
}


# ===========================================================================
# STEP 4 — Effect size distribution
# ===========================================================================

.primed_step4 <- function(df, fence_multiplier) {

  plots <- list()
  tables <- list()

  # --- 4.1 Marginal density with outlier fences ---
  q1 <- stats::quantile(df$.es, 0.25, na.rm = TRUE)
  q3 <- stats::quantile(df$.es, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_fence <- q1 - fence_multiplier * iqr
  upper_fence <- q3 + fence_multiplier * iqr

  n_effects <- nrow(df)

  plots$step4_es_density <- ggplot2::ggplot(df, ggplot2::aes(x = .es)) +
    ggplot2::geom_density(fill = "#c4973b", alpha = 0.5) +
    ggplot2::geom_rug(sides = "b", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = c(q1, q3), linetype = "solid", alpha = 0.6) +
    ggplot2::geom_vline(xintercept = c(lower_fence, upper_fence), linetype = "dashed", alpha = 0.6) +
    ggplot2::labs(
      title = paste0("Effect-Level Distribution (K = ", n_effects, " effects)"),
      x = "Effect Size Estimate",
      y = NULL
    ) +
    ggplot2::theme_minimal()

  # Sample-level average density
  sample_avgs <- df |>
    dplyr::summarise(.es_avg = mean(.es, na.rm = TRUE), .by = c(.study, .sample_id))

  n_samples <- nrow(sample_avgs)
  q1_s <- stats::quantile(sample_avgs$.es_avg, 0.25, na.rm = TRUE)
  q3_s <- stats::quantile(sample_avgs$.es_avg, 0.75, na.rm = TRUE)
  iqr_s <- q3_s - q1_s

  plots$step4_sample_es_density <- ggplot2::ggplot(
    sample_avgs, ggplot2::aes(x = .es_avg)
  ) +
    ggplot2::geom_density(fill = "#2c5f7c", alpha = 0.5) +
    ggplot2::geom_rug(sides = "b", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = c(q1_s, q3_s), linetype = "solid", alpha = 0.6) +
    ggplot2::geom_vline(
      xintercept = c(q1_s - fence_multiplier * iqr_s, q3_s + fence_multiplier * iqr_s),
      linetype = "dashed", alpha = 0.6
    ) +
    ggplot2::labs(
      title = paste0("Sample-Level Average Distribution (M = ", n_samples, " samples)"),
      x = "Sample-Level Average Effect Size Estimate",
      y = NULL
    ) +
    ggplot2::theme_minimal()

  # Outlier table
  outliers_effect <- df |>
    dplyr::filter(.es < lower_fence | .es > upper_fence) |>
    dplyr::select(.study, .sample_id, .es, .se)

  outliers_sample <- sample_avgs |>
    dplyr::filter(
      .es_avg < (q1_s - fence_multiplier * iqr_s) |
        .es_avg > (q3_s + fence_multiplier * iqr_s)
    )

  tables$step4_outliers_effect <- outliers_effect
  tables$step4_outliers_sample <- outliers_sample

  # --- 4.2 Hierarchical forest plot ---
  # Compute 95% CI for each effect
  df_forest <- df |>
    dplyr::mutate(
      .ci_lo = .es - stats::qnorm(0.975) * .se,
      .ci_hi = .es + stats::qnorm(0.975) * .se
    )

  # Sort samples by average ES
  sample_order <- df_forest |>
    dplyr::summarise(.avg = mean(.es, na.rm = TRUE), .by = c(.study, .sample_id)) |>
    dplyr::arrange(dplyr::desc(.avg)) |>
    dplyr::mutate(.row = dplyr::row_number())

  df_forest <- df_forest |>
    dplyr::left_join(
      dplyr::select(sample_order, .sample_id, .row),
      by = ".sample_id"
    )

  # Jitter within row so overlapping points are visible
  set.seed(42)
  df_forest <- df_forest |>
    dplyr::mutate(.row_j = .row + stats::runif(dplyr::n(), -0.2, 0.2))

  # Build labels: study (sample)
  sample_labels <- sample_order |>
    dplyr::mutate(
      .label = as.character(.sample_id)
    )

  if (dplyr::n_distinct(df$.sample_id) <= 80) {
    plots$step4_forest <- ggplot2::ggplot(df_forest) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey50", linetype = "solid") +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .ci_lo, xend = .ci_hi,
          y = .row_j, yend = .row_j,
          colour = factor(.sample_id)
        ),
        alpha = 0.3, show.legend = FALSE
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = .es, y = .row_j, colour = factor(.sample_id)),
        size = 2, alpha = 0.7, show.legend = FALSE
      ) +
      ggplot2::scale_y_continuous(
        breaks = sample_labels$.row,
        labels = sample_labels$.label,
        trans = "reverse"
      ) +
      ggplot2::labs(
        title = "Hierarchical Forest Plot of Dependent Effect Sizes",
        x = "Effect Size Estimate",
        y = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6)
      )
  }

  # --- 4.3 Funnel plot ---
  plots$step4_funnel <- ggplot2::ggplot(df, ggplot2::aes(x = .es, y = .se)) +
    ggplot2::geom_point(
      ggplot2::aes(colour = factor(.sample_id)),
      alpha = 0.5, show.legend = FALSE
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      title = "Funnel Plot",
      x = "Effect Size Estimate",
      y = "Standard Error"
    ) +
    ggplot2::theme_minimal()

  # Sample-level funnel
  sample_funnel <- df |>
    dplyr::summarise(
      .es_avg = mean(.es, na.rm = TRUE),
      .se_avg = mean(.se, na.rm = TRUE),
      .by = c(.study, .sample_id)
    )

  plots$step4_funnel_sample <- ggplot2::ggplot(
    sample_funnel, ggplot2::aes(x = .es_avg, y = .se_avg)
  ) +
    ggplot2::geom_point(alpha = 0.6, colour = "#2c5f7c") +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      title = "Sample-Level Funnel Plot",
      x = "Sample-Level Average Effect Size",
      y = "Average Standard Error"
    ) +
    ggplot2::theme_minimal()

  list(plots = plots, tables = tables)
}
