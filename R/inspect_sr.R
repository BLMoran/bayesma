# ============================================================================
# INSPECT-SR: Main Assessment Function
# ============================================================================
# Single entry point: inspect_sr(). Takes a one-row-per-study data frame,
# runs automated Domain 4 checks, reads manual items, and returns an
# inspect_sr (frequentist) or bayes_inspect_sr (Bayesian) object.
# ============================================================================


#' Run the INSPECT-SR Trustworthiness Assessment
#'
#' Takes a data frame with one row per study and runs the automated Domain 4
#' checks (Carlisle's test, participant-number consistency, GRIM, p-value
#' verification). Manual items (D1, D2, D3, and the non-automated D4 items)
#' are read straight from the input. Domain-level and overall judgements are
#' derived per INSPECT-SR guidance (overall = most severe domain).
#'
#' @param data A data frame or tibble with one row per study. See
#'   **Expected columns** in the package vignette, or the bundled
#'   [inspect_sr_example] dataset for the exact layout.
#' @param studyvar Unquoted column name identifying the study (tidyeval).
#'   Defaults to `study`.
#' @param bayes Logical. If `FALSE` (default) produces frequentist pass/fail
#'   judgements. If `TRUE` produces Bayes factors and a posterior probability
#'   of trustworthiness.
#' @param prior_prob_trustworthy Numeric in (0, 1). Prior probability that
#'   each study is trustworthy, used only when `bayes = TRUE` (default 0.90).
#' @param pvalue_tolerance Numeric. Tolerance for the frequentist p-value
#'   check (default 0.01).
#' @param carlisle_method `"fisher"` (default) or `"ks"` — see
#'   [carlisle_test()].
#' @param verbose Logical. Print a summary to the console (default `TRUE`).
#'
#' @return
#' If `bayes = FALSE`: an object of class `inspect_sr` (a data frame with
#' columns `Study`, `D1`, `D2`, `D3`, `D4`, `Overall`), with per-study
#' details in `attr(x, "details")`.
#'
#' If `bayes = TRUE`: an object of class `bayes_inspect_sr` (a data frame
#' with columns `Study`, `Prior`, `Posterior`, `Combined_BF`,
#' `Interpretation`), with individual Bayes factors in `attr(x, "details")`.
#'
#' @seealso [inspect_sr_table()] for a per-check gt table;
#'   [inspect_plot()] for the traffic-light visualisation;
#'   [filter_trustworthy()] for filtering a meta-analysis dataset.
#'
#' @examples
#' \dontrun{
#' data(inspect_sr_example)
#'
#' # Frequentist
#' res <- inspect_sr(inspect_sr_example, studyvar = study)
#'
#' # Bayesian
#' res_bayes <- inspect_sr(inspect_sr_example, studyvar = study, bayes = TRUE)
#' }
#'
#' @export
inspect_sr <- function(data,
                       studyvar = study,
                       bayes = FALSE,
                       prior_prob_trustworthy = 0.90,
                       pvalue_tolerance = 0.01,
                       carlisle_method = "fisher",
                       verbose = TRUE) {

  studyvar_quo <- rlang::enquo(studyvar)
  studyvar_name <- rlang::as_name(studyvar_quo)

  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame or tibble.")
  }
  if (!studyvar_name %in% names(data)) {
    rlang::abort(paste0("Column `", studyvar_name, "` not found in `data`."))
  }

  study_ids <- data[[studyvar_name]]
  if (any(is.na(study_ids)) || any(duplicated(study_ids))) {
    rlang::abort(paste0("`", studyvar_name, "` must be unique and non-missing."))
  }

  # Run automated checks for each row
  details <- lapply(seq_len(nrow(data)), function(i) {
    run_automated_checks(
      row = data[i, , drop = FALSE],
      studyvar_name = studyvar_name,
      bayes = bayes,
      pvalue_tolerance = pvalue_tolerance,
      carlisle_method = carlisle_method
    )
  })

  if (isTRUE(bayes)) {
    out <- build_bayes_summary(data, details, studyvar_name,
                               prior_prob_trustworthy)
    if (verbose) print_bayes_inspect_sr_summary(out)
  } else {
    out <- build_freq_summary(data, details, studyvar_name)
    if (verbose) print_inspect_sr_summary(out)
  }

  invisible(out)
}


#' @export
print.inspect_sr <- function(x, ...) {
  print_inspect_sr_summary(x)
  invisible(x)
}

#' @export
print.bayes_inspect_sr <- function(x, ...) {
  print_bayes_inspect_sr_summary(x)
  invisible(x)
}


# ============================================================================
# Per-study automated checks
# ============================================================================

#' @noRd
run_automated_checks <- function(row, studyvar_name, bayes,
                                 pvalue_tolerance, carlisle_method) {

  study <- row[[studyvar_name]]
  get_n <- function(x) {
    if (!x %in% names(row)) return(NA_integer_)
    v <- row[[x]]
    if (is.null(v) || length(v) == 0 || is.na(v)) NA_integer_ else as.integer(v)
  }

  baseline   <- first_list_col(row, "baseline")
  statistics <- first_list_col(row, "statistics")

  n_int   <- get_n("n_randomised_int")
  n_ctrl  <- get_n("n_randomised_ctrl")
  n_total <- get_n("n_randomised_total")
  n_ai    <- get_n("n_analysed_int")
  n_ac    <- get_n("n_analysed_ctrl")
  n_li    <- get_n("n_lost_int")
  n_lc    <- get_n("n_lost_ctrl")

  list(
    study         = study,
    grim          = run_grim_block(baseline, n_int, n_ctrl, bayes),
    pvalue        = run_pvalue_block(statistics, bayes, pvalue_tolerance),
    carlisle      = run_carlisle_block(baseline, bayes, carlisle_method),
    n_consistency = run_nconsistency_block(n_int, n_ctrl, n_ai, n_ac,
                                           n_li, n_lc, n_total, bayes),
    outcome_stats = run_outcome_stats_block(row)
  )
}


#' @noRd
first_list_col <- function(row, name) {
  if (!name %in% names(row)) return(NULL)
  val <- row[[name]]
  if (is.list(val) && length(val) == 1) val <- val[[1]]
  if (is.null(val) || (is.data.frame(val) && nrow(val) == 0)) return(NULL)
  val
}


#' @noRd
run_grim_block <- function(baseline, n_int, n_ctrl, bayes) {
  if (is.null(baseline) || !"integer_scale" %in% names(baseline)) {
    return(list(results = list(), n_total = 0L, n_failed = 0L,
                judgement = NA_character_, bf = 1))
  }

  keep <- !is.na(baseline$integer_scale) & baseline$integer_scale == TRUE
  integer_vars <- baseline[keep, , drop = FALSE]
  if (nrow(integer_vars) == 0) {
    return(list(results = list(), n_total = 0L, n_failed = 0L,
                judgement = NA_character_, bf = 1))
  }

  results <- list()
  bfs <- numeric()

  for (i in seq_len(nrow(integer_vars))) {
    r <- integer_vars[i, ]
    observed_max <- max(c(r$mean_int, r$mean_ctrl), na.rm = TRUE)
    max_items <- max(10, ceiling(observed_max * 2))

    for (arm in c("int", "ctrl")) {
      m <- r[[paste0("mean_", arm)]]
      n_arm <- if (arm == "int") n_int else n_ctrl
      if (is.na(m) || is.na(n_arm)) next

      freq_res <- grim_test(m, n_arm)
      bf <- if (isTRUE(bayes)) {
        bayes_grim_test(m, n_arm, max_items = max_items)$bf_inconsistent
      } else NA_real_

      results[[length(results) + 1]] <- list(
        variable = r$variable,
        group = if (arm == "int") "Intervention" else "Control",
        mean_value = m, n = n_arm,
        consistent = freq_res$consistent, bf = bf
      )
      if (isTRUE(bayes)) bfs <- c(bfs, bf)
    }
  }

  n_total <- length(results)
  n_failed <- sum(!vapply(results, `[[`, logical(1), "consistent"))

  list(
    results = results, n_total = n_total, n_failed = n_failed,
    judgement = if (n_total == 0) NA_character_
    else if (n_failed > 0) "Serious concerns"
    else "No concerns",
    bf = if (isTRUE(bayes) && length(bfs) > 0) max(bfs) else 1
  )
}


#' @noRd
run_pvalue_block <- function(statistics, bayes, tolerance) {
  if (is.null(statistics) || nrow(statistics) == 0) {
    return(list(results = list(), n_total = 0L, n_failed = 0L,
                n_major = 0L, judgement = NA_character_, bf = 1))
  }

  results <- list()
  bfs <- numeric()

  for (i in seq_len(nrow(statistics))) {
    r <- statistics[i, ]
    if (any(is.na(c(r$statistic, r$reported_p, r$test_type)))) next

    df_val <- if (!is.na(r$test_type) && r$test_type == "f" &&
                  !is.na(r$df) && !is.na(r$df2)) {
      c(r$df, r$df2)
    } else if (!is.na(r$df)) {
      r$df
    } else if (!is.na(r$test_type) && r$test_type == "z") {
      NULL
    } else next

    freq_res <- tryCatch(
      verify_pvalue(r$test_type, r$statistic, df = df_val,
                    reported_p = r$reported_p, tolerance = tolerance),
      error = function(e) NULL
    )
    if (is.null(freq_res)) next

    bf <- if (isTRUE(bayes)) {
      tryCatch(
        bayes_verify_pvalue(r$test_type, r$statistic, df = df_val,
                            reported_p = r$reported_p)$bf_inconsistent,
        error = function(e) NA_real_
      )
    } else NA_real_

    results[[length(results) + 1]] <- list(
      test_type = r$test_type, statistic = r$statistic,
      reported_p = r$reported_p,
      recalculated_p = freq_res$recalculated_p,
      difference = freq_res$difference,
      consistent = freq_res$consistent,
      context = if ("context" %in% names(r)) r$context else NA_character_,
      bf = bf
    )
    if (isTRUE(bayes) && !is.na(bf)) bfs <- c(bfs, bf)
  }

  n_total <- length(results)
  n_failed <- sum(!vapply(results, `[[`, logical(1), "consistent"))
  n_major <- sum(vapply(results, function(r) {
    !r$consistent && isTRUE(r$difference > 0.05)
  }, logical(1)))

  list(
    results = results, n_total = n_total, n_failed = n_failed,
    n_major = n_major,
    judgement = if (n_total == 0) NA_character_
    else if (n_major > 0) "Serious concerns"
    else if (n_failed > 0) "Some concerns"
    else "No concerns",
    bf = if (isTRUE(bayes) && length(bfs) > 0) prod(bfs) else 1
  )
}


#' @noRd
run_carlisle_block <- function(baseline, bayes, method) {
  if (is.null(baseline) || !"p_value" %in% names(baseline)) {
    return(list(result = NULL, judgement = NA_character_, bf = 1))
  }
  pvals <- baseline$p_value[!is.na(baseline$p_value)]
  if (length(pvals) < 2) {
    return(list(result = NULL, judgement = NA_character_, bf = 1))
  }

  freq_res <- tryCatch(carlisle_test(pvals, method = method),
                       error = function(e) NULL)
  if (is.null(freq_res)) {
    return(list(result = NULL, judgement = NA_character_, bf = 1))
  }

  judgement <- switch(freq_res$interpretation,
                      too_similar = "Serious concerns",
                      too_different = "Some concerns",
                      plausible = "No concerns",
                      NA_character_
  )

  bf <- if (isTRUE(bayes)) {
    res <- tryCatch(bayes_carlisle_test(pvals), error = function(e) NULL)
    if (is.null(res)) 1 else res$bf_nonuniform
  } else 1

  list(result = freq_res, judgement = judgement, bf = bf)
}


#' @noRd
run_nconsistency_block <- function(n_int, n_ctrl, n_ai, n_ac,
                                   n_li, n_lc, n_total, bayes) {
  if (is.na(n_int) || is.na(n_ctrl)) {
    return(list(result = NULL, judgement = NA_character_, bf = 1))
  }

  res <- check_n_consistency(
    n_randomised_int   = n_int,
    n_randomised_ctrl  = n_ctrl,
    n_analysed_int     = if (is.na(n_ai)) NULL else n_ai,
    n_analysed_ctrl    = if (is.na(n_ac)) NULL else n_ac,
    n_lost_int         = if (is.na(n_li)) NULL else n_li,
    n_lost_ctrl        = if (is.na(n_lc)) NULL else n_lc,
    n_randomised_total = if (is.na(n_total)) NULL else n_total
  )

  judgement <- if (res$n_checks == 0) NA_character_
  else if (res$n_failed > 1) "Serious concerns"
  else if (res$n_failed == 1) "Some concerns"
  else "No concerns"

  bf <- if (isTRUE(bayes)) {
    if (res$n_failed == 0) 1 else 20^(res$n_failed)
  } else 1

  list(result = res, judgement = judgement, bf = bf)
}


#' @noRd
run_outcome_stats_block <- function(row) {
  if (!"outcome_estimate" %in% names(row)) return(NULL)
  est <- row$outcome_estimate
  if (is.na(est)) return(NULL)

  check_statistics_consistency(
    estimate  = est,
    ci_lower  = if ("outcome_ci_lower" %in% names(row) &&
                    !is.na(row$outcome_ci_lower)) row$outcome_ci_lower else NULL,
    ci_upper  = if ("outcome_ci_upper" %in% names(row) &&
                    !is.na(row$outcome_ci_upper)) row$outcome_ci_upper else NULL,
    se        = if ("outcome_se" %in% names(row) &&
                    !is.na(row$outcome_se)) row$outcome_se else NULL,
    log_scale = if ("outcome_log_scale" %in% names(row) &&
                    !is.na(row$outcome_log_scale)) isTRUE(row$outcome_log_scale)
    else FALSE
  )
}


# ============================================================================
# Build summary objects
# ============================================================================

#' @noRd
inspect_sr_item_names <- function() {
  c(paste0("d1_", 1:3), paste0("d2_", 1:5),
    paste0("d3_", 1:2), paste0("d4_", 1:11))
}


#' @noRd
most_severe <- function(values) {
  severity <- c("No concerns" = 1, "Some concerns" = 2, "Serious concerns" = 3)
  values <- values[!is.na(values)]
  if (length(values) == 0) return(NA_character_)
  names(severity)[severity == max(severity[values], na.rm = TRUE)]
}


#' @noRd
build_freq_summary <- function(data, details, studyvar_name) {

  items <- inspect_sr_item_names()

  summary_df <- data.frame(
    Study = data[[studyvar_name]],
    D1 = NA_character_, D2 = NA_character_,
    D3 = NA_character_, D4 = NA_character_,
    Overall = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(data))) {
    d <- details[[i]]

    item_vals <- lapply(items, function(item) {
      if (item %in% c("d4_3", "d4_6", "d4_8", "d4_9")) {
        key <- switch(item,
                      d4_3 = "carlisle", d4_6 = "n_consistency",
                      d4_8 = "grim", d4_9 = "pvalue")
        d[[key]]$judgement
      } else if (item %in% names(data)) {
        val <- data[[item]][i]
        if (is.na(val)) NA_character_ else as.character(val)
      } else NA_character_
    })
    names(item_vals) <- items

    summary_df$D1[i] <- most_severe(unlist(item_vals[grep("^d1_", items)]))
    summary_df$D2[i] <- most_severe(unlist(item_vals[grep("^d2_", items)]))
    summary_df$D3[i] <- most_severe(unlist(item_vals[grep("^d3_", items)]))
    summary_df$D4[i] <- most_severe(unlist(item_vals[grep("^d4_", items)]))
    summary_df$Overall[i] <- most_severe(
      c(summary_df$D1[i], summary_df$D2[i],
        summary_df$D3[i], summary_df$D4[i])
    )
  }

  attr(summary_df, "details") <- details
  attr(summary_df, "data") <- data
  attr(summary_df, "studyvar") <- studyvar_name
  class(summary_df) <- c("inspect_sr", "data.frame")
  summary_df
}


#' @noRd
build_bayes_summary <- function(data, details, studyvar_name,
                                prior_prob_trustworthy) {

  prior_odds <- prior_prob_trustworthy / (1 - prior_prob_trustworthy)

  combined_bfs <- vapply(details, function(d) {
    bfs <- c(d$grim$bf, d$pvalue$bf, d$carlisle$bf, d$n_consistency$bf)
    bfs <- bfs[!is.na(bfs)]
    if (length(bfs) == 0) 1 else prod(bfs)
  }, numeric(1))

  posterior <- vapply(combined_bfs, function(bf) {
    if (is.infinite(bf)) return(0)
    post_odds <- prior_odds / bf
    max(0, min(1, post_odds / (1 + post_odds)))
  }, numeric(1))

  interpretation <- vapply(posterior, function(p) {
    if (p > 0.90) "Trustworthy (no concerns)"
    else if (p > 0.50) "Uncertain (some concerns)"
    else "Likely untrustworthy (serious concerns)"
  }, character(1))

  out <- data.frame(
    Study = data[[studyvar_name]],
    Prior = prior_prob_trustworthy,
    Posterior = posterior,
    Combined_BF = combined_bfs,
    Interpretation = interpretation,
    stringsAsFactors = FALSE
  )

  attr(out, "details") <- details
  attr(out, "data") <- data
  attr(out, "studyvar") <- studyvar_name
  class(out) <- c("bayes_inspect_sr", "data.frame")
  out
}


# ============================================================================
# Printing
# ============================================================================

#' @noRd
print_inspect_sr_summary <- function(x) {

  abbrev <- function(j) {
    if (is.na(j)) return("--")
    switch(j, "No concerns" = "OK", "Some concerns" = "SOME",
           "Serious concerns" = "SERIOUS", j)
  }

  cat("\nINSPECT-SR Trustworthiness Assessment\n")
  cat(strrep("=", 50), "\n\n")
  cat(sprintf("%-25s %-5s %-5s %-5s %-8s %-18s\n",
              "Study", "D1", "D2", "D3", "D4", "Overall"))
  cat(strrep("-", 70), "\n")

  for (i in seq_len(nrow(x))) {
    cat(sprintf("%-25s %-5s %-5s %-5s %-8s %-18s\n",
                substr(x$Study[i], 1, 25),
                abbrev(x$D1[i]), abbrev(x$D2[i]),
                abbrev(x$D3[i]), abbrev(x$D4[i]),
                x$Overall[i]))
  }

  cat("\nDomains: D1 post-publication, D2 conduct/governance,\n")
  cat("         D3 text/figures, D4 results (auto-filled for 4.3/4.6/4.8/4.9)\n")
  cat("OK = No concerns, SOME = Some concerns, SERIOUS = Serious concerns\n")
  cat("-- = Not assessed\n\n")
  cat("For a per-check table, call inspect_sr_table().\n\n")
}


#' @noRd
print_bayes_inspect_sr_summary <- function(x) {

  cat("\nBayesian INSPECT-SR Trustworthiness Assessment\n")
  cat(strrep("=", 70), "\n\n")
  cat(sprintf("%-25s %8s %10s %10s  %-30s\n",
              "Study", "Prior", "Posterior", "BF", "Interpretation"))
  cat(strrep("-", 90), "\n")

  for (i in seq_len(nrow(x))) {
    bf_fmt <- if (is.infinite(x$Combined_BF[i])) "Inf"
    else sprintf("%10.3f", x$Combined_BF[i])
    cat(sprintf("%-25s %8.3f %10.3f %10s  %-30s\n",
                substr(x$Study[i], 1, 25),
                x$Prior[i], x$Posterior[i], bf_fmt,
                x$Interpretation[i]))
  }

  cat("\nBF: combined Bayes factor against trustworthiness (>1 = evidence against)\n")
  cat("Posterior: P(trustworthy | data), prior =", unique(x$Prior)[1], "\n")
  cat(">0.90 no concerns, 0.50-0.90 some concerns, <0.50 serious\n\n")
}


# ============================================================================
# Filter helper
# ============================================================================

#' Filter Studies by INSPECT-SR Trustworthiness
#'
#' Filters a meta-analysis dataset based on INSPECT-SR results.
#'
#' @param data A data frame with one row per study.
#' @param inspect_results An object returned by [inspect_sr()].
#' @param studyvar Unquoted column name for the study identifier in `data`
#'   (tidyeval). Must match the Study column in `inspect_results`.
#' @param exclude `"serious"` (default) drops studies with "Serious concerns";
#'   `"some"` drops both "Some concerns" and "Serious concerns".
#' @param domain Which domain: `"Overall"` (default), `"D1"`, `"D2"`, `"D3"`,
#'   or `"D4"`.
#'
#' @return The filtered data frame, with attribute `excluded_studies`.
#'
#' @examples
#' \dontrun{
#' res <- inspect_sr(inspect_sr_example)
#' clean <- filter_trustworthy(my_data, res, studyvar = study)
#' }
#'
#' @export
filter_trustworthy <- function(data,
                               inspect_results,
                               studyvar = study,
                               exclude = c("serious", "some"),
                               domain = "Overall") {

  exclude <- rlang::arg_match(exclude)
  studyvar_name <- rlang::as_name(rlang::enquo(studyvar))

  if (!inherits(inspect_results, "inspect_sr")) {
    rlang::abort("`inspect_results` must come from `inspect_sr()`.")
  }
  if (!studyvar_name %in% names(data)) {
    rlang::abort(paste0("Column `", studyvar_name, "` not found in data."))
  }
  if (!domain %in% c("Overall", "D1", "D2", "D3", "D4")) {
    rlang::abort("`domain` must be 'Overall', 'D1', 'D2', 'D3', or 'D4'.")
  }

  judgements <- inspect_results[[domain]]
  study_names <- inspect_results$Study
  exclude_lev <- if (exclude == "serious") "Serious concerns"
  else c("Some concerns", "Serious concerns")

  excluded <- study_names[judgements %in% exclude_lev]
  filtered <- data[!data[[studyvar_name]] %in% excluded, ]

  attr(filtered, "excluded_studies") <- excluded
  attr(filtered, "exclude_level") <- exclude

  if (length(excluded) > 0) {
    message(sprintf("Excluded %d study/studies: %s",
                    length(excluded), paste(excluded, collapse = ", ")))
  } else {
    message("No studies excluded.")
  }

  filtered
}
