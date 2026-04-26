# ============================================================================
# INSPECT-SR: Per-Check GT Table
# ============================================================================


#' Per-Check INSPECT-SR Results Table
#'
#' Unpacks the automated Domain 4 checks from an [inspect_sr()] result into
#' a publication-ready [gt::gt()] table, grouped by study with visual
#' separators.
#'
#' @param x An object returned by [inspect_sr()].
#' @param study Optional character vector of study names to restrict the table
#'   to (useful for per-study tabsets in Quarto). Default `NULL` = all studies.
#' @param check Character vector. Restrict to specific checks:
#'   `"grim"`, `"pvalue"`, `"carlisle"`, `"n_consistency"`. Default: all.
#' @param only_failed Logical. If `TRUE`, show only flagged rows.
#'   Default `FALSE`.
#'
#' @return A `gt` table object.
#'
#' @examples
#' \dontrun{
#' data(inspect_sr_example)
#' res <- inspect_sr(inspect_sr_example, verbose = FALSE)
#' inspect_sr_table(res)
#' inspect_sr_table(res, only_failed = TRUE)
#' inspect_sr_table(res, study = "Doe (1995)")
#' }
#'
#' @export
inspect_sr_table <- function(x,
                             study = NULL,
                             check = c("grim", "pvalue", "carlisle",
                                       "n_consistency"),
                             only_failed = FALSE) {

  if (!inherits(x, c("inspect_sr", "bayes_inspect_sr"))) {
    rlang::abort("`x` must be an object returned by inspect_sr().")
  }

  details <- attr(x, "details")
  if (is.null(details)) {
    rlang::abort("`x` has no `details` attribute.")
  }

  if (!is.null(study)) {
    all_studies <- vapply(details, `[[`, character(1), "study")
    missing <- setdiff(study, all_studies)
    if (length(missing) > 0) {
      rlang::abort(paste0(
        "Unknown study name(s): ", paste(missing, collapse = ", "),
        ". Available: ", paste(all_studies, collapse = ", ")
      ))
    }
    details <- details[all_studies %in% study]
  }

  check <- rlang::arg_match(check, multiple = TRUE)
  is_bayes <- inherits(x, "bayes_inspect_sr")

  # Build the long-format data frame
  rows <- list()

  add_row <- function(Study, Check, Item, Detail, Result, BF = NA_real_) {
    rows[[length(rows) + 1]] <<- data.frame(
      Study = Study, Check = Check, Item = Item,
      Detail = Detail, Result = Result, BF = BF,
      stringsAsFactors = FALSE
    )
  }

  fmt_bf <- function(bf) if (is_bayes) bf else NA_real_

  for (d in details) {
    study <- d$study

    # GRIM
    if ("grim" %in% check) {
      for (r in d$grim$results) {
        add_row(
          study, "GRIM",
          paste0(r$variable, " (", r$group, ")"),
          sprintf("mean = %.4g, n = %d", r$mean_value, r$n),
          if (isTRUE(r$consistent)) "Pass" else "Fail",
          fmt_bf(r$bf)
        )
      }
    }

    # P-value
    if ("pvalue" %in% check) {
      for (r in d$pvalue$results) {
        item <- if (is.null(r$context) || is.na(r$context)) {
          paste0(toupper(r$test_type), "-test")
        } else r$context
        add_row(
          study, "P-value", item,
          sprintf("reported p = %.4g, recalculated p = %.4g (diff %.4g)",
                  r$reported_p, r$recalculated_p, r$difference),
          if (isTRUE(r$consistent)) "Pass" else "Fail",
          fmt_bf(r$bf)
        )
      }
    }

    # Carlisle
    if ("carlisle" %in% check && !is.null(d$carlisle$result)) {
      r <- d$carlisle$result
      add_row(
        study, "Carlisle", "Baseline p-value distribution",
        sprintf("k = %d, %s combined p = %.4g, %s",
                r$n_comparisons, r$method, r$combined_p, r$interpretation),
        switch(r$interpretation,
               plausible = "Pass", too_similar = "Fail",
               too_different = "Fail", NA_character_),
        fmt_bf(d$carlisle$bf)
      )
    }

    # N consistency
    if ("n_consistency" %in% check && !is.null(d$n_consistency$result)) {
      r <- d$n_consistency$result
      if (is.data.frame(r$checks) && nrow(r$checks) > 0) {
        for (i in seq_len(nrow(r$checks))) {
          add_row(
            study, "N-consistency", r$checks$check[i],
            sprintf("expected = %g, observed = %g",
                    r$checks$expected[i], r$checks$observed[i]),
            if (isTRUE(r$checks$pass[i])) "Pass" else "Fail",
            fmt_bf(d$n_consistency$bf)
          )
        }
      }
    }
  }

  if (length(rows) == 0) {
    out <- data.frame(
      Study = character(), Check = character(), Item = character(),
      Detail = character(), Result = character(), BF = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, rows)
  }

  if (isTRUE(only_failed)) {
    keep <- out$Result == "Fail"
    if (is_bayes) keep <- keep | (!is.na(out$BF) & out$BF > 1)
    out <- out[keep, , drop = FALSE]
  }

  if (!is_bayes) out$BF <- NULL
  rownames(out) <- NULL

  # --- Build gt table ---
  # Two-level grouping:
  #   Study — outer row group, rendered as a stub column on the left.
  #   Check — inner grouping within each Study, rendered via the Check column
  #   with repeated values blanked and horizontal borders at each boundary.
  out <- out[order(out$Study, out$Check), , drop = FALSE]
  rownames(out) <- NULL

  # Row indices where Check changes within a Study (for borders).
  study_check <- paste(out$Study, out$Check, sep = "\r")
  check_boundary <- which(study_check != c("", utils::head(study_check, -1)))

  # Blank out repeated Check values within the same Study.
  display_check <- out$Check
  dup <- c(FALSE, utils::tail(study_check, -1) == utils::head(study_check, -1))
  display_check[dup] <- ""
  out$Check <- display_check

  tbl <- out |>
    gt::gt(groupname_col = "Study") |>
    gt::cols_align(align = "left", columns = "Check") |>
    gt::cols_align(align = "center", columns = "Result") |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = "Check")
    ) |>
    gt::tab_style(
      style = gt::cell_borders(sides = "top", color = "#cccccc",
                               weight = gt::px(1), style = "solid"),
      locations = gt::cells_body(rows = check_boundary)
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = "#e32400", weight = "bold"),
      locations = gt::cells_body(
        columns = "Result",
        rows = out$Result == "Fail"
      )
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = "#77bb41"),
      locations = gt::cells_body(
        columns = "Result",
        rows = out$Result == "Pass"
      )
    ) |>
    gt::tab_header(
      title = gt::md("**INSPECT-SR Automated Check Details**"),
      subtitle = "Domain 4 per-check results"
    ) |>
    gt::tab_options(
      row_group.as_column = TRUE,
      row_group.font.weight = "bold",
      row_group.border.top.color = "#333333",
      row_group.border.top.width = gt::px(2),
      row_group.border.bottom.color = "#cccccc",
      table.font.size = gt::px(13),
      column_labels.font.weight = "bold"
    ) |>
    gt::cols_label(
      Check = "Check",
      Item = "Item",
      Detail = "Detail",
      Result = "Result"
    )

  if (is_bayes) {
    tbl <- tbl |>
      gt::cols_label(BF = "BF") |>
      gt::fmt_number(columns = "BF", decimals = 2) |>
      gt::sub_missing(columns = "BF", missing_text = "\u2014")
  }

  tbl
}
