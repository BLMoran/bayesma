# Helper Functions for sensitivity_plot

#' Create Density Plot for Sensitivity Analysis
#'
#' @description
#' Creates the central density plot showing posterior distributions
#' across different sensitivity analyses.
#'
#' @param df Data frame containing posterior draws
#' @param measure Effect measure type
#' @param split_color_by_null Logical. Color based on null hypothesis
#' @param color_overall_posterior Fill color for densities
#' @param color_overall_posterior_outline Outline color for densities
#' @param color_favours_control Color when favoring control
#' @param color_favours_intervention Color when favoring intervention
#' @param label_control Label for control group
#' @param label_intervention Label for intervention group
#' @param xlim X-axis limits
#' @param x_breaks X-axis break points
#' @param null_value Null hypothesis value
#' @param null_range Range for practical equivalence
#' @param add_null_range Whether to show null range
#' @param color_null_range Color for null range region
#' @param font Font family
#'
#' @return A ggplot object
#'
#' @keywords internal
#' @noRd
sensitivity_density_plot_fn <- function(df,
                                        measure,
                                        split_color_by_null = FALSE,
                                        color_overall_posterior = "dodgerblue",
                                        color_overall_posterior_outline = "blue",
                                        color_favours_control = "firebrick",
                                        color_favours_intervention = "dodgerblue",
                                        label_control = "Control",
                                        label_intervention = "Intervention",
                                        xlim = NULL,
                                        x_breaks = NULL,
                                        null_value = NULL,
                                        null_range = NULL,
                                        add_null_range = NULL,
                                        color_null_range = "#77bb41",
                                        font = NULL) {

  props <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value
  breaks <- x_breaks %||% ggplot2::waiver()

  # For ratio measures (OR, RR, HR, IRR), compute density on log scale
  # This ensures the KDE bandwidth is appropriate for the log-scale display
  use_log_density <- isTRUE(props$log_scale)
  if (use_log_density) {
    df <- df |>
      dplyr::mutate(x_plot = log(x))
    null_value_plot <- log(null_value)
    null_range_plot <- if (!is.null(null_range)) log(null_range) else NULL
    xlim_plot <- if (!is.null(xlim)) log(xlim) else NULL
  } else {
    df <- df |>
      dplyr::mutate(x_plot = x)
    null_value_plot <- null_value
    null_range_plot <- null_range
    xlim_plot <- xlim
  }

  section_levels <- unique(df$section_label)

  df <- df |>
    dplyr::mutate(
      section_label = factor(section_label, levels = section_levels)
    )

  # Create plot_row levels in correct order: section first, then prior within section
  # This ensures the density rows match the table rows
  # Use factor levels if available (respects prior_order), otherwise use unique
  if (is.factor(df$prior_label)) {
    prior_levels <- levels(df$prior_label)
  } else {
    prior_levels <- unique(df$prior_label)
  }

  plot_row_levels <- expand.grid(
    prior_label = prior_levels,
    section_label = section_levels,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(plot_row = paste(section_label, "-", prior_label)) |>
    dplyr::pull(plot_row)

  df <- df |>
    dplyr::mutate(
      plot_row = factor(
        paste(section_label, "-", prior_label),
        levels = plot_row_levels
      ),
      plot_row = forcats::fct_rev(plot_row)
    )


  # Get the numeric positions of all factor levels
  all_levels <- levels(df$plot_row)
  level_positions <- seq_along(all_levels)

  # Find positions where section changes occur
  section_info <- df |>
    dplyr::distinct(plot_row, section_label) |>
    dplyr::arrange(plot_row) |>
    dplyr::mutate(
      level_num = match(plot_row, all_levels),
      section_changes = section_label != dplyr::lag(section_label, default = first(section_label))
    )

  # Get positions for lines between sections
  # We want lines after each section ends (at positions 3.5 and 6.5)
  section_breaks <- section_info |>
    dplyr::summarise(max_level = max(level_num), .by = section_label) |>
    dplyr::filter(max_level < max(level_positions)) |>
    dplyr::pull(max_level) |>
    purrr::map_dbl(~ . + 0.5)

  # Position for top line (above the highest level)
  top_line_pos <- max(level_positions) + 0.5

  # Fix xlim calculation
  calc_xlim <- if (!is.null(xlim_plot)) xlim_plot else range(df$x_plot, na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(y = plot_row)) +
    {if (isTRUE(add_null_range)) {
      ggplot2::annotate("rect",
                        xmin = null_range_plot[1], xmax = null_range_plot[2],
                        ymin = -Inf, ymax = Inf,
                        fill = scales::alpha(color_null_range, 0.3),
                        color = NA)
    }} +
    {if (isTRUE(split_color_by_null)) {
      ggdist::stat_slab(
        ggplot2::aes(x = x_plot, fill = ggplot2::after_stat(x > null_value_plot)),
        normalize = "groups",
        scale  = 0.7,
        colour = color_overall_posterior_outline,
        linewidth = 0.5,
        alpha = 0.7,
        position = ggplot2::position_nudge(y = -0.5)
      )
    } else {
      ggdist::stat_slab(
        ggplot2::aes(x = x_plot),
        normalize = "groups",
        scale  = 0.7,
        fill = color_overall_posterior,
        colour = color_overall_posterior_outline,
        linewidth = 0.5,
        alpha = 0.7,
        position = ggplot2::position_nudge(y = -0.5)
      )
    }} +
    {if (isTRUE(split_color_by_null)) {
      ggplot2::scale_fill_manual(
        values = c(
          "FALSE" = color_favours_intervention,
          "TRUE"  = color_favours_control
        ),
        guide = "none"
      )
    }} +

    # Null line
    ggplot2::geom_vline(
      xintercept = null_value_plot,
      linewidth = 1,
      color = "black"
    ) +

    # Add horizontal lines between sections (thinner and lighter)
    {if (length(section_breaks) > 0) {
      purrr::map(section_breaks, ~ {
        ggplot2::geom_hline(
          yintercept = .x,
          color = "grey85",  # Very light grey
          linewidth = 0.3    # Thin line
        )
      })
    }} +

    # Add grey line at the top of the plot
    ggplot2::geom_hline(
      yintercept = top_line_pos,
      color = "grey60",  # Medium grey for top line
      linewidth = 0.5
    ) +

    # Add grey line at the bottom of the plot
    ggplot2::geom_hline(
      yintercept = 0.5,
      color = "grey60",  # Medium grey for bottom line
      linewidth = 0.5
    ) +

    # Extend y-axis slightly to show top and bottom lines
    ggplot2::scale_y_discrete(expand = c(0, 0)) +

    # Add favours labels
    ggplot2::annotation_custom(
      grid::textGrob(
        label = paste(" Favours\n", label_control),
        x = grid::unit(0.97, "npc"),
        y = grid::unit(1.02, "npc"),
        just = c("right", "bottom"),
        gp = grid::gpar(col = "grey30", fontsize = 9, fontfamily = font)
      ),
      xmin = calc_xlim[1], xmax = calc_xlim[2]- 0.01, ymin = -Inf, ymax = Inf
    ) +

    ggplot2::annotation_custom(
      grid::textGrob(
        label = paste(" Favours\n", label_intervention),
        x = grid::unit(0.05, "npc"),
        y = grid::unit(1.02, "npc"),
        just = c("left", "bottom"),
        gp = grid::gpar(col = "grey30", fontsize = 9, fontfamily = font)
      ),
      xmin = calc_xlim[1] - 0.01, xmax = calc_xlim[2], ymin = -Inf, ymax = Inf
    ) +

    ggplot2::coord_cartesian(
      xlim = calc_xlim,
      ylim = c(0.5, max(level_positions) + 0.5),
      clip = "off"
    ) +


    ggplot2::theme_light() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(vjust = -0.5, family = font),
      panel.grid   = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(0, 0, 0, 0),
      axis.line.x.bottom = ggplot2::element_line(
        linewidth = 0.5, color = "grey60"),
      axis.line.x.top = ggplot2::element_line(
        linewidth = 0.5, color = "grey60"),
    ) +

    ggplot2::labs(x = props$x_label)

  # Apply appropriate scale
  # For ratio measures, data is already log-transformed, so we use
  # a custom transformation to display original scale labels
  if (use_log_density) {
    # Create custom breaks on the original scale
    if (is.null(x_breaks)) {
      # Default breaks for log scale
      orig_breaks <- c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10)
      # Filter to those within the data range
      data_range <- range(df$x, na.rm = TRUE)
      orig_breaks <- orig_breaks[orig_breaks >= data_range[1] * 0.9 &
                                   orig_breaks <= data_range[2] * 1.1]
      if (length(orig_breaks) < 3) {
        orig_breaks <- pretty(data_range, n = 5)
      }
    } else {
      orig_breaks <- x_breaks
    }

    p <- p + ggplot2::scale_x_continuous(
      breaks = log(orig_breaks),
      labels = orig_breaks,
      limits = calc_xlim,
      expand = c(0, 0)
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      breaks = breaks,
      limits = calc_xlim,
      expand = c(0, 0)
    )
  }

  p
}

#' Create Left Table for Sensitivity Plot
#'
#' @description
#' Creates the left-side table showing prior specifications for the
#' sensitivity analysis plot.
#'
#' @param df Data frame containing sensitivity analysis results
#' @param font Optional font family for the table
#' @param math_font Font for mathematical symbols. Default is "STIX Two Math"
#'
#' @return A gt table object
#'
#' @keywords internal
#' @noRd
sensitivity_table_left <- function(
    df,
    font = NULL,
    math_font = "STIX Two Math"
) {

  math_cols <- c("mu_prior_unicode", "tau_prior_unicode")
  math_cols_present <- intersect(math_cols, names(df))

  # Data arrives already sorted correctly from summarise_sensitivity_posteriors
  # Just select columns and add the last_in_group marker
  gt_tbl <- df |>
    dplyr::select(
      dplyr::any_of(c("section_label", "prior_label", "mu_prior_unicode", "tau_prior_unicode"))
    ) |>
    dplyr::mutate(.last_in_group = dplyr::row_number() == dplyr::n(), .by = section_label) |>
    gt::gt(groupname_col = "section_label") |>
    gt::tab_stubhead(label = "Meta-analysis Model") |>
    gt::cols_label(
      prior_label        = "Prior type",
      mu_prior_unicode   = paste0("\u03BC", " Prior"),
      tau_prior_unicode  = paste0("\u03C4", " Prior")
    ) |>
    gt::cols_align(align = "left") |>
    gt::tab_options(
      row_group.as_column = TRUE,
      row_group.font.weight = "bold",

      column_labels.font.weight = "bold",
      column_labels.border.bottom.color = "grey60",
      column_labels.border.bottom.width = gt::px(2),

      table.border.top.color = "white",
      table.border.bottom.color = "grey60",
      table.border.bottom.width = gt::px(2)
    ) |>
    gt::tab_style(
      style = gt::cell_borders(
        sides = "bottom",
        color = "grey60",
        weight = gt::px(2)
      ),
      locations = gt::cells_body(rows = .last_in_group)
    )

  # Apply normal font to non-math columns
  if (!is.null(font)) {
    gt_tbl <- gt_tbl |>
      gt::tab_style(
        style = gt::cell_text(font = font),
        locations = list(
          gt::cells_body(columns = -dplyr::any_of(math_cols_present)),
          gt::cells_column_labels(columns = -dplyr::any_of(math_cols_present)),
          gt::cells_stub(),
          gt::cells_row_groups()
        )
      )
  }

  # Apply math font to math columns (if present)
  if (length(math_cols_present) > 0) {
    gt_tbl <- gt_tbl |>
      gt::tab_style(
        style = gt::cell_text(font = math_font),
        locations = gt::cells_body(columns = dplyr::all_of(math_cols_present))
      )
  }

  gt_tbl |>
    gt::cols_hide(columns = ".last_in_group")
}

#' Create Right Table for Sensitivity Plot
#'
#' @description
#' Creates the right-side table showing effect estimates and probabilities
#' for the sensitivity analysis plot.
#'
#' @param df Data frame containing sensitivity analysis results
#' @param measure Effect measure type
#' @param add_probs Logical. Whether to include probability columns
#' @param font Optional font family for the table
#'
#' @return A gt table object
#'
#' @keywords internal
#' @noRd
sensitivity_table_right <- function(
    df,
    measure,
    add_probs = FALSE,
    add_probs_null_range = FALSE,
    font = NULL
) {

  # ---- Compatibility aliases (so table code doesn't depend on upstream naming) ----
  if (!"pr_no_benefit" %in% names(df) && "pr_harm" %in% names(df)) {
    df <- dplyr::rename(df, pr_no_benefit = pr_harm)
  }
  if (!"pr_benefit_null_range" %in% names(df) && "pr_benefit_gt_delta" %in% names(df)) {
    df <- dplyr::rename(df, pr_benefit_null_range = pr_benefit_gt_delta)
  }
  if (!"pr_no_benefit_null_range" %in% names(df) && "pr_harm_gt_delta" %in% names(df)) {
    df <- dplyr::rename(df, pr_no_benefit_null_range = pr_harm_gt_delta)
  }

  # Base probability columns (always present in the data)
  pr_base_cols <- intersect(
    c("pr_benefit", "pr_no_benefit"),
    names(df)
  )

  # Null range probability columns (only if they exist AND are not all NA)
  pr_null_cols <- character()
  if (isTRUE(add_probs_null_range)) {
    pr_null_candidates <- intersect(
      c("pr_benefit_null_range", "pr_no_benefit_null_range"),
      names(df)
    )
    # Only include if the columns actually have non-NA values
    pr_null_cols <- purrr::keep(pr_null_candidates, function(col) {
      !all(is.na(df[[col]]))
    })
  }

  pr_cols_present <- c(pr_base_cols, pr_null_cols)

  # Build estimate column and keep only what's available
  # Data arrives already sorted correctly from summarise_sensitivity_posteriors
  df <- df |>
    dplyr::mutate(
      estimate = sprintf("%.3f  [%.3f, %.3f]", median, l95, u95)
    ) |>
    dplyr::select(
      dplyr::any_of(c("section_label", "prior_label")),
      estimate,
      dplyr::any_of(pr_cols_present)
    ) |>
    dplyr::mutate(.last_in_group = dplyr::row_number() == dplyr::n(), .by = section_label)

  # Labels only for columns that exist
  labels <- c(
    estimate = paste(measure, "[95% CrI]"),
    pr_benefit = "Pr(Benefit)",
    pr_no_benefit = "Pr(Harm)",
    pr_benefit_null_range = paste("Pr(Benefit>", "\u03B4", ")"),
    pr_no_benefit_null_range = paste("Pr(Harm>", "\u03B4", ")")
  )
  labels <- labels[names(labels) %in% names(df)]

  table_right <- df |>
    gt::gt() |>
    gt::cols_label(.list = labels) |>
    gt::cols_align(align = "center", columns = gt::everything()) |>
    gt::tab_options(
      table.font.names = font,
      column_labels.font.weight = "bold",
      column_labels.border.bottom.color = "grey60",
      column_labels.border.bottom.width = gt::px(2),
      row_group.font.weight = "bold",
      row_group.font.size = gt::px(13),
      table.border.top.color = "white",
      table.border.bottom.color = "grey60",
      table.border.bottom.width = gt::px(2)
    ) |>
    gt::tab_style(
      style = gt::cell_borders(
        sides = "bottom",
        color = "grey60",
        weight = gt::px(2)
      ),
      locations = gt::cells_body(rows = .last_in_group)
    ) |>
    gt::cols_hide(columns = c("section_label", "prior_label", ".last_in_group"))

  # Percent formatting only if there are probability columns present
  if (length(pr_cols_present) > 0) {
    table_right <- table_right |>
      gt::fmt_percent(
        columns = dplyr::all_of(pr_cols_present),
        decimals = 1,
        scale_values = FALSE
      )
  }

  # Hide probability columns if user doesn't want them
  if (isFALSE(add_probs) && length(pr_cols_present) > 0) {
    table_right <- table_right |>
      gt::cols_hide(columns = dplyr::all_of(pr_cols_present))
  }

  table_right
}

# Function for RoBMA draw extraction
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

  # bayesma_robma output: averaged_draws$mu
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

  # Track which draws are from null models (exactly 0 on log scale)
  is_null_draw <- x == 0
  n_total <- length(x)
  n_null <- sum(is_null_draw)
  pct_null <- round(100 * n_null / n_total, 1)

  cli::cli_alert_info("RoBMA {prior}: {n_total} draws ({pct_null}% null model weight)")

  # RoBMA returns mu on log scale for ratio measures
  if (measure %in% c("OR", "RR", "HR", "IRR")) {
    x <- exp(x)
  }

  tibble::tibble(
    x = x,
    prior = prior,
    prior_label = prior_label,
    section_label = section_label,
    is_null_draw = is_null_draw
  )
}


# ============================================================================
# Internal: Posterior summaries for the sensitivity table
# ============================================================================

#' @noRd
summarise_sensitivity_posteriors <- function(draws, null_value, null_range) {

  # Preserve factor levels if they exist, otherwise use unique
  if (is.factor(draws$section_label)) {
    section_levels <- levels(draws$section_label)
  } else {
    section_levels <- unique(draws$section_label)
  }

  if (is.factor(draws$prior_label)) {
    prior_label_levels <- levels(draws$prior_label)
  } else {
    prior_label_levels <- unique(draws$prior_label)
  }

  # For prior (id), use the order from prior_label if available
  prior_levels <- unique(draws$prior)

  has_null_range <- !is.null(null_range) && length(null_range) == 2

  out <- draws |>
    dplyr::mutate(
      section_label = factor(section_label, levels = section_levels),
      prior_label   = factor(prior_label, levels = prior_label_levels),
      prior         = factor(prior, levels = prior_levels)
    ) |>
    dplyr::summarise(
      median = stats::median(x),
      l95    = stats::quantile(x, 0.025),
      u95    = stats::quantile(x, 0.975),

      # Base probability columns (always computed)
      pr_benefit    = round(mean(x < null_value) * 100, 1),
      pr_no_benefit = round(mean(x > null_value) * 100, 1),

      # Null range columns (only meaningful when null_range exists)
      pr_benefit_null_range    = if (has_null_range) {
        round(mean(x < null_range[1]) * 100, 1)
      } else {
        NA_real_
      },
      pr_no_benefit_null_range = if (has_null_range) {
        round(mean(x > null_range[2]) * 100, 1)
      } else {
        NA_real_
      },

      .by = c(section_label, prior, prior_label)
    )

  # Create row order based on factor levels, then convert to character
  # This ensures gt displays rows in the correct order
  out <- out |>
    dplyr::arrange(
      as.integer(section_label),
      as.integer(prior_label)
    ) |>
    dplyr::mutate(
      section_label = as.character(section_label),
      prior_label = as.character(prior_label)
    )

  # Drop null range columns entirely if no null range was specified
  if (!has_null_range) {
    out <- out |>
      dplyr::select(-dplyr::any_of(c(
        "pr_benefit_null_range", "pr_no_benefit_null_range"
      )))
  }

  out
}


# ============================================================================
# Internal: parallel mapping helpers
# ============================================================================

#' @noRd
bayesma_ensure_daemons <- function(workers) {
  daemons_were_set <- tryCatch({
    s <- mirai::status()
    is.data.frame(s$daemons) && nrow(s$daemons) > 0L
  }, error = function(e) FALSE)

  if (!daemons_were_set) {
    n_cores <- workers %||% max(1L, parallel::detectCores() - 1L)
    mirai::daemons(n_cores)
    tryCatch(
      mirai::everywhere(library(bayesma)),
      error = function(e) NULL
    )
    return(TRUE)
  }
  FALSE
}

#' @noRd
bayesma_future_map_dfr <- function(.x, .f, parallel = FALSE, workers = NULL, seed = TRUE, ...) {
  if (!isTRUE(parallel)) {
    return(purrr::map(.x, .f, ...) |> purrr::list_rbind())
  }

  dots <- list(...)

  if (requireNamespace("mirai", quietly = TRUE)) {
    teardown <- bayesma_ensure_daemons(workers)
    if (teardown) on.exit(mirai::daemons(0), add = TRUE)

    results <- mirai::mirai_map(
      .x,
      function(.x_i, .f, .dots) do.call(.f, c(list(.x_i), .dots)),
      .args = list(.f = .f, .dots = dots)
    )
    return(dplyr::bind_rows(results[]))
  }

  if (.Platform$OS.type != "windows") {
    n_cores <- workers %||% max(1L, parallel::detectCores() - 1L)
    results <- parallel::mclapply(.x, function(.x_i) {
      do.call(.f, c(list(.x_i), dots))
    }, mc.cores = n_cores, mc.set.seed = TRUE)
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
