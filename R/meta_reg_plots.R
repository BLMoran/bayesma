#' Plot Method for bayesma_coef_evidence Objects
#'
#' Creates a forest-style plot for meta-regression coefficients using
#' posterior density slabs, consistent with the `bayes_forest` visual style.
#'
#' @param x A `bayesma_coef_evidence` object (created by `coefficient_evidence()`).
#' @param model The `bayesma_metareg` model. Required to extract posterior
#'   draws for density slabs.
#' @param include_intercept Logical. Include intercept in plot (default: FALSE).
#' @param show_null_range Logical. Show null range region if available (default: TRUE).
#' @param null_value Numeric. Value for the null reference line (default: 0).
#' @param output Character. Either "density" (default) for density slabs or
#'   "pointinterval" for point + interval display.
#' @param split_color_by_null Logical. If TRUE, posterior densities are split and
#'   coloured based on whether values fall above or below the null value (default: FALSE).
#' @param color_posterior Color for coefficient posterior densities (default: "dodgerblue").
#' @param color_posterior_outline Color for posterior outlines (
#'   default: "blue").
#' @param color_favours_positive Colour for density regions > null when
#'   `split_color_by_null = TRUE` (default: "dodgerblue").
#' @param color_favours_negative Colour for density regions < null when
#'   `split_color_by_null = TRUE` (default: "firebrick").
#' @param color_null_range Color for null range shading (default: "grey50").
#' @param color_pointinterval Color for point intervals when `output = "pointinterval"`
#'   (default: "blue").
#' @param xlim Numeric vector of length 2 specifying x-axis limits. Default is NULL
#'   (auto-scaled).
#' @param x_breaks Numeric vector specifying custom x-axis break points. Default
#'   is NULL (auto).
#' @param xlab Character. X-axis label (default: "Coefficient Estimate").
#' @param title Character. Plot title (default: "Meta-Regression Coefficients").
#' @param subtitle Character. Plot subtitle. Default shows null range if provided.
#' @param add_table Logical. Add a table with coefficient summaries on the right
#'   (default: TRUE).
#' @param table_width Numeric. Relative width of the table vs plot (default: 0.4).
#' @param font Character. Font family for text elements (default: NULL).
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object (or patchwork object if `add_table = TRUE`).
#'
#' @details
#' This function creates a forest-style visualisation for meta-regression

#' coefficients that matches the aesthetic of `bayes_forest()`. Each coefficient
#' is displayed as a posterior density slab (or point interval), allowing
#' visualisation of the full posterior distribution rather than just point
#' estimates and credible intervals
#'
#' When `split_color_by_null = TRUE`, the density is split at the null value
#' and coloured to show the proportion of the posterior on each side.
#'
#' @examples
#' \dontrun{
#' fit <- meta_reg(data, studyvar = "author", yi = "yi", vi = "vi",
#'                 mods = ~ year + quality)
#'
#' # Get evidence summary
#' ev <- coefficient_evidence(fit, null_range = c(-0.1, 0.1))
#'
#' # Basic forest-style plot
#' metareg_mod_plot(ev, model = fit)
#'
#' # With split coloring
#' metareg_mod_plot(ev, model = fit, split_color_by_null = TRUE)
#'
#' # Point interval style
#' metareg_mod_plot(ev, model = fit, output = "pointinterval")
#' }
#'
#' @export
metareg_mod_plot <- function(x,
                                       model = NULL,
                                       include_intercept = FALSE,
                                       show_null_range = TRUE,
                                       null_value = 0,
                                       output = c("density", "pointinterval"),
                                       split_color_by_null = FALSE,
                                       color_posterior = "dodgerblue",
                                       color_posterior_outline = "blue",
                                       color_favours_positive = "dodgerblue",
                                       color_favours_negative = "firebrick",
                                       color_null_range = "grey50",
                                       color_pointinterval = "blue",
                                       xlim = NULL,
                                       x_breaks = NULL,
                                       xlab = "Coefficient Estimate",
                                       title = "Meta-Regression Coefficients",
                                       subtitle = NULL,
                                       add_table = TRUE,
                                       table_width = 0.4,
                                       font = NULL,
                                       ...) {

  # Match output argument

  output <- rlang::arg_match(output)

  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }

  if (!requireNamespace("ggdist", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggdist} is required for density slabs.")
  }

  # Get null_range from attributes
  null_range <- attr(x, "null_range")

  # Filter data
  plot_data <- x
  if (!include_intercept) {
    plot_data <- dplyr::filter(plot_data, .data$term != "(Intercept)")
  }

  # Order terms by estimate
  plot_data <- plot_data |>
    dplyr::mutate(term = factor(.data$term, levels = .data$term[order(.data$estimate)]))

  # Extract posterior draws if model provided (needed for density slabs)
  has_draws <- !is.null(model) && inherits(model, "bayesma_metareg")

  if (has_draws) {
    # Build draws dataframe for ggdist
    coef_names <- model$meta$coef_names

    draws_list <- purrr::map(plot_data$term, function(term) {
      if (term == "(Intercept)") {
        draws <- as.vector(
          posterior::subset_draws(model$fit$draws("mu"), variable = "mu")
        )
      } else {
        k <- which(coef_names == term)
        vn <- paste0("beta[", k, "]")
        draws <- as.vector(
          posterior::subset_draws(model$fit$draws(vn), variable = vn)
        )
      }
      tibble::tibble(term = term, value = draws)
    })

    draws_df <- dplyr::bind_rows(draws_list) |>
      dplyr::mutate(term = factor(.data$term, levels = levels(plot_data$term)))
  }

  # Calculate xlim if not provided
  if (is.null(xlim)) {
    if (has_draws) {
      xlim <- c(
        min(plot_data$ci_lower, na.rm = TRUE) * 1.1,
        max(plot_data$ci_upper, na.rm = TRUE) * 1.1
      )
    } else {
      xlim <- c(
        min(plot_data$ci_lower, na.rm = TRUE) - 0.1,
        max(plot_data$ci_upper, na.rm = TRUE) + 0.1
      )
    }
    # Ensure null_value is visible
    xlim[1] <- min(xlim[1], null_value - 0.1)
    xlim[2] <- max(xlim[2], null_value + 0.1)
  }

  # Set breaks
  breaks <- if (!is.null(x_breaks)) x_breaks else ggplot2::waiver()

  # Build subtitle
  if (is.null(subtitle) && !is.null(null_range) && show_null_range) {
    subtitle <- sprintf("Shaded region = null range [%.2f, %.2f]",
                        null_range[1], null_range[2])
  }

  # Build the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = .data$term))

  # Add null range shading first (background)
  if (show_null_range && !is.null(null_range)) {
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = null_range[1], xmax = null_range[2],
        ymin = -Inf, ymax = Inf,
        fill = scales::alpha(color_null_range, 0.2),
        color = NA
      )
  }

  # Add density or point interval
  if (has_draws && output == "density") {
    if (split_color_by_null) {
      p <- p +
        ggdist::stat_slab(
          data = draws_df,
          ggplot2::aes(x = .data$value, fill = ggplot2::after_stat(.data$x > null_value)),
          slab_linewidth = 0.5,
          alpha = 0.7,
          height = 0.85,
          normalize = "groups",
          colour = color_posterior_outline
        ) +
        ggplot2::scale_fill_manual(
          values = c("FALSE" = color_favours_negative, "TRUE" = color_favours_positive),
          guide = "none"
        )
    } else {
      p <- p +
        ggdist::stat_slab(
          data = draws_df,
          ggplot2::aes(x = .data$value),
          slab_linewidth = 0.5,
          alpha = 0.7,
          height = 0.85,
          normalize = "groups",
          colour = color_posterior_outline,
          fill = color_posterior
        )
    }
  } else if (has_draws && output == "pointinterval") {
    p <- p +
      ggdist::stat_pointinterval(
        data = draws_df,
        ggplot2::aes(x = .data$value),
        .width = c(0.66, 0.95),
        color = color_pointinterval,
        point_size = 3,
        interval_size_range = c(0.8, 2)
      )
  } else {
    # Fallback: use error bars if no model provided
    p <- p +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
        height = 0.2,
        linewidth = 0.8,
        color = color_posterior_outline
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$estimate),
        size = 3,
        color = color_posterior
      )
  }

  # Add null reference line
  p <- p +
    ggplot2::geom_vline(
      xintercept = null_value,
      color = "black",
      linewidth = 1
    )

  # Add reference lines at pooled estimate median (dashed)
  if (has_draws) {
    mu_median <- stats::median(as.vector(
      posterior::subset_draws(model$fit$draws("mu"), variable = "mu")
    ))
    p <- p +
      ggplot2::geom_vline(
        xintercept = mu_median,
        color = "grey60",
        linewidth = 0.8,
        linetype = "dashed"
      )
  }

  # Apply theme and styling (matching bayes_forest)
  p <- p +
    ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
    ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.5, 0.5))) +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(vjust = -0.5, family = font),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        size = 10,
        hjust = 1,
        family = font,
        color = "black"
      ),
      axis.text.x = ggplot2::element_text(colour = "black", family = font),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      panel.border = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.75),
      plot.title = ggplot2::element_text(family = font, face = "bold"),
      plot.subtitle = ggplot2::element_text(family = font, color = "grey40")
    ) +
    ggplot2::labs(
      x = xlab,
      y = NULL,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::geom_hline(
      yintercept = 0.5,
      color = "black",
      linewidth = 0.75
    )

  # Add coefficient table if requested
  if (add_table && requireNamespace("gt", quietly = TRUE) &&
      requireNamespace("patchwork", quietly = TRUE)) {

    # Create table data
    table_data <- plot_data |>
      dplyr::mutate(
        estimate_fmt = sprintf("%.2f", .data$estimate),
        ci_fmt = sprintf("[%.2f, %.2f]", .data$ci_lower, .data$ci_upper),
        pd_fmt = sprintf("%.1f%%", .data$pd * 100)
      ) |>
      dplyr::select(term, estimate_fmt, ci_fmt, pd_fmt)

    # Add p_outside_null if available
    if ("p_outside_null" %in% names(plot_data) && !all(is.na(plot_data$p_outside_null))) {
      table_data <- table_data |>
        dplyr::left_join(
          plot_data |> dplyr::select(term, p_outside_null),
          by = "term"
        ) |>
        dplyr::mutate(
          p_outside_fmt = dplyr::if_else(
            is.na(.data$p_outside_null),
            "",
            sprintf("%.1f%%", .data$p_outside_null * 100)
          )
        ) |>
        dplyr::select(-p_outside_null)
    }

    # Create gt table
    tbl <- table_data |>
      gt::gt() |>
      gt::cols_label(
        term = "",
        estimate_fmt = gt::md("**Est.**"),
        ci_fmt = gt::md("**95% CI**"),
        pd_fmt = gt::md("**PD**")
      ) |>
      gt::cols_align(align = "right", columns = -term) |>
      gt::cols_align(align = "left", columns = term) |>
      gt::tab_options(
        column_labels.font.weight = "bold",
        table.border.top.color = "white",
        table.border.bottom.color = "white",
        table.font.names = font,
        table.font.size = gt::px(11)
      ) |>
      gt::opt_table_lines(extent = "none")

    # Add p_outside column if present
    if ("p_outside_fmt" %in% names(table_data)) {
      tbl <- tbl |>
        gt::cols_label(p_outside_fmt = gt::md("**P(outside)**"))
    }

    # Combine plot and table using patchwork
    combined <- patchwork::wrap_plots(
      p,
      patchwork::wrap_elements(full = gt::as_gtable(tbl)),
      widths = c(1, table_width)
    )

    return(combined)
  }

  p
}


# Bubble Plot for Continuous Moderators

#' Bubble Plot for Meta-Regression
#'
#' Creates a bubble plot showing the relationship between a continuous
#' moderator and the effect size, with point sizes proportional to
#' study precision (inverse variance).
#'
#' @param object A `bayesma_reg` object.
#' @param mod Character. Name of the moderator variable to plot.
#'   Must be a continuous moderator.
#' @param ci Logical. Show confidence/credible band for the regression
#'   line (default: TRUE).
#' @param ci_level Numeric. Credible interval level (default: 0.95).
#' @param size_scale Numeric. Scaling factor for bubble sizes (default: 1).
#' @param show_studies Logical. Label study points (default: FALSE).
#' @param xlab,ylab Character. Axis labels. If NULL, uses variable names.
#' @param title Character. Plot title.
#' @param color_palette Character vector of length 2. Colors for points
#'   and regression line.
#' @param theme A ggplot2 theme (default: `theme_minimal()`).
#'
#' @return A ggplot object.
#'
#' @details
#' The bubble plot is a standard visualization for meta-regression with
#' continuous moderators. Each study is represented by a circle:
#'
#' - **Position**: x = moderator value, y = effect size
#' - **Size**: Proportional to study weight (1 / variance)
#' - **Line**: Regression line showing the moderator effect
#' - **Band**: Credible interval for the regression line
#'
#' For centered moderators, the x-axis shows the centered values by default.
#' Use `centered = FALSE` in the original `meta_reg()` call if you want
#' the original scale.
#'
#' @examples
#' \dontrun{
#' fit <- meta_reg(data, studyvar = "author", yi = "yi", vi = "vi",
#'                 mods = ~ year + dose)
#'
#' # Basic bubble plot
#' bubble_plot(fit, mod = "year")
#'
#' # Customized plot
#' bubble_plot(fit, mod = "dose",
#'             xlab = "Dose (mg)",
#'             ylab = "Log Odds Ratio",
#'             title = "Dose-Response Relationship")
#' }
#'
#' @export
bubble_plot <- function(object,
                        mod,
                        ci = TRUE,
                        ci_level = 0.95,
                        size_scale = 1,
                        show_studies = FALSE,
                        xlab = NULL,
                        ylab = NULL,
                        title = NULL,
                        color_palette = c("#4292C6", "#08519C"),
                        theme = ggplot2::theme_minimal()) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }

  if (!inherits(object, "bayesma_reg")) {
    cli::cli_abort("{.arg object} must be a {.cls bayesma_reg} object.")
  }

  # Check moderator exists
  coef_names <- object$meta$coef_names
  if (!mod %in% coef_names) {
    cli::cli_abort(c(
      "Moderator {.val {mod}} not found.",
      "i" = "Available moderators: {.val {coef_names}}"
    ))
  }

  # Check it's continuous
  mod_info <- object$meta$mod_info
  mod_idx <- which(colnames(mod_info$X) == mod)

  if (!mod_info$continuous_cols[mod_idx]) {
    cli::cli_abort(c(
      "Bubble plots are for continuous moderators.",
      "i" = "{.val {mod}} appears to be categorical.",
      "i" = "Consider using a forest plot grouped by this variable."
    ))
  }

  # Get data for plotting
  es <- object$meta$es
  X <- mod_info$X

  # Uncenter moderator values for plotting (if centered)
  x_vals <- X[, mod]
  if (mod_info$center) {
    x_vals <- x_vals + mod_info$center_values[mod]
  }
  if (mod_info$scale) {
    x_vals <- x_vals * mod_info$scale_values[mod]
  }

  # Prepare plot data
  plot_data <- tibble::tibble(
    study = object$meta$study_labels,
    x = x_vals,
    y = es$yi,
    se = es$sei,
    weight = 1 / es$vi
  )

  # Normalize weights for plotting
  plot_data$size <- sqrt(plot_data$weight / max(plot_data$weight)) * 10 * size_scale

  # Get regression coefficients
  mu_draws <- as.vector(
    posterior::subset_draws(object$fit$draws("mu"), variable = "mu")
  )

  beta_idx <- which(coef_names == mod)
  beta_var <- paste0("beta[", beta_idx, "]")
  beta_draws <- as.vector(
    posterior::subset_draws(object$fit$draws(beta_var), variable = beta_var)
  )

  # Compute regression line
  x_range <- range(x_vals)
  x_seq <- seq(x_range[1], x_range[2], length.out = 100)

  # Transform x_seq back to model scale (centered/scaled)
  x_seq_model <- x_seq
  if (mod_info$center) {
    x_seq_model <- x_seq_model - mod_info$center_values[mod]
  }
  if (mod_info$scale && mod_info$scale_values[mod] > 0) {
    x_seq_model <- x_seq_model / mod_info$scale_values[mod]
  }

  # Compute predictions
  # For single moderator: y = mu + beta * x
  # For multiple moderators: hold others at 0 (centered) or reference

  # Get other moderator values (hold at 0 for centered continuous, 0 for categorical)
  other_mods <- setdiff(coef_names, mod)
  other_beta_draws <- list()

  if (length(other_mods) > 0) {
    for (om in other_mods) {
      om_idx <- which(coef_names == om)
      om_var <- paste0("beta[", om_idx, "]")
      other_beta_draws[[om]] <- as.vector(
        posterior::subset_draws(object$fit$draws(om_var), variable = om_var)
      )
    }

    # Compute contribution from other moderators at their mean/reference
    # For centered continuous: contribution is 0
    # For categorical: contribution is 0 (reference level)
    other_contrib <- 0
  } else {
    other_contrib <- 0
  }

  # Prediction: mu + beta_mod * x + other_contrib
  pred_matrix <- purrr::map(x_seq_model, function(x) {
    mu_draws + beta_draws * x + other_contrib
  }) |> purrr::list_cbind()

  # Summarize predictions
  ci_probs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)

  line_data <- tibble::tibble(
    x = x_seq,
    y = apply(pred_matrix, 2, stats::median),
    lower = apply(pred_matrix, 2, stats::quantile, probs = ci_probs[1]),
    upper = apply(pred_matrix, 2, stats::quantile, probs = ci_probs[2])
  )

  # Build plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = .data$x, y = .data$y, size = .data$size),
      color = color_palette[1],
      alpha = 0.6
    ) +
    ggplot2::scale_size_identity()

  if (ci) {
    p <- p +
      ggplot2::geom_ribbon(
        data = line_data,
        ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
        fill = color_palette[2],
        alpha = 0.2
      )
  }

  p <- p +
    ggplot2::geom_line(
      data = line_data,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = color_palette[2],
      linewidth = 1
    )

  if (show_studies) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      cli::cli_warn("Install {.pkg ggrepel} for better label placement.")
      p <- p +
        ggplot2::geom_text(
          data = plot_data,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$study),
          size = 2.5,
          vjust = -1
        )
    } else {
      p <- p +
        ggrepel::geom_text_repel(
          data = plot_data,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$study),
          size = 2.5,
          max.overlaps = 15
        )
    }
  }

  # Labels
  if (is.null(xlab)) xlab <- mod
  if (is.null(ylab)) ylab <- object$meta$effect_label
  if (is.null(title)) title <- paste("Meta-Regression:", mod)

  p <- p +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    theme

  # Add coefficient annotation
  beta_est <- object$coefficients |>
    dplyr::filter(.data$term == mod)

  coef_label <- sprintf(
    "\u03B2 = %.3f (95%% CI: %.3f, %.3f)",
    beta_est$estimate, beta_est$q2.5, beta_est$q97.5
  )

  p <- p +
    ggplot2::annotate(
      "text",
      x = x_range[1],
      y = max(plot_data$y + plot_data$se),
      label = coef_label,
      hjust = 0,
      vjust = 1,
      size = 3.5
    )

  p
}


#' Multi-panel Bubble Plots
#'
#' Creates bubble plots for all continuous moderators in a meta-regression.
#'
#' @param object A `bayesma_reg` object.
#' @param ncol Integer. Number of columns in the plot grid (default: 2).
#' @param ... Additional arguments passed to [bubble_plot()].
#'
#' @return A combined ggplot object (using patchwork if available).
#'
#' @export
multi_bubble_plots <- function(object, ncol = 2, ...) {

  if (!inherits(object, "bayesma_reg")) {
    cli::cli_abort("{.arg object} must be a {.cls bayesma_reg} object.")
  }

  mod_info <- object$meta$mod_info
  continuous_mods <- names(mod_info$continuous_cols)[mod_info$continuous_cols]

  if (length(continuous_mods) == 0) {
    cli::cli_abort("No continuous moderators found for bubble plots.")
  }

  plots <- purrr::map(continuous_mods, function(mod) {
    bubble_plot(object, mod = mod, ...)
  })

  patchwork::wrap_plots(plots, ncol = ncol)

}

