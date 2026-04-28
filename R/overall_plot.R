#' Create posterior plots for Bayesian meta-analysis
#'
#' Produces density plots of the pooled effect (\eqn{\mu}), heterogeneity
#' (\eqn{\tau}), and/or prediction distribution from a fitted bayesma model.
#'
#' @param model A \code{bayesma} object.
#' @param estimand Character string for the effect measure or marginal estimand.
#'   Relative-effect: \code{"OR"}, \code{"HR"}, \code{"RR"}, \code{"IRR"},
#'   \code{"MD"}, \code{"SMD"}. Marginal: \code{"RD"}/\code{"ARR"}, \code{"ATE"},
#'   \code{"ATT"}, \code{"CATE"}. For marginal estimands the pooled posterior is
#'   drawn from \code{model$marginal$draws}.
#' @param add_tau Logical, whether to include the heterogeneity (tau) plot.
#' @param add_mu_prior Logical, whether to include the prior distribution for mu.
#' @param add_tau_prior Logical, whether to include the prior distribution for tau.
#' @param add_pred Logical, whether to include the prediction distribution plot.
#' @param color_pred_posterior Color for the prediction posterior distribution.
#' @param color_pred_posterior_outline Outline color for the prediction posterior distribution.
#' @param null_value Numeric, the null value for the effect (default from get_measure_properties).
#' @param null_range Numeric vector of length 2, range for null region (e.g., c(0.9, 1.1)).
#' @param add_null_range Logical, whether to add shaded null range region.
#' @param color_null_range Color for the null range shading.
#' @param label_control Character string for control group label.
#' @param label_intervention Character string for intervention group label.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param title_align Alignment for title ("left", "center", "right").
#' @param mu_xlim Numeric vector of length 2 for mu plot x-axis limits.
#' @param tau_xlim Numeric vector of length 2 for tau plot x-axis limits.
#' @param pred_xlim Numeric vector of length 2 for prediction plot x-axis limits.
#' @param x_breaks Numeric vector for x-axis breaks on mu plot.
#' @param tau_breaks Numeric vector for x-axis breaks on tau plot.
#' @param pred_breaks Numeric vector for x-axis breaks on prediction plot.
#' @param color_overall_posterior Color for the overall posterior distribution.
#' @param color_overall_posterior_outline Outline color for the posterior distribution.
#' @param split_color_by_null Logical, whether to split posterior color by null value.
#' @param color_favours_control Color for region favoring control.
#' @param color_favours_intervention Color for region favoring intervention.
#' @param tau_posterior_color Base color for tau posterior.
#' @param tau_posterior_outline Outline color for tau posterior.
#' @param tau_slab_scale Numeric, scaling factor for tau slab height relative to mu (default 0.7).
#' @param font Font family for text elements.
#' @param plot_arrangement Character, "vertical" (stacked) or "horizontal" (side by side).
#'
#' @return A ggplot object (or combined plot if add_tau = TRUE or add_pred = TRUE).
#' @export
#'
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats dnorm dt dcauchy
#' @importFrom ggplot2 ggplot aes theme_light element_blank element_rect element_line unit annotate geom_line geom_vline scale_fill_manual scale_x_log10 scale_x_continuous coord_cartesian labs theme margin element_text after_stat
#' @importFrom ggdist stat_slab stat_pointinterval scale_fill_ramp_discrete scale_thickness_shared
#' @importFrom patchwork plot_spacer plot_layout plot_annotation
#' @importFrom posterior as_draws_df
#' @importFrom dplyr filter group_by summarise mutate n
#' @importFrom tidyr complete
#' @importFrom scales pretty_breaks
#'
overall_plot <- function(model,
                         estimand,
                         add_tau = FALSE,
                         add_mu_prior = FALSE,
                         add_tau_prior = FALSE,
                         add_pred = FALSE,
                         color_pred_posterior = "forestgreen",
                         color_pred_posterior_outline = "darkgreen",
                         null_value = NULL,
                         null_range = NULL,
                         add_null_range = FALSE,
                         color_null_range = "#77bb41",
                         label_control = "Control",
                         label_intervention = "Intervention",
                         title = NULL,
                         subtitle = NULL,
                         title_align = "left",
                         mu_xlim = NULL,
                         tau_xlim = NULL,
                         pred_xlim = NULL,
                         x_breaks = NULL,
                         tau_breaks = NULL,
                         pred_breaks = NULL,
                         color_overall_posterior = "dodgerblue",
                         color_overall_posterior_outline = "blue",
                         split_color_by_null = FALSE,
                         color_favours_control = "firebrick",
                         color_favours_intervention = "dodgerblue",
                         tau_posterior_color = "#77bb41",
                         tau_posterior_outline = NULL,
                         tau_slab_scale = 0.7,
                         font = NULL,
                         plot_arrangement = "vertical") {

  if (!inherits(model, "bayesma")) {
    cli::cli_abort("{.arg model} must be a {.cls bayesma} object.")
  }

  props        <- get_measure_properties(estimand)
  null_value   <- null_value %||% props$null_value
  is_log_scale <- props$log_scale
  x_label      <- props$x_label

  tau_x_label <- if (estimand %in% c("MD", "SMD") ||
                     (is_marginal_estimand(estimand) && model$meta$likelihood == "gaussian")) {
    "Standard Deviation"
  } else {
    "Standard Deviation (Log Scale)"
  }

  # Set default tau outline color (darker version of base color)
  if (is.null(tau_posterior_outline)) {
    tau_rgb <- grDevices::col2rgb(tau_posterior_color)
    tau_posterior_outline <- grDevices::rgb(
      pmax(0, tau_rgb[1] * 0.5),
      pmax(0, tau_rgb[2] * 0.5),
      pmax(0, tau_rgb[3] * 0.5),
      maxColorValue = 255
    )
  }

  # ========================================================================
  # Extract posterior samples from bayesma object
  # ========================================================================
  is_re <- has_random_effects(model)

  mu_samples <- if (is_marginal_estimand(estimand)) {
    if (is.null(model$marginal)) {
      cli::cli_abort(
        c(
          "Marginal estimand {.val {estimand}} requires {.code $marginal} draws.",
          "i" = "Re-fit the model with {.code estimand = {.val {estimand}}} to attach marginal draws."
        )
      )
    }
    model$marginal$draws
  } else {
    as.numeric(model$draws[["mu"]])
  }
  tau_samples <- if (is_re && "tau" %in% names(model$draws)) {
    as.numeric(model$draws[["tau"]])
  } else {
    NULL
  }

  # --- MU PLOT ---

  # Transform function based on measure type
  transform_fn <- if (is_log_scale) exp else identity


  # Prepare mu data — remove non-finite values to prevent seq.default errors
  # in ggdist::stat_slab when combined with scale_x_log10
  mu_transformed <- transform_fn(mu_samples)
  keep <- is.finite(mu_transformed)
  if (is_log_scale) keep <- keep & mu_transformed > 0
  mu_transformed <- mu_transformed[keep]
  mu_df <- data.frame(mu = mu_transformed)

  # Set default x limits for mu plot with rounding
  mu_xlim <- if (!is.null(mu_xlim)) {
    mu_xlim
  } else {
    if (is_log_scale) {
      c(0.25, 4)
    } else {
      r <- range(mu_transformed, na.rm = TRUE)
      padding <- 0.1 * diff(r)
      c(floor((r[1] - padding) * 10) / 10, ceiling((r[2] + padding) * 10) / 10)
    }
  }

  # Set default x breaks — include lower limit
  if (is.null(x_breaks)) {
    if (is_log_scale) {
      x_breaks <- sort(unique(c(mu_xlim[1], 0.5, 1, 2, mu_xlim[2])))
    } else {
      x_breaks <- scales::pretty_breaks(n = 5)(mu_xlim)
      if (!mu_xlim[1] %in% x_breaks) {
        x_breaks <- sort(unique(c(mu_xlim[1], x_breaks)))
      }
    }
  } else {
    if (!mu_xlim[1] %in% x_breaks) {
      x_breaks <- sort(unique(c(mu_xlim[1], x_breaks)))
    }
  }

  # Determine title alignment
  title_hjust <- switch(title_align,
                        "left" = 0,
                        "center" = 0.5,
                        "right" = 1,
                        0)

  # Define axis line thickness
  axis_line_size <- 0.8

  # Base theme
  base_theme <- ggplot2::theme_light() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white"),
      axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = axis_line_size),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(color = "black", linewidth = axis_line_size),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      plot.title = ggplot2::element_text(hjust = title_hjust),
      plot.subtitle = ggplot2::element_text(hjust = title_hjust),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5, "pt")
    )

  # Start building mu plot
  mu_plot <- ggplot2::ggplot(mu_df, ggplot2::aes(x = mu))

  # Add null range shading if requested
  if (add_null_range && !is.null(null_range)) {
    mu_plot <- mu_plot +
      ggplot2::annotate("rect",
                        xmin = null_range[1], xmax = null_range[2],
                        ymin = -Inf, ymax = Inf,
                        fill = color_null_range, alpha = 0.2)
  }

  # Add prior FIRST (so it renders behind posterior) if requested
  if (add_mu_prior) {
    mu_prior_df <- build_mu_prior_overlay(
      model        = model,
      is_log_scale = is_log_scale,
      mu_xlim      = mu_xlim
    )

    if (!is.null(mu_prior_df)) {
      mu_plot <- mu_plot +
        ggplot2::geom_line(data = mu_prior_df,
                           ggplot2::aes(x = x, y = density),
                           color = "grey50",
                           linewidth = 0.8,
                           linetype = "solid",
                           inherit.aes = FALSE)
    }
  }

  # Add posterior distribution
  if (split_color_by_null) {
    mu_plot <- mu_plot +
      ggdist::stat_slab(ggplot2::aes(fill = ggplot2::after_stat(x < null_value)),
                        color = color_overall_posterior_outline,
                        alpha = 0.8,
                        limits = mu_xlim) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = color_favours_intervention,
                   "FALSE" = color_favours_control),
        guide = "none"
      )
  } else {
    mu_plot <- mu_plot +
      ggdist::stat_slab(fill = color_overall_posterior,
                        color = color_overall_posterior_outline,
                        alpha = 0.8,
                        limits = mu_xlim)
  }

  # Add point interval
  mu_plot <- mu_plot +
    ggdist::stat_pointinterval(.width = c(.66, .80, .95),
                               color = "black",
                               point_size = 2)

  # Add null reference line (black, solid)
  mu_plot <- mu_plot +
    ggplot2::geom_vline(xintercept = null_value, linetype = "solid",
                        color = "black", linewidth = 0.8)

  # Add direction labels (same y-level, bold) — at top
  label_y <- 0.95
  if (is_log_scale) {
    left_x  <- mu_xlim[1] * 1.4
    right_x <- mu_xlim[2] * 0.65
  } else {
    range_width <- diff(mu_xlim)
    left_x  <- mu_xlim[1] + range_width * 0.12
    right_x <- mu_xlim[2] - range_width * 0.12
  }

  mu_plot <- mu_plot +
    ggplot2::annotate("text", x = left_x, y = label_y,
                      label = paste0("Favours\n", label_intervention),
                      fontface = "bold.italic", size = 3, color = "grey30", vjust = 1) +
    ggplot2::annotate("text", x = right_x, y = label_y,
                      label = paste0("Favours\n", label_control),
                      fontface = "bold.italic", size = 3, color = "grey30", vjust = 1) +
    ggplot2::annotate("segment", x = mu_xlim[1], xend = mu_xlim[2],
                      y = Inf, yend = Inf,
                      linewidth = axis_line_size, color = "grey60")

  # Scale and coordinate system — minimal gap at bottom (0.03)
  if (is_log_scale) {
    mu_plot <- mu_plot +
      ggplot2::scale_x_log10(breaks = x_breaks,
                             limits = mu_xlim,
                             labels = function(x) sprintf("%.2g", x),
                             expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = mu_xlim, ylim = c(0.03, 1), clip = "off")
  } else {
    mu_plot <- mu_plot +
      ggplot2::scale_x_continuous(breaks = x_breaks,
                                  limits = mu_xlim,
                                  labels = function(x) sprintf("%.2g", x),
                                  expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = mu_xlim, ylim = c(0.03, 1), clip = "off")
  }

  mu_plot <- mu_plot +
    ggdist::scale_thickness_shared() +
    base_theme +
    ggplot2::labs(x = x_label,
                  y = NULL,
                  title = expression(paste("Overall Effect (", mu, ")")))

  # Apply custom font if specified
  if (!is.null(font)) {
    mu_plot <- mu_plot +
      ggplot2::theme(text = ggplot2::element_text(family = font))
  }

  # --- PREDICTION PLOT ---
  pred_plot <- NULL

  if (isTRUE(add_pred)) {
    # bayesma stores mu_new draws directly (for RE models)
    if ("mu_new" %in% names(model$draws)) {
      pred_raw <- as.numeric(model$draws[["mu_new"]])
    } else {
      cli::cli_warn(
        c("No prediction draws ({.val mu_new}) found in the bayesma object.",
          "i" = "Prediction plots require a random-effects model.")
      )
      pred_raw <- NULL
    }

    if (!is.null(pred_raw)) {
      pred_transformed <- transform_fn(pred_raw)
      pred_df <- data.frame(pred = pred_transformed)

      # Set default prediction x limits
      pred_xlim <- pred_xlim %||% mu_xlim

      # Set default prediction breaks
      if (is.null(pred_breaks)) pred_breaks <- x_breaks

      # Build prediction plot
      pred_plot <- ggplot2::ggplot(pred_df, ggplot2::aes(x = pred))

      # Add null range shading if requested
      if (add_null_range && !is.null(null_range)) {
        pred_plot <- pred_plot +
          ggplot2::annotate("rect",
                            xmin = null_range[1], xmax = null_range[2],
                            ymin = -Inf, ymax = Inf,
                            fill = color_null_range, alpha = 0.2)
      }

      # Add posterior distribution
      if (split_color_by_null) {
        pred_plot <- pred_plot +
          ggdist::stat_slab(
            ggplot2::aes(fill = ggplot2::after_stat(x < null_value)),
            color = color_pred_posterior_outline,
            alpha = 0.8) +
          ggplot2::scale_fill_manual(
            values = c("TRUE" = color_favours_intervention,
                       "FALSE" = color_favours_control),
            guide = "none")
      } else {
        pred_plot <- pred_plot +
          ggdist::stat_slab(fill = color_pred_posterior,
                            color = color_pred_posterior_outline,
                            alpha = 0.8)
      }

      # Add point interval
      pred_plot <- pred_plot +
        ggdist::stat_pointinterval(.width = c(.66, .80, .95),
                                   color = "black",
                                   point_size = 2)

      # Add null reference line
      pred_plot <- pred_plot +
        ggplot2::geom_vline(xintercept = null_value, linetype = "solid",
                            color = "black", linewidth = 0.8)

      # Add direction labels
      if (is_log_scale) {
        pred_left_x  <- pred_xlim[1] * 1.4
        pred_right_x <- pred_xlim[2] * 0.65
      } else {
        pred_range_width <- diff(pred_xlim)
        pred_left_x  <- pred_xlim[1] + pred_range_width * 0.12
        pred_right_x <- pred_xlim[2] - pred_range_width * 0.12
      }

      pred_plot <- pred_plot +
        ggplot2::annotate("text", x = pred_left_x, y = label_y,
                          label = paste0("Favours\n", label_intervention),
                          fontface = "bold.italic", size = 3,
                          color = "grey30", vjust = 1) +
        ggplot2::annotate("text", x = pred_right_x, y = label_y,
                          label = paste0("Favours\n", label_control),
                          fontface = "bold.italic", size = 3,
                          color = "grey30", vjust = 1) +
        ggplot2::annotate("segment",
                          x = pred_xlim[1], xend = pred_xlim[2],
                          y = Inf, yend = Inf,
                          linewidth = axis_line_size, color = "grey60")

      # Scale and coordinate system
      if (is_log_scale) {
        pred_plot <- pred_plot +
          ggplot2::scale_x_log10(breaks = pred_breaks,
                                 limits = pred_xlim,
                                 labels = function(x) sprintf("%.2g", x),
                                 expand = c(0, 0)) +
          ggplot2::coord_cartesian(xlim = pred_xlim, ylim = c(0.03, 1),
                                   clip = "off")
      } else {
        pred_plot <- pred_plot +
          ggplot2::scale_x_continuous(breaks = pred_breaks,
                                      limits = pred_xlim,
                                      labels = function(x) sprintf("%.2g", x),
                                      expand = c(0, 0)) +
          ggplot2::coord_cartesian(xlim = pred_xlim, ylim = c(0.03, 1),
                                   clip = "off")
      }

      pred_plot <- pred_plot +
        ggdist::scale_thickness_shared() +
        base_theme +
        ggplot2::labs(x = x_label, y = NULL, title = "Prediction")

      # Apply custom font if specified
      if (!is.null(font)) {
        pred_plot <- pred_plot +
          ggplot2::theme(text = ggplot2::element_text(family = font))
      }
    }
  }

  # --- TAU PLOT ---
  tau_plot <- NULL

  if (add_tau) {
    if (is.null(tau_samples)) {
      cli::cli_warn(
        c("No heterogeneity parameter found.",
          "i" = "Tau plot requires a random-effects model.")
      )
    } else {
      tau_df <- data.frame(tau = tau_samples)

      # Define heterogeneity cutpoints
      tau_cuts   <- c(0, 0.1, 0.5, 1, Inf)
      tau_labels_vec <- c("Low", "Reasonable", "Fairly high", "Fairly extreme")

      # Categorize tau samples and calculate percentages
      tau_df$category <- cut(tau_df$tau,
                             breaks = tau_cuts,
                             labels = tau_labels_vec,
                             include.lowest = TRUE,
                             right = FALSE)

      # Calculate percentage in each category
      tau_percentages <- tau_df |>
        dplyr::summarise(n = dplyr::n(), .by = category) |>
        dplyr::mutate(pct = round(n / sum(n) * 100, 1)) |>
        tidyr::complete(
          category = factor(tau_labels_vec, levels = tau_labels_vec),
          fill = list(n = 0, pct = 0)
        )

      pct_low         <- tau_percentages$pct[tau_percentages$category == "Low"]
      pct_reasonable  <- tau_percentages$pct[tau_percentages$category == "Reasonable"]
      pct_fairly_high <- tau_percentages$pct[tau_percentages$category == "Fairly high"]
      pct_extreme     <- tau_percentages$pct[tau_percentages$category == "Fairly extreme"]

      # Set default tau x limits with rounding
      tau_xlim <- if (!is.null(tau_xlim)) {
        tau_xlim
      } else {
        r <- range(tau_df$tau, na.rm = TRUE)
        padding <- 0.1 * diff(r)
        c(max(0, floor((r[1] - padding) * 10) / 10),
          ceiling((r[2] + padding) * 10) / 10)
      }

      # Set default tau breaks
      if (is.null(tau_breaks)) {
        reference_points <- c(0, 0.10, 0.25, 0.50, 0.75, 1.00)
        if (tau_xlim[2] > 1.5) {
          extra_breaks <- seq(1.5, tau_xlim[2], by = 0.5)
          extra_breaks <- extra_breaks[extra_breaks < tau_xlim[2]]
          reference_points <- c(reference_points, extra_breaks)
        }
        tau_breaks <- reference_points[reference_points >= tau_xlim[1] &
                                         reference_points <= tau_xlim[2]]
        tau_breaks <- sort(unique(c(tau_xlim[1], tau_breaks, tau_xlim[2])))
      } else {
        if (!tau_xlim[1] %in% tau_breaks)
          tau_breaks <- sort(unique(c(tau_xlim[1], tau_breaks)))
        if (!tau_xlim[2] %in% tau_breaks)
          tau_breaks <- sort(unique(c(tau_breaks, tau_xlim[2])))
      }

      # Build tau plot with scaled slab height
      tau_plot <- ggplot2::ggplot(tau_df, ggplot2::aes(x = tau)) +
        ggdist::stat_slab(
          ggplot2::aes(
            fill_ramp = ggplot2::after_stat(
              cut(x,
                  breaks = c(-Inf, 0.1, 0.5, 1, Inf),
                  labels = c("Low", "Reasonable",
                             "Fairly high", "Fairly extreme"))
            )),
          fill = tau_posterior_color,
          color = tau_posterior_outline,
          linewidth = 0.8,
          scale = tau_slab_scale
        ) +
        ggdist::scale_fill_ramp_discrete(
          range = c(0.8, 0.2),
          guide = "none"
        ) +
        ggdist::stat_pointinterval(.width = c(.66, .80, .95),
                                   color = "black",
                                   point_size = 2)

      # Add tau prior if requested
      if (add_tau_prior) {
        tau_prior_df <- build_tau_prior_overlay(
          model          = model,
          tau_xlim       = tau_xlim,
          tau_slab_scale = tau_slab_scale
        )

        if (!is.null(tau_prior_df)) {
          tau_plot <- tau_plot +
            ggplot2::geom_line(data = tau_prior_df,
                               ggplot2::aes(x = x, y = density),
                               color = "grey50",
                               linewidth = 0.8,
                               linetype = "solid",
                               inherit.aes = FALSE)
        }
      }

      # Add heterogeneity category reference lines (dashed, grey)
      tau_plot <- tau_plot +
        ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed",
                            color = "grey40", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                            color = "grey40", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                            color = "grey40", linewidth = 0.5)

      very_high_x <- (1 + tau_xlim[2]) / 2

      # Add labels with percentages — only if region within tau_xlim
      if (tau_xlim[1] < 0.1 && tau_xlim[2] > 0) {
        tau_plot <- tau_plot +
          ggplot2::annotate("text", x = 0.05, y = 0.95,
                            label = paste0("Low\n(", pct_low, "%)"),
                            fontface = "bold.italic", size = 2.6,
                            color = "grey30", vjust = 1)
      }
      if (tau_xlim[1] < 0.5 && tau_xlim[2] > 0.1) {
        tau_plot <- tau_plot +
          ggplot2::annotate("text", x = 0.3, y = 0.95,
                            label = paste0("Moderate\n(", pct_reasonable, "%)"),
                            fontface = "bold.italic", size = 2.6, vjust = 1)
      }
      if (tau_xlim[1] < 1.0 && tau_xlim[2] > 0.5) {
        tau_plot <- tau_plot +
          ggplot2::annotate("text", x = 0.75, y = 0.95,
                            label = paste0("High\n(", pct_fairly_high, "%)"),
                            fontface = "bold.italic", size = 2.6, vjust = 1)
      }
      if (tau_xlim[2] > 1.0) {
        tau_plot <- tau_plot +
          ggplot2::annotate("text", x = very_high_x, y = 0.95,
                            label = paste0("Very High\n(", pct_extreme, "%)"),
                            fontface = "bold.italic", size = 2.6, vjust = 1)
      }

      tau_plot <- tau_plot +
        ggplot2::annotate("segment",
                          x = tau_xlim[1], xend = tau_xlim[2],
                          y = Inf, yend = Inf,
                          linewidth = axis_line_size, color = "grey60")

      # Scale and theme for tau plot
      tau_plot <- tau_plot +
        ggplot2::scale_x_continuous(breaks = tau_breaks,
                                    limits = tau_xlim,
                                    labels = function(x) sprintf("%.2g", x),
                                    expand = c(0, 0)) +
        ggplot2::coord_cartesian(xlim = tau_xlim, ylim = c(0.03, 1),
                                 clip = "off") +
        ggdist::scale_thickness_shared() +
        base_theme +
        ggplot2::labs(x = tau_x_label,
                      y = NULL,
                      title = expression(paste("Heterogeneity (", tau, ")")))

      # Apply custom font if specified
      if (!is.null(font)) {
        tau_plot <- tau_plot +
          ggplot2::theme(text = ggplot2::element_text(family = font))
      }
    }
  }

  # Collect all plots that exist
  plot_list <- list(mu_plot)
  if (isTRUE(add_pred) && !is.null(pred_plot)) {
    plot_list <- c(plot_list, list(pred_plot))
  }
  if (isTRUE(add_tau) && !is.null(tau_plot)) {
    plot_list <- c(plot_list, list(tau_plot))
  }

  # If only mu plot, return it directly
  if (length(plot_list) == 1) {
    final_plot <- mu_plot
    if (!is.null(title) || !is.null(subtitle)) {
      final_plot <- final_plot +
        ggplot2::labs(
          title = if (!is.null(title)) title
          else expression(paste("Overall Effect (", mu, ")")),
          subtitle = subtitle
        )
    }
    return(final_plot)
  }

  # Multiple plots — combine with patchwork
  if (plot_arrangement == "horizontal") {
    spaced_list <- list(plot_list[[1]])
    for (i in seq_along(plot_list)[-1]) {
      spaced_list <- c(spaced_list,
                       list(patchwork::plot_spacer()),
                       list(plot_list[[i]]))
    }
    n_plots   <- length(plot_list)
    n_spacers <- n_plots - 1
    widths <- as.list(rep(c(1, 0.05), n_plots))
    widths <- widths[seq_len(n_plots + n_spacers)]
    combined_plot <- Reduce("+", spaced_list) +
      patchwork::plot_layout(ncol = n_plots + n_spacers,
                             widths = unlist(widths))
  } else {
    combined_plot <- Reduce("/", plot_list) +
      patchwork::plot_layout(ncol = 1)
  }

  # Add title and subtitle using patchwork
  if (!is.null(title) || !is.null(subtitle)) {
    annotation_theme <- ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = title_hjust,
                                         face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = title_hjust, size = 11)
    )

    if (!is.null(font)) {
      annotation_theme <- annotation_theme +
        ggplot2::theme(
          text = ggplot2::element_text(family = font),
          plot.title = ggplot2::element_text(
            hjust = title_hjust, face = "bold", size = 14, family = font),
          plot.subtitle = ggplot2::element_text(
            hjust = title_hjust, size = 11, family = font)
        )
    }

    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title    = title,
        subtitle = subtitle,
        theme    = annotation_theme
      )
  }

  combined_plot
}



# Internal helpers for prior density overlays

#' Build mu prior density curve for overlay
#'
#' Computes the prior density on the appropriate scale (log or natural)
#' for bayesma models.
#'
#' @param model A bayesma model object.
#' @param is_log_scale Logical.
#' @param mu_xlim Numeric(2).
#'
#' @return data.frame(x, density) or NULL.
#'
#' @noRd
build_mu_prior_overlay <- function(model, is_log_scale, mu_xlim) {

  prior_obj <- model$meta$priors$mu
  if (is.null(prior_obj)) return(NULL)

  # Compute density using the bayesma_prior object directly
  bayesma_prior_density(prior_obj,
                        is_log_scale = is_log_scale,
                        xlim         = mu_xlim,
                        half         = FALSE)
}


#' Build tau prior density curve for overlay
#'
#' @param model A bayesma model object.
#' @param tau_xlim Numeric(2).
#' @param tau_slab_scale Numeric scaling factor.
#'
#' @return data.frame(x, density) or NULL.
#'
#' @noRd
build_tau_prior_overlay <- function(model, tau_xlim, tau_slab_scale = 0.7) {

  prior_obj <- model$meta$priors$tau
  if (is.null(prior_obj)) return(NULL)

  bayesma_prior_density(prior_obj,
                        is_log_scale = FALSE,
                        xlim         = tau_xlim,
                        half         = TRUE,
                        scale_factor = tau_slab_scale * 0.5)
}


#' Compute prior density from a bayesma_prior object
#'
#' @param prior_obj A bayesma_prior object.
#' @param is_log_scale Logical — evaluate on log scale then exponentiate x?
#' @param xlim Numeric(2) axis limits.
#' @param half Logical — double density for half-distributions?
#' @param scale_factor Numeric normalisation ceiling (default 0.5).
#'
#' @return data.frame(x, density) or NULL.
#'
#' @noRd
bayesma_prior_density <- function(prior_obj, is_log_scale, xlim,
                                  half = FALSE, scale_factor = 0.5) {
  fam <- prior_obj$family

  if (is_log_scale) {
    x_seq <- seq(log(xlim[1] * 0.1), log(xlim[2] * 5), length.out = 500)
  } else if (half) {
    x_seq <- seq(0.001, xlim[2] * 1.5, length.out = 500)
  } else {
    x_seq <- seq(xlim[1] - diff(xlim) * 2,
                 xlim[2] + diff(xlim) * 2,
                 length.out = 500)
  }

  d <- switch(
    fam,
    normal = stats::dnorm(x_seq, mean = prior_obj$mean, sd = prior_obj$sd),
    half_normal = {
      raw <- stats::dnorm(x_seq, mean = prior_obj$mean, sd = prior_obj$sd)
      raw * 2
    },
    half_cauchy = {
      raw <- stats::dcauchy(x_seq,
                            location = prior_obj$location,
                            scale    = prior_obj$scale)
      raw * 2
    },
    half_student_t = {
      raw <- stats::dt((x_seq - prior_obj$location) / prior_obj$scale,
                       df = prior_obj$df) / prior_obj$scale
      raw * 2
    },
    exponential = stats::dexp(x_seq, rate = prior_obj$rate),
    NULL
  )

  if (is.null(d) || max(d) == 0) return(NULL)

  x_out <- if (is_log_scale) exp(x_seq) else x_seq

  data.frame(
    x       = x_out,
    density = d / max(d) * scale_factor
  )
}
