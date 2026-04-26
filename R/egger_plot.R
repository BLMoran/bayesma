#' Plot method for bayesma_egger
#'
#' @description
#' Creates diagnostic plots for the Bayesian Egger's test, including
#' a funnel plot with posterior latent SEs and a posterior distribution
#' plot for the slope parameter.
#'
#' @param x A bayesma_egger object
#' @param type Character. Type of plot: "funnel" (default), "beta_posterior",
#'   or "both".
#' @param show_observed Logical. Show observed SEs alongside latent SEs?
#'   Default TRUE.
#' @param contour_lines Logical. Add significance contour regions? Default TRUE.
#' @param contour_alpha Numeric. Alpha transparency for contour shading. Default 0.15.
#' @param color_observed Color for observed SE points. Default "grey50".
#' @param color_latent Color for latent SE points. Default "#D55E00".
#' @param color_pooled Color for pooled estimate line. Default "blue".
#' @param color_contour Color for contour fills. Default "grey50".
#' @param point_size Numeric. Size of study points. Default 3.
#' @param ... Additional arguments (currently unused)
#'
#' @export
egger_plot <- function(x, type = c("funnel", "beta_posterior", "both"),
                       show_observed = TRUE,
                       contour_lines = TRUE,
                       contour_alpha = 0.15,
                       color_observed = "grey50",
                       color_latent = "#D55E00",
                       color_pooled = "blue",
                       color_contour = "grey50",
                       point_size = 3,
                       ...) {
  type <- rlang::arg_match(type)

  if (type == "funnel" || type == "both") {
    # Use sigma_summary if available (binomial models), else fall back to observed
    if (!is.null(x$sigma_summary)) {
      plot_data <- tibble::tibble(
        yi = x$meta$yi,
        sei = x$meta$sei,
        sigma_median = x$sigma_summary$sigma_median,
        study = x$meta$study_labels
      )
      has_latent <- TRUE
    } else {
      plot_data <- tibble::tibble(
        yi = x$meta$yi,
        sei = x$meta$sei,
        sigma_median = x$meta$sei,
        study = x$meta$study_labels
      )
      has_latent <- FALSE
    }

    # Get pooled estimate (d for binomial, alpha for generic)
    if ("d" %in% names(x$draws)) {
      pooled_estimate <- stats::median(x$draws$d)
    } else {
      pooled_estimate <- stats::median(x$draws$alpha)
    }

    # Determine SE range for contours
    se_max <- max(c(plot_data$sei, plot_data$sigma_median), na.rm = TRUE) * 1.1

    # ---- Build funnel contours ----
    contour_data <- NULL
    border_data <- NULL
    ci_labels <- c("99% CI", "95% CI", "90% CI")

    if (isTRUE(contour_lines)) {
      z_vals <- c(stats::qnorm(0.995), stats::qnorm(0.975), stats::qnorm(0.95))
      se_seq <- seq(0, se_max, length.out = 300)

      contour_data <- purrr::map(seq_along(z_vals), function(i) {
        lower <- pooled_estimate - z_vals[i] * se_seq
        upper <- pooled_estimate + z_vals[i] * se_seq

        tibble::tibble(
          x = c(lower, rev(upper)),
          y = c(se_seq, rev(se_seq)),
          level = ci_labels[i]
        )
      }) |>
        purrr::list_rbind() |>
        dplyr::mutate(level = factor(level, levels = ci_labels))

      border_data <- purrr::map(seq_along(z_vals), function(i) {
        lower <- pooled_estimate - z_vals[i] * se_seq
        upper <- pooled_estimate + z_vals[i] * se_seq

        dplyr::bind_rows(
          tibble::tibble(x = lower, y = se_seq, side = "lower", level = ci_labels[i]),
          tibble::tibble(x = upper, y = se_seq, side = "upper", level = ci_labels[i])
        )
      }) |> purrr::list_rbind()
    }

    # ---- Build the plot ----
    p_funnel <- ggplot2::ggplot(plot_data, ggplot2::aes(x = yi, y = sei))

    # Contour shading (widest first)
    if (isTRUE(contour_lines) && !is.null(contour_data)) {
      alpha_steps <- c(contour_alpha, contour_alpha * 1.5, contour_alpha * 2.5)

      for (i in seq_along(ci_labels)) {
        cd <- contour_data |> dplyr::filter(.data$level == ci_labels[i])
        p_funnel <- p_funnel +
          ggplot2::geom_polygon(
            data = cd,
            ggplot2::aes(x = x, y = y),
            fill = color_contour,
            alpha = alpha_steps[i],
            color = NA,
            inherit.aes = FALSE
          )
      }

      # Dashed contour border lines
      p_funnel <- p_funnel +
        ggplot2::geom_line(
          data = border_data,
          ggplot2::aes(x = x, y = y, group = interaction(level, side)),
          color = color_contour, linetype = "dashed", linewidth = 0.4,
          inherit.aes = FALSE
        )
    }

    # Pooled estimate line
    p_funnel <- p_funnel +
      ggplot2::geom_vline(
        xintercept = pooled_estimate,
        color = color_pooled, linewidth = 0.75, linetype = "solid"
      )

    # Add observed SEs if requested and we have latent
    if (show_observed && has_latent) {
      p_funnel <- p_funnel +
        ggplot2::geom_point(
          ggplot2::aes(y = sei, shape = "Observed SE"),
          fill = color_observed, color = "grey30",
          size = point_size, stroke = 0.6
        )
    }

    # Add latent SEs (or observed if no latent available)
    if (has_latent) {
      p_funnel <- p_funnel +
        ggplot2::geom_point(
          ggplot2::aes(y = sigma_median, shape = "Latent SE"),
          fill = color_latent, color = "black",
          size = point_size, stroke = 0.6
        )
    } else {
      p_funnel <- p_funnel +
        ggplot2::geom_point(
          ggplot2::aes(y = sei),
          fill = color_observed, color = "grey30",
          shape = 21, size = point_size, stroke = 0.6
        )
    }

    # Shape scale for legend
    if (has_latent) {
      if (show_observed) {
        p_funnel <- p_funnel +
          ggplot2::scale_shape_manual(
            values = c("Observed SE" = 21, "Latent SE" = 24),
            name = "SE Type"
          )
      } else {
        p_funnel <- p_funnel +
          ggplot2::scale_shape_manual(
            values = c("Latent SE" = 24),
            name = "SE Type"
          )
      }
    }

    # Axis and theme
    p_funnel <- p_funnel +
      ggplot2::scale_y_reverse(expand = c(0, 0), limits = c(se_max, 0)) +
      ggplot2::scale_x_continuous(expand = c(0.02, 0)) +
      ggplot2::theme_light() +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "grey92"),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.75),
        axis.line.y.left = ggplot2::element_line(color = "black", linewidth = 0.75),
        axis.text = ggplot2::element_text(color = "black"),
        plot.margin = ggplot2::margin(10, 10, 10, 10),
        legend.position = "bottom"
      ) +
      ggplot2::labs(
        title = "Funnel Plot with Latent Standard Errors",
        subtitle = "Shi et al. (2020) Bayesian Egger's Test",
        x = "Effect Size (log OR)",
        y = "Standard Error"
      )

    if (type == "funnel") return(p_funnel)
  }

  if (type == "beta_posterior" || type == "both") {
    beta_draws <- x$draws$beta

    p_beta <- ggplot2::ggplot(tibble::tibble(beta = beta_draws), ggplot2::aes(x = beta)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                              bins = 50, fill = "steelblue", alpha = 0.7) +
      ggplot2::geom_density(color = "darkblue", linewidth = 1) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_vline(xintercept = x$beta_summary$lower, linetype = "dotted") +
      ggplot2::geom_vline(xintercept = x$beta_summary$upper, linetype = "dotted") +
      ggplot2::labs(
        title = "Posterior Distribution of Small-Study Effect (beta)",
        subtitle = sprintf("%.0f%% CrI: [%.3f, %.3f]",
                           x$meta$credible_level * 100,
                           x$beta_summary$lower,
                           x$beta_summary$upper),
        x = expression(beta),
        y = "Density"
      ) +
      ggplot2::theme_minimal()

    if (type == "beta_posterior") return(p_beta)
  }

  if (type == "both") {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(patchwork::wrap_plots(p_funnel, p_beta, ncol = 2))
    } else {
      cli::cli_warn("Install {.pkg patchwork} to combine plots. Returning funnel plot only.")
      return(p_funnel)
    }
  }
}
