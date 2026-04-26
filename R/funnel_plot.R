#' Create a Funnel Plot for Bayesian Meta-Analysis
#'
#' @description
#' Creates a publication-ready funnel plot for meta-analysis, displaying
#' study-level effect sizes against their precision (standard error).
#' Supports both brmsfit and bayesma model objects. The plot aesthetic
#' matches the bayesfoRest package style.
#'
#' @param model A fitted model object. Either a brmsfit object (class
#'   'brmsfit') or a bayesma object (class 'bayesma').
#' @param data A data frame containing the study data used for the meta-analysis.
#' @param measure Character string specifying the effect measure. Must be one of:
#'   "OR" (Odds Ratio), "HR" (Hazard Ratio), "RR" (Risk Ratio),
#'   "IRR" (Incidence Rate Ratio), "MD" (Mean Difference), or "SMD" (Standardized Mean Difference).
#' @param studyvar Column name containing study identifiers/authors. Default is NULL.
#' @param year Column name containing publication years. Default is NULL.
#' @param c_n Column name containing control group sample sizes.
#' @param i_n Column name containing intervention group sample sizes.
#' @param c_event Column name containing control group event counts.
#' @param i_event Column name containing intervention group event counts.
#' @param c_mean Column name containing control group means.
#' @param i_mean Column name containing intervention group means.
#' @param c_sd Column name containing control group standard deviations.
#' @param i_sd Column name containing intervention group standard deviations.
#' @param c_time Column name containing control group time periods.
#' @param i_time Column name containing intervention group time periods.
#' @param subgroup Logical indicating whether to colour points by subgroup. Default is FALSE.
#' @param subgroup_var Character string. Name of the variable in data to use for subgroup colouring.
#' @param contour_lines Logical indicating whether to add significance contour
#'   lines (at p = 0.10, 0.05, 0.01). Default is TRUE.
#' @param contour_alpha Numeric. Alpha transparency for contour shading. Default is 0.15.
#' @param pooled_line Character string specifying the pooled estimate line style.
#'   Options: "posterior" (default, uses posterior median), "fixed" (uses the
#'   fixed/common effect), or "none" (no pooled line).
#' @param color_points Color for study points. Default is "dodgerblue".
#' @param color_points_outline Color for study point outlines. Default is "blue".
#' @param color_pooled Color for the pooled effect line. Default is "blue".
#' @param color_null Color for the null effect line. Default is "black".
#' @param color_contour Color for significance contour fills. Default is "grey50".
#' @param point_size Numeric. Size of the study points. Default is 3.
#' @param xlim Numeric vector of length 2 specifying x-axis limits. Default is NULL (auto-scaled).
#' @param x_breaks Numeric vector specifying custom x-axis break points. Default is NULL.
#' @param null_value Numeric. X-axis value of the null line. Default is NULL (measure-specific).
#' @param add_null_line Logical. Whether to draw a vertical line at the
#'   null value. Default is FALSE.
#' @param label_studies Logical indicating whether to label study points. Default is FALSE.
#' @param label_size Numeric. Size of study labels. Default is 3.
#' @param title Character string for the plot title. Default is NULL.
#' @param subtitle Character string for the plot subtitle. Default is NULL.
#' @param title_align Character string specifying title alignment.
#'   Options: "left" (default), "center"/"centre", "right".
#' @param font Character string specifying the font family. Default is NULL.
#'
#' @return A ggplot object containing the funnel plot.
#'
#' @export
funnel_plot <- function(model,
                        data,
                        measure,
                        studyvar = NULL,
                        year = NULL,
                        c_n = NULL,
                        i_n = NULL,
                        c_event = NULL,
                        i_event = NULL,
                        c_mean = NULL,
                        i_mean = NULL,
                        c_sd = NULL,
                        i_sd = NULL,
                        c_time = NULL,
                        i_time = NULL,
                        subgroup = FALSE,
                        subgroup_var = NULL,
                        contour_lines = TRUE,
                        contour_alpha = 0.15,
                        pooled_line = c("posterior", "fixed", "none"),
                        color_points = "dodgerblue",
                        color_points_outline = "blue",
                        color_pooled = "blue",
                        color_null = "black",
                        color_contour = "grey50",
                        point_size = 3,
                        xlim = NULL,
                        x_breaks = NULL,
                        null_value = NULL,
                        add_null_line = FALSE,
                        label_studies = FALSE,
                        label_size = 3,
                        title = NULL,
                        subtitle = NULL,
                        title_align = "left",
                        font = NULL) {

  pooled_line <- rlang::arg_match(pooled_line)

  # ---- Detect model class ----
  is_bayesma <- inherits(model, "bayesma")
  is_brms    <- inherits(model, "brmsfit")

  if (!is_bayesma && !is_brms) {
    rlang::abort("`model` must be a {.cls brmsfit} or {.cls bayesma} object.")
  }

  # ---- Get measure properties ----
  props <- get_measure_properties(measure)
  null_value <- if (!is.null(null_value)) null_value else props$null_value

  # ---- Column renaming ----
  if (measure %in% c("OR", "RR")) {
    data <- data |>
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{year}},
        N_Control = {{c_n}},
        N_Intervention = {{i_n}},
        Event_Control = {{c_event}},
        Event_Intervention = {{i_event}}
      )
  } else if (measure %in% c("MD", "SMD")) {
    data <- data |>
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{year}},
        N_Control = {{c_n}},
        N_Intervention = {{i_n}},
        Mean_Control = {{c_mean}},
        Mean_Intervention = {{i_mean}},
        SD_Control = {{c_sd}},
        SD_Intervention = {{i_sd}}
      )
  } else if (measure == "IRR") {
    data <- data |>
      dplyr::rename(
        Author = {{studyvar}},
        Year = {{year}},
        Time_Control = {{c_time}},
        Time_Intervention = {{i_time}},
        Event_Control = {{c_event}},
        Event_Intervention = {{i_event}}
      )
  }

  # Handle subgroup variable renaming
  if (isTRUE(subgroup)) {
    if (!rlang::quo_is_null(rlang::enquo(subgroup_var))) {
      subgroup_var_name <- rlang::as_name(rlang::ensym(subgroup_var))
      if (subgroup_var_name != "Subgroup") {
        data <- data |>
          dplyr::rename(Subgroup = {{subgroup_var}})
      }
    }
  }

  # ---- Compute yi and sei ----
  if (is_bayesma) {
    es <- model$meta$es
    if (!is.null(es)) {
      data$yi  <- es$yi
      data$sei <- es$sei
    } else {
      fdf <- model$forest_df |> dplyr::filter(.data$type == "study")
      data$yi  <- fdf$estimate
      data$sei <- (fdf$upper - fdf$lower) / (2 * 1.96)
    }
  } else {
    # brmsfit: yi and sei should already be in data, or compute from raw columns
    if (!"yi" %in% names(data)) {
      if (measure %in% c("OR", "RR")) {
        a <- data$Event_Intervention
        b <- data$N_Intervention - data$Event_Intervention
        c <- data$Event_Control
        d <- data$N_Control - data$Event_Control
        cc <- dplyr::if_else(a == 0 | b == 0 | c == 0 | d == 0, 0.5, 0)
        data$yi  <- log(((a + cc) * (d + cc)) / ((b + cc) * (c + cc)))
        data$sei <- sqrt(1 / (a + cc) + 1 / (b + cc) + 1 / (c + cc) + 1 / (d + cc))
      } else if (measure %in% c("MD", "SMD")) {
        data$yi  <- data$Mean_Intervention - data$Mean_Control
        sp <- sqrt(((data$N_Control - 1) * data$SD_Control^2 +
                      (data$N_Intervention - 1) * data$SD_Intervention^2) /
                     (data$N_Control + data$N_Intervention - 2))
        data$sei <- sp * sqrt(1 / data$N_Control + 1 / data$N_Intervention)
      } else if (measure == "IRR") {
        data$yi  <- log((data$Event_Intervention / data$Time_Intervention) /
                          (data$Event_Control / data$Time_Control))
        data$sei <- sqrt(1 / data$Event_Intervention + 1 / data$Event_Control)
      }
    }
    if (!"sei" %in% names(data) && "vi" %in% names(data)) {
      data$sei <- sqrt(data$vi)
    }
  }

  # ---- Extract pooled estimate ----
  fixef_summary <- extract_fixef(model)
  pooled_estimate <- fixef_summary[1, 1]  # Posterior median of mu

  # ---- Prepare plot data ----
  plot_data <- data |>
    dplyr::mutate(
      x_val = if (isTRUE(props$log_scale)) exp(yi) else yi,
      Author_label = dplyr::if_else(
        is.na(Year),
        as.character(Author),
        paste0(Author, " (", Year, ")")
      )
    )

  # ---- Axis limits ----
  if (isTRUE(props$log_scale)) {
    if (is.null(xlim)) {
      x_min <- min(exp(data$yi - 2 * data$sei), na.rm = TRUE) * 0.8
      x_max <- max(exp(data$yi + 2 * data$sei), na.rm = TRUE) * 1.2
      xlim <- c(max(0.05, x_min), x_max)
    }
  } else {
    if (is.null(xlim)) {
      x_range <- max(abs(data$yi) + 2 * data$sei, na.rm = TRUE)
      xlim <- c(-x_range, x_range)
    }
  }

  breaks <- if (!is.null(x_breaks)) x_breaks else ggplot2::waiver()

  se_max <- max(data$sei, na.rm = TRUE) * 1.1

  # ---- Build funnel contours ----
  # Confidence interval contour regions centred on the pooled estimate
  contour_data <- NULL
  if (isTRUE(contour_lines)) {
    z_vals <- c(qnorm(0.995), qnorm(0.975), qnorm(0.95))  # 2.576, 1.96, 1.645
    ci_labels <- c("99% CI", "95% CI", "90% CI")

    se_seq <- seq(0, se_max, length.out = 300)

    contour_data <- purrr::map(seq_along(z_vals), function(i) {
      if (isTRUE(props$log_scale)) {
        lower <- exp(pooled_estimate - z_vals[i] * se_seq)
        upper <- exp(pooled_estimate + z_vals[i] * se_seq)
      } else {
        lower <- pooled_estimate - z_vals[i] * se_seq
        upper <- pooled_estimate + z_vals[i] * se_seq
      }

      # Build a polygon: go down the left edge, then back up the right edge
      tibble::tibble(
        x = c(lower, rev(upper)),
        y = c(se_seq, rev(se_seq)),
        level = ci_labels[i]
      )
    }) |>
      purrr::list_rbind() |>
      dplyr::mutate(level = factor(level, levels = ci_labels))
  }

  # ---- Build the plot ----
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_val, y = sei))

  # Contour shading (widest first so narrower regions overlay on top)
  if (isTRUE(contour_lines) && !is.null(contour_data)) {
    alpha_steps <- c(contour_alpha, contour_alpha * 1.5, contour_alpha * 2.5)

    for (i in seq_along(ci_labels)) {
      cd <- contour_data |> dplyr::filter(level == ci_labels[i])
      p <- p +
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
    border_data <- purrr::map(seq_along(z_vals), function(i) {
      if (isTRUE(props$log_scale)) {
        lower <- exp(pooled_estimate - z_vals[i] * se_seq)
        upper <- exp(pooled_estimate + z_vals[i] * se_seq)
      } else {
        lower <- pooled_estimate - z_vals[i] * se_seq
        upper <- pooled_estimate + z_vals[i] * se_seq
      }
      dplyr::bind_rows(
        tibble::tibble(x = lower, y = se_seq, side = "lower", level = ci_labels[i]),
        tibble::tibble(x = upper, y = se_seq, side = "upper", level = ci_labels[i])
      )
    }) |> purrr::list_rbind()

    p <- p +
      ggplot2::geom_line(
        data = border_data,
        ggplot2::aes(x = x, y = y, group = interaction(level, side)),
        color = color_contour, linetype = "dashed", linewidth = 0.4,
        inherit.aes = FALSE
      )
  }

  # Null line
  if (isTRUE(add_null_line)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = null_value,
        color = color_null, linewidth = 0.75
      )
  }

  # Pooled estimate line
  if (pooled_line != "none") {
    pooled_x <- if (isTRUE(props$log_scale)) exp(pooled_estimate) else pooled_estimate
    p <- p +
      ggplot2::geom_vline(
        xintercept = pooled_x,
        color = color_pooled, linewidth = 0.75, linetype = "solid"
      )
  }

  # Study points
  if (isTRUE(subgroup) && "Subgroup" %in% names(plot_data)) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(fill = Subgroup),
        shape = 21, size = point_size,
        color = color_points_outline, stroke = 0.6
      ) +
      ggplot2::scale_fill_brewer(palette = "Set2", name = "Subgroup")
  } else {
    p <- p +
      ggplot2::geom_point(
        fill = color_points, color = color_points_outline,
        shape = 21, size = point_size, stroke = 0.6
      )
  }

  # Study labels
  if (isTRUE(label_studies)) {
    p <- p +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = Author_label),
        size = label_size, color = "grey30",
        max.overlaps = 20, family = font
      )
  }

  # ---- Reverse y-axis (SE: 0 at top, max at bottom) ----
  p <- p +
    ggplot2::scale_y_reverse(
      expand = c(0, 0),
      limits = c(se_max, 0)
    )

  # ---- X-axis scale ----
  if (isTRUE(props$log_scale)) {
    p <- p +
      ggplot2::scale_x_log10(
        breaks = breaks,
        expand = c(0, 0),
        limits = xlim
      )
  } else {
    p <- p +
      ggplot2::scale_x_continuous(
        breaks = breaks,
        expand = c(0, 0),
        limits = xlim
      )
  }

  # ---- Theme ----
  p <- p +
    ggplot2::theme_light() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.75),
      axis.line.y.left = ggplot2::element_line(color = "black", linewidth = 0.75),
      axis.text = ggplot2::element_text(color = "black", family = font),
      axis.title = ggplot2::element_text(family = font),
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", family = font),
      legend.text = ggplot2::element_text(family = font)
    ) +
    ggplot2::labs(
      x = props$x_label,
      y = "Standard Error"
    )

  # ---- Title / subtitle ----
  if (!is.null(title) || !is.null(subtitle)) {
    hjust_val <- switch(title_align,
                        "left" = 0,
                        "center" = 0.5,
                        "centre" = 0.5,
                        "right" = 1,
                        0)

    if (!is.null(title)) {
      p <- p + ggplot2::labs(title = title) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 16, face = "bold", hjust = hjust_val,
            margin = ggplot2::margin(b = if (is.null(subtitle)) 10 else 5),
            family = font
          )
        )
    }
    if (!is.null(subtitle)) {
      p <- p + ggplot2::labs(subtitle = subtitle) +
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(
            size = 14, hjust = hjust_val,
            margin = ggplot2::margin(b = 10),
            color = "gray30", family = font
          )
        )
    }
  }

  return(p)
}
