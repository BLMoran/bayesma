#' ECDF Plot Comparing Prior Sensitivity
#'
#' @description
#' Creates an empirical cumulative distribution function (ECDF) plot comparing
#' posterior distributions across different prior specifications for a single
#' model type (or small set of model types).
#'
#' Use this plot to assess how sensitive conclusions are to the choice of
#' prior distributions.
#'
#' @param model A fitted \code{bayesma} object.
#' @param data A data frame containing the study data used to fit the model.
#' @param priors A named list of prior specifications. Each element must be a
#'   list with at least \code{mu_prior} and optionally \code{tau_prior}, and
#'   may include \code{name} for display labels.
#' @param estimand Effect estimand string (e.g., "OR", "RR", "HR", "IRR", "MD", "SMD").
#' @param model_types Character vector specifying which model type(s) to include.
#'   Maximum of 2 model types to avoid clutter. Valid values: \code{"common_effect"},
#'   \code{"random_effect"}, \code{"bias_corrected"}, \code{"selection_copas"},
#'   \code{"selection_weight"}, \code{"pet_peese"}, \code{"robust"}, \code{"robma"}.
#'   Default is \code{"random_effect"}.
#' @param prior_order Optional character vector specifying the display order of
#'   priors. Should contain the names (IDs) of the priors in the desired order.
#'   If NULL, priors are displayed in the order they appear in \code{priors}.
#' @param prob_reference Character string. What reference to use for
#'   probability axis labels. Either \code{"null"} (compare to
#'   \code{null_value}) or \code{"null_range"} (compare to
#'   \code{null_range} boundaries). Default is \code{"null"}.
#' @param null_value Null hypothesis value. If NULL, uses estimand default.
#' @param null_range Numeric vector of length 2 giving null range bounds.
#' @param add_null_range Logical. If TRUE and \code{null_range} is NULL,
#'   uses estimand-appropriate defaults.
#' @param color_null_range Fill colour for the null range band.
#'   Default \code{"#77bb41"}.
#' @param label_control Label for control group. Default \code{"Control"}.
#' @param label_intervention Label for intervention group. Default
#'   \code{"Intervention"}.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle. If NULL and a single model type is
#'   selected, displays the model type.
#' @param xlim Optional x-axis limits.
#' @param x_breaks Optional x-axis break points.
#' @param color_palette Optional named vector of colours for priors.
#' @param linetype_by_model Logical. If TRUE and multiple model types are
#'   selected, uses different linetypes for each model type. Default TRUE.
#' @param show_density Logical. If TRUE, includes a density plot
#'   below the ECDF. Default is TRUE.
#' @param font Optional font family.
#'
#' @return A ggplot object (or patchwork object if \code{show_density = TRUE}).
#'
#' @details
#' The plot shows one ECDF line per prior specification (and per model type if
#' multiple are selected). This allows direct comparison of how different
#' prior choices affect the posterior distribution.
#'
#' When two model types are selected, colour represents prior and linetype
#' represents model type, allowing comparison of both dimensions simultaneously.
#'
#' The left y-axis shows P(effect < x) and the right y-axis shows P(effect > x).
#'
#' @seealso \code{\link{ecdf_model_plot}} for comparing models within a prior.
#'
#' @export
#'
#' @importFrom rlang %||%
#' @importFrom dplyr filter mutate bind_rows select
#' @importFrom purrr map map_dfr map_chr imap_dfr
#' @importFrom ggplot2 ggplot aes annotate geom_vline stat_ecdf
#'   scale_y_continuous sec_axis scale_color_manual scale_linetype_manual
#'   theme_light theme element_rect element_text element_blank element_line
#'   unit margin labs guides guide_legend scale_x_log10 scale_x_continuous
#'   coord_cartesian expansion waiver annotation_custom
#' @importFrom ggdist stat_slab
#' @importFrom scales alpha percent_format
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#' @importFrom grid textGrob gpar
#' @importFrom patchwork plot_spacer plot_layout
#' @importFrom tibble tibble
#'
ecdf_prior_plot <- function(model,
                            data,
                            priors,
                            estimand,
                            model_types = "random_effect",
                            prior_order = NULL,
                            prob_reference = "null",
                            null_value = NULL,
                            null_range = NULL,
                            add_null_range = FALSE,
                            color_null_range = "#77bb41",
                            label_control = "Control",
                            label_intervention = "Intervention",
                            title = NULL,
                            subtitle = NULL,
                            xlim = NULL,
                            x_breaks = NULL,
                            color_palette = NULL,
                            linetype_by_model = TRUE,
                            show_density = TRUE,
                            font = NULL) {

  # ---------------------------
  # 1. Validation
  # ---------------------------
  if (!inherits(model, "bayesma")) {
    cli::cli_abort("This function requires a {.cls bayesma} model.")
  }


  validate_inputs_sens_plot(
    model   = model,
    data    = data,
    estimand = estimand,
    priors  = priors
  )

  if (!prob_reference %in% c("null", "null_range")) {
    cli::cli_abort("{.arg prob_reference} must be one of {.val null} or {.val null_range}.")
  }

  # Validate priors structure
  if (!is.list(priors) || length(priors) == 0) {
    cli::cli_abort("{.arg priors} must be a non-empty list.")
  }

  prior_ids <- names(priors)
  if (is.null(prior_ids) || any(prior_ids == "")) {
    cli::cli_abort("{.arg priors} must be a named list (each prior must have an ID).")
  }

  # Handle prior_order
  if (!is.null(prior_order)) {
    prior_order <- as.character(prior_order)
    missing_priors <- setdiff(prior_order, prior_ids)
    if (length(missing_priors) > 0) {
      cli::cli_abort(c(
        "Some priors in {.arg prior_order} not found in {.arg priors}:",
        "x" = paste(missing_priors, collapse = ", ")
      ))
    }
    extra_priors <- setdiff(prior_ids, prior_order)
    if (length(extra_priors) > 0) {
      cli::cli_warn(c(
        "Some priors not specified in {.arg prior_order} will be appended:",
        "i" = paste(extra_priors, collapse = ", ")
      ))
      prior_order <- c(prior_order, extra_priors)
    }
    prior_ids <- prior_order
  }

  # Build prior label map
  prior_label_map <- purrr::map_chr(prior_ids, function(pid) {
    priors[[pid]]$name %||% pid
  })
  names(prior_label_map) <- prior_ids

  # Validate model_types
  valid_model_types <- c(
    "common_effect", "random_effect", "bias_corrected",
    "selection_copas", "selection_weight", "pet_peese", "robust", "robma"
  )

  invalid_types <- setdiff(model_types, valid_model_types)
  if (length(invalid_types) > 0) {
    cli::cli_abort(c(
      "Invalid {.arg model_types}: {.val {invalid_types}}",
      "i" = "Valid values: {.val {valid_model_types}}"
    ))
  }

  if (length(model_types) > 2) {
    cli::cli_abort(c(
      "Maximum of 2 model types allowed to avoid clutter.",
      "i" = "You specified {length(model_types)}: {.val {model_types}}"
    ))
  }

  # Check RoBMA if requested
  if ("robma" %in% model_types && !has_robma_sensitivity(model)) {
    cli::cli_abort(c(
      "RoBMA requested but sensitivity fits not found on model.",
      "i" = "Use {.fn run_robma_sensitivity} and {.fn attach_robma_sensitivity} first."
    ))
  }

  props      <- get_measure_properties(estimand)
  null_value <- null_value %||% props$null_value

  # ---------------------------
  # 2. Null range handling
  # ---------------------------
  if (is.null(null_range) && isTRUE(add_null_range)) {
    null_range <- switch(
      estimand,
      OR  = c(0.9, 1.1),
      RR  = c(0.9, 1.1),
      HR  = c(0.9, 1.1),
      IRR = c(0.9, 1.1),
      SMD = c(-0.1, 0.1),
      cli::cli_abort("For MD, {.arg null_range} must be supplied.")
    )
  } else if (!is.null(null_range)) {
    null_range <- if (length(null_range) == 1) {
      c(null_value - null_range, null_value + null_range)
    } else {
      sort(null_range)
    }
  }

  # ---------------------------
  # 3. Model type labels and overrides
  # ---------------------------
  stage <- model$meta$stage %||% model$meta$call_args$stage
  stage_label <- switch(stage, one_stage = "One-Stage", two_stage = "Two-Stage", stage)
  base_model_type <- model$meta$call_args$model_type %||% "random_effect"

  model_type_config <- list(
    common_effect = list(
      label = paste0("Common Effect (", stage_label, ")"),
      overrides = list(model_type = "common_effect")
    ),
    random_effect = list(
      label = paste0("Random Effects (", stage_label, ")"),
      overrides = list(model_type = "random_effect")
    ),
    bias_corrected = list(
      label = "Bias Adjusted (Jung)",
      overrides = list(model_type = "bias_corrected", stage = "two_stage")
    ),
    selection_copas = list(
      label = "Selection Model (Copas)",
      overrides = list(model_type = "selection_copas", stage = "two_stage")
    ),
    selection_weight = list(
      label = "Selection Model (Weight)",
      overrides = list(model_type = "selection_weight", stage = "two_stage")
    ),
    pet_peese = list(
      label = "PET-PEESE",
      overrides = list(model_type = "pet_peese", stage = "two_stage")
    ),
    robust = list(
      label = "Robust Mixture Model",
      overrides = list(model_type = base_model_type, robust = TRUE)
    ),
    robma = list(
      label = "RoBMA (Model-Averaged)",
      overrides = NULL
    )
  )

  # ---------------------------
  # 4. Extract draws for each prior × model_type
  # ---------------------------
  transform_mu <- function(mu_raw) {
    if (estimand %in% c("OR", "RR", "HR", "IRR")) exp(mu_raw) else mu_raw
  }

  extract_mu_draws <- function(m) {
    if (is_marginal_estimand(estimand) && !is.null(m$marginal)) {
      m$marginal$draws
    } else {
      as.numeric(m$draws[["mu"]])
    }
  }

  priors_match <- function(p1, p2) {
    if (is.null(p1) && is.null(p2)) return(TRUE)
    if (is.null(p1) || is.null(p2)) return(FALSE)
    if (!inherits(p1, "bayesma_prior") || !inherits(p2, "bayesma_prior")) return(FALSE)
    identical(unclass(p1), unclass(p2))
  }

  orig_mu_prior   <- model$meta$priors$mu %||% model$meta$call_args$mu_prior
  orig_tau_prior  <- model$meta$priors$tau %||% model$meta$call_args$tau_prior
  orig_model_type <- model$meta$model_type %||% model$meta$call_args$model_type %||% "random_effect"
  orig_stage      <- model$meta$stage %||% model$meta$call_args$stage %||% "one_stage"
  orig_estimand   <- model$meta$call_args$estimand %||% "OR"

  cli::cli_h3("Fitting models for ECDF plot")
  cli::cli_alert_info("Model types: {.val {model_types}}")
  cli::cli_alert_info("Priors: {.val {prior_ids}}")

  draws_list <- list()

  for (mt in model_types) {
    config <- model_type_config[[mt]]
    model_label <- config$label

    if (mt == "robma") {
      # Handle RoBMA separately
      robma_sens <- model$robma_sensitivity

      for (pid in prior_ids) {
        if (pid %in% names(robma_sens) && !is.null(robma_sens[[pid]])) {
          cli::cli_alert_info("{model_label} + {pid}: extracting draws")

          all_draws <- robma_to_sensitivity_draws(
            robma_fit     = robma_sens[[pid]],
            estimand       = estimand,
            prior         = pid,
            prior_label   = prior_label_map[[pid]],
            section_label = model_label
          )

          draws_list[[length(draws_list) + 1]] <- all_draws |>
            dplyr::select(x, section_label) |>
            dplyr::mutate(
              prior       = pid,
              prior_label = prior_label_map[[pid]],
              model_type  = model_label
            )
        } else {
          cli::cli_warn("Prior {.val {pid}} not found in RoBMA fits.")
        }
      }
    } else {
      # Standard bayesma model types
      for (pid in prior_ids) {
        ps <- priors[[pid]]
        mu_prior <- ps$mu_prior
        tau_prior <- ps$tau_prior

        if (is.null(mu_prior)) {
          cli::cli_abort("Prior {.val {pid}} must include {.val mu_prior}.")
        }

        # Check reuse
        mu_matches  <- priors_match(mu_prior, orig_mu_prior)
        tau_matches <- is.null(tau_prior) || priors_match(tau_prior, orig_tau_prior)

        sec_model_type <- config$overrides$model_type %||% orig_model_type
        sec_stage <- config$overrides$stage %||% orig_stage

        model_matches    <- (sec_model_type == orig_model_type) && (sec_stage == orig_stage)
        estimand_matches <- identical(estimand, orig_estimand)
        can_reuse <- mu_matches && tau_matches && model_matches && estimand_matches

        if (can_reuse) {
          cli::cli_alert_success("{model_label} + {pid}: reusing original draws")
          mu_raw <- extract_mu_draws(model)
          x <- if (is_marginal_estimand(estimand)) mu_raw else transform_mu(mu_raw)

          draws_list[[length(draws_list) + 1]] <- tibble::tibble(
            x           = x,
            prior       = pid,
            prior_label = prior_label_map[[pid]],
            model_type  = model_label
          )
        } else {
          cli::cli_alert_info("{model_label} + {pid}: fitting model...")

          fit <- tryCatch(
            {
              refit_args <- c(
                list(
                  model     = model,
                  data      = data,
                  mu_prior  = mu_prior,
                  tau_prior = tau_prior,
                  estimand  = estimand
                ),
                config$overrides
              )
              do.call(refit_bayesma_update, refit_args)
            },
            error = function(e) {
              cli::cli_warn(c(
                "Failed to refit {.val {model_label}} + {.val {pid}}.",
                "i" = e$message
              ))
              NULL
            }
          )

          if (!is.null(fit)) {
            mu_raw <- extract_mu_draws(fit)
            x <- if (is_marginal_estimand(estimand)) mu_raw else transform_mu(mu_raw)

            draws_list[[length(draws_list) + 1]] <- tibble::tibble(
              x           = x,
              prior       = pid,
              prior_label = prior_label_map[[pid]],
              model_type  = model_label
            )
          }
        }
      }
    }
  }

  draws <- dplyr::bind_rows(draws_list)

  if (nrow(draws) == 0) {
    cli::cli_abort(c(
      "No draws available for ECDF plot.",
      "i" = "Check that model refits succeeded."
    ))
  }

  # ---------------------------
  # 5. Factor ordering
  # ---------------------------
  # Prior labels in specified order
  prior_labels_ordered <- prior_label_map[prior_ids]
  draws <- draws |>
    dplyr::mutate(prior_label = factor(prior_label, levels = prior_labels_ordered))

  # Model types in specified order
  model_labels_ordered <- purrr::map_chr(model_types, ~ model_type_config[[.x]]$label)
  draws <- draws |>
    dplyr::mutate(model_type = factor(model_type, levels = model_labels_ordered))

  # ---------------------------
  # 6. Colors and linetypes
  # ---------------------------
  n_priors <- length(prior_ids)

  if (is.null(color_palette)) {
    n_colors <- max(3, min(9, n_priors))
    color_values <- RColorBrewer::brewer.pal(n_colors, "Set1")[seq_len(n_priors)]
    color_palette <- stats::setNames(color_values, prior_labels_ordered)
  }

  use_linetype <- length(model_types) > 1 && isTRUE(linetype_by_model)
  linetype_values <- stats::setNames(c("solid", "dashed")[seq_along(model_labels_ordered)], model_labels_ordered)

  # ---------------------------
  # 7. Axis limits and labels
  # ---------------------------
  calc_xlim <- if (!is.null(xlim)) xlim else range(draws$x, na.rm = TRUE)
  breaks    <- x_breaks %||% ggplot2::waiver()

  if (prob_reference == "null") {
    y_left_label  <- paste0("Probability ", estimand, " < ", null_value)
    y_right_label <- paste0("Probability ", estimand, " > ", null_value)
  } else {
    y_left_label  <- paste0("Probability ", estimand, " < ", null_range[1])
    y_right_label <- paste0("Probability ", estimand, " > ", null_range[2])
  }

  # Default subtitle shows model type(s)
  if (is.null(subtitle)) {
    subtitle <- paste0("Model: ", paste(model_labels_ordered, collapse = " & "))
  }

  # ---------------------------
  # 8. Build ECDF plot
  # ---------------------------
  if (use_linetype) {
    p <- ggplot2::ggplot(draws, ggplot2::aes(x = x, color = prior_label, linetype = model_type))
  } else {
    p <- ggplot2::ggplot(draws, ggplot2::aes(x = x, color = prior_label))
  }

  p <- p +
    {if (isTRUE(add_null_range) && !is.null(null_range)) {
      ggplot2::annotate("rect",
                        xmin = null_range[1], xmax = null_range[2],
                        ymin = -Inf, ymax = Inf,
                        fill = scales::alpha(color_null_range, 0.2),
                        color = NA)
    }} +
    {if (isTRUE(add_null_range) && !is.null(null_range)) {
      ggplot2::geom_vline(xintercept = null_range, linetype = "dotted", color = "grey30")
    }} +
    ggplot2::stat_ecdf(geom = "step", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = null_value, linewidth = 0.8, color = "black") +
    ggplot2::scale_y_continuous(
      name   = y_left_label,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = scales::percent_format(),
      sec.axis = ggplot2::sec_axis(
        ~ 1 - .,
        name   = y_right_label,
        breaks = seq(0, 1, 0.1),
        labels = scales::percent_format()
      )
    ) +
    ggplot2::scale_color_manual(
      values = color_palette,
      limits = prior_labels_ordered,
      name   = "Prior"
    ) +
    {if (use_linetype) {
      ggplot2::scale_linetype_manual(
        values = linetype_values,
        limits = model_labels_ordered,
        name   = "Model"
      )
    }} +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position      = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.background    = ggplot2::element_rect(fill = "white", color = "grey60", linewidth = 0.3),
      legend.title         = ggplot2::element_text(family = font, face = "bold"),
      legend.text          = ggplot2::element_text(family = font),
      axis.title           = ggplot2::element_text(family = font),
      axis.text            = ggplot2::element_text(family = font),
      panel.grid.minor     = ggplot2::element_blank(),
      plot.title           = ggplot2::element_text(family = font, face = "bold"),
      plot.subtitle        = ggplot2::element_text(family = font),
      axis.line.x.bottom   = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.x         = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length.x  = ggplot2::unit(0.15, "cm")
    ) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linetype = "solid", shape = NA)),
      linetype = ggplot2::guide_legend()
    )

  # X-axis scale
  if (isTRUE(props$log_scale)) {
    p <- p +
      ggplot2::scale_x_log10(breaks = breaks, expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off") +
      ggplot2::labs(x = props$x_label)
  } else {
    p <- p +
      ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0)) +
      ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off") +
      ggplot2::labs(x = props$x_label)
  }

  # ---------------------------
  # 9. Density panel (optional)
  # ---------------------------
  if (isTRUE(show_density)) {
    p <- p +
      ggplot2::theme(
        axis.title.x       = ggplot2::element_blank(),
        axis.text.x        = ggplot2::element_blank(),
        axis.ticks.x       = ggplot2::element_blank(),
        axis.line.x.bottom = ggplot2::element_blank(),
        plot.margin        = ggplot2::margin(t = 5.5, r = 5.5, b = 2, l = 5.5, unit = "pt")
      )

    if (use_linetype) {
      density_p <- ggplot2::ggplot(draws, ggplot2::aes(x = x, color = prior_label, linetype = model_type))
    } else {
      density_p <- ggplot2::ggplot(draws, ggplot2::aes(x = x, color = prior_label))
    }

    density_p <- density_p +
      {if (isTRUE(add_null_range) && !is.null(null_range)) {
        ggplot2::annotate("rect",
                          xmin = null_range[1], xmax = null_range[2],
                          ymin = -Inf, ymax = Inf,
                          fill = scales::alpha(color_null_range, 0.2),
                          color = NA)
      }} +
      {if (isTRUE(add_null_range) && !is.null(null_range)) {
        ggplot2::geom_vline(xintercept = null_range, linetype = "dotted", color = "grey30")
      }} +
      ggdist::stat_slab(fill = NA, slab_linewidth = 0.9, normalize = "all") +
      ggplot2::geom_vline(xintercept = null_value, linewidth = 0.8, color = "black") +
      ggplot2::scale_color_manual(values = color_palette, limits = prior_labels_ordered, name = "Prior") +
      {if (use_linetype) {
        ggplot2::scale_linetype_manual(values = linetype_values, limits = model_labels_ordered, name = "Model")
      }} +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::theme_light() +
      ggplot2::theme(
        legend.position     = "none",
        axis.title.y        = ggplot2::element_blank(),
        axis.text.y         = ggplot2::element_blank(),
        axis.ticks.y        = ggplot2::element_blank(),
        axis.title.x        = ggplot2::element_text(family = font),
        axis.text.x         = ggplot2::element_text(family = font),
        panel.grid.minor    = ggplot2::element_blank(),
        axis.line.x.bottom  = ggplot2::element_line(linewidth = 0.8, color = "black"),
        axis.ticks.x        = ggplot2::element_line(linewidth = 0.8, color = "black"),
        axis.ticks.length.x = ggplot2::unit(0.15, "cm"),
        plot.margin         = ggplot2::margin(t = 2, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
      ) +
      ggplot2::labs(x = props$x_label)

    # X-axis scale for density
    if (isTRUE(props$log_scale)) {
      density_p <- density_p +
        ggplot2::scale_x_log10(breaks = breaks, expand = c(0, 0)) +
        ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off")
    } else {
      density_p <- density_p +
        ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0)) +
        ggplot2::coord_cartesian(xlim = calc_xlim, clip = "off")
    }

    # Favours labels
    density_p <- density_p +
      ggplot2::annotation_custom(
        grid::textGrob(
          label = paste0("Favours\n", label_intervention),
          x = grid::unit(0.05, "npc"), y = grid::unit(0.6, "npc"),
          just = c("left", "bottom"),
          gp = grid::gpar(col = "grey30", fontsize = 9, fontface = "bold.italic", fontfamily = font)
        ),
        xmin = calc_xlim[1] - 0.01, xmax = calc_xlim[2], ymin = -Inf, ymax = Inf
      ) +
      ggplot2::annotation_custom(
        grid::textGrob(
          label = paste0("Favours\n", label_control),
          x = grid::unit(0.97, "npc"), y = grid::unit(0.6, "npc"),
          just = c("right", "bottom"),
          gp = grid::gpar(col = "grey30", fontsize = 9, fontface = "bold.italic", fontfamily = font)
        ),
        xmin = calc_xlim[1], xmax = calc_xlim[2] - 0.01, ymin = -Inf, ymax = Inf
      )

    # Combine
    final_plot <- (p / patchwork::plot_spacer() / density_p) +
      patchwork::plot_layout(ncol = 1, heights = c(1, 0.0, 0.35), axes = "collect_x", axis_titles = "collect_x") &
      ggplot2::coord_cartesian(xlim = calc_xlim)

    return(final_plot)
  }

  p
}
