#' ECDF Plot Comparing Model Strategies
#'
#' @description
#' Creates an empirical cumulative distribution function (ECDF) plot comparing
#' posterior distributions across different model strategies (e.g., random effects,
#' selection models, PET-PEESE) under a single prior specification.
#'
#' Use this plot to assess how robust conclusions are across different
#' modelling assumptions.
#'
#' @param model A fitted \code{bayesma} object.
#' @param data A data frame containing the study data used to fit the model.
#' @param prior A list containing the prior specification with elements:
#'   \describe{
#'     \item{\code{mu_prior}}{(Required) Prior specification for mu.}
#'     \item{\code{tau_prior}}{(Optional) Prior specification for tau.}
#'     \item{\code{name}}{Display label for the prior. Required when
#'       \code{incl_robma = TRUE} as it is used to match against priors in
#'       \code{model$robma_sensitivity}. If not provided and
#'       \code{incl_robma = FALSE}, defaults to "User-specified prior".}
#'   }
#' @param measure Effect measure string (e.g., "OR", "RR", "HR", "IRR", "MD", "SMD").
#' @param rob_var Optional. Name of the risk-of-bias variable (unquoted).
#' @param exclude_high_rob Logical. If TRUE, includes an "Excluding High RoB" section.
#' @param incl_common_effect Logical. Include common effect model. Default FALSE.
#' @param incl_random_effect Logical. Include random effects model. Default TRUE.
#' @param incl_bias_corrected Logical. Include bias-corrected model. Default FALSE.
#' @param incl_selection_copas Logical. Include Copas selection model. Default FALSE.
#' @param incl_selection_weight Logical. Include weight-function selection model.
#'   Default FALSE.
#' @param incl_pet_peese Logical. Include PET-PEESE model. Default FALSE.
#' @param incl_robust Logical. Include robust mixture model. Default FALSE.
#' @param incl_robma Logical. Include RoBMA model. Requires
#'   \code{model$robma_sensitivity} to be present and \code{prior$name} to match
#'   one of the priors used when fitting. Default FALSE.
#' @param prob_reference Character string. What reference to use for
#'   probability axis labels. Either \code{"null"} (compare to
#'   \code{null_value}) or \code{"null_range"} (compare to
#'   \code{null_range} boundaries). Default is \code{"null"}.
#' @param null_value Null hypothesis value. If NULL, uses measure default.
#' @param null_range Numeric vector of length 2 giving null range bounds.
#' @param add_null_range Logical. If TRUE and \code{null_range} is NULL,
#'   uses measure-appropriate defaults.
#' @param color_null_range Fill colour for the null range band.
#'   Default \code{"#77bb41"}.
#' @param label_control Label for control group. Default \code{"Control"}.
#' @param label_intervention Label for intervention group. Default
#'   \code{"Intervention"}.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle. If NULL and a prior name is available,
#'   displays the prior name.
#' @param xlim Optional x-axis limits.
#' @param x_breaks Optional x-axis break points.
#' @param color_palette Optional named vector of colours for model sections.
#' @param show_density Logical. If TRUE, includes a density plot
#'   below the ECDF. Default is TRUE.
#' @param font Optional font family.
#'
#' @return A ggplot object (or patchwork object if \code{show_density = TRUE}).
#'
#' @details
#' The plot shows one ECDF line per model strategy, all using the same prior.
#' This allows direct comparison of how different modelling assumptions
#' (random effects, selection models, bias correction, etc.) affect the
#' posterior distribution.
#'
#' The left y-axis shows P(effect < x) and the right y-axis shows P(effect > x).
#'
#' @seealso \code{\link{ecdf_prior_plot}} for comparing priors within a model type.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' priors <- list(
#'   vague = list(
#'     name = "Vague",
#'     mu_prior = normal(0, 10),
#'     tau_prior = half_cauchy(0, 1)
#'   ),
#'   weak_reg = list(
#'     name = "Weakly Regularising",
#'     mu_prior = normal(0, 1),
#'     tau_prior = half_cauchy(0, 0.5)
#'   )
#' )
#'
#' ecdf_model_plot(
#'   model = model,
#'   data = dat,
#'   measure = "OR",
#'   prior = priors$weak_reg,
#'   incl_random_effect = TRUE,
#'   incl_pet_peese = TRUE
#' )
#'
#' ecdf_model_plot(
#'   model = model,
#'   data = dat,
#'   measure = "OR",
#'   prior = priors$weak_reg,
#'   incl_robma = TRUE
#' )
#' }
#'
#' @importFrom rlang enquo quo_is_null %||%
#' @importFrom dplyr filter mutate bind_rows select
#' @importFrom purrr keep map_dfr map_chr
#' @importFrom ggplot2 ggplot aes annotate geom_vline stat_ecdf
#'   scale_y_continuous sec_axis scale_color_manual theme_light theme
#'   element_rect element_text element_blank element_line unit margin labs
#'   guides guide_legend scale_x_log10 scale_x_continuous coord_cartesian
#'   expansion waiver annotation_custom
#' @importFrom ggdist stat_slab
#' @importFrom scales alpha percent_format
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#' @importFrom grid textGrob gpar
#' @importFrom patchwork plot_spacer plot_layout
#' @importFrom tibble tibble
#'
ecdf_model_plot <- function(model,
                            data,
                            prior,
                            measure,
                            rob_var = NULL,
                            exclude_high_rob = FALSE,
                            incl_common_effect = FALSE,
                            incl_random_effect = TRUE,
                            incl_bias_corrected = FALSE,
                            incl_selection_copas = FALSE,
                            incl_selection_weight = FALSE,
                            incl_pet_peese = FALSE,
                            incl_robust = FALSE,
                            incl_robma = FALSE,
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
    measure = measure,
    priors  = prior
  )

  if (!prob_reference %in% c("null", "null_range")) {
    cli::cli_abort("{.arg prob_reference} must be one of {.val null} or {.val null_range}.")
  }

  # Validate prior structure
  if (!is.list(prior)) {
    cli::cli_abort("{.arg prior} must be a list.")
  }

  mu_prior <- prior$mu_prior
  tau_prior <- prior$tau_prior
  prior_label <- prior$name %||% "User-specified prior"

  if (is.null(mu_prior)) {
    cli::cli_abort("{.arg prior} must include {.val mu_prior}.")
  }

  # Validate RoBMA requirements
  if (isTRUE(incl_robma)) {
    if (!has_robma_sensitivity(model)) {
      cli::cli_abort(c(
        "RoBMA sensitivity fits not found on model.",
        "i" = "Use {.fn run_robma_sensitivity} to pre-compute RoBMA fits, then",
        "i" = "attach them with {.fn attach_robma_sensitivity}."
      ))
    }

    prior_name <- prior$name

    if (is.null(prior_name)) {
      cli::cli_abort("{.arg prior} must include {.val name} when {.arg incl_robma} is TRUE.")
    }

    robma_sens <- model$robma_sensitivity

    prior_id <- purrr::detect(
      names(robma_sens),
      ~ identical(robma_sens[[.x]]$name, prior_name)
    )

    if (is.null(prior_id)) {
      available_names <- purrr::map_chr(robma_sens, "name")
      cli::cli_abort(c(
        "Prior {.val {prior_name}} not found in {.code model$robma_sensitivity}.",
        "i" = "Available priors: {.val {available_names}}"
      ))
    }

    robma_fit <- robma_sens[[prior_id]]
  }

  props      <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value

  # ---------------------------
  # 2. Null range handling
  # ---------------------------
  if (is.null(null_range) && isTRUE(add_null_range)) {
    null_range <- switch(
      measure,
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
  # 3. Define sections (model strategies)
  # ---------------------------
  stage <- model$meta$stage %||% model$meta$call_args$stage

  stage_label <- switch(stage, one_stage = "One-Stage", two_stage = "Two-Stage", stage)

  rob_var_q <- rlang::enquo(rob_var)
  rob_data <- NULL
  if (isTRUE(exclude_high_rob)) {
    if (rlang::quo_is_null(rob_var_q)) {
      cli::cli_warn(c(
        "{.arg exclude_high_rob} is TRUE but {.arg rob_var} is not specified.",
        "i" = "The 'Excluding High RoB' section will not be included."
      ))
    } else {
      rob_data <- dplyr::filter(data, !!rob_var_q != "High" | is.na(!!rob_var_q))
    }
  }

  base_model_type <- model$meta$call_args$model_type %||% "random_effect"

  sections <- list(
    list(
      id            = "common_effect",
      include       = isTRUE(incl_common_effect),
      label         = paste0("Common Effect (", stage_label, ")"),
      overrides     = list(model_type = "common_effect"),
      data_override = NULL
    ),
    list(
      id            = "random_effect",
      include       = isTRUE(incl_random_effect),
      label         = paste0("Random Effects (", stage_label, ")"),
      overrides     = list(model_type = "random_effect"),
      data_override = NULL
    ),
    list(
      id            = "excluding_high_rob",
      include       = isTRUE(exclude_high_rob) && !is.null(rob_data),
      label         = "Excluding High RoB",
      overrides     = list(model_type = "random_effect"),
      data_override = rob_data
    ),
    list(
      id            = "bias_adjusted",
      include       = isTRUE(incl_bias_corrected),
      label         = "Bias Adjusted (Jung)",
      overrides     = list(model_type = "bias_corrected", stage = "two_stage"),
      data_override = NULL
    ),
    list(
      id            = "selection_copas",
      include       = isTRUE(incl_selection_copas),
      label         = "Selection Model (Copas)",
      overrides     = list(model_type = "selection_copas", stage = "two_stage"),
      data_override = NULL
    ),
    list(
      id            = "selection_weight",
      include       = isTRUE(incl_selection_weight),
      label         = "Selection Model (Weight)",
      overrides     = list(model_type = "selection_weight", stage = "two_stage"),
      data_override = NULL
    ),
    list(
      id            = "pet_peese",
      include       = isTRUE(incl_pet_peese),
      label         = "PET-PEESE",
      overrides     = list(model_type = "pet_peese", stage = "two_stage"),
      data_override = NULL
    ),
    list(
      id            = "robust",
      include       = isTRUE(incl_robust),
      label         = "Robust Mixture Model",
      overrides     = list(model_type = base_model_type, robust = TRUE),
      data_override = NULL
    )
  )

  active_sections <- purrr::keep(sections, ~ isTRUE(.x$include))

  if (!length(active_sections) && !isTRUE(incl_robma)) {
    cli::cli_abort("No model sections selected. Set at least one {.arg incl_*} argument to TRUE.")
  }

  # ---------------------------
  # 4. Extract draws for each section

  # ---------------------------
  transform_mu <- function(mu_raw) {
    if (measure %in% c("OR", "RR", "HR", "IRR")) exp(mu_raw) else mu_raw
  }

  priors_match <- function(p1, p2) {
    if (is.null(p1) && is.null(p2)) return(TRUE)
    if (is.null(p1) || is.null(p2)) return(FALSE)
    if (!inherits(p1, "bayesma_prior") || !inherits(p2, "bayesma_prior")) return(FALSE)
    identical(unclass(p1), unclass(p2))
  }

  orig_mu_prior  <- model$meta$priors$mu %||% model$meta$call_args$mu_prior
  orig_tau_prior <- model$meta$priors$tau %||% model$meta$call_args$tau_prior
  orig_model_type <- model$meta$model_type %||% model$meta$call_args$model_type %||% "random_effect"
  orig_stage <- model$meta$stage %||% model$meta$call_args$stage %||% "one_stage"

  cli::cli_h3("Fitting models for ECDF plot")
  cli::cli_alert_info("Prior: {prior_label}")

  draws_list <- purrr::map(active_sections, function(sec) {
    # Check if we can reuse original draws
    mu_matches  <- priors_match(mu_prior, orig_mu_prior)
    tau_matches <- is.null(tau_prior) || priors_match(tau_prior, orig_tau_prior)

    sec_model_type <- sec$overrides$model_type %||% orig_model_type
    sec_stage <- sec$overrides$stage %||% orig_stage

    model_matches <- (sec_model_type == orig_model_type) && (sec_stage == orig_stage)
    can_reuse <- mu_matches && tau_matches && model_matches && is.null(sec$data_override)

    if (can_reuse) {
      cli::cli_alert_success("{sec$label}: reusing original draws")
      mu_raw <- as.numeric(model$draws[["mu"]])
      x <- transform_mu(mu_raw)
      return(tibble::tibble(
        section_label = sec$label,
        x             = x
      ))
    }

    cli::cli_alert_info("{sec$label}: fitting model...")

    fit <- tryCatch(
      {
        refit_args <- c(
          list(
            model     = model,
            data      = sec$data_override %||% data,
            mu_prior  = mu_prior,
            tau_prior = tau_prior
          ),
          sec$overrides
        )
        do.call(refit_bayesma_update, refit_args)
      },
      error = function(e) {
        cli::cli_warn(c(
          "Failed to refit {.val {sec$label}}.",
          "i" = e$message
        ))
        NULL
      }
    )

    if (is.null(fit)) {
      return(tibble::tibble(section_label = character(), x = numeric()))
    }

    mu_raw <- as.numeric(fit$draws[["mu"]])
    x <- transform_mu(mu_raw)

    tibble::tibble(
      section_label = sec$label,
      x             = x
    )
  })

  draws_ma <- dplyr::bind_rows(draws_list)

  # ---------------------------
  # 5. RoBMA draws (if requested)
  # ---------------------------
  draws_robma <- NULL

  if (isTRUE(incl_robma)) {
    robma_sens <- model$robma_sensitivity

    if (prior %in% names(robma_sens)) {
      robma_fit <- robma_sens[[prior]]

      if (!is.null(robma_fit)) {
        cli::cli_alert_info("RoBMA (Model-Averaged): extracting draws")

        all_draws <- robma_to_sensitivity_draws(
          robma_fit     = robma_fit,
          measure       = measure,
          prior         = prior,
          prior_label   = prior_label,
          section_label = "RoBMA (Model-Averaged)"
        )

        draws_robma <- all_draws |>
          dplyr::select(x, section_label)
      }
    } else {
      cli::cli_warn("Prior {.val {prior}} not found in pre-computed RoBMA fits.")
    }
  }

  draws <- dplyr::bind_rows(draws_ma, draws_robma)

  if (nrow(draws) == 0) {
    cli::cli_abort(c(
      "No draws available for ECDF plot.",
      "i" = "Check that model refits succeeded."
    ))
  }

  # ---------------------------
  # 6. Section ordering
  # ---------------------------
  section_order <- c(
    "Common Effect",
    "Random Effects",
    "Excluding High RoB",
    "Bias Adjusted (Jung)",
    "Selection Model (Copas)",
    "Selection Model (Weight)",
    "PET-PEESE",
    "Robust Mixture Model",
    "RoBMA (Conditional)",
    "RoBMA (Model-Averaged)"
  )

  .match_section_order <- function(labels, order_vec) {
    ordered_out <- character()
    for (template in order_vec) {
      matching <- labels[startsWith(labels, template) | labels == template]
      ordered_out <- c(ordered_out, matching)
    }
    remaining <- setdiff(labels, ordered_out)
    unique(c(ordered_out, remaining))
  }

  actual_sections <- unique(draws$section_label)
  present_sections <- .match_section_order(actual_sections, section_order)

  draws <- draws |>
    dplyr::mutate(section_label = factor(section_label, levels = present_sections))

  # ---------------------------
  # 7. Colors
  # ---------------------------
  if (is.null(color_palette)) {
    n_sections <- length(present_sections)
    n_colors   <- max(3, min(9, n_sections))
    color_values <- RColorBrewer::brewer.pal(n_colors, "Set1")[seq_len(n_sections)]
    color_palette <- stats::setNames(color_values, present_sections)
  } else {
    missing_sections <- setdiff(present_sections, names(color_palette))
    if (length(missing_sections) > 0) {
      cli::cli_warn(c(
        "Some sections missing from {.arg color_palette}: {.val {missing_sections}}.",
        "i" = "Using default colours for these."
      ))
      extra_colors <- RColorBrewer::brewer.pal(max(3, length(missing_sections)), "Set2")
      color_palette <- c(
        color_palette,
        stats::setNames(extra_colors[seq_along(missing_sections)], missing_sections)
      )
    }
  }

  # ---------------------------
  # 8. Axis limits and labels
  # ---------------------------
  calc_xlim <- if (!is.null(xlim)) xlim else range(draws$x, na.rm = TRUE)
  breaks    <- x_breaks %||% ggplot2::waiver()

  if (prob_reference == "null") {
    y_left_label  <- paste0("Probability ", measure, " < ", null_value)
    y_right_label <- paste0("Probability ", measure, " > ", null_value)
  } else {
    y_left_label  <- paste0("Probability ", measure, " < ", null_range[1])
    y_right_label <- paste0("Probability ", measure, " > ", null_range[2])
  }

  # Default subtitle shows prior
  if (is.null(subtitle)) {
    subtitle <- paste0("Prior: ", prior_label)
  }

  # ---------------------------
  # 9. Build ECDF plot
  # ---------------------------
  p <- ggplot2::ggplot(draws, ggplot2::aes(x = x, color = section_label)) +
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
      limits = present_sections,
      name   = "Model"
    ) +
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
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = "solid", shape = NA)))

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
  # 10. Density panel (optional)
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

    density_p <- ggplot2::ggplot(draws, ggplot2::aes(x = x, color = section_label)) +
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
      ggplot2::scale_color_manual(values = color_palette, limits = present_sections, name = "Model") +
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
