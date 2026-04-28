#' Create Bayesian Forest Plot for Meta-Analysis
#'
#' @description
#' This function creates a Bayesian forest plot for meta-analysis from a
#' bayesma model object.
#'
#' @param model A fitted bayesma object (class 'bayesma').
#' @param data A data frame containing the study data used for the meta-analysis.
#' @param estimand Character string specifying the effect measure or marginal estimand. Must be one of:
#'   \itemize{
#'     \item Relative-effect: \code{"OR"} (Odds Ratio), \code{"HR"} (Hazard Ratio),
#'       \code{"RR"} (Risk Ratio), \code{"IRR"} (Incidence Rate Ratio),
#'       \code{"MD"} (Mean Difference), \code{"SMD"} (Standardised Mean Difference).
#'     \item Marginal: \code{"RD"} / \code{"ARR"} (Risk Difference), \code{"ATE"}
#'       (Average Treatment Effect), \code{"ATT"} (Average Treatment Effect on the Treated),
#'       \code{"CATE"} (Conditional Average Treatment Effect).
#'   }
#' @param studyvar Column name containing study identifiers/authors. Default is NULL.
#'   Extracted from the model object if not provided.
#' @param year Column name containing publication years. Default is NULL.
#' @param c_n Column name containing control group sample sizes. Required for OR, RR, MD, SMD.
#' @param i_n Column name containing intervention group sample sizes. Required for OR, RR, MD, SMD.
#' @param c_event Column name containing control group event counts. Required for OR, RR, IRR.
#' @param i_event Column name containing intervention group event counts. Required for OR, RR, IRR.
#' @param c_mean Column name containing control group means. Required for MD, SMD.
#' @param i_mean Column name containing intervention group means. Required for MD, SMD.
#' @param c_sd Column name containing control group standard deviations. Required for MD, SMD.
#' @param i_sd Column name containing intervention group standard deviations. Required for MD, SMD.
#' @param c_time Column name containing control group time periods. Required for IRR.
#' @param i_time Column name containing intervention group time periods. Required for IRR.
#' @param sort_studies_by Character string specifying how to sort studies.
#'   Options: "author" (default), "year", or "effect".
#' @param subgroup Logical indicating whether to create subgroup analysis. Default is FALSE.
#' @param subgroup_var Character string. Name of the variable in data to use for subgroup analysis.
#' @param sort_subgroup_by Character string or vector specifying subgroup ordering.
#'   Options: "alphabetical" (default), "effect", or custom character vector of subgroup names.
#' @param label_outcome Character string for outcome label. Default is "Outcome".
#' @param label_control Character string for control group label. Default is "Control".
#' @param label_intervention Character string for intervention group label. Default is "Intervention".
#' @param title Character string for the plot title. Default is NULL (no title).
#' @param subtitle Character string for the plot subtitle. Default is NULL (no subtitle).
#' @param title_align Character string specifying title alignment. Options: "left" (default), "center"/"centre", "right".
#' @param xlim Numeric vector of length 2 specifying x-axis limits. Default is NULL (auto-scaled).
#' @param x_breaks Numeric vector specifying custom x-axis break points. Default is NULL (uses measure-specific defaults).
#' @param shrinkage_output Character string specifying shrinkage visualization.
#'   Options: "density" (default) or "pointinterval".
#' @param null_value Numeric value specifying x-axis value of the null line. Default is NULL (uses measure-specific defaults).
#' @param add_rope Logical indicating whether to add ROPE (Region of Practical Equivalence). Default is FALSE.
#' @param rope_value Numeric vector of length 2 specifying ROPE range, or single value for symmetric range around null.
#'   Default is NULL (uses Kruschke's recommendations: OR/HR/RR/IRR = c(0.9, 1.1), SMD = c(-0.1, 0.1), MD requires specification).
#' @param rope_color Color for ROPE shading. Default is transparent grey.
#' @param color_palette Character vector of colors for the plot. Default is NULL.
#' @param color_study_posterior_null_left Color for left side of null study posteriors. Default is "deepskyblue".
#' @param color_study_posterior_null_right Color for right side of null study posteriors. Default is "violet".
#' @param color_study_posterior Color for study posterior densities. Default is "dodgerblue".
#' @param color_study_posterior_outline Color for study posterior outlines. Default is "blue".
#' @param color_overall_posterior Color for overall posterior. Default is "blue".
#' @param color_shrinkage_pointinterval Color for shrinkage point intervals (used when
#'   \code{shrinkage_output = "pointinterval"}). Default is "purple".
#' @param color_shrinkage_outline Color for shrinkage plot outlines. Default is "purple".
#' @param color_shrinkage_fill Color for shrinkage plot fill. Default is NULL.
#' @param split_color_by_null Logical. If TRUE, posterior densities are split and
#' coloured based on whether values fall above or below the null value.
#' @param color_favours_control Colour used for density regions favouring the control group when \code{split_color_by_null = TRUE}.
#' @param color_favours_intervention Colour used for density regions favouring the
#' intervention group when \code{split_color_by_null = TRUE}.
#' @param add_arm_labels Logical indicating whether to display "Favours Control" /
#'   "Favours Intervention" labels above the density plot. Default is TRUE.
#' @param reverse_arm_labels Logical indicating whether to swap the positions of
#'   the "Favours Control" and "Favours Intervention" labels. Default is FALSE.
#' @param add_pred Logical indicating whether to add a prediction interval row
#'   beneath the Pooled Effect. Default is FALSE.
#' @param add_pred_subgroup Logical indicating whether to add prediction
#'   interval rows for each subgroup when \code{subgroup = TRUE}. If
#'   \code{FALSE} (the default), a prediction row is only added for the
#'   overall pooled effect. Ignored when \code{subgroup = FALSE}.
#' @param pred_output Character string specifying the visualisation for the
#'   prediction interval row. Options: "density" (default) or "pointinterval".
#' @param color_pred_posterior Color for the prediction interval density fill.
#'   Default is "forestgreen".
#' @param color_pred_outline Color for the prediction interval density outline.
#'   Default is "darkgreen".
#' @param color_pred_pointinterval Color for the prediction interval point
#'   interval. Default is "forestgreen".
#' @param plot_width Numeric value specifying the relative width of the plot component. Default is 4.
#' @param add_rob Logical indicating whether to add Risk of Bias assessment. Default is FALSE.
#' @param rob_tool Character string specifying RoB tool. Options: "rob2" (default).
#' @param add_rob_legend Logical indicating whether to add RoB legend. Default is FALSE.
#' @param exclude_high_rob Logical indicating whether to exclude high risk of bias studies
#'   and refit the model. Default is FALSE.
#' @param font Character string specifying the font family to use throughout the plot.
#'   Default is NULL (uses system defaults).
#'
#' @return A patchwork object containing the complete forest plot with study information table,
#'   density plots, and effect size table.
#'
#' @export
forest <- function(model,
                         data,
                         estimand,
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
                         sort_studies_by = "author",
                         subgroup = FALSE,
                         subgroup_var = NULL,
                         sort_subgroup_by = "alphabetical",
                         label_outcome = "Outcome",
                         label_control = "Control",
                         label_intervention = "Intervention",
                         title = NULL,
                         subtitle = NULL,
                         title_align = "left",
                         xlim = NULL,
                         x_breaks = NULL,
                         add_rope = FALSE,
                         rope_value = NULL,
                         rope_color = "grey50",
                         shrinkage_output = "density",
                         null_value = NULL,
                         color_palette = NULL,
                         color_study_posterior_null_left = "deepskyblue",
                         color_study_posterior_null_right = "violet",
                         color_study_posterior = "dodgerblue",
                         color_study_posterior_outline = "blue",
                         color_overall_posterior = "blue",
                         color_shrinkage_pointinterval = "purple",
                         color_shrinkage_outline = "purple",
                         color_shrinkage_fill = NULL,
                         split_color_by_null = FALSE,
                         color_favours_control = "firebrick",
                         color_favours_intervention = "dodgerblue",
                         add_arm_labels = TRUE,
                         reverse_arm_labels = FALSE,
                         add_pred = FALSE,
                         add_pred_subgroup = FALSE,
                         pred_output = c("density", "pointinterval"),
                         color_pred_posterior = "forestgreen",
                         color_pred_outline = "darkgreen",
                         color_pred_pointinterval = "forestgreen",
                         plot_width = 4,
                         add_rob = FALSE,
                         rob_tool = c("rob2", "robins_i", "quadas2", "robins_e"),
                         add_rob_legend = FALSE,
                         exclude_high_rob = FALSE,
                         font = NULL) {

  # Match pred_output argument
  pred_output <- rlang::arg_match(pred_output)

  if (!inherits(model, "bayesma")) {
    cli::cli_abort("{.arg model} must be a {.cls bayesma} object.")
  }
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  # ---- Input validation ----
  validate_inputs_bayesma(
    model           = model,
    data            = data,
    estimand         = estimand,
    studyvar        = {{studyvar}},
    year            = {{year}},
    subgroup        = subgroup,
    subgroup_var    = {{subgroup_var}},
    sort_studies_by = sort_studies_by,
    shrinkage_output = shrinkage_output
  )

  # ---- Column renaming ----
  # Marginal estimands use the same columns as their underlying likelihood model
  col_estimand <- if (is_marginal_estimand(estimand)) {
    switch(model$meta$likelihood, binomial = "OR", poisson = "IRR", gaussian = "MD")
  } else {
    estimand
  }

  if (col_estimand %in% c("OR", "RR", "HR")) {
    data <- data |>
      dplyr::rename(
        Author             = {{studyvar}},
        Year               = {{year}},
        N_Control          = {{c_n}},
        N_Intervention     = {{i_n}},
        Event_Control      = {{c_event}},
        Event_Intervention = {{i_event}}
      )
  } else if (col_estimand %in% c("MD", "SMD")) {
    data <- data |>
      dplyr::rename(
        Author             = {{studyvar}},
        Year               = {{year}},
        N_Control          = {{c_n}},
        N_Intervention     = {{i_n}},
        Mean_Control       = {{c_mean}},
        Mean_Intervention  = {{i_mean}},
        SD_Control         = {{c_sd}},
        SD_Intervention    = {{i_sd}}
      )
  } else if (col_estimand == "IRR") {
    data <- data |>
      dplyr::rename(
        Author             = {{studyvar}},
        Year               = {{year}},
        Time_Control       = {{c_time}},
        Time_Intervention  = {{i_time}},
        Event_Control      = {{c_event}},
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

  # Disambiguate duplicate Author names
  data <- make_authors_unique(data)

  # ---- Compute yi/vi ----
  # For marginal binary estimands: show observed risk differences (p_int - p_ctrl)
  # rather than log-ORs so the forest is consistent with the estimand scale.
  is_marginal_binary <- is_marginal_estimand(estimand) &&
    identical(model$meta$likelihood, "binomial")

  es <- model$meta$es
  if (is_marginal_binary &&
      all(c("Event_Intervention", "N_Intervention", "Event_Control", "N_Control") %in% names(data))) {
    p_int  <- data$Event_Intervention / data$N_Intervention
    p_ctrl <- data$Event_Control      / data$N_Control
    data$yi <- p_int - p_ctrl
    data$vi <- p_int * (1 - p_int) / data$N_Intervention +
               p_ctrl * (1 - p_ctrl) / data$N_Control
  } else if (!is.null(es)) {
    data$yi <- es$yi
    data$vi <- es$sei^2
  } else {
    fdf <- model$forest_df |> dplyr::filter(.data$type == "study")
    data$yi <- fdf$estimate
    data$vi <- ((fdf$upper - fdf$lower) / (2 * 1.96))^2
  }

  # ---- RoB handling (shared) ----
  rob_vars <- c("D1", "D2", "D3", "Overall")
  missing_vars <- setdiff(rob_vars, names(data))

  if ((length(missing_vars) > 0) && isTRUE(add_rob)) {
    cli::cli_abort("Risk of Bias columns must be provided for addition to the forest plot.")
  }

  if (length(missing_vars) > 0) {
    data <- data |> dplyr::mutate(
      D1 = NA_character_,
      Overall = NA_character_)
  }

  if (!shrinkage_output %in% c("density", "pointinterval")) {
    cli::cli_abort("{.arg shrinkage_output} must be either {.val density} or {.val pointinterval}.")
  }

  if (is.null(color_shrinkage_fill)) {
    color_shrinkage_fill <- NA
  }

  # ---- Exclude high RoB (shared logic, but refit differs) ----
  if (isTRUE(exclude_high_rob)) {
    if ("Overall" %in% names(data)) {
      data <- data |> dplyr::filter(Overall != "High" | is.na(Overall))

      message("Re-fitting model after excluding high risk of bias studies...")

        model <- refit_bayesma(model, data)
      }
    }

  # Count distinct studies after filtering
  num_studies <- data |> dplyr::distinct(Author) |> nrow()

  if (!subgroup) {
    if (num_studies < 2) {
      cli::cli_abort("Cannot run meta-analysis with fewer than 2 studies.")
    }
  } else {
    subgroup_counts <- data |>
      dplyr::filter(!is.na(Subgroup)) |>
      dplyr::summarise(n = dplyr::n_distinct(Author), .by = Subgroup)

    if (length(subgroup_counts) == 0) {
      warning("Subgroups with <2 studies will only show individual results and not pooled results for that subgroup. The study will still contribute to the overall effect.")
    }
  }

  # Create subgroup order if needed
  subgroup_order <- if (isTRUE(subgroup)) {
    get_subgroup_order(data, sort_subgroup_by)
  } else {
    NULL
  }

  # ---- Detect random effects ----
  has_re <- has_random_effects(model)

  # Prediction intervals are not meaningful for common-effect models
  if (isFALSE(has_re) && isTRUE(add_pred)) {
    warning("Prediction intervals are not available for common-effect models ",
            "(no between-study heterogeneity). `add_pred` will be ignored.",
            call. = FALSE)
    add_pred <- FALSE
  }

  if (isFALSE(has_re) && isTRUE(add_pred_subgroup)) {
    add_pred_subgroup <- FALSE
  }

  # ---- Extract fixef-equivalent for density plot ----
  # For marginal estimands use the marginal distribution so the reference lines
  # appear on the correct (ARR) scale rather than the underlying log-OR scale.
  fixef_summary <- if (is_marginal_estimand(estimand) && !is.null(model$marginal)) {
    marg_s <- model$marginal$summary
    matrix(
      c(marg_s$median, NA_real_, marg_s$lower, marg_s$upper),
      nrow = 1,
      dimnames = list("Intercept", c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
    )
  } else {
    extract_fixef(model)
  }

  # ---- Forest data extraction (dispatch by model class) ----
  if (!subgroup) {
    # Single Forest Plot Workflow
    study.effect.draws <- forest_data_fn(
      data            = data,
      model           = model,
      estimand        = estimand,
      subgroup        = FALSE,
      sort_studies_by = sort_studies_by,
      subgroup_order  = NULL,
      add_pred        = add_pred,
      add_pred_subgroup = FALSE,
      has_re          = has_re
    )

    # Create the density plot
    study.plot <- study.density.plot_fn(
      df = study.effect.draws,
      fixef_summary = fixef_summary,
      estimand = estimand,
      subgroup = FALSE,
      subgroup_order = NULL,
      has_re = has_re,
      color_palette = color_palette,
      color_study_posterior = color_study_posterior,
      color_study_posterior_outline = color_study_posterior_outline,
      color_overall_posterior = color_overall_posterior,
      color_shrinkage_outline = color_shrinkage_outline,
      color_shrinkage_pointinterval = color_shrinkage_pointinterval,
      color_shrinkage_fill = color_shrinkage_fill,
      split_color_by_null = split_color_by_null,
      color_favours_control = color_favours_control,
      color_favours_intervention = color_favours_intervention,
      label_control = label_control,
      label_intervention = label_intervention,
      add_arm_labels = add_arm_labels,
      reverse_arm_labels = reverse_arm_labels,
      shrinkage_output = shrinkage_output,
      xlim = xlim,
      x_breaks = x_breaks,
      null_value = null_value,
      add_rope = add_rope,
      rope_value = rope_value,
      rope_color = rope_color,
      add_pred = add_pred,
      add_pred_subgroup = add_pred_subgroup,
      pred_output = pred_output,
      color_pred_posterior = color_pred_posterior,
      color_pred_outline = color_pred_outline,
      color_pred_pointinterval = color_pred_pointinterval,
      font = font
    )

    # Create summary data for tables
    forest.data.summary <- forest.data.summary_fn(
      spread_df = study.effect.draws,
      data = data,
      estimand = estimand,
      sort_studies_by = sort_studies_by,
      subgroup = FALSE,
      add_pred = add_pred,
      add_pred_subgroup = add_pred_subgroup,
      has_re = has_re
    )

  } else {
    # Subgroup Forest Plot Workflow
    subgroup.effect.draws <- forest_data_fn(
      data              = data,
      model             = model,
      estimand          = estimand,
      subgroup          = TRUE,
      sort_studies_by   = sort_studies_by,
      subgroup_order    = subgroup_order,
      add_pred          = add_pred,
      add_pred_subgroup = add_pred_subgroup,
      has_re            = has_re
    )

    # Create subgroup summary
    forest.data.summary <- subgroup.effect.draws |>
      dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
      tidyr::nest(.by = Subgroup) |>
      dplyr::rename(spread_df = data) |>
      dplyr::mutate(
        spread_df = purrr::map2(spread_df, Subgroup, ~ dplyr::mutate(.x, Subgroup = .y)),
        subgroup.forest.summary = purrr::map(spread_df, ~ forest.data.summary_fn(
          spread_df = .x,
          data = data,
          estimand = estimand,
          sort_studies_by = sort_studies_by,
          add_pred = add_pred,
          add_pred_subgroup = add_pred_subgroup,
          has_re = has_re))) |>
      tidyr::unnest(subgroup.forest.summary) |>
      dplyr::select(-spread_df)

    if (isTRUE(subgroup) && estimand %in% c("MD", "SMD")){
      forest.data.summary <- forest.data.summary |>
        dplyr::mutate(
          N_int = dplyr::case_when(
            Author == "Overall Effect" ~ as.character(sum(N_Intervention, na.rm = TRUE)),
            TRUE ~ N_int),
          N_ctrl = dplyr::case_when(
            Author == "Overall Effect" ~ as.character(sum(N_Control, na.rm = TRUE)),
            TRUE ~ N_ctrl),
          int_mean_sd = dplyr::case_when(
            Author == "Overall Effect" ~ NA,
            TRUE ~ int_mean_sd),
          ctrl_mean_sd = dplyr::case_when(
            Author == "Overall Effect" ~ NA,
            TRUE ~ ctrl_mean_sd))
    }

    if (isTRUE(subgroup) && estimand %in% c("OR", "HR", "RR", "IRR")) {
      forest.data.summary <- forest.data.summary |>
        dplyr::mutate(
          Subgroup = dplyr::case_when(
            Author == "Overall Effect" ~ "Overall",
            TRUE ~ Subgroup),
          control_outcome_frac = dplyr::case_when(
            Author == "Overall Effect" & control_outcome_frac == "NA/NA" ~
              paste0(sum(Event_Control, na.rm = TRUE), "/", sum(N_Control, na.rm = TRUE)),
            TRUE ~ control_outcome_frac),
          int_outcome_frac = dplyr::case_when(
            Author == "Overall Effect" & int_outcome_frac == "NA/NA" ~
              paste0(sum(Event_Intervention, na.rm = TRUE), "/", sum(N_Intervention, na.rm = TRUE)),
            TRUE ~ int_outcome_frac))
    }

    if (isTRUE(exclude_high_rob)) {
      forest.data.summary <- forest.data.summary |>
        dplyr::left_join(subgroup_counts, by = dplyr::join_by(Subgroup)) |>
        dplyr::mutate(
          Author = dplyr::case_when(
            n == 1 & Author %in% c("Pooled Effect", "Overall Effect") ~ "No Pooled Effect",
            TRUE ~ Author))
    } else {
      forest.data.summary
    }

    # Create plot data with subgroup (inserting spacers)
    subgroup.plot.data <- purrr::map(unique(subgroup.effect.draws$Subgroup), ~ {
      subgroup_data <- subgroup.effect.draws |> dplyr::filter(Subgroup == .x)
      spacer_row <- subgroup_data[1, ]
      spacer_row[1, ] <- NA
      spacer_row$Subgroup[1] <- .x
      spacer_row$Author[1] <- paste0("--- ", .x, " ---")
      dplyr::bind_rows(spacer_row, subgroup_data)
    }) |> purrr::list_rbind()

    # Create lookup for author ordering
    seen_combos <- character(0)
    counter <- 1
    author_lookup <- purrr::map(1:nrow(subgroup.plot.data), ~ {
      current_combo <- paste(subgroup.plot.data$Subgroup[.x], subgroup.plot.data$Author[.x], sep = "_")

      if(!current_combo %in% seen_combos) {
        seen_combos <<- c(seen_combos, current_combo)
        result <- data.frame(
          Subgroup = subgroup.plot.data$Subgroup[.x],
          Author = subgroup.plot.data$Author[.x],
          Author_ordered = counter
        )
        counter <<- counter + 1
        return(result)
      } else {
        return(NULL)
      }
    }) |>
      purrr::compact() |>
      purrr::list_rbind()

    subgroup.plot.data <- subgroup.plot.data |>
      dplyr::left_join(author_lookup, by = dplyr::join_by(Subgroup, Author))

    study.plot <- study.density.plot_fn(
      df = subgroup.plot.data,
      fixef_summary = fixef_summary,
      estimand = estimand,
      subgroup = subgroup,
      subgroup_order = subgroup_order,
      has_re = has_re,
      color_palette = color_palette,
      color_study_posterior = color_study_posterior,
      color_study_posterior_outline = color_study_posterior_outline,
      color_overall_posterior = color_overall_posterior,
      color_shrinkage_outline = color_shrinkage_outline,
      color_shrinkage_pointinterval = color_shrinkage_pointinterval,
      color_shrinkage_fill = color_shrinkage_fill,
      split_color_by_null = split_color_by_null,
      color_favours_control = color_favours_control,
      color_favours_intervention = color_favours_intervention,
      label_control = label_control,
      label_intervention = label_intervention,
      add_arm_labels = add_arm_labels,
      reverse_arm_labels = reverse_arm_labels,
      shrinkage_output = shrinkage_output,
      xlim = xlim,
      x_breaks = x_breaks,
      null_value = null_value,
      add_rope = add_rope,
      rope_value = rope_value,
      rope_color = rope_color,
      add_pred = add_pred,
      pred_output = pred_output,
      color_pred_posterior = color_pred_posterior,
      color_pred_outline = color_pred_outline,
      color_pred_pointinterval = color_pred_pointinterval,
      font = font
    )
  }

  # Restore original Author names for table display
  if ("Author_original" %in% names(forest.data.summary)) {
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        Author = dplyr::if_else(
          is.na(Author_original),
          as.character(Author),
          as.character(Author_original)
        )
      )
  }

  # Create left table (study information)
  forest.table.left <- forest_table_left_fn(
    forest.data.summary = forest.data.summary,
    subgroup            = subgroup,
    label_control       = label_control,
    label_intervention  = label_intervention,
    estimand            = estimand,
    font                = font
  )

  # Create right table (effect sizes and optionally RoB)
  forest.table.right <- forest_table_right_fn(
    df       = forest.data.summary,
    subgroup = subgroup,
    add_rob  = add_rob,
    estimand = estimand,
    has_re   = has_re,
    font     = font
  )

  # Validate RoB legend parameters
  if (add_rob_legend == TRUE && add_rob == FALSE) {
    warning("Risk of Bias legend cannot be added when add_rob = FALSE. ",
            "Set add_rob = TRUE to include the RoB table and legend, ",
            "or set add_rob_legend = FALSE to suppress this warning.",
            call. = FALSE)
    add_rob_legend <- FALSE
  }

  # Combine everything using patchwork
  forest_plot <- patchwork_fn(
    table.left = forest.table.left,
    study.density.plot = study.plot,
    table.right = forest.table.right,
    plot_width = plot_width,
    title = title,
    subtitle = subtitle,
    title_align = title_align,
    add_rob_legend = add_rob_legend,
    rob_tool = rob_tool,
    font = font
  )

  # Attach recommended figure height based on number of rows
  n_rows <- nrow(forest.data.summary)
  has_title <- !is.null(title) || !is.null(subtitle)
  recommended_height <- n_rows * 0.55 + 1.5 + if (has_title) 0.8 else 0
  attr(forest_plot, "recommended_height") <- recommended_height
  attr(forest_plot, "recommended_width")  <- 14

  return(forest_plot)
}
