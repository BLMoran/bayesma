#' Generate Sensitivity Analysis Plot for Bayesian Meta-Analysis
#'
#' @description
#' Creates a sensitivity analysis visualisation showing how meta-analytic
#' estimates vary across different model strategies and user-specified priors.
#'
#' For RoBMA results, models must be pre-computed using `run_robma_sensitivity()`
#' and attached to the bayesma object before calling this function.
#'
#' @param model A fitted `bayesma` object. For RoBMA results, must have
#'   `$robma_sensitivity` attached via `attach_robma_sensitivity()`.
#' @param data A data frame containing the study data used to fit the model.
#' @param priors A named list of prior specifications. Each element must be a list
#'   with at least `mu_prior` and (optionally) `tau_prior`, and may include
#'   `name` used for display.
#' @param measure Effect measure string (e.g., "OR", "RR", "HR", "IRR", "MD", "SMD").
#' @param prior_order Optional character vector specifying the display order of priors.
#'   Should contain the names (IDs) of the priors in the desired order,
#'   e.g., `c("vague", "weak_reg", "informative")`. If NULL (default), priors
#'   are displayed in the order they appear in the `priors` list.
#' @param model_order Optional. Specifies the display order of model sections.
#'   Can be provided as unquoted names (using rlang) or as a character vector.
#'
#'   Valid values match the `incl_*` argument names (without the `incl_` prefix):
#'   `common_effect`, `random_effect`, `bias_corrected`, `selection_copas`,
#'   `selection_weight`, `pet_peese`, `robust`, `robma`.
#'
#'   For RoBMA, use `robma` to include both conditional and model-averaged estimates
#'   (conditional first), or `robma_conditional` / `robma_model_averaged` separately.
#'
#'   Example: `model_order = c(common_effect, random_effect, pet_peese, robma)`
#'
#'   If NULL (default), sections are displayed in the default order.
#' @param rob_var Optional. Name of the risk-of-bias variable (unquoted).
#' @param exclude_high_rob Logical. If TRUE, runs an "Excluding High RoB" section.
#' @param incl_common_effect,incl_random_effect,incl_bias_corrected,incl_selection_copas,incl_selection_weight,incl_pet_peese,incl_robust
#'   Logical. Which model strategy sections to include.
#' @param incl_robma Logical. If TRUE, include RoBMA section. Requires
#'   `model$robma_sensitivity` to be present (see Details).
#' @param parallel Logical. If TRUE, uses parallel processing for non-RoBMA refits.
#' @param workers Optional integer. Number of parallel workers.
#' @param seed Logical. If TRUE (default), uses parallel-safe seeding.
#' @param add_probs Logical. Add probability columns to the results table.
#' @param null_value,null_range,add_null_range,color_null_range ROPE settings.
#' @param label_control,label_intervention Group labels for the plot.
#' @param title,subtitle,title_align Title settings.
#' @param xlim,x_breaks Density axis settings.
#' @param color_palette,color_overall_posterior,color_overall_posterior_outline Colour settings.
#' @param split_color_by_null Logical. If TRUE, colour the posterior split at
#'   the null value using `color_favours_control` and `color_favours_intervention`.
#' @param color_favours_control Colour for the side of the posterior favouring control.
#' @param color_favours_intervention Colour for the side favouring intervention.
#' @param plot_width Width ratio for the density plot section.
#' @param font Optional font family.
#'
#' @return A `bayesma_sensitivity_plot` object (patchwork combining tables and plots).
#'
#' @details
#' ## RoBMA Results
#'
#' To include RoBMA in the sensitivity plot, you must pre-compute the RoBMA fits:
#'
#' ```r
#' # Step 1: Run RoBMA sensitivity analysis
#' robma_sens <- run_robma_sensitivity(
#'   data = my_data,
#'   priors = my_priors,
#'   robma_template = my_robma_fit,
#'   parallel = TRUE
#' )
#'
#' # Step 2: Attach to bayesma model
#' model <- attach_robma_sensitivity(model, robma_sens)
#'
#' # Step 3: Create sensitivity plot with RoBMA
#' sensitivity_plot(model, data, priors, measure = "OR", incl_robma = TRUE)
#' ```
#'
#' This separation allows the computationally expensive RoBMA fitting to be done
#' once and cached, rather than being re-run every time the plot is generated.
#'
#' @export
sensitivity_plot <- function(
    model,
    data,
    priors,
    measure,
    prior_order = NULL,
    model_order = NULL,
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
    parallel = FALSE,
    workers = NULL,
    seed = TRUE,
    add_probs = FALSE,
    null_value = NULL,
    null_range = NULL,
    add_null_range = FALSE,
    color_null_range = "#77bb41",
    label_control = "Control",
    label_intervention = "Intervention",
    title = NULL,
    subtitle = NULL,
    title_align = "left",
    xlim = NULL,
    x_breaks = NULL,
    color_palette = NULL,
    color_overall_posterior = "dodgerblue",
    color_overall_posterior_outline = "blue",
    split_color_by_null = FALSE,
    color_favours_control = "firebrick",
    color_favours_intervention = "dodgerblue",
    plot_width = 4,
    font = NULL
) {

  # ---------------------------
  # 1) Validation
  # ---------------------------
  validate_inputs_sens_plot(
    model   = model,
    data    = data,
    measure = measure,
    priors  = priors
  )

  if (!inherits(model, "bayesma")) {
    cli::cli_abort("This version of sensitivity_plot() currently expects a {.cls bayesma} model.")
  }

  # Validate RoBMA requirements
  if (isTRUE(incl_robma)) {
    if (!has_robma_sensitivity(model)) {
      cli::cli_abort(c(
        "RoBMA sensitivity fits not found on model.",
        "i" = "Use {.fn run_robma_sensitivity} to pre-compute RoBMA fits, then",
        "i" = "attach them with {.fn attach_robma_sensitivity} before calling {.fn sensitivity_plot}.",
        "",
        "Example:",
        "  robma_sens <- run_robma_sensitivity(data, priors, robma_template = robma_fit)",
        "
  model <- attach_robma_sensitivity(model, robma_sens)",
        "  sensitivity_plot(model, data, priors, measure, incl_robma = TRUE)"
      ))
    }
  }

  # Validate priors structure + build labels
  if (!is.list(priors) || length(priors) == 0) {
    cli::cli_abort("{.arg priors} must be a non-empty list.")
  }

  prior_ids <- names(priors)
  if (is.null(prior_ids) || any(prior_ids == "")) {
    cli::cli_abort("{.arg priors} must be a named list (each prior must have an ID).")
  }

  # Handle prior_order
  if (!is.null(prior_order)) {
    # Convert to character if needed
    prior_order <- as.character(prior_order)

    # Validate that all specified priors exist
    missing_priors <- setdiff(prior_order, prior_ids)
    if (length(missing_priors) > 0) {
      cli::cli_abort(c(
        "Some priors in {.arg prior_order} not found in {.arg priors}:",
        "x" = paste(missing_priors, collapse = ", ")
      ))
    }

    # Warn about priors not in order (they will be appended)
    extra_priors <- setdiff(prior_ids, prior_order)
    if (length(extra_priors) > 0) {
      cli::cli_warn(c(
        "Some priors not specified in {.arg prior_order} will be appended:",
        "i" = paste(extra_priors, collapse = ", ")
      ))
      prior_order <- c(prior_order, extra_priors)
    }

    # Reorder prior_ids
    prior_ids <- prior_order
  }

  prior_name_map <- purrr::map_chr(priors, ~ .x$name %||% "")
  prior_label_map <- purrr::map_chr(prior_ids, function(pid) {
    nm <- prior_name_map[[pid]]
    if (!is.null(nm) && nm != "") nm else pid
  })
  names(prior_label_map) <- prior_ids

  # ---------------------------
  # 2) Null range handling
  # ---------------------------
  props      <- get_measure_properties(measure)
  null_value <- null_value %||% props$null_value

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
  # 3) Define sections (model strategies)
  # ---------------------------
  stage <- model$meta$stage %||% model$meta$call_args$stage
  stage_label <- switch(stage, one_stage = "One-Stage", two_stage = "Two-Stage", stage)

  rob_var_q <- rlang::enquo(rob_var)
  rob_data <- NULL
  if (isTRUE(exclude_high_rob)) {
    if (rlang::quo_is_null(rob_var_q)) {
      cli::cli_warn(c(
        "{.arg exclude_high_rob} is TRUE but {.arg rob_var} is not specified.",
        "i" = "The 'Excluding High RoB' section will not be included.",
        "i" = "Specify {.arg rob_var} with the column name containing risk of bias ratings."
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

  # Apply model_order if specified
  if (!is.null(model_order)) {
    # Map section ids to their list elements
    section_ids <- purrr::map_chr(active_sections, ~ .x$id)

    # Filter model_order to only include active sections (excluding robma for now)
    non_robma_order <- setdiff(model_order, c("robma", "robma_conditional", "robma_model_averaged"))
    valid_order <- intersect(non_robma_order, section_ids)

    # Sections not in model_order get appended at the end
    remaining <- setdiff(section_ids, valid_order)
    final_order <- c(valid_order, remaining)

    # Reorder active_sections
    active_sections <- active_sections[match(final_order, section_ids)]
    active_sections <- purrr::compact(active_sections)  # Remove NULLs
  }

  if (!length(active_sections) && !isTRUE(incl_robma)) {
    cli::cli_abort("No model sections selected. Set at least one {.arg incl_*} argument to TRUE.")
  }

  # ---------------------------
  # 4) Build task grid (section × prior) for non-RoBMA models
  # ---------------------------
  draws_ma <- NULL

  if (length(active_sections) > 0) {
    task_grid <- tidyr::crossing(
      tibble::tibble(section_idx = seq_along(active_sections)),
      tibble::tibble(prior_id = prior_ids)
    )

    # helper: effect-scale transform for plotting
    transform_mu <- function(mu_raw) {
      if (measure %in% c("OR", "RR", "HR", "IRR")) exp(mu_raw) else mu_raw
    }

    # Helper: check if two bayesma_prior objects are identical
    priors_match <- function(p1, p2) {
      if (is.null(p1) && is.null(p2)) return(TRUE)
      if (is.null(p1) || is.null(p2)) return(FALSE)
      if (!inherits(p1, "bayesma_prior") || !inherits(p2, "bayesma_prior")) return(FALSE)
      identical(unclass(p1), unclass(p2))
    }

    # Extract original model priors and model_type for comparison
    orig_mu_prior  <- model$meta$priors$mu %||% model$meta$call_args$mu_prior
    orig_tau_prior <- model$meta$priors$tau %||% model$meta$call_args$tau_prior
    orig_model_type <- model$meta$model_type %||% model$meta$call_args$model_type %||% "random_effect"
    orig_stage <- model$meta$stage %||% model$meta$call_args$stage %||% "one_stage"

    cli::cli_h3("Fitting sensitivity models")
    cli::cli_alert_info("Original model: {orig_model_type} ({orig_stage})")

    n_tasks <- nrow(task_grid)
    n_reused <- 0
    n_refitted <- 0

    # task runner: one refit -> draws tibble
    run_bayesma_task <- function(section_idx, prior_id) {
      sec <- active_sections[[section_idx]]
      ps  <- priors[[prior_id]]

      mu_prior  <- ps$mu_prior
      tau_prior <- ps$tau_prior

      if (is.null(mu_prior)) {
        cli::cli_abort(c(
          "Each prior must include {.val mu_prior}.",
          "i" = paste0("Missing mu_prior for prior id: ", prior_id)
        ))
      }

      # Check if we can reuse the original model's draws
      mu_matches  <- priors_match(mu_prior, orig_mu_prior)
      tau_matches <- is.null(tau_prior) || priors_match(tau_prior, orig_tau_prior)

      # Get section's target model_type and stage
      sec_model_type <- sec$overrides$model_type %||% orig_model_type
      sec_stage <- sec$overrides$stage %||% orig_stage

      # Can reuse if: priors match AND model_type matches AND stage matches AND no data override
      model_matches <- (sec_model_type == orig_model_type) && (sec_stage == orig_stage)
      can_reuse <- mu_matches && tau_matches && model_matches && is.null(sec$data_override)

      if (can_reuse) {
        cli::cli_alert_success("{sec$label} + {prior_id}: reusing original draws")
        n_reused <<- n_reused + 1
        mu_raw <- as.numeric(model$draws[["mu"]])
        x <- transform_mu(mu_raw)
        return(tibble::tibble(
          section_label = sec$label,
          prior         = prior_id,
          prior_label   = prior_label_map[[prior_id]],
          x             = x
        ))
      }

      cli::cli_alert_info("{sec$label} + {prior_id}: fitting model...")
      n_refitted <<- n_refitted + 1

      fit <- tryCatch(
        {
          refit_args <- c(
            list(
              model    = model,
              data     = sec$data_override %||% data,
              mu_prior = mu_prior,
              tau_prior = tau_prior
            ),
            sec$overrides
          )
          do.call(refit_bayesma_update, refit_args)
        },
        error = function(e) {
          cli::cli_warn(c(
            "Failed to refit section {.val {sec$label}} with prior {.val {prior_id}}.",
            "i" = e$message
          ))
          NULL
        }
      )

      if (is.null(fit)) {
        return(tibble::tibble(
          section_label = character(),
          prior = character(),
          prior_label = character(),
          x = numeric()
        ))
      }

      # Check for problematic estimates (e.g., Copas model failures)
      mu_raw <- as.numeric(fit$draws[["mu"]])
      x <- transform_mu(mu_raw)

      # Warn and optionally exclude if estimates are unreasonable
      if (measure %in% c("OR", "RR", "HR", "IRR")) {
        median_x <- stats::median(x)
        q_low <- stats::quantile(x, 0.025)
        q_high <- stats::quantile(x, 0.975)

        # Flag extreme estimates (likely model failure)
        if (median_x > 100 || median_x < 0.01 || q_high > 10000) {
          cli::cli_warn(c(
            "!" = "Extreme estimates for {sec$label} + {prior_id}",
            "i" = "Median {measure} = {format(round(median_x, 1), big.mark=',')}",
            "i" = "95% CrI: [{format(round(q_low, 1), big.mark=',')}, {format(round(q_high, 1), big.mark=',')}]",
            "i" = "This typically indicates model convergence failure.",
            "i" = "Consider excluding this model type from sensitivity analysis."
          ))
        }
      }

      tibble::tibble(
        section_label = sec$label,
        prior         = prior_id,
        prior_label   = prior_label_map[[prior_id]],
        x             = x
      )
    }

    # run bayesma tasks
    draws_ma <- bayesma_future_pmap_dfr(
      .l = list(task_grid$section_idx, task_grid$prior_id),
      .f = run_bayesma_task,
      parallel = parallel,
      workers  = workers,
      seed     = seed
    )

    cli::cli_alert_info("Summary: {n_reused} reused, {n_refitted} refitted")
  }

  # ---------------------------
  # 5) RoBMA — extract draws from pre-computed fits
  # ---------------------------
  draws_robma <- NULL

  if (isTRUE(incl_robma)) {
    robma_sens <- model$robma_sensitivity
    robma_meta <- attr(robma_sens, "meta")

    # Validate that priors match
    missing_priors <- setdiff(prior_ids, names(robma_sens))
    if (length(missing_priors) > 0) {
      cli::cli_warn(c(
        "Some priors not found in pre-computed RoBMA fits:",
        "x" = paste(missing_priors, collapse = ", "),
        "i" = "These will be excluded from the RoBMA section."
      ))
    }

    available_priors <- intersect(prior_ids, names(robma_sens))

    if (length(available_priors) > 0) {
      # Create TWO sets of draws for each prior:
      # 1. Conditional estimate (H1 only) - for "RoBMA (Conditional)"
      # 2. Model-averaged estimate (including null) - for "RoBMA (Model-Averaged)"
      draws_robma <- purrr::map(available_priors, function(pid) {
        robma_fit <- robma_sens[[pid]]

        if (is.null(robma_fit)) {
          return(tibble::tibble(
            x = numeric(),
            prior = character(),
            prior_label = character(),
            section_label = character(),
            is_null_draw = logical()
          ))
        }

        # Get all draws with null tracking
        all_draws <- robma_to_sensitivity_draws(
          robma_fit     = robma_fit,
          measure       = measure,
          prior         = pid,
          prior_label   = prior_label_map[[pid]] %||% pid,
          section_label = "RoBMA"  # temporary, will be updated
        )

        # Create conditional draws (H1 only)
        conditional_draws <- all_draws |>
          dplyr::filter(!is_null_draw) |>
          dplyr::mutate(section_label = "RoBMA (Conditional)")

        # Create model-averaged draws (all draws)
        # Add small jitter to null draws so the density shows a spike properly
        model_avg_draws <- all_draws |>
          dplyr::mutate(section_label = "RoBMA (Model-Averaged)")

        # For model-averaged density: add jitter to null draws
        # This creates a narrow spike at null rather than a point mass
        if (any(model_avg_draws$is_null_draw)) {
          h1_draws <- model_avg_draws$x[!model_avg_draws$is_null_draw]
          if (length(h1_draws) > 0) {
            # Get the null value on the transformed scale
            null_val <- if (measure %in% c("OR", "RR", "HR", "IRR")) 1 else 0
            # Jitter width: 2% of H1 range on log scale for ratio measures
            if (measure %in% c("OR", "RR", "HR", "IRR")) {
              log_range <- diff(range(log(h1_draws)))
              jitter_sd <- exp(log_range * 0.02) - 1
              n_null <- sum(model_avg_draws$is_null_draw)
              model_avg_draws$x[model_avg_draws$is_null_draw] <-
                null_val * exp(stats::rnorm(n_null, 0, log_range * 0.02))
            } else {
              h1_sd <- sd(h1_draws)
              jitter_sd <- h1_sd * 0.02
              n_null <- sum(model_avg_draws$is_null_draw)
              model_avg_draws$x[model_avg_draws$is_null_draw] <-
                null_val + stats::rnorm(n_null, 0, jitter_sd)
            }
          }
        }

        # Return both
        dplyr::bind_rows(conditional_draws, model_avg_draws)
      }) |> purrr::list_rbind()
    }
  }

  draws <- dplyr::bind_rows(draws_ma, draws_robma)

  # ---------------------------
  # 5b) Apply model_order to section ordering
  # ---------------------------
  # Build the final section order based on model_order argument
  all_section_labels <- unique(draws$section_label)

  # Process model_order - allow unquoted names via substitute
  model_order_expr <- substitute(model_order)
  if (!is.null(model_order_expr) && !identical(model_order_expr, as.name("NULL"))) {
    # Convert expression to character vector
    if (is.call(model_order_expr) && identical(model_order_expr[[1]], as.name("c"))) {
      # c(common_effect, random_effect, ...) - extract symbols
      model_order <- vapply(as.list(model_order_expr)[-1], as.character, character(1))
    } else if (is.character(model_order)) {
      # Already a character vector (e.g., c("common_effect", "random_effect"))
      model_order <- model_order
    } else if (is.symbol(model_order_expr)) {
      # Single unquoted name
      model_order <- as.character(model_order_expr)
    } else {
      # Fallback - try to coerce
      model_order <- as.character(model_order)
    }

    # Create mapping from model_order ids to section labels
    # Names match incl_* arguments (without incl_ prefix) plus common aliases
    id_to_label <- list(
      common_effect = paste0("Common Effect (", stage_label, ")"),
      random_effect = paste0("Random Effects (", stage_label, ")"),
      exclude_high_rob = "Excluding High RoB",
      excluding_high_rob = "Excluding High RoB",  # alias
      bias_corrected = "Bias Adjusted (Jung)",
      bias_adjusted = "Bias Adjusted (Jung)",     # alias (matches section id)
      selection_copas = "Selection Model (Copas)",
      selection_weight = "Selection Model (Weight)",
      pet_peese = "PET-PEESE",
      robust = "Robust Mixture Model",
      robma = c("RoBMA (Conditional)", "RoBMA (Model-Averaged)"),
      robma_conditional = "RoBMA (Conditional)",
      robma_model_averaged = "RoBMA (Model-Averaged)"
    )

    # Validate model_order values
    invalid_ids <- setdiff(model_order, names(id_to_label))
    if (length(invalid_ids) > 0) {
      cli::cli_warn(c(
        "Unknown values in {.arg model_order} will be ignored:",
        "x" = paste(invalid_ids, collapse = ", "),
        "i" = "Valid values: common_effect, random_effect, exclude_high_rob, bias_corrected, selection_copas, selection_weight, pet_peese, robust, robma, robma_conditional, robma_model_averaged"
      ))
    }

    # Warn about sections in model_order that aren't actually included
    requested_labels <- unlist(id_to_label[intersect(model_order, names(id_to_label))])
    not_included <- setdiff(requested_labels, all_section_labels)
    if (length(not_included) > 0) {
      # Map back to the model_order ids for clearer message
      not_included_ids <- model_order[purrr::map_lgl(model_order, function(mo) {
        any(id_to_label[[mo]] %in% not_included)
      })]
      cli::cli_warn(c(
        "Some sections in {.arg model_order} are not included in the plot:",
        "x" = paste(not_included_ids, collapse = ", "),
        "i" = "Set the corresponding {.arg incl_*} argument to TRUE to include them."
      ))
    }

    # Build ordered section labels
    ordered_labels <- character(0)
    for (mo in model_order) {
      labels <- id_to_label[[mo]]
      if (!is.null(labels)) {
        # Only include labels that actually exist in the data
        labels <- intersect(labels, all_section_labels)
        ordered_labels <- c(ordered_labels, labels)
      }
    }

    # Append any sections not in model_order
    remaining <- setdiff(all_section_labels, ordered_labels)
    ordered_labels <- c(ordered_labels, remaining)

    # Convert section_label to factor with correct order
    draws <- draws |>
      dplyr::mutate(section_label = factor(section_label, levels = ordered_labels))
  } else {
    # Default order: keep as-is (order from bind_rows)
    draws <- draws |>
      dplyr::mutate(section_label = factor(section_label, levels = all_section_labels))
  }

  # Set prior_label as factor with correct order (respecting prior_order)
  # prior_label_map is ordered by prior_ids which respects prior_order
  ordered_prior_labels <- unname(prior_label_map[prior_ids])
  draws <- draws |>
    dplyr::mutate(prior_label = factor(prior_label, levels = ordered_prior_labels))

  # Debug info

  cli::cli_h3("Sensitivity Plot Data Summary")
  if (!is.null(draws_ma) && nrow(draws_ma) > 0) {
    n_ma_sections <- length(unique(draws_ma$section_label))
    n_ma_priors <- length(unique(draws_ma$prior))
    cli::cli_alert_success("Meta-analysis draws: {nrow(draws_ma)} rows ({n_ma_sections} sections x {n_ma_priors} priors)")
  } else {
    cli::cli_alert_warning("Meta-analysis draws: 0 rows (check if model refits succeeded)")
  }

  if (!is.null(draws_robma) && nrow(draws_robma) > 0) {
    n_robma_priors <- length(unique(draws_robma$prior))
    cli::cli_alert_success("RoBMA draws: {nrow(draws_robma)} rows ({n_robma_priors} priors)")
  } else if (isTRUE(incl_robma)) {
    cli::cli_alert_warning("RoBMA draws: 0 rows (check if RoBMA fits succeeded)")
  }

  if (nrow(draws) == 0) {
    cli::cli_abort(c(
      "No draws available for sensitivity plot.",
      "i" = "Check that model refits succeeded.",
      "i" = "For RoBMA, ensure robma_sensitivity fits are attached and not all NULL.",
      "i" = "Run run_robma_sensitivity() with parallel=FALSE to see error details."
    ))
  }

  cli::cli_alert_info("Total draws: {nrow(draws)} rows")

  # ---------------------------
  # 6) Prior table (mu/tau printing)
  # ---------------------------
  priors_mu_tau <- purrr::map(priors, ~ list(mu_prior = .x$mu_prior, tau_prior = .x$tau_prior))
  prior_table <- extract_mu_tau_priors(priors = priors_mu_tau, model = model) |>
    dplyr::mutate(prior_label = dplyr::coalesce(
      prior_label_map[prior], prior_label
    ))

  # ---------------------------
  # 7) Posterior summary + join priors
  # ---------------------------
  # Note: For RoBMA, we now have separate section_labels for
  # "RoBMA (Conditional)" and "RoBMA (Model-Averaged)", so no filtering needed.
  summary_df <- summarise_sensitivity_posteriors(
    draws      = draws,
    null_value = null_value,
    null_range = null_range
  ) |>
    dplyr::mutate(prior = as.character(prior)) |>
    dplyr::left_join(prior_table, by = dplyr::join_by(prior, prior_label))

  # ---------------------------
  # 8) Tables
  # ---------------------------
  has_null_range <- !is.null(null_range) && length(null_range) == 2

  table_left <- sensitivity_table_left(summary_df, font = font)
  table_right <- sensitivity_table_right(
    summary_df,
    measure              = measure,
    add_probs            = add_probs,
    add_probs_null_range = has_null_range && isTRUE(add_probs),
    font                 = font
  )

  # ---------------------------
  # 9) Density plot
  # ---------------------------
  sens_plot <- sensitivity_density_plot_fn(
    df                              = draws,
    measure                         = measure,
    split_color_by_null             = split_color_by_null,
    color_overall_posterior         = color_overall_posterior,
    color_overall_posterior_outline = color_overall_posterior_outline,
    color_favours_control           = color_favours_control,
    color_favours_intervention      = color_favours_intervention,
    label_control                   = label_control,
    label_intervention              = label_intervention,
    xlim                            = xlim,
    x_breaks                        = x_breaks,
    null_value                      = null_value,
    null_range                      = null_range,
    color_null_range                = color_null_range,
    add_null_range                  = add_null_range,
    font                            = font
  )

  # ---------------------------
  # 10) Build output object
  # ---------------------------
  density_args <- list(
    df                              = draws,
    measure                         = measure,
    split_color_by_null             = split_color_by_null,
    color_overall_posterior         = color_overall_posterior,
    color_overall_posterior_outline = color_overall_posterior_outline,
    color_favours_control           = color_favours_control,
    color_favours_intervention      = color_favours_intervention,
    label_control                   = label_control,
    label_intervention              = label_intervention,
    xlim                            = xlim,
    x_breaks                        = x_breaks,
    null_value                      = null_value,
    null_range                      = null_range,
    color_null_range                = color_null_range,
    add_null_range                  = add_null_range,
    font                            = font
  )

  out <- structure(
    list(
      table_left    = table_left,
      density_plot  = sens_plot,
      table_right   = table_right,
      density_args  = density_args,
      plot_width    = plot_width,
      title         = title,
      subtitle      = subtitle,
      title_align   = title_align,
      font          = font,
      # State flags for add_* functions
      has_probs     = isTRUE(add_probs),
      has_null_range = has_null_range,
      color_null_range = color_null_range,
      # Store everything needed for add_probs / add_null / re-rendering
      summary_df    = summary_df,
      prior_table   = prior_table,
      draws         = draws,
      measure       = measure,
      null_value    = null_value,
      null_range    = null_range
    ),
    class = "bayesma_sensitivity_plot"
  )

  # Attach recommended figure dimensions based on number of rows
  n_rows <- nrow(summary_df)
  has_title_sub <- !is.null(title) || !is.null(subtitle)
  recommended_height <- n_rows * 0.55 + 1.5 + if (has_title_sub) 0.8 else 0
  attr(out, "recommended_height") <- recommended_height
  attr(out, "recommended_width")  <- 14

  out
}


# S3: print / plot — auto-render the patchwork

#' @export
#' @noRd
print.bayesma_sensitivity_plot <- function(x, ...) {
  pw <- render_sensitivity_patchwork(x)
  print(pw)
  invisible(x)
}

#' @export
plot.bayesma_sensitivity_plot <- function(x, ...) {
  pw <- render_sensitivity_patchwork(x)
  print(pw)
  invisible(x)
}

#' Render the stored components into a patchwork object
#'
#' @param x A `bayesma_sensitivity_plot` object.
#' @return A patchwork object.
#' @export
render_sensitivity_patchwork <- function(x) {
  # Use patchwork_fn (defined in patchwork_fn_bayesma.R)
  patchwork_fn(
    table.left         = x$table_left,
    study.density.plot = x$density_plot,
    table.right        = x$table_right,
    plot_width         = x$plot_width,
    title              = x$title,
    subtitle           = x$subtitle,
    title_align        = x$title_align,
    add_rob_legend     = FALSE,
    rob_tool           = "rob2",
    font               = x$font
  )
}

# Post-Rendering additions (to avoid having to re-run the whole plot)

#' Post-render modifications to a sensitivity plot
#'
#' Helpers for adjusting an existing `bayesma_sensitivity_plot` object without
#' re-running the full pipeline.
#'
#' @param x A `bayesma_sensitivity_plot` object.
#' @param show Logical. Whether to display posterior probability columns.
#' @param range Numeric vector of length 2 giving the null/ROPE range, or
#'   `NULL` to use the default for the model's measure.
#' @param color_null_range Colour used to shade the null range.
#' @param title,subtitle Character. Plot titles.
#' @param align Title alignment: `"left"`, `"center"`, or `"right"`.
#' @param xlim Numeric vector of length 2. New x-axis limits.
#' @param x_breaks Optional numeric vector of x-axis break locations.
#' @param plot_width Positive numeric. Relative width of the density panel.
#'
#' @return The updated `bayesma_sensitivity_plot` object.
#'
#' @name sens_add
#' @export
sens_add_probs <- function(x, show = TRUE) {
  if (!inherits(x, "bayesma_sensitivity_plot")) {
    cli::cli_abort("{.fn sens_add_probs} requires a {.cls bayesma_sensitivity_plot} object.")
  }
  x$has_probs <- isTRUE(show)
  x$table_right <- sensitivity_table_right(
    x$summary_df,
    measure              = x$measure,
    add_probs            = x$has_probs,
    add_probs_null_range = x$has_null_range && x$has_probs,
    font                 = x$font
  )
  x
}

#' @rdname sens_add
#' @export
sens_add_null <- function(x, range = NULL, color_null_range = "#77bb41") {
  if (!inherits(x, "bayesma_sensitivity_plot")) {
    cli::cli_abort("{.fn sens_add_null} requires a {.cls bayesma_sensitivity_plot} object.")
  }

  # Resolve default null range based on measure
  if (is.null(range)) {
    range <- switch(
      x$measure,
      OR  = c(0.9, 1.1),
      RR  = c(0.9, 1.1),
      HR  = c(0.9, 1.1),
      IRR = c(0.9, 1.1),
      SMD = c(-0.1, 0.1),
      cli::cli_abort(
        "No default null range for measure {.val {x$measure}}. Please supply {.arg range} explicitly."
      )
    )
  } else if (length(range) == 1 && is.numeric(range)) {
    # Single value: symmetric around null_value
    range <- c(x$null_value - range, x$null_value + range)
  } else if (!is.numeric(range) || length(range) != 2) {
    cli::cli_abort("{.arg range} must be NULL, a single number, or a numeric vector of length 2.")
  }

  range <- sort(range)
  x$null_range       <- range
  x$has_null_range   <- TRUE
  x$color_null_range <- color_null_range

  # Recalculate summary_df with the new null_range
  x$summary_df <- summarise_sensitivity_posteriors(
    draws      = x$draws,
    null_value = x$null_value,
    null_range = x$null_range
  ) |>
    dplyr::mutate(prior = as.character(prior)) |>
    dplyr::left_join(x$prior_table, by = dplyr::join_by(prior, prior_label))

  # Rebuild the density plot with null range shading
  da <- x$density_args
  da$null_range      <- x$null_range
  da$add_null_range  <- TRUE
  da$color_null_range <- x$color_null_range
  x$density_plot <- do.call(sensitivity_density_plot_fn, da)
  x$density_args <- da

  # Rebuild the right table (shows δ columns if probs are active)
  x$table_right <- sensitivity_table_right(
    x$summary_df,
    measure              = x$measure,
    add_probs            = x$has_probs,
    add_probs_null_range = x$has_probs,
    font                 = x$font
  )

  x
}

#' @rdname sens_add
#' @export
sens_add_titles <- function(x, title = NULL, subtitle = NULL, align = NULL) {
  if (!inherits(x, "bayesma_sensitivity_plot")) {
    cli::cli_abort("{.fn sens_add_titles} requires a {.cls bayesma_sensitivity_plot} object.")
  }
  if (!is.null(title))    x$title       <- title
  if (!is.null(subtitle)) x$subtitle    <- subtitle
  if (!is.null(align))    x$title_align <- align
  x
}

#' @rdname sens_add
#' @export
sens_add_x_lim <- function(x, xlim, x_breaks = NULL) {
  if (!inherits(x, "bayesma_sensitivity_plot")) {
    cli::cli_abort("{.fn sens_add_x_lim} requires a {.cls bayesma_sensitivity_plot} object.")
  }
  if (!is.numeric(xlim) || length(xlim) != 2) {
    cli::cli_abort("{.arg xlim} must be a numeric vector of length 2.")
  }

  # Rebuild the density plot from scratch with updated limits
  da       <- x$density_args
  da$xlim  <- xlim
  if (!is.null(x_breaks)) da$x_breaks <- x_breaks

  x$density_plot <- do.call(sensitivity_density_plot_fn, da)
  x$density_args <- da
  x
}

#' @rdname sens_add
#' @export
sens_add_plot_width <- function(x, plot_width) {
  if (!inherits(x, "bayesma_sensitivity_plot")) {
    cli::cli_abort("{.fn sens_add_plot_width} requires a {.cls bayesma_sensitivity_plot} object.")
  }
  if (!is.numeric(plot_width) || length(plot_width) != 1 || plot_width <= 0) {
    cli::cli_abort("{.arg plot_width} must be a single positive number.")
  }
  x$plot_width <- plot_width
  x
}


# Register S3 methods
register_sensitivity_plot_methods <- function(envir = parent.frame()) {
  registerS3method("print", "bayesma_sensitivity_plot",
                   print.bayesma_sensitivity_plot,
                   envir = envir)
  registerS3method("plot", "bayesma_sensitivity_plot",
                   plot.bayesma_sensitivity_plot,
                   envir = envir)
}

# Run S3 registration when this file is source()
register_sensitivity_plot_methods()
