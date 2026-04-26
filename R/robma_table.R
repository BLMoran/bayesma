#' Create a gt Table for RoBMA Results
#'
#' @description
#' Produces a publication-ready \code{gt} table summarising the Robust Bayesian
#' Meta-Analysis results. Includes component summary (posterior probabilities
#' and Bayes factors), model-averaged estimates, and direction probabilities.
#'
#' @param x A \code{bayesma_robma} object from \code{robma()}.
#' @param digits Integer. Number of decimal places for numeric values.
#'   Default is 3.
#' @param include_components Logical. If \code{TRUE} (default), includes the
#'   components summary table (Effect, Heterogeneity, Bias).
#' @param include_estimates Logical. If \code{TRUE} (default), includes the
#'   model-averaged parameter estimates table.
#' @param include_direction Logical. If \code{TRUE} (default), includes
#'   direction probabilities.
#' @param include_models Logical. If \code{TRUE}, includes individual model
#'   posterior probabilities (bridge method only). Default is \code{FALSE}.
#' @param exponentiate Logical. If \code{TRUE} and the effect is on a log
#'
#' @return A \code{gt} table object.
#'
#' @examples
#' \dontrun{
#' fit <- robma(data, studyvar = "study", ...)
#' robma_table(fit)
#'
#' # Options
#' robma_table(fit, include_models = TRUE)
#' robma_table(fit, digits = 4)
#'
#' # Save
#' robma_table(fit) |> gt::gtsave("robma_results.html")
#' robma_table(fit) |> gt::gtsave("robma_results.docx")
#' }
#'
#' @export
robma_table <- function(
    x,
    digits = 3,
    include_components = TRUE,
    include_estimates = TRUE,

    include_direction = TRUE,
    include_models = FALSE,
    exponentiate = FALSE
) {

  if (!inherits(x, "bayesma_robma")) {
    cli::cli_abort("{.arg x} must be a {.cls bayesma_robma} object.")
  }

  # ---- Extract metadata ----

  method     <- x$meta$method %||% "bridge"
  eff_label  <- x$meta$effect_label
  is_ratio   <- eff_label %in% c("log_or", "log_rr")
  pp         <- x$posterior_probs
  bf         <- x$inclusion_bf
  nrp        <- x$meta$null_range_probs
  nr         <- x$meta$null_range
  pooled     <- dplyr::filter(x$forest_df, .data$type == "pooled")

  # Format helpers
  fmt_num <- function(val, d = digits) {
    sprintf(paste0("%.", d, "f"), val)
  }

  fmt_bf <- function(val) {
    if (is.na(val)) return("\u2014")
    if (is.infinite(val)) return("\u221E")
    if (val < 0.001) return(formatC(val, format = "e", digits = 1))
    if (val > 1000) return(formatC(val, format = "e", digits = 1))
    fmt_num(val, digits)
  }

  fmt_prob <- function(val, d = digits) {
    if (val < 0.001 && val > 0) {
      return(formatC(val, format = "e", digits = 1))
    }
    fmt_num(val, d)
  }

  # ---- Build Components Summary Table ----
  components_tbl <- NULL
  if (include_components) {
    # Calculate prior probs (assuming 0.5 each for H0/H1)
    n_models <- x$meta$n_models %||% length(x$component_fits)

    # Count models with each component
    if (!is.null(x$model_table)) {
      n_effect <- sum(!x$model_table$null_effect)
      n_hetero <- sum(x$model_table$has_heterogeneity)
      n_bias   <- sum(x$model_table$has_bias)
    } else {
      # Fallback for spike-slab
      n_effect <- n_models / 2
      n_hetero <- n_models / 2
      n_bias   <- n_models / 2
    }

    components_df <- tibble::tibble(
      component = c("Effect", "Heterogeneity", "Publication Bias"),
      models = c(
        sprintf("%d/%d", n_effect, n_models),
        sprintf("%d/%d", n_hetero, n_models),
        sprintf("%d/%d", n_bias, n_models)
      ),
      prior_prob = rep("0.500", 3),
      post_prob = c(
        fmt_prob(pp$effect),
        fmt_prob(pp$heterogeneity),
        fmt_prob(pp$bias)
      ),
      inclusion_bf = c(
        fmt_bf(bf$effect),
        fmt_bf(bf$heterogeneity),
        fmt_bf(bf$bias)
      )
    )

    components_tbl <- gt::gt(components_df) |>
      gt::cols_label(
        component = "Component",
        models = "Models",
        prior_prob = "Prior Prob.",
        post_prob = "Post. Prob.",
        inclusion_bf = "Inclusion BF"
      ) |>
      gt::tab_header(
        title = "Robust Bayesian Meta-Analysis",
        subtitle = if (method == "bridge") {
          sprintf("Bridge sampling (%d component models)", n_models)
        } else {
          bi <- x$meta$bias_indicator %||% "bias_corrected"
          sprintf("Spike-and-slab (bias indicator: %s)", bi)
        }
      ) |>
      gt::tab_spanner(
        label = "Components Summary",
        columns = c("component", "models", "prior_prob", "post_prob", "inclusion_bf")
      )
  }

  # ---- Build Model-Averaged Estimates Table ----
  estimates_tbl <- NULL
  if (include_estimates) {
    # Get mu and tau draws
    mu_draws  <- x$averaged_draws$mu
    tau_draws <- if (!is.null(x$averaged_draws$tau)) {
      x$averaged_draws$tau
    } else {
      # Try to extract from component fits
      NULL
    }

    estimates_list <- list(
      list(
        param = "\u03BC (effect)",
        mean = mean(mu_draws),
        median = stats::median(mu_draws),
        lower = stats
        ::quantile(mu_draws, 0.025),
        upper = stats::quantile(mu_draws, 0.975)
      )
    )

    # Add exponentiated if ratio scale
    if (is_ratio && exponentiate) {
      ratio_name <- toupper(gsub("log_", "", eff_label))
      estimates_list <- c(estimates_list, list(
        list(
          param = ratio_name,
          mean = mean(exp(mu_draws)),
          median = stats::median(exp(mu_draws)),
          lower = stats::quantile(exp(mu_draws), 0.025),
          upper = stats::quantile(exp(mu_draws), 0.975)
        )
      ))
    }

    # Add tau if available
    if (!is.null(tau_draws)) {
      estimates_list <- c(estimates_list, list(
        list(
          param = "\u03C4 (heterogeneity)",
          mean = mean(tau_draws),
          median = stats::median(tau_draws),
          lower = stats::quantile(tau_draws, 0.025),
          upper = stats::quantile(tau_draws, 0.975)
        )
      ))
    }

    estimates_df <- purrr::map(estimates_list, function(est) {
      tibble::tibble(
        parameter = est$param,
        mean = fmt_num(est$mean),
        median = fmt_num(est$median),
        ci_95 = sprintf("[%s, %s]", fmt_num(est$lower), fmt_num(est$upper))
      )
    }) |> purrr::list_rbind()

    # Add pooled estimate row from forest_df
    pooled_row <- tibble::tibble(
      parameter = sprintf("Pooled (%s)", eff_label),
      mean = "\u2014",
      median = fmt_num(pooled$estimate),
      ci_95 = sprintf("[%s, %s]", fmt_num(pooled$lower), fmt_num(pooled$upper))
    )

    estimates_df <- dplyr::bind_rows(pooled_row, estimates_df)

    estimates_tbl <- gt::gt(estimates_df) |>
      gt::cols_label(
        parameter = "Parameter",
        mean = "Mean",
        median = "Median",
        ci_95 = "95% CrI"
      ) |>
      gt::tab_spanner(
        label = "Model-Averaged Estimates",
        columns = c("parameter", "mean", "median", "ci_95")
      )
  }

  # ---- Build Direction Probabilities Table ----
  direction_tbl <- NULL
  if (include_direction && !is.null(nrp)) {
    if (is.null(nr)) {
      # Point null
      direction_df <- tibble::tibble(
        direction = c(
          "P(\u03BC < 0)",
          "P(\u03BC > 0)",
          "P(\u03BC = 0)"
        ),
        probability = c(
          fmt_prob(nrp$p_negative),
          fmt_prob(nrp$p_positive),
          if (nrp$p_null > 0.001) fmt_prob(nrp$p_null) else "\u2248 0"
        )
      )
    } else {
      # ROPE / null range
      if (is_ratio) {
        nr_nat <- nrp$null_range_natural
        ratio_name <- toupper(gsub("log_", "", eff_label))
        direction_df <- tibble::tibble(
          direction = c(
            sprintf("P(%s < %s)", ratio_name, fmt_num(nr_nat[1], 3)),
            sprintf("P(practically null)"),
            sprintf("P(%s > %s)", ratio_name, fmt_num(nr_nat[2], 3))
          ),
          probability = c(
            fmt_prob(nrp$p_negative),
            fmt_prob(nrp$p_null),
            fmt_prob(nrp$p_positive)
          )
        )
      } else {
        direction_df <- tibble::tibble(
          direction = c(
            sprintf("P(\u03BC < %s)", fmt_num(nr[1], 3)),
            "P(practically null)",
            sprintf("P(\u03BC > %s)", fmt_num(nr[2], 3))
          ),
          probability = c(
            fmt_prob(nrp$p_negative),
            fmt_prob(nrp$p_null),
            fmt_prob(nrp$p_positive)
          )
        )
      }
    }

    direction_tbl <- gt::gt(direction_df) |>
      gt::cols_label(
        direction = "Direction",
        probability = "Probability"
      ) |>
      gt::tab_spanner(
        label = "Direction Probabilities (Model-Averaged Posterior)",
        columns = c("direction", "probability")
      )
  }

  # ---- Build Individual Models Table (optional) ----
  models_tbl <- NULL
  if (include_models && method == "bridge" && !is.null(x$model_table)) {
    mt <- x$model_table |>
      dplyr::mutate(
        post_prob_fmt = purrr::map_chr(.data$post_prob, function(p) {
          if (p < 0.001 && p > 0) formatC(p, format = "e", digits = 1)
          else fmt_num(p, digits)
        }),
        log_ml_fmt = fmt_num(.data$log_ml, 1),
        effect_fmt = dplyr::if_else(.data$null_effect, "\u2717", "\u2713"),
        hetero_fmt = dplyr::if_else(.data$has_heterogeneity, "\u2713", "\u2717"),
        bias_fmt = dplyr::if_else(.data$has_bias, "\u2713", "\u2717")
      ) |>
      dplyr::select(
        "model", "effect_fmt", "hetero_fmt", "bias_fmt",
        "post_prob_fmt", "log_ml_fmt"
      )

    models_tbl <- gt::gt(mt) |>
      gt::cols_label(
        model = "Model",
        effect_fmt = "Effect",
        hetero_fmt = "Hetero.",
        bias_fmt = "Bias",
        post_prob_fmt = "Post. Prob.",
        log_ml_fmt = "log(ML)"
      ) |>
      gt::tab_spanner(
        label = "Component Models",
        columns = c("model", "effect_fmt", "hetero_fmt", "bias_fmt",
                    "post_prob_fmt", "log_ml_fmt")
      )
  }

  # ---- Combine tables ----
  # Start with components table as base
  gt_tbl <- if (!is.null(components_tbl)) {
    components_tbl
  } else if (!is.null(estimates_tbl)) {
    estimates_tbl
  } else if (!is.null(direction_tbl)) {
    direction_tbl
  } else {
    cli::cli_abort("At least one table section must be included.")
  }

  # ---- Apply consistent bayesma styling ----
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_spanners()
    ) |>
    gt::cols_align(align = "center", columns = dplyr::everything()) |>
    gt::cols_align(align = "left", columns = 1) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(14),
      heading.title.font.weight = "bold",
      heading.subtitle.font.size = gt::px(11),
      column_labels.font.size = gt::px(11),
      column_labels.border.bottom.color = "grey40",
      column_labels.border.bottom.width = gt::px(2),
      table_body.border.bottom.color = "grey40",
      table.border.top.color = "grey40",
      table.border.top.width = gt::px(2),
      table.border.bottom.width = gt::px(2),
      data_row.padding = gt::px(6)
    )

  # ---- Highlight strong evidence ----
  # Highlight BF > 10 in green, BF < 0.1 in red
  if (include_components) {
    strong_bf_rows <- which(
      purrr
      ::map_lgl(c(bf$effect, bf$heterogeneity, bf$bias), ~ !is.na(.x) && .x > 10)
    )
    weak_bf_rows <- which(
      purrr::map_lgl(c(bf$effect, bf$heterogeneity, bf$bias), ~ !is.na(.x) && .x < 0.1)
    )

    if (length(strong_bf_rows) > 0) {
      gt_tbl <- gt_tbl |>
        gt::tab_style(
          style = gt::cell_fill(color = "#E8F4E8"),
          locations = gt::cells_body(rows = strong_bf_rows)
        )
    }
    if (length(weak_bf_rows) > 0) {
      gt_tbl <- gt_tbl |>
        gt::tab_style(
          style = gt::cell_fill(color = "#FCE4E4"),
          locations = gt::cells_body(rows = weak_bf_rows)
        )
    }
  }

  # ---- Add source note ----
  gt_tbl <- gt_tbl |>
    gt::tab_source_note(
      source_note = gt::md(
        "*BF*: Inclusion Bayes Factor (BF > 10 = strong evidence, BF < 0.1 = strong evidence against). *CrI*: Credible Interval."
      )
    )

  # ---- Return single table or list of tables ----
  if (sum(c(!is.null(estimates_tbl), !is.null(direction_tbl), !is.null(models_tbl))) > 0 &&
      !is.null(components_tbl)) {
    # Return a list of tables for flexible use
    result <- list(
      components = gt_tbl
    )
    if (!is.null(estimates_tbl)) {
      result$estimates <- estimates_tbl |>
        apply_bayesma_style()
    }
    if (!is.null(direction_tbl)) {
      result$direction <- direction_tbl |>
        apply_bayesma_style()
    }
    if (!is.null(models_tbl)) {
      result$models <- models_tbl |>
        apply_bayesma_style()
    }

    # Return combined gt table instead
    return(build_combined_robma_table(x, digits, include_components,
                                      include_estimates, include_direction,
                                      include_models, exponentiate))
  }

  gt_tbl
}


#' Apply consistent bayesma styling to gt table
#' @noRd
apply_bayesma_style <- function(gt_tbl) {
  gt_tbl |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_spanners()
    ) |>
    gt::cols_align(align = "center", columns = dplyr::everything()) |>
    gt::cols_align(align = "left", columns = 1) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      column_labels.font.size = gt::px(11),
      column_labels.border.bottom.color = "grey40",
      column_labels.border.bottom.width = gt::px(2),
      table_body.border.bottom.color = "grey40",
      table.border.top.color = "grey40",
      table.border.top.width = gt::px(2),
      table.border.bottom.width = gt::px(2),
      data_row.padding = gt::px(6)
    )
}


#' Build a single combined gt table for RoBMA results
#' @noRd
build_combined_robma_table <- function(
    x, digits, include_components, include_estimates,
    include_direction, include_models, exponentiate
) {

  method     <- x$meta$method %||% "bridge"
  eff_label  <- x$meta$effect_label
  is_ratio   <- eff_label %in% c("log_or", "log_rr")
  pp         <- x$posterior_probs
  bf         <- x$inclusion_bf
  nrp        <- x$meta$null_range_probs
  nr         <- x$meta$null_range
  pooled     <- dplyr::filter(x$forest_df, .data$type == "pooled")
  n_models   <- x$meta$n_models %||% length(x$component_fits)

  # Format helpers
  fmt_num <- function(val, d = digits) sprintf(paste0("%.", d, "f"), val)

  fmt_bf <- function(val) {
    if (is.na(val)) return("\u2014")
    if (is.infinite(val)) return("\u221E")
    if (val < 0.001) return(formatC(val, format = "e", digits = 1))
    if (val > 1000) return(formatC(val, format = "e", digits = 1))
    fmt_num(val, digits)
  }

  fmt_prob <- function(val, d = digits) {
    if (val < 0.001 && val > 0) return(formatC(val, format = "e", digits = 1))
    fmt_num(val, d)
  }

  # ---- Build unified data frame ----
  rows <- list()

  # Components section
  if (include_components) {
    if (!is.null(x$model_table)) {
      n_effect <- sum(!x$model_table$null_effect)
      n_hetero <- sum(x$model_table$has_heterogeneity)
      n_bias   <- sum(x$model_table$has_bias)
    } else {
      n_effect <- n_hetero <- n_bias <- n_models / 2
    }

    rows <- c(rows, list(
      tibble::tibble(
        section = "Components",
        item = "Effect",
        col1 = sprintf("%d/%d", n_effect, n_models),
        col2 = "0.500",
        col3 = fmt_prob(pp$effect),
        col4 = fmt_bf(bf$effect)
      ),
      tibble::tibble(
        section = "Components",
        item = "Heterogeneity",
        col1 = sprintf("%d/%d", n_hetero, n_models),
        col2 = "0.500",
        col3 = fmt_prob(pp$heterogeneity),
        col4 = fmt_bf(bf$heterogeneity)
      ),
      tibble::tibble(
        section = "Components",
        item = "Publication Bias",
        col1 = sprintf("%d/%d", n_bias, n_models),
        col2 = "0.500",
        col3 = fmt_prob(pp$bias),
        col4 = fmt_bf(bf$bias)
      )
    ))
  }

  # Estimates section
  if (include_estimates) {
    mu_draws <- x$averaged_draws$mu

    rows <- c(rows, list(
      tibble::tibble(
        section = "Estimates",
        item = sprintf("Pooled (%s)", eff_label),
        col1 = fmt_num(mean(mu_draws)),
        col2 = fmt_num(stats::median(mu_draws)),
        col3 = fmt_num(stats::quantile(mu_draws, 0.025)),
        col4 = fmt_num(stats::quantile(mu_draws, 0.975))
      )
    ))

    if (is_ratio && exponentiate) {
      ratio_name <- toupper(gsub("log_", "", eff_label))
      rows <- c(rows, list(
        tibble::tibble(
          section = "Estimates",
          item = ratio_name,
          col1 = fmt_num(mean(exp(mu_draws))),
          col2 = fmt_num(stats::median(exp(mu_draws))),
          col3 = fmt_num(stats::quantile(exp(mu_draws), 0.025)),
          col4 = fmt_num(stats::quantile(exp(mu_draws), 0.975))
        )
      ))
    }
  }

  # Direction section
  if (include_direction && !is.null(nrp)) {
    if (is.null(nr)) {
      rows <- c(rows, list(
        tibble::tibble(
          section = "Direction",
          item = "P(\u03BC < 0)",
          col1 = fmt_prob(nrp$p_negative),
          col2 = "", col3 = "", col4 = ""
        ),
        tibble::tibble(
          section = "Direction",
          item = "P(\u03BC > 0)",
          col1 = fmt_prob(nrp$p_positive),
          col2 = "", col3 = "", col4 = ""
        )
      ))
    } else {
      rows <- c(rows, list(
        tibble::tibble(
          section = "Direction",
          item = "P(harmful)",
          col1 = fmt_prob(nrp$p_negative),
          col2 = "", col3 = "", col4 = ""
        ),
        tibble::tibble(
          section = "Direction",
          item = "P(practically null)",
          col1 = fmt_prob(nrp$p_null),
          col2 = "", col3 = "", col4 = ""
        ),
        tibble::tibble(
          section = "Direction",
          item = "P(beneficial)",
          col1 = fmt_prob(nrp$p_positive),
          col2 = "", col3 = "", col4 = ""
        )
      ))
    }
  }

  combined_df <- dplyr::bind_rows(rows)

  # ---- Build gt table with row groups ----
  gt_tbl <- gt::gt(combined_df, groupname_col = "section") |>
    gt::cols_label(
      item = "",
      col1 = gt::md("**Value 1**"),
      col2 = gt::md("**Value 2**"),
      col3 = gt::md("**Value 3**"),
      col4 = gt::md("**Value 4**")
    ) |>
    gt::tab_header(
      title = "Robust Bayesian Meta-Analysis",
      subtitle = if (method == "bridge") {
        sprintf("Bridge sampling (%d component models)", n_models)
      } else {
        bi <- x$meta$bias_indicator %||% "bias_corrected"
        sprintf("Spike-and-slab (bias: %s)", bi)
      }
    )

  # Update column labels based on section content
  if (include_components && include_estimates) {
    gt_tbl <- gt_tbl |>
      gt::cols_label(
        item = "",
        col1 = "Models / Mean",
        col2 = "Prior / Median",
        col3 = "Post. / 2.5%",
        col4 = "BF / 97.5%"
      )
  } else if (include_components) {
    gt_tbl <- gt_tbl |>
      gt::cols_label(
        item = "",
        col1 = "Models",
        col2 = "Prior Prob.",
        col3 = "Post. Prob.",
        col4 = "Inclusion BF"
      )
  } else if (include_estimates) {
    gt_tbl <- gt_tbl |>
      gt::cols_label(
        item = "",
        col1 = "Mean",
        col2 = "Median",
        col3 = "2.5%",
        col4 = "97.5%"
      )
  }

  # Apply bayesma styling
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_row_groups()
    ) |>
    gt::cols_align(align = "center", columns = c("col1", "col2", "col3", "col4")) |>
    gt::cols_align(align = "left", columns = "item") |>
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(14),
      heading.title.font.weight = "bold",
      heading.subtitle.font.size = gt::px(11),
      column_labels.font.size = gt::px(11),
      column_labels.border.bottom.color = "grey40",
      column_labels.border.bottom.width = gt::px(2),
      table_body.border.bottom.color = "grey40",
      table.border.top.color = "grey40",
      table.border.top.width = gt::px(2),
      table.border.bottom.width = gt::px(2),
      row_group.border.top.width = gt::px(1),
      row_group.border.top.color = "grey70",
      row_group.padding = gt::px(8),
      data_row.padding = gt::px(6)
    ) |>
    gt::tab_source_note(
      source_note = gt::md(
        "*BF*: Inclusion Bayes Factor (> 10 strong evidence for, < 0.1 strong against). Model-averaged estimates from posterior draws."
      )
    )

  gt_tbl
}
