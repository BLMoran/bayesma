#' Internal function to resolve default parameters
#'
#' @noRd
default_pars <- function(object, draws_array = NULL) {

  stage      <- object$meta$stage %||% "unknown"
  model_type <- object$meta$model_type %||% "random_effect"

  # Get available parameters from draws if provided, otherwise from summary
  if (!is.null(draws_array)) {
    all_pars <- dimnames(draws_array)[[3]]
  } else {
    all_pars <- object$summary$variable
  }

  # Filter to non-internal parameters
  model_pars <- all_pars[!grepl("^(lp__|lf__|log_lik|__|y_rep\\[)", all_pars)]

  # Prefer exact scalar parameters first, then indexed versions as fallback
  result <- character()

  # mu — always want this

  if ("mu" %in% model_pars) {
    result <- c(result, "mu")
  }

  # tau — only for random effects / bias-corrected models
  if ("tau" %in% model_pars && model_type != "common_effect") {
    result <- c(result, "tau")
  }

  # Model-specific extras
  if ("nu" %in% model_pars)          result <- c(result, "nu")
  if ("alpha_skew" %in% model_pars)  result <- c(result, "alpha_skew")
  if ("pi_main" %in% model_pars)     result <- c(result, "pi_main")
  if ("B" %in% model_pars)           result <- c(result, "B")
  if ("p_bias" %in% model_pars)      result <- c(result, "p_bias")
  if ("rho" %in% model_pars)         result <- c(result, "rho")
  if ("beta_bias" %in% model_pars)   result <- c(result, "beta_bias")
  if ("gamma0" %in% model_pars)      result <- c(result, "gamma0")
  if ("gamma1" %in% model_pars)      result <- c(result, "gamma1")

  # Cap at 4 for readability
  if (length(result) > 0) {
    return(utils::head(result, 4))
  }

  # Fallback: first few available (excluding indexed arrays with many elements)
  scalar_pars <- model_pars[!grepl("\\[", model_pars)]
  if (length(scalar_pars) > 0) {
    return(utils::head(scalar_pars, 4))
  }

  utils::head(model_pars, 4)
}

#' Internal function to extract numeric diagnostics from CmdStanMCMC
#'
#' @noRd
extract_diagnostics <- function(object) {

  summary_df <- object$summary

  # Rhat
  rhat_vals   <- summary_df$rhat
  max_rhat    <- max(rhat_vals, na.rm = TRUE)
  n_high_rhat <- sum(rhat_vals > 1.01, na.rm = TRUE)
  rhat_ok     <- max_rhat <= 1.01

  # ESS (use bulk; compute ratio from total draws)
  total_draws <- nrow(object$draws)
  ess_bulk    <- summary_df$ess_bulk
  ess_ratio   <- ess_bulk / total_draws
  min_ess     <- min(ess_ratio, na.rm = TRUE)
  n_low_ess   <- sum(ess_ratio < 0.1, na.rm = TRUE)
  ess_ok      <- min_ess >= 0.1

  # CmdStan diagnostic summary (per chain)
  diag_summary <- object$fit$diagnostic_summary()

  n_divergent  <- sum(diag_summary$num_divergent)
  divergent_ok <- n_divergent == 0

  max_td       <- max(diag_summary$num_max_treedepth)
  treedepth_ok <- max_td == 0

  ebfmi        <- diag_summary$ebfmi
  min_ebfmi    <- min(ebfmi)
  ebfmi_ok     <- min_ebfmi > 0.3

  all_ok <- rhat_ok && ess_ok && divergent_ok && treedepth_ok && ebfmi_ok

  list(
    max_rhat           = max_rhat,
    n_high_rhat        = n_high_rhat,
    rhat_ok            = rhat_ok,
    min_ess_ratio      = min_ess,
    n_low_ess          = n_low_ess,
    ess_ok             = ess_ok,
    n_divergent        = n_divergent,
    divergent_ok       = divergent_ok,
    n_max_treedepth    = max_td,
    treedepth_ok       = treedepth_ok,
    min_ebfmi          = min_ebfmi,
    ebfmi_ok           = ebfmi_ok,
    ebfmi_per_chain    = ebfmi,
    all_ok             = all_ok
  )
}

#' Internal function to construct a consistent base theme
#'
#' @noRd
diag_theme <- function() {
  ggplot2::theme_minimal(base_family = "") +
    ggplot2::theme(
      text          = ggplot2::element_text(family = ""),
      plot.title    = ggplot2::element_text(face = "bold", size = 10),
      axis.text.y   = ggplot2::element_text(size = 7),
      axis.text.x   = ggplot2::element_text(size = 7),
      strip.text    = ggplot2::element_text(size = 8),
      legend.text   = ggplot2::element_text(size = 8)
    )
}

#' Internal functions to construct panel builders
#' Trace Plot
#'
#' @noRd
panel_trace <- function(draws_array, pars) {

  available <- dimnames(draws_array)[[3]]
  pars_use  <- intersect(pars, available)

  if (length(pars_use) == 0) {
    pars_use <- available[!grepl("^(lp__|lf__|log_lik|__|y_rep\\[)", available)]
    pars_use <- utils::head(pars_use, 4)
  }

  if (length(pars_use) == 0) {
    return(.empty_panel("Trace Plots", "No matching parameters found"))
  }

  tryCatch(
    bayesplot::mcmc_trace(draws_array, pars = pars_use) +
      ggplot2::ggtitle("Trace Plots") +
      diag_theme(),
    error = function(e) {
      empty_panel("Trace Plots", conditionMessage(e))
    }
  )
}

#' PP_Check
#'
#' @noRd
panel_ppc <- function(object, ndraws = 100) {
  tryCatch({
    stage      <- object$meta$stage %||% "two_stage"
    likelihood <- object$meta$likelihood %||% "gaussian"

    y <- extract_observed(object, stage, likelihood)
    if (is.null(y)) {
      return(empty_panel("Posterior Predictive Check", "No observed data found"))
    }

    all_vars   <- posterior::variables(object$fit$draws())
    y_rep_vars <- grep("^y_rep\\[", all_vars, value = TRUE)

    if (length(y_rep_vars) == 0) {
      return(empty_panel("Posterior Predictive Check",
                          "No y_rep in draws.\nAdd generated quantities block."))
    }

    y_rep <- as.matrix(
      posterior::as_draws_matrix(object$fit$draws(variables = y_rep_vars))
    )

    if (nrow(y_rep) > ndraws) {
      idx   <- sample.int(nrow(y_rep), ndraws)
      y_rep <- y_rep[idx, , drop = FALSE]
    }

    if (ncol(y_rep) != length(y)) {
      return(empty_panel("Posterior Predictive Check",
                          sprintf("Dimension mismatch: y=%d, y_rep=%d",
                                  length(y), ncol(y_rep))))
    }

    # Always use density overlay
    y_dens <- stats::density(y)
    y_df <- data.frame(x = y_dens$x, y = y_dens$y)

    rep_list <- lapply(seq_len(nrow(y_rep)), function(i) {
      d <- stats::density(y_rep[i, ], from = min(y_dens$x), to = max(y_dens$x))
      data.frame(x = d$x, y = d$y, draw = paste0("rep_", i))
    })
    rep_df <- do.call(rbind, rep_list)

    ggplot2::ggplot() +
      ggplot2::geom_line(
        data = rep_df,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$draw),
        colour = "skyblue", alpha = 0.3, linewidth = 0.3
      ) +
      ggplot2::geom_line(
        data = y_df,
        ggplot2::aes(x = .data$x, y = .data$y),
        colour = "black", linewidth = 1
      ) +
      ggplot2::labs(
        title = "Posterior Predictive Check",
        x = NULL, y = "Density"
      ) +
      diag_theme()
  },
  error = function(e) {
    empty_panel("Posterior Predictive Check", conditionMessage(e))
  })
}

#' Rhat
#'
#' @noRd
panel_rhat <- function(summary_df) {

  plot_df <- summary_df |>
    dplyr::filter(!grepl("^(lp__|lf__|log_lik)", .data$variable)) |>
    dplyr::mutate(
      variable = forcats::fct_reorder(.data$variable, .data$rhat),
      flag     = .data$rhat > 1.01
    )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$rhat, y = .data$variable)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$flag), size = 2) +
    ggplot2::geom_vline(xintercept = 1.01, linetype = "dashed", colour = "red") +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "red"),
      guide  = "none"
    ) +
    ggplot2::labs(
      title = "Rhat",
      x     = expression(hat(R)),
      y     = NULL
    ) +
    diag_theme()
}

#' ESS
#'
#' @noRd
panel_neff <- function(summary_df) {

  plot_df <- summary_df |>
    dplyr::filter(!grepl("^(lp__|lf__|log_lik)", .data$variable)) |>
    dplyr::mutate(
      variable = forcats::fct_reorder(.data$variable, .data$ess_bulk),
      flag     = .data$ess_bulk < 400
    )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$ess_bulk, y = .data$variable)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$flag), size = 2) +
    ggplot2::geom_vline(xintercept = 400, linetype = "dashed", colour = "red") +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "red"),
      guide  = "none"
    ) +
    ggplot2::labs(
      title = "Effective Sample Size (Bulk)",
      x     = "ESS bulk",
      y     = NULL
    ) +
    diag_theme()
}

#' ACF Plot
#'
#' @noRd
panel_acf <- function(draws_array, pars) {

  available <- dimnames(draws_array)[[3]]
  pars_use  <- intersect(pars, available)

  if (length(pars_use) == 0) {
    pars_use <- available[!grepl("^(lp__|lf__|log_lik|__|y_rep\\[)", available)]
    pars_use <- utils::head(pars_use, 4)
  }

  if (length(pars_use) == 0) {
    return(.empty_panel("Autocorrelation", "No matching parameters found"))
  }

  tryCatch(
    bayesplot::mcmc_acf(draws_array, pars = pars_use) +
      ggplot2::ggtitle("Autocorrelation") +
      diag_theme(),
    error = function(e) {
      .empty_panel("Autocorrelation", conditionMessage(e))
    }
  )
}

#' Summary Table
#'
#' @noRd
panel_diagnostics_table <- function(diagnostics) {

  ebfmi_str <- paste0(
    "min = ", sprintf("%.3f", diagnostics$min_ebfmi),
    " [", paste(sprintf("%.3f", diagnostics$ebfmi_per_chain), collapse = ", "), "]"
  )

  ok_flags <- c(
    diagnostics$rhat_ok, diagnostics$ess_ok, diagnostics$divergent_ok,
    diagnostics$treedepth_ok, diagnostics$ebfmi_ok, diagnostics$all_ok
  )

  tbl_df <- tibble::tibble(
    Diagnostic = c(
      "Rhat (all < 1.01)",
      "ESS ratio (all >= 0.1)",
      "Divergent transitions",
      "Max treedepth hits",
      "E-BFMI (all > 0.3)",
      "Overall"
    ),
    Value = c(
      sprintf("max = %.4f (%d flagged)", diagnostics$max_rhat, diagnostics$n_high_rhat),
      sprintf("min = %.4f (%d flagged)", diagnostics$min_ess_ratio, diagnostics$n_low_ess),
      as.character(diagnostics$n_divergent),
      as.character(diagnostics$n_max_treedepth),
      ebfmi_str,
      dplyr::if_else(diagnostics$all_ok, "All checks passed", "Issues detected")
    ),
    Status = dplyr::if_else(ok_flags, "Pass", "Fail")
  )

  gt_tbl <- gt::gt(tbl_df) |>
    gt::tab_header(title = "MCMC Diagnostics Summary") |>
    gt::tab_style(
      style = gt::cell_text(color = "forestgreen", weight = "bold"),
      locations = gt::cells_body(
        columns = "Status",
        rows = ok_flags
      )
    ) |>
    gt::tab_style(
      style = gt::cell_text(color = "red", weight = "bold"),
      locations = gt::cells_body(
        columns = "Status",
        rows = !ok_flags
      )
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        columns = "Diagnostic",
        rows = Diagnostic == "Overall"
      )
    ) |>
    gt::tab_style(
      style = gt::cell_fill(color = "grey95"),
      locations = gt::cells_body(rows = c(1, 3, 5))
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    ) |>
    gt::tab_options(
      table.font.size = gt::px(11),
      heading.title.font.size = gt::px(13),
      heading.title.font.weight = "bold",
      column_labels.border.bottom.color = "grey40",
      table_body.border.bottom.color = "grey40",
      table.border.top.color = "grey40"
    )

  patchwork::wrap_table(gt_tbl, panel = "full", space = "fixed")
}


#' Empty Panel
#'
#' @noRd
empty_panel <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = message,
                      size = 4, colour = "grey50") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_void(base_family = "") +
    ggplot2::theme(plot.title = ggplot2::element_text(family = "", face = "bold", size = 10))
}
