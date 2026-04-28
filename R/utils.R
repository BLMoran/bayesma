utils::globalVariables(c(
  ".avg", ".ci_hi", ".ci_lo", ".df", ".empty_panel", ".es", ".es_avg",
  ".last_in_group", ".lower", ".mod", ".n", ".row_j", ".sample",
  ".sample_id", ".sample_mean", ".scaled_se", ".se", ".se_avg",
  ".study", ".upper", ".width", ":=",
  "Author", "Author_label", "Author_ordered", "Author_original",
  "Author_pooled", "Code", "D1", "Description", "Diagnostic", "Overall",
  "Subgroup", "Year", "a", "b", "b_Intercept", "category", "ci_fmt",
  "ctrl_mean_sd", "density", "domain", "effect_draws", "estimate",
  "estimate_fmt", "int_mean_sd", "is_null_draw", "judgement", "k", "l95",
  "level", "level_num", "m", "max_level", "mean_effect", "mean_w", "med",
  "med_se", "model_type", "mu", "n_missing", "null_line_to_use",
  "p_outside_null", "pd_fmt", "pi_lower", "pi_median", "pi_upper",
  "plot_row", "pr_benefit_gt_delta", "pr_harm", "pr_harm_gt_delta",
  "pred", "prior", "prior_label", "qnorm", "r_Author", "rho",
  "sd_Author__Intercept", "section_label", "sei", "side", "sigma",
  "sigma_median", "spread_df", "study", "study_count",
  "subgroup.forest.summary", "subgroup_model", "tau", "term", "u95",
  "unweighted_effect", "vi", "w", "w_norm", "x", "x_plot", "x_studies",
  "x_val", "xdist", "y", "yi"
))

# ---- CmdStan summary helper ----

#' @noRd
stan_summary <- function(fit, variables = NULL) {
  fit$summary(
    variables = variables,
    "median", "mean", "mad", "sd",
    ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
    "ess_bulk", "ess_tail", "rhat"
  ) |> tibble::as_tibble()
}

# ---- CmdStan model cache (compile once per Stan code hash) ----
.bayesma_cmdstan_cache <- new.env(parent = emptyenv())

#' @noRd
get_cmdstan_model_cached <- function(stan_code, compile = TRUE,
                                     force_recompile = FALSE,
                                     compile_model_methods = FALSE, ...) {
  if (!compile) return(NULL)

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cli::cli_abort("{.pkg cmdstanr} is required to compile Stan models.")
  }

  # Use digest if available; otherwise fall back to no caching
  has_digest <- requireNamespace("digest", quietly = TRUE)
  # Include compile_model_methods in cache key so bridge-sampling-ready

  # models are cached separately from standard ones
  key <- if (has_digest) {
    digest::digest(list(stan_code, compile_model_methods = compile_model_methods))
  } else {
    NULL
  }

  if (!force_recompile && !is.null(key) && exists(key, envir = .bayesma_cmdstan_cache, inherits = FALSE)) {
    return(get(key, envir = .bayesma_cmdstan_cache, inherits = FALSE))
  }

  stan_file <- cmdstanr::write_stan_file(stan_code)
  mod <- cmdstanr::cmdstan_model(stan_file,
                                 compile_model_methods = compile_model_methods, ...)

  if (!is.null(key)) {
    assign(key, mod, envir = .bayesma_cmdstan_cache)
  }

  mod
}


#' Internal function to disambiguate duplicate Author names
#'
#' When multiple studies share the same Author name, this function appends a
#' letter suffix to make each row unique. Two layers of disambiguation are
#' applied:
#'
#' 1. **Display-level (Year suffix):** When the same Author appears with the
#'    same Year (e.g. multi-arm trials), a letter is appended to the Year
#'    column (e.g. "1999" -> "1999a", "1999b"). This ensures the tables show
#'    "Smith (1999a)" and "Smith (1999b)".
#'
#' 2. **Internal-level (Author suffix):** When the same Author name appears
#'    more than once (regardless of Year), the Author column itself is made
#'    unique by appending a letter (e.g. "Smith" -> "Smith_a", "Smith_b").
#'    This is needed so that ggplot2 maps each study to its own row on the
#'    y-axis.
#'
#' The original Author name is preserved in a new column `Author_original` so
#' it can be restored for display purposes later.
#'
#' @param data A data frame containing at least `Author` and `Year` columns.
#'
#' @return The data frame with unique `Author` values, disambiguated `Year`
#'   values where needed, and a new `Author_original` column containing the
#'   original (possibly duplicated) names.
#'
#' @noRd
make_authors_unique <- function(data) {
  # Preserve the original name for display in tables
  data$Author_original <- data$Author

  # --- Step 1: Disambiguate Year for same Author + same Year combos ---
  data$Year <- as.character(data$Year)

  author_year_combos <- paste(data$Author, data$Year, sep = "|||")
  dup_combos <- author_year_combos[duplicated(author_year_combos)]

  if (length(dup_combos) > 0) {
    for (combo in unique(dup_combos)) {
      idx <- which(author_year_combos == combo)
      suffixes <- letters[seq_along(idx)]
      original_year <- strsplit(combo, "|||", fixed = TRUE)[[1]][2]
      data$Year[idx] <- paste0(original_year, suffixes)
    }
  }

  # --- Step 2: Disambiguate Author names that appear more than once ---
  dup_authors <- data$Author[duplicated(data$Author)]

  if (length(dup_authors) > 0) {
    for (author_name in unique(dup_authors)) {
      idx <- which(data$Author == author_name)
      suffixes <- letters[seq_along(idx)]
      data$Author[idx] <- paste0(author_name, "_", suffixes)
    }
  }

  return(data)
}

#' Internal function to sort studies
#'
#' @noRd
sort_studies_fn <- function(.data, sort_studies_by = NULL) {
  sort_by <- if (is.null(sort_studies_by)) "author" else sort_studies_by
  has_subgroup <- "Subgroup" %in% names(.data) &&
    !all(is.na(.data$Subgroup)) &&
    dplyr::n_distinct(.data$Subgroup) > 1

  # Separate rows (Pooled Effect & Prediction) from actual studies.
  pooled <- .data |> dplyr::filter(Author == "Pooled Effect")
  prediction <- .data |> dplyr::filter(Author == "Prediction")
  studies <- .data |> dplyr::filter(!Author %in% c("Pooled Effect", "Prediction"))

  # Arrange studies
  studies <- switch(
    sort_by,
    "author" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), Author),
    "year" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), Year, Author),
    "effect" = studies |> dplyr::arrange(dplyr::across(any_of("Subgroup")), yi, Author),
    cli::cli_abort("Invalid value for {.arg sort_studies_by}. Must be one of {.val author}, {.val year}, or {.val effect}.")
  )

  # Recombine: studies, then Pooled Effect, then Prediction
  full <- dplyr::bind_rows(studies, pooled, prediction)

  full <- full |>
    dplyr::mutate(
      Author = factor(Author, levels = unique(Author))
    )

  return(full)
}

#' Get subgroup order based on sorting preference
#'
#' @noRd
get_subgroup_order <- function(data, sort_subgroup_by) {
  if (!"Subgroup" %in% names(data)) {
    return(NULL)
  }

  unique_subgroups <- unique(data$Subgroup[!is.na(data$Subgroup)])

  if (is.character(sort_subgroup_by)) {
    if (length(sort_subgroup_by) == 1) {
      subgroup_order <- switch(
        sort_subgroup_by,
        "alphabetical" = sort(unique_subgroups),
        "effect" = {
          data |>
            dplyr::summarise(mean_effect = mean(yi, na.rm = TRUE), .by = Subgroup) |>
            dplyr::arrange(mean_effect) |>
            dplyr::pull(Subgroup)
        },
        cli::cli_abort(
          "Invalid {.arg sort_subgroup_by} value: {.val {sort_subgroup_by}}. Must be {.val alphabetical}, {.val effect}, or a character vector of subgroup names."
        )
      )
    } else {
      missing_groups <- setdiff(sort_subgroup_by, unique_subgroups)
      if (length(missing_groups) > 0) {
        rlang::warn(
          paste0("These subgroups in sort_subgroup_by were not found in data: ",
                 paste(missing_groups, collapse = ", "))
        )
      }
      subgroup_order <- sort_subgroup_by
    }
  } else {
    cli::cli_abort(
      "{.arg sort_subgroup_by} must be a character vector. Use {.val alphabetical}, {.val effect}, or provide a vector of subgroup names."
    )
  }

  c(subgroup_order, "Overall")
}


#' Detect whether a model has random effects (i.e., is not a common/fixed effect model)
#'
#' Checks `model$meta$model_type` for bayesma objects.
#'
#' @param model A bayesma object.
#'
#' @return Logical. TRUE if the model has random effects, FALSE otherwise.
#'
#' @noRd
has_random_effects <- function(model) {
  if (!inherits(model, "bayesma")) {
    cli::cli_abort("{.arg model} must be a {.cls bayesma} object.")
  }
  model$meta$model_type != "common_effect"
}


#' Extract fixed effect summary from a bayesma model
#'
#' Returns a 1x4 matrix with columns Estimate, Est.Error, Q2.5, Q97.5.
#'
#' @param model A bayesma object.
#'
#' @return A 1x4 matrix.
#'
#' @noRd
extract_fixef <- function(model) {
  if (!inherits(model, "bayesma")) {
    cli::cli_abort("{.arg model} must be a {.cls bayesma} object.")
  }

  # Extract mu draws from the cmdstanr fit
  mu_draws <- as.numeric(
    posterior::subset_draws(model$fit$draws("mu"), variable = "mu")
  )

  mat <- matrix(
    c(
      stats::median(mu_draws),
      stats::mad(mu_draws),
      stats::quantile(mu_draws, 0.025),
      stats::quantile(mu_draws, 0.975)
    ),
    nrow = 1,
    dimnames = list("Intercept", c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  )

  mat
}


# ============================================================================
# refit_bayesma: refit a bayesma model on new data
# ============================================================================

#' Refit a bayesma model on new (subset) data
#'
#' Uses the stored call_args from the original fit to re-run bayesma()
#' with different data. Used internally for subgroup analysis and
#' RoB exclusion.
#'
#' @param model A bayesma object with stored call_args.
#' @param newdata A data frame (subset of original).
#'
#' @return A new bayesma object.
#'
#' @noRd
refit_bayesma <- function(model, newdata) {
  if (is.null(model$meta$call_args)) {
    cli::cli_abort(
      c(
        "Cannot refit: bayesma object does not contain stored call arguments.",
        "i" = "Refit your model with the latest version of bayesma to store call_args."
      )
    )
  }

  args <- model$meta$call_args
  args$data <- newdata
  do.call(bayesma::bayesma, args)
}

#' Refit a bayesma model with updated arguments
#'
#' Starts from model$meta$call_args and overrides any provided arguments.
#' Ensures data is passed in (since call_args intentionally doesn't store it).
#'
#' @param model A bayesma object with stored call_args.
#' @param data Data frame to refit on (required unless stored elsewhere).
#' @param ... Named overrides (e.g., model_type, mu_prior, tau_prior, robust, stage).
#'
#' @return A bayesma object.
#'
#' @noRd
refit_bayesma_update <- function(model, data, ...) {
  if (!inherits(model, "bayesma")) {
    cli::cli_abort("{.fn refit_bayesma_update} requires a {.cls bayesma} object.")
  }
  if (is.null(model$meta$call_args)) {
    cli::cli_abort(c(
      "Cannot refit: bayesma object does not contain stored call arguments.",
      "i" = "Ensure {.fn bayesma} stores {.code meta$call_args}."
    ))
  }
  if (missing(data) || is.null(data)) {
    cli::cli_abort("{.arg data} must be supplied for refitting.")
  }

  args <- model$meta$call_args
  args$data <- data

  dots <- list(...)
  if (length(dots)) {
    for (nm in names(dots)) args[[nm]] <- dots[[nm]]
  }

  do.call(bayesma::bayesma, args)
}

#' Internal function to join sections of forest plot
#'
#' @noRd
patchwork_fn <- function(table.left,
                         study.density.plot,
                         table.right,
                         add_rob_legend = FALSE,
                         rob_tool = "rob2",
                         plot_width = NULL,
                         title = NULL,
                         subtitle = NULL,
                         title_align = "center",
                         font = NULL){
  if (is.null(plot_width)){
    plot_width <- 4
  }

  # Create the base plot layout
  if (isFALSE(add_rob_legend)) {
    base_plot <- patchwork::wrap_table(table.left, space = "fixed") +
      study.density.plot +
      patchwork::wrap_table(table.right, space = "fixed") +
      patchwork::plot_layout(widths = grid::unit(c(-1, plot_width, -1),
                                                 c("null", "cm", "null")))
  } else {
    rob_legend <- rob_legend_fn(rob_tool, font)
    # Include rob_legend to the right of table.right
    base_plot <- patchwork::wrap_table(table.left, space = "fixed") +
      study.density.plot +
      patchwork::wrap_table(table.right, space = "fixed") +
      patchwork::wrap_table(rob_legend, space = "fixed") +
      patchwork::plot_layout(widths = grid::unit(c(-1, plot_width, -1, -1),
                                                 c("null", "cm", "null", "null")))
  }

  # Add title and/or subtitle if provided
  if (!is.null(title) || !is.null(subtitle)) {
    # Set alignment value based on title_align parameter
    hjust_val <- switch(title_align,
                        "left" = 0,
                        "center" = 0.5,
                        "centre" = 0.5,
                        "right" = 1,
                        0.5) # default to center if invalid input

    # Create theme elements
    title_theme <- ggplot2::theme()

    if (!is.null(title)) {
      title_theme <- title_theme +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 16,
            face = "bold",
            hjust = hjust_val,
            margin = ggplot2::margin(b = if (is.null(subtitle)) 10 else 5),
            family = font))
    }

    if (!is.null(subtitle)) {
      title_theme <- title_theme +
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(
            size = 14,
            hjust = hjust_val,
            margin = ggplot2::margin(b = 10),
            color = "gray30",
            family = font))
    }

    base_plot <- base_plot +
      patchwork::plot_annotation(
        title = title,
        subtitle = subtitle,
        theme = title_theme
      )
  }

  return(base_plot)
}


#' Get recommended figure height for a bayesma plot
#'
#' Extracts the recommended figure height (in inches) stored by
#' \code{\link{bayesma::forest}} or \code{\link{sensitivity_plot}}.
#' Useful for setting \code{fig-height} dynamically in Quarto/Rmd chunks.
#'
#' @param plot A patchwork object returned by \code{bayesma::forest} or
#'   \code{sensitivity_plot}.
#' @param default Fallback height if the attribute is not found.
#'
#' @return A numeric scalar (inches).
#'
#' @examples
#' \dontrun{
#' p <- bayesma::forest(model, data, measure = "OR", studyvar = Author)
#' # In a Quarto chunk header:
#' # fig-height: !expr bayesma::fig_height(p)
#' }
#'
#' @noRd
#' @keywords internal
fig_height <- function(plot, default = 8) {
  attr(plot, "recommended_height") %||% default
}


#' Get recommended figure width for a bayesma plot
#'
#' Extracts the recommended figure width (in inches) stored by
#' \code{\link{bayesma::forest}} or \code{\link{sensitivity_plot}}.
#'
#' @inheritParams fig_height
#'
#' @return A numeric scalar (inches).
#' @noRd
#' @keywords internal
fig_width <- function(plot, default = 14) {
  attr(plot, "recommended_width") %||% default
}


#' Calculate recommended figure height from number of rows
#'
#' Computes an appropriate figure height (in inches) for forest plots and
#' sensitivity plots based on the number of study rows. Use this in Quarto
#' chunk headers when the plot object is not yet available.
#'
#' @param n_rows Number of rows (studies + pooled effect + subgroup headers).
#' @param has_title Logical; does the plot have a title/subtitle?
#' @param per_row Height per row in inches. Default 0.55.
#' @param base Baseline height for axes, labels, etc. Default 1.5.
#'
#' @return A numeric scalar (inches).
#'
#' @examples
#' \dontrun{
#' # In a Quarto chunk header:
#' #| fig-height: !expr bayesma::calc_fig_height(nrow(my_data) + 1)
#' }
#'
#' @noRd
#' @keywords internal
calc_fig_height <- function(n_rows, has_title = TRUE, per_row = 0.55, base = 1.5) {
  base + n_rows * per_row + if (has_title) 0.8 else 0
}


#' Internal function to get certain properties dependent on measure
#'
#' @noRd
get_measure_properties <- function(measure) {
  switch(measure,
         "OR" = list(
           null_value = 1,
           log_scale = TRUE,
           x_label = "Odds Ratio (log scale)",
           data_cols = c("Event_Control", "N_Control", "Event_Intervention", "N_Intervention")
         ),
         "HR" = list(
           null_value = 1,
           log_scale = TRUE,
           x_label = "Hazard Ratio (log scale)",
           data_cols = c("Event_Control", "N_Control", "Event_Intervention", "N_Intervention")
         ),
         "RR" = list(
           null_value = 1,
           log_scale = TRUE,
           x_label = "Risk Ratio (log scale)",
           data_cols = c("Event_Control", "N_Control", "Event_Intervention", "N_Intervention")
         ),
         "IRR" = list(
           null_value = 1,
           log_scale = TRUE,
           x_label = "Incident Rate Ratio (log scale)",
           data_cols = c("Event_Control", "N_Control", "Time_Control", "Event_Intervention", "N_Intervention", "Time_Intervention")
         ),
         "MD" = list(
           null_value = 0,
           log_scale = FALSE,
           x_label = "Mean Difference",
           data_cols = c("Mean_Control", "SD_Control", "N_Control", "Mean_Intervention", "SD_Intervention", "N_Intervention")
         ),
         "SMD" = list(
           null_value = 0,
           log_scale = FALSE,
           x_label = "Standardised Mean Difference",
           data_cols = c("Mean_Control", "SD_Control", "N_Control", "Mean_Intervention", "SD_Intervention", "N_Intervention")
         ),
         "RD" = ,
         "ARR" = list(
           null_value = 0,
           log_scale = FALSE,
           x_label = "Risk Difference",
           data_cols = NULL
         ),
         "ATE" = list(
           null_value = 0,
           log_scale = FALSE,
           x_label = "Average Treatment Effect",
           data_cols = NULL
         ),
         "ATT" = list(
           null_value = 0,
           log_scale = FALSE,
           x_label = "Average Treatment Effect on the Treated",
           data_cols = NULL
         ),
         "CATE" = list(
           null_value = 0,
           log_scale = FALSE,
           x_label = "Conditional Average Treatment Effect",
           data_cols = NULL
         ),
         cli::cli_abort("Effect size must be one of: {.val OR}, {.val HR}, {.val RR}, {.val IRR}, {.val MD}, {.val SMD}, {.val RD}, {.val ARR}, {.val ATE}, {.val ATT}, {.val CATE}.")
  )
}

#' Map a marginal estimand to the underlying model-scale measure
#'
#' Used by forest/funnel plots to display individual study effects on the
#' native model scale (log-OR, MD, IRR) rather than the marginal scale.
#'
#' @noRd
underlying_measure <- function(estimand, likelihood) {
  if (!is_marginal_estimand(estimand)) return(estimand)
  switch(likelihood,
    binomial = "OR",
    poisson  = "IRR",
    gaussian = "MD"
  )
}

#' Extract Tau Prior from Model
#'
#' @description
#' Extracts the prior specification for the heterogeneity parameter (tau)
#' from a bayesma model.
#'
#' @param model A fitted bayesma object.
#'
#' @return Character string of the tau prior, or NA if not found.
#'
#' @noRd
extract_tau_from_model <- function(model) {

  if (!inherits(model, "bayesma")) {
    return(NA_character_)
  }

  # Extract tau prior from stored priors in meta
  tau_prior <- model$meta$priors$tau
  if (is.null(tau_prior)) return(NA_character_)

  # Convert bayesma_prior to a string representation
  format(tau_prior)
}

#' Convert Prior String to Unicode Format
#'
#' @description
#' Converts prior distribution strings to Unicode mathematical notation
#' for pretty printing in tables. Works with bayesma prior format strings.
#'
#' @param prior_string Character string representing a prior distribution.
#'
#' @return Unicode-formatted prior string.
#'
#' @details
#' Supports conversion of:
#' \itemize{
#'   \item N() to 𝒩() (normal)
#'   \item HN() to 𝒩⁺() (half-normal)
#'   \item C() to 𝒞() (cauchy)
#'   \item HC() to 𝒞⁺() (half-cauchy)
#'   \item t() to 𝓉() (student-t)
#'   \item Exp() to ℰ() (exponential)
#' }
#'
#' @noRd
prior_to_unicode <- function(prior_string) {
  if (is.null(prior_string) || is.na(prior_string)) {
    return(NA_character_)
  }

  # Unicode math-script letters (STIX Two Math)
  script_N <- "\U0001D4A9"
  script_C <- "\U0001D49E"
  script_t <- "\U0001D4C9"
  script_E <- "\U2130"

  prior_string <- trimws(prior_string)

  # ---- bayesma format strings ----
  # HN(mean, sd) -> half-normal
  if (grepl("^HN\\s*\\(", prior_string)) {
    params <- sub("^HN\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_N, "\u207A(", params, ")\u2009"))
  }

  # HC(location, scale) -> half-cauchy
  if (grepl("^HC\\s*\\(", prior_string)) {
    params <- sub("^HC\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_C, "\u207A(", params, ")\u2009"))
  }

  # N(mean, sd) -> normal
  if (grepl("^N\\s*\\(", prior_string)) {
    params <- sub("^N\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_N, "(", params, ")\u2009"))
  }

  # C(location, scale) -> cauchy
  if (grepl("^C\\s*\\(", prior_string)) {
    params <- sub("^C\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_C, "(", params, ")\u2009"))
  }

  # t(df, location, scale) -> half-student-t
  if (grepl("^t\\s*\\(", prior_string)) {
    params <- sub("^t\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    parts <- strsplit(params, ",")[[1]]
    parts <- trimws(parts)
    if (length(parts) == 3) {
      df  <- gsub("\\s+", "", parts[1])
      loc <- gsub("\\s+", "", parts[2])
      scl <- gsub("\\s+", "", parts[3])
      if (nchar(df) == 1 && grepl("^[0-9]$", df)) {
        df <- c(
          "0" = "\u2080", "1" = "\u2081", "2" = "\u2082", "3" = "\u2083",
          "4" = "\u2084", "5" = "\u2085", "6" = "\u2086", "7" = "\u2087",
          "8" = "\u2088", "9" = "\u2089"
        )[df]
      }
      return(paste0(script_t, df, "(", loc, ", ", scl, ")\u2009"))
    }
  }

  # Exp(rate) -> exponential
  if (grepl("^Exp\\s*\\(", prior_string)) {
    params <- sub("^Exp\\s*\\(([^)]+)\\).*$", "\\1", prior_string)
    params <- gsub("\\s+", "", params)
    return(paste0(script_E, "(", params, ")\u2009"))
  }

  prior_string
}

#' Extract Mu and Tau Priors
#'
#' @description
#' Extracts and formats the prior specifications for both the mean (mu) and
#' heterogeneity (tau) parameters from a list of priors. Supports bayesma
#' priors (named lists of bayesma_prior objects).
#'
#' @param priors Named list of prior specifications.
#' @param model A fitted bayesma object (for extracting default tau prior).
#'
#' @return A tibble with columns: prior, prior_label, mu_prior_unicode, tau_prior_unicode.
#'
#' @noRd
extract_mu_tau_priors <- function(priors, model) {

  model_tau <- extract_tau_from_model(model)

  purrr::imap(priors, function(prior_obj, prior_name) {

    # bayesma priors: named list with bayesma_prior objects
    mu_prior_str <- if (!is.null(prior_obj$mu_prior)) {
      format(prior_obj$mu_prior)
    } else {
      NA_character_
    }

    tau_prior_str <- if (!is.null(prior_obj$tau_prior)) {
      format(prior_obj$tau_prior)
    } else {
      NA_character_
    }

    tibble::tibble(
      prior = prior_name,
      prior_label = dplyr::recode(
        prior_name,
        vague = "Vague",
        weakreg = "Weakly Regularising",
        informative = "Informative",
        .default = stringr::str_to_sentence(prior_name)
      ),
      mu_prior_unicode =
        if (!is.na(mu_prior_str)) prior_to_unicode(mu_prior_str) else NA_character_,
      tau_prior_unicode =
        if (!is.na(tau_prior_str)) {
          prior_to_unicode(tau_prior_str)
        } else if (!is.na(model_tau)) {
          prior_to_unicode(model_tau)
        } else {
          NA_character_
        }
    )
  }) |> purrr::list_rbind()
}

#' Apply Math Font to GT Table
#'
#' @description
#' Applies a mathematical font (like STIX Two Math) to specific columns
#' in a gt table for proper rendering of mathematical symbols.
#'
#' @param gt_tbl A gt table object.
#' @param columns Column names to apply the math font to.
#' @param math_font Character string. Name of the math font. Default is "STIX Two Math".
#'
#' @return Modified gt table object.
#'
#' @noRd
apply_math_font <- function(
    gt_tbl,
    columns,
    math_font = "STIX Two Math"
) {
  gt_tbl |>
    gt::tab_style(
      style = gt::cell_text(font = math_font),
      locations = gt::cells_body(columns = columns))
}

#' Extract Priors from RoBMA Fit
#'
#' @description
#' Extracts and formats prior specifications from a RoBMA fit object
#' for display in tables.
#'
#' @param robma_fit A RoBMA fit object.
#'
#' @return List with mu_prior_unicode and tau_prior_unicode.
#'
#' @noRd
extract_priors_from_robma_fit <- function(robma_fit) {

  prior_to_label <- function(p) {
    fam <- tolower(p$family %||% "")
    if (fam %in% c("normal", "gaussian")) {
      sprintf("N(%.2f, %.2f)", p$mean %||% 0, p$sd %||% NA_real_)
    } else if (fam %in% c("student_t", "student-t", "t", "student")) {
      sprintf("t(df=%.1f, loc=%.2f, scale=%.2f)", p$df, p$location, p$scale)
    } else {
      fam
    }
  }

  if (!is.null(robma_fit$meta$priors)) {

    eff <- robma_fit$meta$priors$effect %||% list()
    het <- robma_fit$meta$priors$heterogeneity %||% list()

    mu_labels  <- unique(purrr::map_chr(eff, prior_to_label))
    tau_labels <- unique(purrr::map_chr(het, prior_to_label))

    mu_prior_unicode <- paste(mu_labels, collapse = " / ")
    tau_prior_unicode <- paste(tau_labels, collapse = " / ")

    return(list(
      mu_prior_unicode = mu_prior_unicode,
      tau_prior_unicode = tau_prior_unicode
    ))
  }

  cli::cli_abort(c(
    "Cannot extract priors from this RoBMA object.",
    "i" = "Expected robma_fit$meta$priors$effect and $heterogeneity."
  ))
}

