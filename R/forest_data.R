#' Remap Author names in draws to match disambiguated data
#'
#' When `make_authors_unique()` has renamed duplicate Author names in the data
#' (e.g. "Wildes" -> "Wildes_a", "Wildes_b"), the draws still use the original
#' names from `model$meta$study_labels`. This function positionally remaps the
#' draws' Author column to match the disambiguated data.
#'
#' @param draws A tibble of posterior draws with an `Author` column.
#' @param model A bayesma model object with `meta$study_labels`.
#' @param data The data frame after `make_authors_unique()` has been applied.
#' @return The draws tibble with Author names remapped to match the data.
#' @noRd
remap_draws_authors <- function(draws, model, data) {
  original_labels <- model$meta$study_labels
  disambiguated_labels <- data$Author[seq_along(original_labels)]

  # Only remap if disambiguation changed any names
  if (identical(as.character(original_labels), as.character(disambiguated_labels))) {
    return(draws)
  }

  # Separate study draws from special rows (Pooled Effect, Prediction, etc.)
  special_authors <- c("Pooled Effect", "Prediction", "Overall Effect", "No Pooled Effect")
  study_draws <- draws |> dplyr::filter(!Author %in% special_authors)
  other_draws <- draws |> dplyr::filter(Author %in% special_authors)

  # Each study has the same number of MCMC draws, and extract_forest_draws()
  # produces them in positional order (study 1 first, study 2, ..., study S).
  n_studies <- length(original_labels)
  n_total   <- nrow(study_draws)
  n_per     <- n_total / n_studies

  # Positionally assign disambiguated names
  study_draws$Author <- rep(disambiguated_labels, each = n_per)

  dplyr::bind_rows(study_draws, other_draws)
}


#' Internal function to extract draws from the posterior
#'
#' Dispatches to brms-specific or bayesma-specific extraction depending
#' on the class of `model`.
#'
#' @noRd
forest_data_fn <- function(data,
                           model,
                           subgroup = FALSE,
                           sort_studies_by = "author",
                           subgroup_order = NULL,
                           add_pred = FALSE,
                           add_pred_subgroup = FALSE,
                           has_re = TRUE) {

  is_bayesma <- inherits(model, "bayesma")

  if (is_bayesma) {
    forest_data_fn_bayesma(
      data = data,
      model = model,
      subgroup = subgroup,
      sort_studies_by = sort_studies_by,
      subgroup_order = subgroup_order,
      add_pred = add_pred,
      add_pred_subgroup = add_pred_subgroup,
      has_re = has_re
    )
  } else {
    forest_data_fn_brms(
      data = data,
      model = model,
      subgroup = subgroup,
      sort_studies_by = sort_studies_by,
      subgroup_order = subgroup_order,
      add_pred = add_pred,
      add_pred_subgroup = add_pred_subgroup,
      has_re = has_re
    )
  }
}


# ============================================================================
# bayesma pathway
# ============================================================================

#' @noRd
forest_data_fn_bayesma <- function(data,
                                   model,
                                   subgroup = FALSE,
                                   sort_studies_by = "author",
                                   subgroup_order = NULL,
                                   add_pred = FALSE,
                                   add_pred_subgroup = FALSE,
                                   has_re = TRUE) {

  if (isFALSE(subgroup) && "Subgroup" %in% names(data)) {
    data <- data |> dplyr::select(-Subgroup)
  }

  if (isFALSE(subgroup)) {
    # ---- No subgroup: extract draws from bayesma object directly ----
    effect.draws <- extract_forest_draws(model)

    # Remap Author names to match disambiguated data (e.g. "Wildes" -> "Wildes_a")
    effect.draws <- remap_draws_authors(effect.draws, model, data)

    # If prediction not requested, drop Prediction rows
    if (isFALSE(add_pred)) {
      effect.draws <- effect.draws |>
        dplyr::filter(Author != "Prediction")
    }

    # Replace dots with spaces in Author names (matching brms convention)
    effect.draws <- effect.draws |>
      dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |>
      dplyr::ungroup() |>
      dplyr::left_join(
        dplyr::select(data, Author, Author_original, Year, yi, vi),
        by = dplyr::join_by(Author)
      ) |>
      sort_studies_fn(sort_studies_by)

  } else {
    # ---- Subgroup: refit bayesma on each subgroup ----
    subgroup_df <- data |>
      tidyr::nest(.by = Subgroup) |>
      dplyr::mutate(
        subgroup_model = purrr::map(data, ~ refit_bayesma(model, .x)),
        study_count = purrr::map_int(data, nrow)
      )

    # Overall effect draws from the main (full-data) model
    overall_draws <- extract_forest_draws(model)
    overall.effect.draws <- overall_draws |>
      dplyr::filter(Author == "Pooled Effect") |>
      dplyr::mutate(
        Author = "Overall Effect",
        Author_original = NA_character_,
        Subgroup = "Overall",
        Year = NA_character_,
        yi = NA_real_,
        vi = NA_real_
      )

    # Extract draws for each subgroup model
    study.effect.draws <- subgroup_df |>
      dplyr::mutate(
        effect_draws = purrr::pmap(
          list(subgroup_model, data, study_count),
          function(sub_model, sub_data, n_studies) {
            draws <- extract_forest_draws(sub_model)

            # Remap Author names to match disambiguated subgroup data
            draws <- remap_draws_authors(draws, sub_model, sub_data)

            # Rename pooled label for small subgroups
            pooled_label <- if (n_studies == 1) "No Pooled Effect" else "Pooled Effect"
            draws <- draws |>
              dplyr::mutate(
                Author = dplyr::if_else(
                  Author == "Pooled Effect", pooled_label, Author
                )
              )

            # Remove prediction draws if not requested for subgroups
            if (isFALSE(add_pred) || isFALSE(add_pred_subgroup)) {
              draws <- draws |> dplyr::filter(Author != "Prediction")
            }

            draws |>
              dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |>
              dplyr::ungroup() |>
              dplyr::left_join(
                dplyr::select(sub_data, Author, Author_original, Year, yi, vi),
                by = dplyr::join_by(Author)
              ) |>
              sort_studies_fn(sort_studies_by)
          }
        )
      ) |>
      tidyr::unnest(effect_draws) |>
      dplyr::select(-data, -subgroup_model, -study_count)

    effect.draws <- dplyr::bind_rows(study.effect.draws, overall.effect.draws)

    # Add overall prediction draws when add_pred = TRUE
    if (isTRUE(add_pred)) {
      pred_from_overall <- overall_draws |>
        dplyr::filter(Author == "Prediction") |>
        dplyr::mutate(
          Author_original = NA_character_,
          Subgroup = "Overall",
          Year = NA_character_,
          yi = NA_real_,
          vi = NA_real_
        )
      if (nrow(pred_from_overall) > 0) {
        effect.draws <- dplyr::bind_rows(effect.draws, pred_from_overall)
      }
    }

    # Custom group order for Subgroup column
    if (!is.null(subgroup_order)) {
      effect.draws <- effect.draws |>
        dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
        dplyr::arrange(Subgroup) |>
        dplyr::mutate(Subgroup = dplyr::case_when(
          is.na(Subgroup) & Author == "Overall Effect" ~ "Overall",
          .default = Subgroup
        ))
    }
  }

  return(effect.draws)
}


# ============================================================================
# brms pathway (original logic, extracted to its own function)
# ============================================================================

#' @noRd
forest_data_fn_brms <- function(data,
                                model,
                                subgroup = FALSE,
                                sort_studies_by = "author",
                                subgroup_order = NULL,
                                add_pred = FALSE,
                                add_pred_subgroup = FALSE,
                                has_re = TRUE) {

  if (subgroup == FALSE && "Subgroup" %in% names(data)) {
    data <- data |> dplyr::select(-Subgroup)
  }
  if (subgroup == FALSE) {
    if (isTRUE(has_re)) {
      # Random effects model — extract study-level and pooled draws
      study.draws <- tidybayes::spread_draws(model, r_Author[Author, ], b_Intercept) |>
        dplyr::mutate(b_Intercept = r_Author + b_Intercept)
      pooled.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
        dplyr::mutate(Author = "Pooled Effect")
    } else {
      # Common effect model — no random effects, study draws use observed yi
      pooled.draws <- tidybayes::spread_draws(model, b_Intercept) |>
        dplyr::mutate(
          Author = "Pooled Effect",
          sd_Author__Intercept = NA_real_
        )

      # For common effect, study-level draws are just the pooled draws
      # (no shrinkage). We still need them in the same structure for plotting
      # the observed likelihood-based densities via yi/vi.
      study_authors <- unique(data$Author)
      study.draws <- purrr::map(study_authors, function(auth) {
        tidybayes::spread_draws(model, b_Intercept) |>
          dplyr::mutate(
            Author = auth,
            r_Author = 0,
            sd_Author__Intercept = NA_real_
          )
      }) |> purrr::list_rbind()
    }

    effect.draws <- dplyr::bind_rows(study.draws, pooled.draws)

    # Generate prediction draws if requested (only meaningful for RE models)
    if (isTRUE(add_pred) && isTRUE(has_re)) {
      nd <- data.frame(Author = "new", sei = 0)
      pred_samples <- brms::posterior_predict(
        object = model,
        newdata = nd,
        re_formula = NULL,
        allow_new_levels = TRUE,
        sample_new_levels = "gaussian"
      )

      pred.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
        dplyr::mutate(
          Author = "Prediction",
          b_Intercept = as.vector(pred_samples)
        )
      effect.draws <- dplyr::bind_rows(effect.draws, pred.draws)
    }

    effect.draws <- effect.draws |>
      dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |>
      dplyr::ungroup() |>
      dplyr::left_join(dplyr::select(data, Author, Author_original, Year, yi, vi), by = dplyr::join_by(Author)) |>
      sort_studies_fn(sort_studies_by)
  } else {
    # With subgroup
    subgroup_df <- data |>
      tidyr::nest(.by = Subgroup) |>
      dplyr::mutate(
        subgroup_model = purrr::map(data, ~ stats::update(model, newdata = .x)),
        study_count = purrr::map_int(data, nrow)
      )

    if (isTRUE(has_re)) {
      overall.effect.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
        dplyr::mutate(
          Author = "Overall Effect",
          Author_original = NA_character_,
          Subgroup = "Overall",
          r_Author = 0,
          Year = NA_character_,
          yi = NA_real_,
          vi = NA_real_
        )
    } else {
      overall.effect.draws <- tidybayes::spread_draws(model, b_Intercept) |>
        dplyr::mutate(
          Author = "Overall Effect",
          Author_original = NA_character_,
          Subgroup = "Overall",
          r_Author = 0,
          sd_Author__Intercept = NA_real_,
          Year = NA_character_,
          yi = NA_real_,
          vi = NA_real_
        )
    }

    study.effect.draws <- subgroup_df |>
      dplyr::mutate(
        effect_draws = purrr::pmap(list(subgroup_model, data, study_count), function(sub_mod, sub_data, n_studies) {
          if (isTRUE(has_re)) {
            study <- tidybayes::spread_draws(sub_mod, r_Author[Author, ], b_Intercept) |>
              dplyr::mutate(b_Intercept = r_Author + b_Intercept)
          } else {
            study_authors <- unique(sub_data$Author)
            study <- purrr::map(study_authors, function(auth) {
              tidybayes::spread_draws(sub_mod, b_Intercept) |>
                dplyr::mutate(Author = auth, r_Author = 0, sd_Author__Intercept = NA_real_)
            }) |> purrr::list_rbind()
          }

          pooled_label <- if (n_studies == 1) "No Pooled Effect" else "Pooled Effect"

          if (isTRUE(has_re)) {
            pooled <- tidybayes::spread_draws(sub_mod, b_Intercept, sd_Author__Intercept) |>
              dplyr::mutate(Author = pooled_label)
          } else {
            pooled <- tidybayes::spread_draws(sub_mod, b_Intercept) |>
              dplyr::mutate(Author = pooled_label, sd_Author__Intercept = NA_real_)
          }

          combined <- dplyr::bind_rows(study, pooled)

          if (isTRUE(add_pred) && isTRUE(add_pred_subgroup) && n_studies > 1 && isTRUE(has_re)) {
            nd <- data.frame(Author = "new", sei = 0)
            pred_samples <- brms::posterior_predict(
              object = sub_mod,
              newdata = nd,
              re_formula = NULL,
              allow_new_levels = TRUE,
              sample_new_levels = "gaussian"
            )

            pred <- tidybayes::spread_draws(sub_mod, b_Intercept, sd_Author__Intercept) |>
              dplyr::mutate(
                Author = "Prediction",
                b_Intercept = as.vector(pred_samples)
              )
            combined <- dplyr::bind_rows(combined, pred)
          }

          combined <- combined |>
            dplyr::mutate(Author = stringr::str_replace_all(Author, "\\.", " ")) |>
            dplyr::ungroup() |>
            dplyr::left_join(dplyr::select(sub_data, Author, Author_original, Year, yi, vi), by = dplyr::join_by(Author)) |>
            sort_studies_fn(sort_studies_by)
          return(combined)
        })
      ) |>
      tidyr::unnest(effect_draws) |>
      dplyr::select(-data, -subgroup_model, -study_count)

    effect.draws <- dplyr::bind_rows(study.effect.draws, overall.effect.draws)

    if (isTRUE(add_pred) && isTRUE(has_re)) {
      nd <- data.frame(Author = "new", sei = 0)
      pred_samples <- brms::posterior_predict(
        object = model,
        newdata = nd,
        re_formula = NULL,
        allow_new_levels = TRUE,
        sample_new_levels = "gaussian"
      )

      overall.pred.draws <- tidybayes::spread_draws(model, b_Intercept, sd_Author__Intercept) |>
        dplyr::mutate(
          Author = "Prediction",
          Author_original = NA_character_,
          Subgroup = "Overall",
          r_Author = 0,
          Year = NA_character_,
          yi = NA_real_,
          vi = NA_real_,
          b_Intercept = as.vector(pred_samples)
        )
      effect.draws <- dplyr::bind_rows(effect.draws, overall.pred.draws)
    }

    if (!is.null(subgroup_order)) {
      effect.draws <- effect.draws |>
        dplyr::mutate(Subgroup = factor(Subgroup, levels = subgroup_order)) |>
        dplyr::arrange(Subgroup) |>
        dplyr::mutate(Subgroup = dplyr::case_when(
          is.na(Subgroup) & Author == "Overall Effect" ~ "Overall",
          .default = Subgroup
        ))
    }
  }
  return(effect.draws)
}


# ============================================================================
# forest.data.summary_fn — unchanged, works with both pathways
# ============================================================================

#' Internal function to summarise data for forest plot
#'
#' @noRd
forest.data.summary_fn <- function(spread_df,
                                   data,
                                   measure,
                                   sort_studies_by = "author",
                                   subgroup = FALSE,
                                   add_pred = FALSE,
                                   add_pred_subgroup = FALSE,
                                   has_re = TRUE) {
  # Get effect size properties
  props <- get_measure_properties(measure)

  if (isFALSE(subgroup)){
    # Study Summaries
    forest.data <- spread_df |>
      dplyr::group_by(Author) |>
      tidybayes::median_qi(b_Intercept)

    # Tau Summary — only meaningful for random effects models
    if (isTRUE(has_re)) {
      tau.summary <- spread_df |>
        dplyr::group_by(Author) |>
        tidybayes::median_qi(sd_Author__Intercept)
    } else {
      # Create a placeholder tau summary with NAs for common effect
      tau.summary <- forest.data |>
        dplyr::select(Author) |>
        dplyr::mutate(
          sd_Author__Intercept = NA_real_,
          .lower = NA_real_,
          .upper = NA_real_,
          .width = 0.95,
          .point = "median",
          .interval = "qi"
        )
    }

  } else if (subgroup == TRUE) {
    forest.data <- spread_df |>
      dplyr::group_by(Subgroup, Author) |>
      tidybayes::median_qi(b_Intercept)

    if (isTRUE(has_re)) {
      tau.summary <- spread_df |>
        dplyr::group_by(Subgroup, Author) |>
        tidybayes::median_qi(sd_Author__Intercept)
    } else {
      tau.summary <- forest.data |>
        dplyr::select(Subgroup, Author) |>
        dplyr::mutate(
          sd_Author__Intercept = NA_real_,
          .lower = NA_real_,
          .upper = NA_real_,
          .width = 0.95,
          .point = "median",
          .interval = "qi"
        )
    }
  }

  # Select only needed join columns — include Author_original so it survives
  join_vars <- unique(c("Author", "Author_original", "Year", props$data_cols))
  forest.data.summary <- forest.data |>
    dplyr::left_join(data |> dplyr::select(dplyr::any_of(join_vars), yi, vi, D1:Overall), by = dplyr::join_by(Author)) |>
    dplyr::left_join(tau.summary, by = dplyr::join_by(Author), suffix = c("", "_sd")) |>
    sort_studies_fn(sort_studies_by)

  # Add formatted effect estimates
  forest.data.summary <- forest.data.summary |>
    dplyr::mutate(
      weighted_effect = if (measure %in% c("MD", "SMD")) {
        paste0(sprintf("%.2f", b_Intercept), " [", sprintf("%.2f", .lower), ", ", sprintf("%.2f", .upper), "]")
      } else {
        paste0(sprintf("%.2f", exp(b_Intercept)), " [", sprintf("%.2f", exp(.lower)), ", ", sprintf("%.2f", exp(.upper)), "]")
      },
      unweighted_effect = if (measure %in% c("MD", "SMD")) {
        paste0(sprintf("%.2f", yi), " [", sprintf("%.2f", yi - 1.96 * sqrt(vi)), ", ", sprintf("%.2f", yi + 1.96 * sqrt(vi)), "]")
      } else {
        paste0(sprintf("%.2f", exp(yi)), " [", sprintf("%.2f", exp(yi - 1.96 * sqrt(vi))), ", ", sprintf("%.2f", exp(yi + 1.96 * sqrt(vi))), "]")
      })

  # Handle the Pooled/Overall/Prediction rows where yi/vi are NA
  forest.data.summary <- forest.data.summary |>
    dplyr::mutate(
      unweighted_effect = dplyr::case_when(
        # Prediction row: always blank
        as.character(Author) == "Prediction" ~ "",
        # Rows with valid observed effects: keep as-is
        unweighted_effect != "NA [NA, NA]" ~ unweighted_effect,
        # Common-effect pooled/overall rows: show posterior estimate here
        # (since weighted_effect column will be hidden)
        isFALSE(has_re) & Author %in% c("Pooled Effect", "Overall Effect") ~
          weighted_effect,
        # Random-effect pooled/overall rows: show tau
        !is.na(sd_Author__Intercept) ~
          paste0("\u03c4 = ", sprintf("%.2f", sd_Author__Intercept),
                 " [", sprintf("%.2f", .lower_sd), ", ",
                 sprintf("%.2f", .upper_sd), "]"),
        # Fallback
        TRUE ~ ""
      )
    )

  # Add study group summary columns depending on measure type
  if (measure %in% c("MD", "SMD")) {
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        N_int = dplyr::case_when(
          Author == "Pooled Effect" ~ as.character(sum(N_Intervention, na.rm = TRUE)),
          TRUE ~ as.character(N_Intervention)),
        int_mean_sd = dplyr::case_when(
          Author %in% c("Pooled Effect", "Prediction") ~ NA_character_,
          TRUE ~ paste0(sprintf("%.2f", Mean_Intervention), " (", sprintf("%.2f", SD_Intervention), ")")),
        N_ctrl = dplyr::case_when(
          Author == "Pooled Effect" ~ as.character(sum(N_Control, na.rm = TRUE)),
          TRUE ~ as.character(N_Control)),
        ctrl_mean_sd = dplyr::case_when(
          Author %in% c("Pooled Effect", "Prediction") ~ NA_character_,
          TRUE ~ paste0(sprintf("%.2f", Mean_Control), " (", sprintf("%.2f", SD_Control), ")"))
      )
  } else {
    forest.data.summary <- forest.data.summary |>
      dplyr::mutate(
        control_outcome_frac = dplyr::case_when(
          Author == "Pooled Effect" ~ paste0(
            sum(Event_Control[Author != "Pooled Effect"], na.rm = TRUE), "/",
            sum(N_Control[Author != "Pooled Effect"], na.rm = TRUE)),
          Author == "Prediction" ~ "",
          TRUE ~ paste0(Event_Control, "/", N_Control)),
        int_outcome_frac = dplyr::case_when(
          Author == "Pooled Effect" ~ paste0(
            sum(Event_Intervention[Author != "Pooled Effect"], na.rm = TRUE), "/",
            sum(N_Intervention[Author != "Pooled Effect"], na.rm = TRUE)),
          Author == "Prediction" ~ "",
          TRUE ~ paste0(Event_Intervention, "/", N_Intervention)))
  }

  # Remove rows where Author == "No Pooled Effect"
  forest.data.summary <- forest.data.summary |>
    dplyr::filter(Author != "No Pooled Effect")

  return(forest.data.summary)
}
