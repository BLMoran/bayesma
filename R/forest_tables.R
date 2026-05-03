#' Internal function to create table for left side of forest plot
#'
#' @noRd
forest_table_left_fn <- function(forest.data.summary,
                                 subgroup = FALSE,
                                 label_control = "Control",
                                 label_intervention = "Intervention",
                                 estimand = "OR",
                                 font = NULL) {

  is_continuous <- estimand %in% c("MD", "SMD") ||
    (is_marginal_estimand(estimand) && "int_mean_sd" %in% names(forest.data.summary))

  # Define column labels based on effect type
  if (is_continuous) {
    control_label <- paste(label_control, "\nMean (SD)")
    int_label <- paste(label_intervention, "\nMean (SD)")
  } else {
    control_label <- paste(label_control, "\n(Events/Total)")
    int_label <- paste(label_intervention, "\n(Events/Total)")
  }

  # Choose correct columns for data
  if (is_continuous) {
    selected_cols <- c("Author", "Year", "N_int", "int_mean_sd", "N_ctrl", "ctrl_mean_sd")
  } else {
    selected_cols <- c("Author", "Year", "int_outcome_frac", "control_outcome_frac")
  }
  if (isTRUE(subgroup)) {
    selected_cols <- c("Author", "Year", "Subgroup", selected_cols[!(selected_cols %in% c("Author", "Year"))])
  }

  df <- forest.data.summary |> dplyr::select(dplyr::any_of(selected_cols))

  df <- df |> dplyr::mutate(
    Author = dplyr::if_else(
      is.na(Year),
      Author,
      paste0(Author, " (", Year, ")")
    )
  )

  # Replace NA values with blanks for the 'Pooled Effect' and 'Prediction' rows
  if (is_continuous) {
    df <- df |> dplyr::mutate(
      int_mean_sd = dplyr::if_else(Author %in% c("Pooled Effect", "Prediction"), NA_character_, int_mean_sd),
      ctrl_mean_sd = dplyr::if_else(Author %in% c("Pooled Effect", "Prediction"), NA_character_, ctrl_mean_sd)
    )
  }

  # Base gt table
  if (isFALSE(subgroup)) {
    forest.table.left <- df |>
      gt::gt() |>
      gt::cols_label(Author = "Study") |>
      gt::cols_align(align = "left") |>
      gt::tab_style(
        style = gt::cell_text(style = "italic", weight = "bold"),
        locations = gt::cells_body(rows = Author %in% c("Pooled Effect", "Prediction")))
  } else {
    forest.table.left <- df |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_options(row_group.font.weight = "bold") |>
      gt::cols_label(Author = gt::md("Subgroup/ \nStudy")) |>
      gt::cols_align(align = "left") |>
      gt::tab_style(
        gt::cell_text(color = "white"),
        locations = gt::cells_row_groups(groups = "Overall")) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold", style = "italic"),
        locations = gt::cells_body(rows = Author == "Overall Effect")) |>
      gt::tab_style(
        style = gt::cell_text(style = "italic",  weight = "bold", color = "grey60"),
        locations = gt::cells_body(rows = Author == "Pooled Effect")) |>
      # Subgroup-level Prediction: same style as Pooled Effect (grey, italic, bold)
      gt::tab_style(
        style = gt::cell_text(style = "italic", weight = "bold", color = "grey60"),
        locations = gt::cells_body(rows = Author == "Prediction" & Subgroup != "Overall")) |>
      # Overall Prediction: same style as Overall Effect (bold, italic, no grey)
      gt::tab_style(
        style = gt::cell_text(style = "italic", weight = "bold"),
        locations = gt::cells_body(rows = Author == "Prediction" & Subgroup == "Overall")) |>
      gt::tab_style(
        style = gt::cell_text(indent = gt::px(10)),
        locations = gt::cells_body(rows = !Author %in% c("Overall Effect") &
                                     !(Author == "Prediction" & Subgroup == "Overall"),
                                   columns = Author))
  }

  # Column label mapping
  col_labels <- if (is_continuous) {
    c(N_int = "N", int_mean_sd = gt::md(int_label),
      N_ctrl = "N", ctrl_mean_sd = gt::md(control_label))
  } else {
    c(control_outcome_frac = gt::md(control_label),
      int_outcome_frac = gt::md(int_label))
  }

  # Final styling
  forest.table.left <- forest.table.left |>
    gt::cols_label(!!!col_labels) |>
    gt::tab_options(
      column_labels.font.weight = "bold",
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font) |>
    gt::cols_hide(columns = Year) |>
    gt::opt_table_lines(extent = "none") |>
    gt::tab_style(
      style = gt::cell_fill(color = "grey95"),
      locations = gt::cells_body(rows = Author %in% c("Pooled Effect", "Overall Effect", "Prediction"))
    )

  if (is_continuous) {
    forest.table.left <- forest.table.left |>
      gt::sub_missing(columns = c(int_mean_sd, ctrl_mean_sd), missing_text = "")
  }

  return(forest.table.left)
}


#' Internal function to create table on right of forest plot (including RoB)
#'
#' @noRd
forest_table_right_fn <- function(df,
                                  subgroup = FALSE,
                                  estimand = "OR",
                                  add_rob = FALSE,
                                  has_re = TRUE,
                                  incl_shrinkage = TRUE,
                                  font = NULL) {

  # Determine which effect columns to include
  if (isTRUE(has_re) && isTRUE(incl_shrinkage)) {
    effect_cols <- c("shrinkage_display", "unweighted_effect")
  } else {
    effect_cols <- c("unweighted_effect")
  }

  # Identify RoB columns present in data
  rob_cols <- names(df)[grepl("^D\\d+$|^Overall$", names(df))]

  if (isFALSE(subgroup)) {
    select_cols <- c("Author", effect_cols, rob_cols)

    forest.table.right <- df |>
      dplyr::select(dplyr::all_of(select_cols)) |>
      gt::gt() |>
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold")),
        locations = gt::cells_body(
          rows = Author %in% c("Pooled Effect", "Prediction")))

  } else if (isTRUE(subgroup)) {
    select_cols <- c("Author", "Subgroup", effect_cols, rob_cols)

    forest.table.right <- df |>
      dplyr::select(dplyr::all_of(select_cols)) |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_style(
        gt::cell_text(color = "white"),
        locations = gt::cells_row_groups(groups = gt::everything())) |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold", style = "italic"),
        locations = gt::cells_body(rows = Author == "Overall Effect")) |>
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold",
            color = "grey60")),
        locations = gt::cells_body(
          rows = Author == "Pooled Effect")) |>
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold",
            color = "grey60")),
        locations = gt::cells_body(
          rows = Author == "Prediction" & Subgroup != "Overall")) |>
      gt::tab_style(
        style = list(
          gt::cell_text(
            style = "italic",
            weight = "bold")),
        locations = gt::cells_body(
          rows = Author == "Prediction" & Subgroup == "Overall"))
  }

  forest.table.right <- forest.table.right |>
    gt::tab_options(column_labels.font.weight = "bold") |>
    gt::cols_align(align = "right") |>
    gt::cols_label(.list = create_rob_labels(df, estimand, has_re = has_re, incl_shrinkage = incl_shrinkage)) |>
    gt::tab_style(
      style = gt::cell_text(align = "center"),
      locations = gt::cells_body(
        columns = dplyr::any_of("shrinkage_display"),
        rows = shrinkage_display == "—")) |>
    gt::tab_options(
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font) |>
    gt::opt_table_lines(extent = "none") |>
    gt::tab_style(
      style = gt::cell_fill(color = "grey95"),
      locations = gt::cells_body(
        rows = Author %in% c("Pooled Effect", "Overall Effect", "Prediction"),
        columns = dplyr::all_of(effect_cols)))

  if (add_rob == FALSE) {
    # Hide RoB columns, but keep Author as a 1px invisible column
    # so patchwork always has >= 2 columns (avoids zero-length unit error)
    forest.table.right <- forest.table.right |>
      gt::cols_hide(columns = dplyr::all_of(rob_cols)) |>
      gt::cols_label(Author = "") |>
      gt::cols_width(Author ~ gt::px(1)) |>
      gt::tab_style(
        style = gt::cell_text(color = "white", size = gt::px(1)),
        locations = gt::cells_body(columns = Author)) |>
      gt::tab_style(
        style = gt::cell_borders(sides = "all", color = "white", weight = gt::px(0)),
        locations = gt::cells_body(columns = Author))

  } else if (add_rob == TRUE) {
    forest.table.right <- forest.table.right |>
      gt::opt_table_lines(extent = "none") |>
      gt::sub_missing(columns = dplyr::all_of(rob_cols), missing_text = "") |>
      gt::cols_align(
        align = "center",
        columns = dplyr::all_of(rob_cols)) |>
      gt::tab_style(
        style = "padding-left:0px; padding-right:0px;",
        locations = gt::cells_body(columns = get_rob_text_columns(df))) |>
      gt::tab_style(
        style = "padding-left:0px; padding-right:0px;",
        locations = gt::cells_column_labels(columns = get_rob_text_columns(df))) |>
      gt::text_case_match(
        "High" ~ fontawesome::fa(name = "circle-plus", fill = "#e32400", height = "1.1em"),
        "Low" ~ fontawesome::fa(name = "circle-minus", fill = "#77bb41", height = "1.1em"),
        "Some concerns" ~ fontawesome::fa(name = "circle-question", fill = "#f5ec00", height = "1.1em")) |>
      gt::tab_style(
        style = "padding-top:2px; padding-bottom:2px;",
        locations = gt::cells_body(columns = dplyr::all_of(rob_cols))) |>
      gt::tab_style(
        style = "padding-right:2px;",
        locations = gt::cells_body(
          columns = unweighted_effect,
          rows = Author %in% c("Pooled Effect", "Prediction"))) |>
      gt::cols_label(Author = "") |>
      gt::cols_width(Author ~ gt::px(1)) |>
      gt::tab_style(
        style = gt::cell_text(color = "white", size = gt::px(1)),
        locations = gt::cells_body(columns = Author)) |>
      gt::tab_style(
        style = gt::cell_borders(sides = "all", color = "white", weight = gt::px(0)),
        locations = gt::cells_body(columns = Author))
  }

  return(forest.table.right)
}
