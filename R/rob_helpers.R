

#' Create Risk of Bias Table
#'
#' Internal function to create the gt table for risk of bias visualization.
#'
#' @param df Data frame with risk of bias data
#' @param subgroup Logical indicating whether this is a subgroup analysis
#' @param font Font family to use
#'
#' @return A gt table object
#' @noRd
rob_table_fn <- function(df,
                         subgroup = FALSE,
                         font = NULL) {

  df <- df |> dplyr::mutate(
    Author = dplyr::if_else(
      is.na(Year),
      Author,
      paste0(Author, " (", Year, ")")
    )
  )

  # Base gt table
  if (isFALSE(subgroup)) {
    rob.table <- df |>
      gt::gt() |>
      gt::cols_label(Author = "Study")

  } else {
    rob.table <- df |>
      gt::gt(groupname_col = "Subgroup") |>
      gt::tab_options(row_group.font.weight = "bold") |>
      gt::tab_style(
        style = gt::cell_fill(color = "grey95"),
        locations = gt::cells_row_groups()) |>
      gt::cols_label(Author = gt::md("Subgroup/ \nStudy")) |>
      gt::cols_align(align = "left") |>
      gt::tab_style(
        style = gt::cell_text(indent = gt::px(10)),
        locations = gt::cells_body(columns = Author))
  }

  # Final styling
  rob.table <- rob.table |>
    gt::cols_align(align = "left",
                   columns = Author) |>
    gt::cols_align(align = "center",
                   columns = c(D1:Overall)) |>
    gt::tab_options(
      column_labels.font.weight = "bold",
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font) |>
    gt::cols_hide(columns = Year) |>
    gt::opt_table_lines(extent = "none") |>
    gt::text_case_match(
      "High" ~ fontawesome::fa(name = "circle-plus", fill = "#e32400", height = "1.1em"),
      "Low" ~ fontawesome::fa(name = "circle-minus", fill = "#77bb41", height = "1.1em"),
      "Some concerns" ~ fontawesome::fa(name = "circle-question", fill = "#f5ec00", height = "1.1em"))

  return(rob.table)
}


#' Create Risk of Bias Legend
#'
#' Internal function to create a legend for risk of bias symbols.
#'
#' @param rob_tool Character string specifying the risk of bias tool
#' @param font Font family to use
#'
#' @return A gt table object containing the legend
#' @noRd
create_rob_legend <- function(rob_tool, font = NULL) {
  # This is a placeholder - you'll need to implement the actual legend creation
  # based on your specific requirements for each tool type

  legend_data <- data.frame(
    Symbol = c(
      fontawesome::fa(name = "circle-minus", fill = "#77bb41", height = "1.1em"),
      fontawesome::fa(name = "circle-question", fill = "#f5ec00", height = "1.1em"),
      fontawesome::fa(name = "circle-plus", fill = "#e32400", height = "1.1em")
    ),
    Meaning = c("Low risk", "Some concerns", "High risk")
  )

  legend_table <- legend_data |>
    gt::gt() |>
    gt::tab_options(
      column_labels.hidden = TRUE,
      table.font.names = font,
      table.border.top.color = "white",
      table.border.bottom.color = "white"
    ) |>
    gt::opt_table_lines(extent = "none")

  return(legend_table)
}

#' Internal function to create risk of bias legend table
#'
#' @noRd
rob_legend_fn <- function(rob_tool = c("rob2", "robins_i", "quadas2", "robins_e", "nos"),
                          font = NULL) {

  # Get domains for the specified tool
  domains <- get_rob_domains(rob_tool, full_text = TRUE)

  # Remove overall domain
  domains <- domains[names(domains) != "Overall"]

  # Create the data for the table
  legend_data <- data.frame(
    Code = paste0("(", names(domains), ")"),
    Description = unname(domains),
    stringsAsFactors = FALSE)

  # Add empty row for spacing
  legend_data <- rbind(
    legend_data,
    data.frame(Code = "", Description = "", stringsAsFactors = FALSE))

  # Add risk level rows with placeholder text
  risk_rows <- data.frame(
    Code = c(" ", "Low", "Some concerns", "High"),
    Description = c(" ", "Low risk", "Some concerns", "High risk"),
    stringsAsFactors = FALSE)

  legend_data <- rbind(legend_data, risk_rows)

  # Create the gt table
  rob_legend_table <- legend_data |>
    gt::gt() |>
    # Remove column headers
    gt::cols_label(
      Code = "",
      Description = "") |>
    # Add title
    gt::tab_header(
      title = gt::html("<b>Risk of bias legend<b>")) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", size = gt::px(16))),
      locations = gt::cells_title(groups = "title")) |>
    # Style the domain codes (first n rows)
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", size = gt::px(12))),
      locations = gt::cells_body(
        columns = Code,
        rows = 1:length(domains))) |>
    # Transform the risk level text to icons
    gt::text_transform(
      locations = gt::cells_body(
        columns = Code,
        rows = (length(domains) + 2):nrow(legend_data)),
      fn = function(x) {
        dplyr::case_when(
          x == "High" ~ as.character(fontawesome::fa(name = "circle-plus", fill = "#e32400", height = "1.5em")),
          x == "Low" ~ as.character(fontawesome::fa(name = "circle-minus", fill = "#77bb41", height = "1.5em")),
          x == "Some concerns" ~ as.character(fontawesome::fa(name = "circle-question", fill = "#f5ec00", height = "1.5em")),
          TRUE ~ x)
      }) |>
    # Style the risk level descriptions to be bold and italic
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", style = "italic")),
      locations = gt::cells_body(
        columns = Description,
        rows = (length(domains) + 2):nrow(legend_data))) |>
    # Reduce space in between RoB and legend
    gt::tab_style(
      style = "padding-left:2px;",
      locations = gt::cells_body(
        columns = Code,
        rows = gt::everything())) |>
    # Set column widths
    gt::cols_width(
      Code ~ gt::px(40),
      Description ~ gt::px(300)) |>
    # Table options
    gt::tab_options(
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      heading.align = "left",
      heading.padding = gt::px(6),
      table.font.size = gt::px(12),
      table.font.names = font) |>
    gt::opt_table_lines(extent = "none")

  return(rob_legend_table)
}


#' Internal function to get RoB domains for RoB legend
#'
#' @noRd
get_rob_domains <- function(rob_tool = c("rob2", "robins_i", "quadas2", "robins_e"),
                            full_text = TRUE) {
  rob_tool <- rlang::arg_match(rob_tool)

  domains <- switch(rob_tool,
                    rob2 = c(
                      D1 = "Bias arising from the randomization process",
                      D2 = "Bias due to deviations from intended interventions",
                      D3 = "Bias due to missing outcome data",
                      D4 = "Bias in measurement of the outcome",
                      D5 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    ),
                    robins_i = c(
                      D1 = "Bias due to confounding",
                      D2 = "Bias in selection of participants into the study",
                      D3 = "Bias in classification of interventions",
                      D4 = "Bias due to deviations from intended interventions",
                      D5 = "Bias due to missing data",
                      D6 = "Bias in measurement of outcomes",
                      D7 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    ),
                    quadas2 = c(
                      D1 = "Patient selection",
                      D2 = "Index test(s)",
                      D3 = "Reference standard",
                      D4 = "Flow and timing",
                      Overall = "Overall risk of bias"
                    ),
                    robins_e = c(
                      D1 = "Bias due to confounding",
                      D2 = "Bias in selection of participants into the study",
                      D3 = "Bias in classification of exposures",
                      D4 = "Bias due to departures from intended exposures",
                      D5 = "Bias due to missing data",
                      D6 = "Bias in measurement of outcomes",
                      D7 = "Bias in selection of the reported result",
                      Overall = "Overall risk of bias"
                    )
  )

  if (isTRUE(full_text)){
    return(domains)
  } else {
    return(names(domains))
  }
}

#' Internal function to get columns that contain Risk/of/Bias text
#'
#' @noRd
get_rob_text_columns <- function(df) {

  # Get the column names that start with D or are "Overall"
  rob_cols <- names(df)[grepl("^D\\d+$|^Overall$", names(df))]

  # Split "Risk of Bias" across the middle columns
  rob_text <- c("Risk", "of", "Bias")
  n_rob_cols <- length(rob_cols)

  # Calculate which columns get the text
  if (n_rob_cols >= 3) {
    start_pos <- ceiling((n_rob_cols - 3) / 2) + 1
    text_positions <- start_pos:(start_pos + 2)
  } else {
    text_positions <- 1:min(n_rob_cols, 3)
  }

  # Return only the columns that have the Risk/of/Bias text
  return(rob_cols[text_positions])
}

#' Internal function to create RoB labels based on how many domains are present
#'
#' @noRd
create_rob_labels <- function(df, measure, has_re = TRUE) {

  # Get the column names that start with D or are "Overall"
  rob_cols <- names(df)[grepl("^D\\d+$|^Overall$", names(df))]

  # Create base labels — adjust depending on whether random effects are present
  if (isTRUE(has_re)) {
    base_labels <- list(
      weighted_effect = gt::md(paste("Shrinkage", measure, "\n", "[95% CrI]")),
      unweighted_effect = gt::md(paste("Observed", measure, "\n", "[95% CI]"))
    )
  } else {
    # Common-effect: no shrinkage column; single column with both CI and CrI
    base_labels <- list(
      unweighted_effect = gt::md(paste("Observed", measure, "\n", "[95% CI]"))
    )
  }

  # Split "Risk of Bias" across the middle columns
  rob_text <- c("Risk", "of", "Bias")
  n_rob_cols <- length(rob_cols)

  # Calculate which columns get the text
  if (n_rob_cols >= 3) {
    start_pos <- ceiling((n_rob_cols - 3) / 2) + 1
    text_positions <- start_pos:(start_pos + 2)
  } else {
    text_positions <- 1:min(n_rob_cols, 3)
    rob_text <- rob_text[1:length(text_positions)]
  }

  # Create labels for ROB columns using purrr
  rob_labels <- purrr::imap(rob_cols, ~ {
    col_name <- .x
    position <- .y

    if (position %in% text_positions) {
      text_idx <- which(text_positions == position)
      gt::md(paste0("**", rob_text[text_idx], "**\n", col_name))
    } else {
      gt::md(paste0("&nbsp;\n", col_name))
    }
  }) |>
    purrr::set_names(rob_cols)

  # Combine base labels with ROB labels
  return(c(base_labels, rob_labels))
}
