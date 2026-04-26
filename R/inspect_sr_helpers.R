# ============================================================================
# INSPECT-SR: Helper Functions for Visualisation
# ============================================================================


#' Get INSPECT-SR Domain Names
#'
#' Returns the full domain names for the INSPECT-SR tool.
#'
#' @param full_text Logical. If TRUE (default), returns full descriptions.
#'   If FALSE, returns short codes.
#'
#' @return A named character vector.
#' @noRd
get_inspect_sr_domains <- function(full_text = TRUE) {
  domains <- c(
    D1 = "Inspecting post-publication notices",
    D2 = "Inspecting conduct, governance, and transparency",
    D3 = "Inspecting text and publication details",
    D4 = "Inspecting results in the study",
    Overall = "Overall trustworthiness judgement"
  )

  if (isTRUE(full_text)) {
    return(domains)
  } else {
    return(names(domains))
  }
}


#' Create INSPECT-SR Table
#'
#' Internal function to create the gt table for INSPECT-SR visualisation.
#'
#' @param df Data frame with INSPECT-SR judgements (columns: Study, D1-D4,
#'   Overall).
#' @param font Font family to use
#'
#' @return A gt table object.
#' @noRd
inspect_sr_table_fn <- function(df, font = NULL) {

  judgement_cols <- c("D1", "D2", "D3", "D4", "Overall")

  trust_text <- c("Trust-", "worthi-", "ness", "")
  col_labels <- list(
    Study   = "Study",
    D1      = gt::md(paste0("**", trust_text[1], "**\n", "D1")),
    D2      = gt::md(paste0("**", trust_text[2], "**\n", "D2")),
    D3      = gt::md(paste0("**", trust_text[3], "**\n", "D3")),
    D4      = gt::md(paste0("**", trust_text[4], "**\n", "D4")),
    Overall = gt::md("&nbsp;\nOverall")
  )

  df |>
    gt::gt() |>
    gt::cols_label(.list = col_labels) |>
    gt::cols_align(align = "left", columns = "Study") |>
    gt::cols_align(align = "center", columns = judgement_cols) |>
    gt::tab_options(
      column_labels.font.weight = "bold",
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font
    ) |>
    gt::opt_table_lines(extent = "none") |>
    gt::text_case_match(
      "Serious concerns" ~ fontawesome::fa(
        name = "circle-xmark", fill = "#e32400", height = "1.1em"
      ),
      "No concerns" ~ fontawesome::fa(
        name = "circle-check", fill = "#77bb41", height = "1.1em"
      ),
      "Some concerns" ~ fontawesome::fa(
        name = "circle-exclamation", fill = "#f5ec00", height = "1.1em"
      ),
      .locations = gt::cells_body(columns = judgement_cols)
    ) |>
    gt::sub_missing(columns = judgement_cols, missing_text = "\u2014")
}


#' Create Expanded INSPECT-SR Table (per-item, domain spanners)
#'
#' Internal. Builds a gt table with one column per INSPECT-SR item, spanner
#' headers grouping columns by domain, each domain's "Overall" column, and a
#' final study-level "Overall".
#'
#' @param df Data frame with Study, d1_1..d4_11 item columns, domain-level
#'   D1..D4 columns, and Overall.
#' @param font Font family.
#'
#' @return A gt table object.
#' @noRd
inspect_sr_item_table_fn <- function(df, font = NULL) {

  item_cols <- list(
    D1 = paste0("d1_", 1:3),
    D2 = paste0("d2_", 1:5),
    D3 = paste0("d3_", 1:2),
    D4 = paste0("d4_", 1:11)
  )

  judgement_cols <- unlist(unname(c(
    lapply(names(item_cols), function(dom) c(item_cols[[dom]], dom)),
    list("Overall")
  )))

  df <- df[, c("Study", judgement_cols), drop = FALSE]

  # Column labels: item columns become "1.1", "1.2", ...; domain columns
  # labelled "Overall"; the study-level Overall column has its label blanked
  # so the spanner "Overall" (at the Domain-1..Domain-4 level) serves as the
  # visible header.
  col_labels <- list(Study = "Study", Overall = gt::html("&nbsp;"))
  for (dom in names(item_cols)) {
    for (col in item_cols[[dom]]) {
      parts <- strsplit(sub("^d", "", col), "_")[[1]]
      col_labels[[col]] <- paste0(parts[1], ".", parts[2])
    }
    col_labels[[dom]] <- gt::md("**Overall**")
  }

  inspect_table <- df |>
    gt::gt() |>
    gt::cols_label(.list = col_labels) |>
    gt::cols_align(align = "left", columns = "Study") |>
    gt::cols_align(align = "center", columns = judgement_cols) |>
    gt::tab_spanner(label = "Domain 1", id = "Domain 1",
                    columns = c(item_cols$D1, "D1")) |>
    gt::tab_spanner(label = "Domain 2", id = "Domain 2",
                    columns = c(item_cols$D2, "D2")) |>
    gt::tab_spanner(label = "Domain 3", id = "Domain 3",
                    columns = c(item_cols$D3, "D3")) |>
    gt::tab_spanner(label = "Domain 4", id = "Domain 4",
                    columns = c(item_cols$D4, "D4")) |>
    gt::tab_spanner(label = "Overall", id = "OverallSpan",
                    columns = "Overall") |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_spanners()
    ) |>
    gt::tab_style(
      style = gt::cell_borders(sides = "left", color = "#888888",
                               weight = gt::px(1.5), style = "solid"),
      locations = list(
        gt::cells_body(columns = c("d1_1", "d2_1", "d3_1", "d4_1", "Overall")),
        gt::cells_column_labels(columns = c("d1_1", "d2_1", "d3_1", "d4_1",
                                            "Overall")),
        gt::cells_column_spanners(
          spanners = c("Domain 1", "Domain 2", "Domain 3", "Domain 4")
        )
      )
    ) |>
    gt::tab_style(
      style = gt::cell_borders(sides = "left", color = "#cccccc",
                               weight = gt::px(1), style = "solid"),
      locations = list(
        gt::cells_body(columns = names(item_cols)),
        gt::cells_column_labels(columns = names(item_cols))
      )
    ) |>
    gt::tab_options(
      column_labels.font.weight = "bold",
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      table.font.names = font,
      table.font.size = gt::px(12)
    ) |>
    gt::opt_table_lines(extent = "none") |>
    gt::text_case_match(
      "Serious concerns" ~ fontawesome::fa(
        name = "circle-xmark", fill = "#e32400", height = "1.1em"
      ),
      "No concerns" ~ fontawesome::fa(
        name = "circle-check", fill = "#77bb41", height = "1.1em"
      ),
      "Some concerns" ~ fontawesome::fa(
        name = "circle-exclamation", fill = "#f5ec00", height = "1.1em"
      ),
      .locations = gt::cells_body(columns = judgement_cols)
    ) |>
    gt::sub_missing(columns = judgement_cols, missing_text = "\u2014")

  inspect_table
}


#' Create INSPECT-SR Legend
#'
#' Internal function to create a legend for INSPECT-SR symbols and domains.
#'
#' @param font Font family to use
#'
#' @return A gt table object containing the legend.
#' @noRd
inspect_sr_legend_fn <- function(font = NULL) {

  domains <- get_inspect_sr_domains()
  domains <- domains[names(domains) != "Overall"]

  # Domain descriptions
  legend_data <- data.frame(
    Code = paste0("(", names(domains), ")"),
    Description = unname(domains),
    stringsAsFactors = FALSE
  )

  # Spacer row
  legend_data <- rbind(
    legend_data,
    data.frame(Code = "", Description = "", stringsAsFactors = FALSE)
  )

  # Judgement level rows
  judgement_rows <- data.frame(
    Code = c(" ", "No concerns", "Some concerns", "Serious concerns", "Not assessed"),
    Description = c(" ", "No concerns", "Some concerns", "Serious concerns",
                    "Not assessed (\u2014)"),
    stringsAsFactors = FALSE
  )

  legend_data <- rbind(legend_data, judgement_rows)

  # Build table
  legend_table <- legend_data |>
    gt::gt() |>
    gt::cols_label(Code = "", Description = "") |>
    gt::tab_header(title = gt::html("<b>INSPECT-SR legend</b>")) |>
    gt::tab_style(
      style = list(gt::cell_text(weight = "bold", size = gt::px(16))),
      locations = gt::cells_title(groups = "title")
    ) |>
    # Bold domain codes
    gt::tab_style(
      style = list(gt::cell_text(weight = "bold", size = gt::px(12))),
      locations = gt::cells_body(
        columns = "Code",
        rows = seq_len(length(domains))
      )
    ) |>
    # Transform judgement text to icons
    gt::text_transform(
      locations = gt::cells_body(
        columns = "Code",
        rows = (length(domains) + 2):nrow(legend_data)
      ),
      fn = function(x) {
        dplyr::case_when(
          x == "Serious concerns" ~ as.character(
            fontawesome::fa(name = "circle-xmark", fill = "#e32400", height = "1.5em")
          ),
          x == "No concerns" ~ as.character(
            fontawesome::fa(name = "circle-check", fill = "#77bb41", height = "1.5em")
          ),
          x == "Some concerns" ~ as.character(
            fontawesome::fa(name = "circle-exclamation", fill = "#f5ec00", height = "1.5em")
          ),
          x == "Not assessed" ~ "\u2014",
          TRUE ~ x
        )
      }
    ) |>
    # Bold/italic judgement descriptions
    gt::tab_style(
      style = list(gt::cell_text(weight = "bold", style = "italic")),
      locations = gt::cells_body(
        columns = "Description",
        rows = (length(domains) + 2):nrow(legend_data)
      )
    ) |>
    gt::tab_style(
      style = "padding-left:2px;",
      locations = gt::cells_body(columns = "Code", rows = gt::everything())
    ) |>
    gt::cols_width(
      Code ~ gt::px(40),
      Description ~ gt::px(340)
    ) |>
    gt::tab_options(
      table.border.top.color = "white",
      table.border.bottom.color = "white",
      heading.align = "left",
      heading.padding = gt::px(6),
      table.font.size = gt::px(12),
      table.font.names = font
    ) |>
    gt::opt_table_lines(extent = "none")

  return(legend_table)
}
