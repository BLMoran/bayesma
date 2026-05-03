# ============================================================================
# Study Table, Template Download, and App Launcher
# ============================================================================

# RoB domain columns by tool
.rob_cols <- list(
  rob2             = c("D1", "D2", "D3", "D4", "D5", "Overall"),
  robins_i         = c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "Overall"),
  newcastle_ottawa = c("Selection", "Comparability", "Outcome", "Total"),
  quadas2          = c("PD1", "PD2", "PD3", "PD4", "AC1", "AC2", "AC3", "Overall")
)

# Detect RoB tool from column names
.detect_rob_tool <- function(data) {
  nms <- names(data)
  if (all(c("D1", "D2", "D3", "D4", "D5", "D6", "D7") %in% nms)) return("robins_i")
  if (all(c("D1", "D2", "D3", "D4", "D5") %in% nms))              return("rob2")
  if (all(c("PD1", "PD2", "PD3", "PD4") %in% nms))               return("quadas2")
  if (all(c("Selection", "Comparability") %in% nms))               return("newcastle_ottawa")
  "none"
}

# Detect outcome type from column names
.detect_outcome_type <- function(data) {
  nms <- names(data)
  has_cont <- all(c("Mean_Control", "Mean_Intervention") %in% nms)
  has_bin  <- all(c("Event_Control", "Event_Intervention") %in% nms)
  if (has_cont && !has_bin) return("continuous")
  if (has_bin  && !has_cont) return("binary")
  if (has_cont && has_bin) {
    cli::cli_inform(
      "Both continuous and binary columns found; defaulting to {.val continuous}.",
      "Set {.arg outcome_type} explicitly to override."
    )
    return("continuous")
  }
  cli::cli_abort(
    "Could not detect outcome type. Ensure the data contains either \\
    {.val Mean_Control} / {.val Mean_Intervention} (continuous) or \\
    {.val Event_Control} / {.val Event_Intervention} (binary)."
  )
}

# Apply traffic-light RoB formatting
# text_case_match() matches by exact value, so non-RoB cells are unaffected.
.apply_rob_formatting <- function(tbl, rob_tool, present_cols) {
  if (length(present_cols) == 0) return(tbl)

  locs <- gt::cells_body(columns = present_cols)

  if (rob_tool %in% c("rob2", "robins_i")) {
    tbl <- tbl |>
      gt::text_case_match(
        "Low"           ~ fontawesome::fa("circle-minus",       fill = "#77bb41", height = "1.1em"),
        "Some concerns" ~ fontawesome::fa("circle-question",    fill = "#f5ec00", height = "1.1em"),
        "Moderate"      ~ fontawesome::fa("circle-question",    fill = "#f5a623", height = "1.1em"),
        "Serious"       ~ fontawesome::fa("circle-exclamation", fill = "#e07000", height = "1.1em"),
        "High"          ~ fontawesome::fa("circle-plus",        fill = "#e32400", height = "1.1em"),
        "Critical"      ~ fontawesome::fa("circle-plus",        fill = "#8b0000", height = "1.1em"),
        .locations = locs
      )
  } else if (rob_tool == "quadas2") {
    tbl <- tbl |>
      gt::text_case_match(
        "Low"     ~ fontawesome::fa("circle-minus",    fill = "#77bb41", height = "1.1em"),
        "High"    ~ fontawesome::fa("circle-plus",     fill = "#e32400", height = "1.1em"),
        "Unclear" ~ fontawesome::fa("circle-question", fill = "#f5ec00", height = "1.1em"),
        .locations = locs
      )
  }
  tbl
}


#' Create a Table of Included Studies
#'
#' Produces a formatted [gt::gt()] table from a study-level data frame. Handles
#' continuous and binary outcomes, subgroup grouping, and traffic-light
#' risk-of-bias formatting for RoB 2, ROBINS-I, NOS, and QUADAS-2.
#'
#' @param data A data frame or tibble with one row per study. Expected columns
#'   depend on `outcome_type` and `rob_tool` — see the bundled templates
#'   ([download_template()]) for the exact layout.
#' @param rob_tool Risk-of-bias tool. One of `"auto"` (detect from column
#'   names), `"none"`, `"rob2"`, `"robins_i"`, `"newcastle_ottawa"`, or
#'   `"quadas2"`.
#' @param outcome_type Whether outcomes are `"continuous"`, `"binary"`, or
#'   `"auto"` (detect from column names).
#' @param subgroup_col Unquoted name of the subgroup column, or `NULL` (no
#'   grouping).
#' @param font Font family for the table. `NULL` uses the gt default.
#'
#' @return A [gt::gt()] table object.
#'
#' @examples
#' \dontrun{
#' data(binary_outcome)
#' study_table(binary_outcome)
#'
#' data(cont_outcome)
#' study_table(cont_outcome)
#' }
#'
#' @export
study_table <- function(data,
                        rob_tool     = c("auto", "none", "rob2", "robins_i",
                                         "newcastle_ottawa", "quadas2"),
                        outcome_type = c("auto", "continuous", "binary"),
                        subgroup_col = NULL,
                        font         = NULL) {

  rob_tool     <- rlang::arg_match(rob_tool)
  outcome_type <- rlang::arg_match(outcome_type)

  if (rob_tool     == "auto") rob_tool     <- .detect_rob_tool(data)
  if (outcome_type == "auto") outcome_type <- .detect_outcome_type(data)

  # ensym() captures the bare name without evaluating — handles NULL gracefully
  subgroup_nm <- tryCatch(
    rlang::as_name(rlang::ensym(subgroup_col)),
    error = function(e) NULL
  )

  # Format study label
  data <- data |>
    dplyr::mutate(
      Study = dplyr::if_else(
        is.na(Year),
        as.character(Author),
        paste0(Author, " (", Year, ")")
      ),
      .before = 1
    ) |>
    dplyr::select(-dplyr::any_of(c("Author", "Year")))

  # Select outcome columns
  outcome_cols <- if (outcome_type == "continuous") {
    c("N_Control", "N_Intervention",
      "Mean_Control", "SD_Control",
      "Mean_Intervention", "SD_Intervention")
  } else {
    c("N_Control", "N_Intervention",
      "Event_Control", "Event_Intervention")
  }

  rob_domain_cols <- .rob_cols[[rob_tool]] %||% character(0)

  keep <- c("Study", subgroup_nm, outcome_cols, rob_domain_cols)
  df   <- data |> dplyr::select(dplyr::any_of(keep))

  # Build gt table
  tbl <- if (!is.null(subgroup_nm) && subgroup_nm %in% names(df)) {
    df |>
      gt::gt(groupname_col = subgroup_nm) |>
      gt::tab_options(row_group.font.weight = "bold") |>
      gt::tab_style(
        style     = gt::cell_fill(color = "grey95"),
        locations = gt::cells_row_groups()
      )
  } else {
    df |> gt::gt()
  }

  # Column labels
  label_map <- list(
    N_Control           = gt::md("N\n(Control)"),
    N_Intervention      = gt::md("N\n(Intervention)"),
    Mean_Control        = gt::md("Mean\n(Control)"),
    SD_Control          = gt::md("SD\n(Control)"),
    Mean_Intervention   = gt::md("Mean\n(Intervention)"),
    SD_Intervention     = gt::md("SD\n(Intervention)"),
    Event_Control       = gt::md("Events\n(Control)"),
    Event_Intervention  = gt::md("Events\n(Intervention)")
  )
  active_labels <- label_map[intersect(names(label_map), names(df))]

  # RoB domain labels
  domain_labels <- switch(rob_tool,
    rob2 = list(
      D1 = "D1\nRandomisation",
      D2 = "D2\nDeviations",
      D3 = "D3\nMissing data",
      D4 = "D4\nMeasurement",
      D5 = "D5\nReporting",
      Overall = "Overall"
    ),
    robins_i = list(
      D1 = "D1\nConfounding",
      D2 = "D2\nSelection",
      D3 = "D3\nClassification",
      D4 = "D4\nDeviations",
      D5 = "D5\nMissing data",
      D6 = "D6\nMeasurement",
      D7 = "D7\nReporting",
      Overall = "Overall"
    ),
    quadas2 = list(
      PD1 = "PD1\nPt. Selection",
      PD2 = "PD2\nIndex Test",
      PD3 = "PD3\nRef. Standard",
      PD4 = "PD4\nFlow & Timing",
      AC1 = "AC1\nPt. Selection",
      AC2 = "AC2\nIndex Test",
      AC3 = "AC3\nRef. Standard",
      Overall = "Overall"
    ),
    list()
  )
  active_domain_labels <- domain_labels[intersect(names(domain_labels), names(df))]
  all_labels <- c(active_labels, active_domain_labels)

  tbl <- tbl |>
    gt::cols_label(.list = all_labels) |>
    gt::cols_align(align = "left",   columns = Study) |>
    gt::cols_align(align = "center", columns = -Study) |>
    gt::tab_options(
      column_labels.font.weight = "bold",
      table.border.top.color    = "white",
      table.border.bottom.color = "white",
      table.font.names          = font
    ) |>
    gt::opt_table_lines(extent = "none")

  # Traffic-light icons for RoB / QUADAS-2
  present_rob <- intersect(rob_domain_cols, names(df))
  tbl <- .apply_rob_formatting(tbl, rob_tool, present_rob)

  # NOS: shade Total column by score
  if (rob_tool == "newcastle_ottawa" && "Total" %in% names(df)) {
    tbl <- tbl |>
      gt::data_color(
        columns = Total,
        method  = "numeric",
        palette = c("#e32400", "#f5ec00", "#77bb41"),
        domain  = c(0, 9)
      )
  }

  tbl
}


#' Download a Data Entry Template
#'
#' Writes a CSV template to disk so you know exactly which columns are expected
#' for each input data frame. Templates have header rows only; fill in one row
#' per study.
#'
#' @param type Template type. One of `"continuous"`, `"binary"`, `"rob2"`,
#'   `"robins_i"`, `"newcastle_ottawa"`, `"quadas2"`, or `"inspect_sr"`.
#' @param path Directory to write the template into. Defaults to the current
#'   working directory.
#'
#' @return Invisibly returns the path to the written file.
#'
#' @examples
#' \dontrun{
#' download_template("binary")
#' download_template("rob2", path = "~/Desktop")
#' }
#'
#' @export
download_template <- function(
    type = c("continuous", "binary", "rob2", "robins_i",
             "newcastle_ottawa", "quadas2", "inspect_sr"),
    path = NULL) {

  type <- rlang::arg_match(type)
  path <- path %||% getwd()

  src  <- system.file("templates", paste0(type, ".csv"), package = "bayesma")
  if (!nzchar(src)) {
    cli::cli_abort("Template file for {.val {type}} not found in package installation.")
  }
  dest <- file.path(path, paste0(type, "_template.csv"))
  file.copy(src, dest, overwrite = TRUE)
  cli::cli_inform("Template written to {.path {dest}}")
  invisible(dest)
}


#' Launch the bayesma Data Entry App
#'
#' Opens a local Shiny app for entering study details, risk-of-bias
#' assessments, and INSPECT-SR judgements. Data can be exported as CSV files
#' ready for use with [study_table()] and [inspect_sr()].
#'
#' @return Called for its side effect; returns `NULL` invisibly.
#'
#' @examples
#' \dontrun{
#' launch_data_entry()
#' }
#'
#' @export
launch_data_entry <- function() {
  app_path <- system.file("shiny", "data_entry", package = "bayesma")
  if (!nzchar(app_path)) {
    cli::cli_abort("Shiny app not found. Re-install {.pkg bayesma} to restore it.")
  }
  shiny::runApp(app_path, launch.browser = TRUE)
  invisible(NULL)
}
