library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(readxl)

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ── Constants ──────────────────────────────────────────────────────────────────

BASE_STUDY_COLS <- c(
  "Author", "Year", "Study_Design", "Subgroup",
  "Control", "Intervention", "Follow_Up", "Outcomes"
)

OUTCOME_TYPE_COLS <- list(
  binary     = c("N_Control", "N_Intervention",
                 "Event_Control", "Event_Intervention"),
  continuous = c("N_Control", "N_Intervention",
                 "Mean_Control", "SD_Control", "Mean_Intervention", "SD_Intervention"),
  count      = c("N_Control", "N_Intervention",
                 "Count_Control", "Time_Control",
                 "Count_Intervention", "Time_Intervention"),
  tte        = c("N_Control", "N_Intervention", "HR", "CI_Lower", "CI_Upper"),
  proportion = c("N", "Events")
)

OUTCOME_TYPE_LABELS <- c(
  binary     = "Binary",
  continuous = "Continuous",
  count      = "Count",
  tte        = "Time-to-event",
  proportion = "Proportion"
)

ROB_DOMAIN_COLS <- list(
  rob2             = c("D1", "D2", "D3", "D4", "D5", "Overall"),
  rob2_crt         = c("D1", "D2", "D3", "D4", "D5", "D6", "Overall"),
  rob2_xo          = c("D1", "D2", "D3", "D4", "D5", "Overall"),
  robins_i         = c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "Overall"),
  robins_ii        = c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "Overall"),
  robins_e         = c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "Overall"),
  newcastle_ottawa = c("Selection", "Comparability", "Outcome", "Total"),
  quadas2          = c("PD1", "PD2", "PD3", "PD4", "AC1", "AC2", "AC3", "Overall")
)

ROB_DOMAIN_LABELS <- list(
  rob2 = c(
    D1 = "D1: Randomisation process",
    D2 = "D2: Deviations from intended interventions",
    D3 = "D3: Missing outcome data",
    D4 = "D4: Measurement of the outcome",
    D5 = "D5: Selection of reported results",
    Overall = "Overall"
  ),
  rob2_crt = c(
    D1 = "D1: Randomisation process",
    D2 = "D2: Identification and recruitment of participants",
    D3 = "D3: Deviations from intended interventions",
    D4 = "D4: Missing outcome data",
    D5 = "D5: Measurement of the outcome",
    D6 = "D6: Selection of reported results",
    Overall = "Overall"
  ),
  rob2_xo = c(
    D1 = "D1: Randomisation process (incl. carry-over)",
    D2 = "D2: Deviations from intended interventions",
    D3 = "D3: Missing outcome data",
    D4 = "D4: Measurement of the outcome",
    D5 = "D5: Selection of reported results",
    Overall = "Overall"
  ),
  robins_i = c(
    D1 = "D1: Confounding",
    D2 = "D2: Selection of participants",
    D3 = "D3: Classification of interventions",
    D4 = "D4: Deviations from intended interventions",
    D5 = "D5: Missing data",
    D6 = "D6: Measurement of outcomes",
    D7 = "D7: Selection of reported results",
    Overall = "Overall"
  ),
  robins_ii = c(
    D1 = "D1: Confounding",
    D2 = "D2: Selection of participants",
    D3 = "D3: Classification of interventions",
    D4 = "D4: Deviations from intended interventions",
    D5 = "D5: Missing data",
    D6 = "D6: Measurement of outcomes",
    D7 = "D7: Selection of reported results",
    Overall = "Overall"
  ),
  robins_e = c(
    D1 = "D1: Confounding",
    D2 = "D2: Selection of participants",
    D3 = "D3: Classification of exposures",
    D4 = "D4: Departures from intended exposures",
    D5 = "D5: Missing data",
    D6 = "D6: Measurement of outcomes",
    D7 = "D7: Selection of reported results",
    Overall = "Overall"
  ),
  newcastle_ottawa = c(
    Selection = "Selection (0-4)",
    Comparability = "Comparability (0-2)",
    Outcome = "Outcome (0-3)",
    Total = "Total (0-9)"
  ),
  quadas2 = c(
    PD1 = "PD1: Patient selection",
    PD2 = "PD2: Index test",
    PD3 = "PD3: Reference standard",
    PD4 = "PD4: Flow & timing",
    AC1 = "AC1: Patient selection",
    AC2 = "AC2: Index test",
    AC3 = "AC3: Reference standard",
    Overall = "Overall"
  )
)

INSPECT_DOMAIN_COLS <- c(
  "d1_1", "d1_2", "d1_3",
  "d2_1", "d2_2", "d2_3", "d2_4", "d2_5",
  "d3_1", "d3_2",
  "d4_1", "d4_2", "d4_3", "d4_4", "d4_5",
  "d4_6", "d4_7", "d4_8", "d4_9", "d4_10", "d4_11"
)

INSPECT_CALC_COLS <- c(
  "d4_3_pvalues",
  "d4_6_n_control", "d4_6_n_intervention", "d4_6_n_total",
  "d4_8_mean", "d4_8_n", "d4_8_decimals",
  "d4_9_statistic", "d4_9_df", "d4_9_pvalue"
)

INSPECT_LABELS <- c(
  d1_1  = "D1.1 Registration prospective",
  d1_2  = "D1.2 Registration consistent",
  d1_3  = "D1.3 SAP consistent",
  d2_1  = "D2.1 Ethical approval",
  d2_2  = "D2.2 Primary registry",
  d2_3  = "D2.3 Registered before enrolment",
  d2_4  = "D2.4 Protocol vs conduct",
  d2_5  = "D2.5 Other conduct concerns",
  d3_1  = "D3.1 All outcomes reported",
  d3_2  = "D3.2 Selective reporting",
  d4_1  = "D4.1 Flow consistent",
  d4_2  = "D4.2 Analysis populations",
  d4_3  = "D4.3 Carlisle test",
  d4_4  = "D4.4 Data errors",
  d4_5  = "D4.5 Data duplication",
  d4_6  = "D4.6 N-consistency",
  d4_7  = "D4.7 Per guidelines",
  d4_8  = "D4.8 GRIM test",
  d4_9  = "D4.9 p-value verify",
  d4_10 = "D4.10 Retraction",
  d4_11 = "D4.11 Other concerns"
)

INSPECT_CALC_LABELS <- c(
  d4_3_pvalues         = "D4.3: Baseline p-values (;-separated)",
  d4_6_n_control       = "D4.6: N (Control)",
  d4_6_n_intervention  = "D4.6: N (Intervention)",
  d4_6_n_total         = "D4.6: N (Total)",
  d4_8_mean            = "D4.8: Reported mean",
  d4_8_n               = "D4.8: Sample size (n)",
  d4_8_decimals        = "D4.8: Decimal places",
  d4_9_statistic       = "D4.9: Test statistic",
  d4_9_df              = "D4.9: Degrees of freedom",
  d4_9_pvalue          = "D4.9: Reported p-value"
)

ALL_INSPECT_COLS   <- c("Author", "Year", INSPECT_DOMAIN_COLS, INSPECT_CALC_COLS)
ALL_INSPECT_LABELS <- c(INSPECT_LABELS, INSPECT_CALC_LABELS)

# ── Helpers ────────────────────────────────────────────────────────────────────

COL_TYPE_CHOICES <- c("Text" = "text", "Integer" = "integer",
                      "Numeric" = "numeric", "Date" = "date")

coerce_col_type <- function(x, type) {
  if (is.null(type) || identical(type, "text")) return(as.character(x))
  v <- as.character(x)
  v[v == ""] <- NA
  out <- switch(type,
    integer = suppressWarnings(as.integer(v)),
    numeric = suppressWarnings(as.numeric(v)),
    date    = suppressWarnings(as.Date(v)),
    as.character(v)
  )
  out
}

apply_col_types <- function(df, types) {
  if (is.null(types) || length(types) == 0 || nrow(df) == 0) return(df)
  for (col in intersect(names(df), names(types))) {
    df[[col]] <- coerce_col_type(df[[col]], types[[col]])
  }
  df
}

empty_df <- function(cols) {
  structure(
    replicate(length(cols), character(0), simplify = FALSE),
    names = cols, class = "data.frame", row.names = integer(0)
  )
}

blank_row <- function(cols, vals = list()) {
  row <- setNames(as.list(rep("", length(cols))), cols)
  for (nm in names(vals)) row[[nm]] <- vals[[nm]]
  as.data.frame(row, stringsAsFactors = FALSE)
}

study_label <- function(author, year) paste0(author, " (", year, ")")

outcome_prefix_cols <- function(label, type) {
  data_cols <- OUTCOME_TYPE_COLS[[type]]
  if (nzchar(label))
    c(paste0(label, "_Name"), paste0(label, "_Type"), paste0(label, "_", data_cols))
  else
    data_cols
}

sync_rows <- function(target, study_df, id_cols, extra_cols) {
  all_cols <- c(id_cols, extra_cols)
  if (nrow(study_df) == 0) return(target[integer(0), , drop = FALSE])

  target_keys <- if (nrow(target) > 0)
    paste(target[[id_cols[1]]], target[[id_cols[2]]], sep = "|||")
  else character(0)

  keep <- logical(nrow(target))
  seen <- character(0)
  need <- empty_df(all_cols)

  for (i in seq_len(nrow(study_df))) {
    if (!nzchar(study_df$Author[i])) next
    k <- paste(study_df$Author[i], study_df$Year[i], sep = "|||")
    if (k %in% seen) next
    seen <- c(seen, k)

    idx <- which(target_keys == k & !keep)
    if (length(idx) > 0) {
      keep[idx[1]] <- TRUE
    } else {
      need <- bind_rows(need, blank_row(all_cols, list(
        Author = study_df$Author[i], Year = study_df$Year[i]
      )))
    }
  }
  bind_rows(target[keep, , drop = FALSE], need)
}

sync_outcome_rows <- function(target, study_df) {
  if (nrow(study_df) == 0) return(target[integer(0), , drop = FALSE])
  has_sg <- "Study_Group" %in% names(study_df)

  # Add Study_Group column to target if study has it but target doesn't
  if (has_sg && !"Study_Group" %in% names(target)) {
    idx_yr <- which(names(target) == "Year")
    if (length(idx_yr) && nrow(target) > 0) {
      left   <- target[, seq_len(idx_yr), drop = FALSE]
      right  <- target[, seq(idx_yr + 1, ncol(target)), drop = FALSE]
      target <- cbind(left,
                      data.frame(Study_Group = rep("", nrow(target)), stringsAsFactors = FALSE),
                      right)
    } else if (!nrow(target)) {
      target$Study_Group <- character(0)
    }
  }

  all_cols <- if (length(names(target))) names(target) else {
    if (has_sg) c("Author", "Year", "Study_Group") else c("Author", "Year")
  }

  # Without Study_Group: plain Author+Year sync
  if (!has_sg || !"Study_Group" %in% names(target)) {
    target_keys <- if (nrow(target) > 0)
      paste(target$Author, target$Year, sep = "|||")
    else character(0)
    keep <- logical(nrow(target))
    seen <- character(0)
    need <- empty_df(all_cols)
    for (i in seq_len(nrow(study_df))) {
      if (!nzchar(study_df$Author[i])) next
      k <- paste(study_df$Author[i], study_df$Year[i], sep = "|||")
      if (k %in% seen) next
      seen <- c(seen, k)
      idx <- which(target_keys == k & !keep)
      if (length(idx) > 0) {
        keep[idx[1]] <- TRUE
      } else {
        need <- bind_rows(need, blank_row(all_cols,
                          list(Author = study_df$Author[i], Year = study_df$Year[i])))
      }
    }
    return(bind_rows(target[keep, , drop = FALSE], need))
  }

  # With Study_Group: 1:1 mapping from each study row to a target row.
  # Promote-unassigned handles transition from blank to filled Study_Group.
  # Unmatched target rows are removed.
  target_full <- paste(target$Author, target$Year, target$Study_Group, sep = "|||")
  target_ay   <- paste(target$Author, target$Year, sep = "|||")
  keep <- logical(nrow(target))
  need <- empty_df(all_cols)

  for (i in seq_len(nrow(study_df))) {
    if (!nzchar(study_df$Author[i])) next

    sg <- study_df$Study_Group[i]
    ay <- paste(study_df$Author[i], study_df$Year[i], sep = "|||")
    k  <- paste(study_df$Author[i], study_df$Year[i], sg, sep = "|||")

    # Exact match — claim the first unclaimed target row with this key
    idx <- which(target_full == k & !keep)
    if (length(idx) > 0) {
      keep[idx[1]] <- TRUE
      next
    }

    # Promote: target row with same Author+Year but blank Study_Group
    if (nzchar(sg)) {
      idx <- which(target_ay == ay & !nzchar(target$Study_Group) & !keep)
      if (length(idx) > 0) {
        target$Study_Group[idx[1]] <- sg
        target_full[idx[1]] <- k
        keep[idx[1]] <- TRUE
        next
      }
    }

    vals <- list(Author = study_df$Author[i], Year = study_df$Year[i], Study_Group = sg)
    need <- bind_rows(need, blank_row(all_cols, vals))
  }

  bind_rows(target[keep, , drop = FALSE], need)
}

sync_inspect_rows <- function(target, study_df) {
  if (nrow(study_df) == 0) return(target[integer(0), , drop = FALSE])

  target_keys <- if (nrow(target) > 0)
    paste(target$Author, target$Year, sep = "|||")
  else character(0)

  keep <- logical(nrow(target))
  seen <- character(0)
  need <- empty_df(ALL_INSPECT_COLS)

  for (i in seq_len(nrow(study_df))) {
    if (!nzchar(study_df$Author[i])) next
    k <- paste(study_df$Author[i], study_df$Year[i], sep = "|||")
    if (k %in% seen) next
    seen <- c(seen, k)

    idx <- which(target_keys == k & !keep)
    if (length(idx) > 0) {
      keep[idx[1]] <- TRUE
    } else {
      need <- bind_rows(need, blank_row(ALL_INSPECT_COLS, list(
        Author = study_df$Author[i], Year = study_df$Year[i]
      )))
    }
  }
  bind_rows(target[keep, , drop = FALSE], need)
}

add_study_col <- function(df) df

merge_all_sections <- function(study, outcome, rob, inspect) {
  dfs <- list(study = study, outcome = outcome, rob = rob, inspect = inspect)

  base_cols <- c("Author", "Year")
  use_sg    <- any(vapply(dfs, function(d) "Study_Group" %in% names(d), logical(1)))

  # Sections WITH Study_Group: tag arm rows with .arm_n to disambiguate duplicates.
  # Sections WITHOUT Study_Group: dedupe to one row per (Author, Year) so they
  # broadcast across all arms via a plain (Author, Year) left-join.
  prepared <- lapply(dfs, function(d) {
    if (nrow(d) == 0 || !all(base_cols %in% names(d))) return(d)
    if ("Study_Group" %in% names(d)) {
      d |> dplyr::mutate(.arm_n = dplyr::row_number(),
                         .by = dplyr::all_of(c("Author", "Year", "Study_Group")))
    } else {
      d[!duplicated(d[, base_cols, drop = FALSE]), , drop = FALSE]
    }
  })

  # Spine: arm-level keys from SG-aware sections + any (A, Y) found only in
  # non-SG sections (those broadcast as a single empty-Study_Group entry).
  if (use_sg) {
    sg_parts <- purrr::map(prepared, function(d) {
      if (nrow(d) == 0 || !"Study_Group" %in% names(d)) return(NULL)
      dplyr::distinct(d[, c(base_cols, "Study_Group", ".arm_n"), drop = FALSE])
    }) |> purrr::compact()

    spine <- if (length(sg_parts))
      purrr::list_rbind(sg_parts) |> dplyr::distinct()
    else
      data.frame(Author = character(), Year = character(),
                 Study_Group = character(), .arm_n = integer(),
                 stringsAsFactors = FALSE)

    nonsg_ay <- purrr::map(prepared, function(d) {
      if (nrow(d) == 0 || "Study_Group" %in% names(d)) return(NULL)
      dplyr::distinct(d[, base_cols, drop = FALSE])
    }) |> purrr::compact()

    if (length(nonsg_ay) > 0) {
      nonsg_ay <- purrr::list_rbind(nonsg_ay) |> dplyr::distinct()
      spine_ay <- dplyr::distinct(spine[, base_cols, drop = FALSE])
      missing  <- dplyr::anti_join(nonsg_ay, spine_ay, by = base_cols)
      if (nrow(missing) > 0) {
        missing$Study_Group <- ""
        missing$.arm_n      <- 1L
        spine <- dplyr::bind_rows(spine, missing)
      }
    }
  } else {
    spine <- purrr::map(prepared, function(d) {
      if (nrow(d) == 0 || !all(base_cols %in% names(d))) return(NULL)
      dplyr::distinct(d[, base_cols, drop = FALSE])
    }) |> purrr::compact() |> purrr::list_rbind() |> dplyr::distinct()
  }

  result <- spine
  for (d in prepared) {
    if (nrow(d) == 0) next
    has_d_sg <- "Study_Group" %in% names(d)
    join_keys <- if (use_sg && has_d_sg) c(base_cols, "Study_Group", ".arm_n") else base_cols
    extra <- setdiff(names(d), c(join_keys, ".arm_n"))
    if (length(extra) == 0) next
    if (!has_d_sg && ".arm_n" %in% names(d)) d$.arm_n <- NULL
    result <- dplyr::left_join(result, d, by = join_keys)
  }
  result$.arm_n <- NULL
  result
}

read_uploaded <- function(file) {
  ext <- tolower(tools::file_ext(file$name))
  if (ext == "csv") {
    list(type = "flat",
         df   = read.csv(file$datapath, stringsAsFactors = FALSE, check.names = FALSE))
  } else if (ext %in% c("xls", "xlsx")) {
    list(type = "flat",
         df   = as.data.frame(readxl::read_excel(file$datapath), stringsAsFactors = FALSE))
  } else if (ext %in% c("rda", "rdata")) {
    env  <- new.env(parent = emptyenv())
    load(file$datapath, envir = env)
    dfs  <- Filter(is.data.frame, as.list(env))
    if (length(dfs) == 0) stop("No data frames found in the .rda file.")
    if (length(dfs) == 1) {
      list(type = "flat", df = dfs[[1]])
    } else {
      list(type = "rda", dfs = dfs)
    }
  } else if (ext == "rds") {
    obj <- readRDS(file$datapath)
    if (is.data.frame(obj)) {
      list(type = "flat", df = obj)
    } else if (is.list(obj) && all(vapply(obj, is.data.frame, logical(1)))) {
      list(type = "rda", dfs = obj)
    } else {
      stop("RDS file must contain a data frame or a named list of data frames.")
    }
  } else {
    stop("Unsupported file type. Please use CSV, XLSX, RDA, or RDS.")
  }
}

best_col_match <- function(expected, file_cols) {
  idx <- which(tolower(file_cols) == tolower(expected))
  if (length(idx)) return(file_cols[idx[1]])
  "(Skip)"
}

dt_base <- function(df, id, sel_mode = "single", editable = TRUE) {
  datatable(
    df, editable = editable, selection = sel_mode,
    rownames = FALSE,
    class    = "cell-border stripe",
    options  = list(
      pageLength = 50, scrollX = TRUE, dom = "tip",
      autoWidth = FALSE, scrollCollapse = TRUE,
      initComplete = JS(
        "function(settings, json) {",
        "  var api = this.api();",
        "  $(api.table().body()).on('click', 'td', function() {",
        "    var td = this;",
        "    if ($(td).find('input, select, textarea').length) return;",
        "    setTimeout(function() {",
        "      if (!$(td).find('input, select, textarea').length) {",
        "        $(td).trigger('dblclick');",
        "      }",
        "    }, 50);",
        "  });",
        "}"
      )
    )
  )
}

# ── CSS ────────────────────────────────────────────────────────────────────────

app_css <- "
html, body {
  height: auto !important;
  min-height: 0 !important;
  overflow: visible !important;
}
.container-fluid, .container, .bslib-page-fill, .html-fill-container {
  height: auto !important;
  min-height: 0 !important;
  overflow: visible !important;
}
body { font-size: 0.9rem; color: #2E3A4A; }
.section-card { margin-bottom: 1.5rem; }
.card-header-custom {
  background: #3191bf; color: white;
  padding: 0.6rem 1rem;
  display: flex; align-items: center; justify-content: space-between;
  border-radius: 4px 4px 0 0;
}
.card-header-custom h5 { margin: 0; font-size: 1rem; font-weight: 600; }
.app-blurb {
  background: #EAF4FB; border-left: 4px solid #3191bf;
  padding: 10px 16px; margin-bottom: 20px;
  border-radius: 0 4px 4px 0; font-size: 0.88rem; color: #2E3A4A;
}
/* RoB buttons */
.btn-rob-low      { background:#77bb41; color:white; border:none; }
.btn-rob-low:hover { background:#5e9430; color:white; }
.btn-rob-some     { background:#f5c518; color:#2E3A4A; border:none; }
.btn-rob-some:hover { background:#d4a814; color:#2E3A4A; }
.btn-rob-moderate { background:#f5a623; color:white; border:none; }
.btn-rob-moderate:hover { background:#d48c1c; color:white; }
.btn-rob-serious  { background:#e07000; color:white; border:none; }
.btn-rob-serious:hover  { background:#c06000; color:white; }
.btn-rob-high     { background:#e32400; color:white; border:none; }
.btn-rob-high:hover { background:#c21e00; color:white; }
.btn-rob-critical { background:#8b0000; color:white; border:none; }
.btn-rob-critical:hover { background:#6a0000; color:white; }
.btn-rob-unclear  { background:#adb5bd; color:white; border:none; }
.btn-rob-ni       { background:#6c8ebf; color:white; border:none; }
.btn-rob-ni:hover { background:#5a7aad; color:white; }
/* INSPECT buttons */
.btn-inspect-none    { background:#77bb41; color:white; border:none; }
.btn-inspect-none:hover { background:#5e9430; color:white; }
.btn-inspect-some    { background:#f5c518; color:#2E3A4A; border:none; }
.btn-inspect-some:hover { background:#d4a814; color:#2E3A4A; }
.btn-inspect-serious { background:#e32400; color:white; border:none; }
.btn-inspect-serious:hover { background:#c21e00; color:white; }
/* Layout helpers */
.selection-info {
  padding: 6px 12px; background:#EAF4FB;
  border-left: 3px solid #3191bf;
  font-size: 0.85rem; margin-bottom: 8px;
  border-radius: 0 4px 4px 0;
}
.btn-bar { display:flex; flex-wrap:wrap; gap:6px; align-items:center; margin-bottom:8px; }
.btn-bar .label { font-weight:600; font-size:0.85rem; margin-right:4px; color:#2E3A4A; }
.dl-row { display:flex; flex-wrap:wrap; gap:10px; align-items:flex-start; }
.map-row { display:flex; align-items:center; gap:8px; margin-bottom:4px; }
.map-label { min-width:220px; font-size:0.85rem; font-weight:500; color:#2E3A4A; }
/* Excel-like tables — no hover highlight, soft cell borders */
table.dataTable { border-collapse: collapse !important; }
table.dataTable thead tr th {
  background: #f0f4f8 !important;
  border: 1px solid #c8d6e0 !important;
  font-weight: 600; color: #2E3A4A; white-space: nowrap;
}
table.dataTable.stripe tbody tr.odd > * {
  box-shadow: none; background-color: #f7fafd;
}
table.dataTable.stripe tbody tr.even > * {
  box-shadow: none; background-color: #ffffff;
}
table.dataTable.cell-border tbody td {
  border-right: 1px solid #dce8ef !important;
  border-bottom: 1px solid #dce8ef !important;
  padding: 10px 10px !important;
  line-height: 1.5;
  cursor: cell;
}
table.dataTable tbody tr { height: 44px; }
table.dataTable tbody tr:hover > * { box-shadow: none !important; }
table.dataTable tbody tr.selected > * {
  box-shadow: inset 0 0 0 9999px rgba(49,145,191,0.08) !important;
  color: inherit !important;
}
/* Cell edit input — clear, full-width, obvious edit state */
table.dataTable tbody td input[type='text'],
table.dataTable tbody td input:not([type]),
table.dataTable tbody td textarea,
table.dataTable tbody td select {
  outline: none !important;
  border: 2px solid #3191bf !important;
  border-radius: 3px;
  background: #ffffff !important;
  padding: 3px 6px !important;
  width: 100% !important;
  box-sizing: border-box;
  font-size: 0.9rem;
  color: #2E3A4A;
  box-shadow: 0 0 0 3px rgba(49,145,191,0.15) !important;
}
/* Import box — resizable dropzone */
.import-box {
  resize: vertical;
  overflow: auto;
  min-height: 90px;
  border: 2px dashed #c8d6e0;
  border-radius: 6px;
  padding: 10px 12px;
  background: #fafcfe;
}
/* Import outcome blocks */
.imp-outcome-block {
  border: 1px solid #c8d6e0; border-radius: 5px;
  padding: 10px 12px; margin-bottom: 10px; background: #fafcfe;
}
.imp-outcome-header { display:flex; align-items:flex-end; gap:8px; margin-bottom:4px; }
/* Help icon in card headers */
.btn-help {
  background: none !important; border: none !important; box-shadow: none !important;
  color: rgba(255,255,255,0.8) !important; padding: 0 2px !important;
  font-size: 1rem; line-height: 1; cursor: pointer;
}
.btn-help:hover { color: white !important; }
.btn-help:focus { outline: none !important; box-shadow: none !important; }
"

# ── UI ─────────────────────────────────────────────────────────────────────────

otype_btn <- function(id, label, active = FALSE) {
  cls <- if (active) "btn-sm btn-primary" else "btn-sm btn-outline-secondary"
  actionButton(id, label, class = cls)
}

ui <- fluidPage(
  theme = bs_theme(bootswatch = "cosmo", primary = "#3191bf", secondary = "#3191bf"),
  tags$head(
    tags$style(HTML(app_css)),
    tags$script(src = "https://cdn.jsdelivr.net/npm/@iframe-resizer/child@5")
  ),

  div(style = "max-width:1400px; margin:0 auto; padding:16px;",

    div(class = "app-blurb",
      "Enter study data in each section below. ",
      strong("Double-click"), " any cell to edit it. Use ", strong("Sync studies"),
      " to propagate Author and Year from Study Details into the other sections. ",
      "Run ", code("bayesma::launch_data_entry()"), " locally for an offline version."
    ),

    # ── Upload ────────────────────────────────────────────────────────────────
    div(class = "section-card",
      div(class = "card-header-custom", h5("Upload Data")),
      div(class = "card-body border border-top-0 rounded-bottom p-3",
        div(class = "import-box",
          fluidRow(
            column(5,
              fileInput("upload_file", "Choose file",
                        accept = c(".csv", ".xls", ".xlsx", ".rda", ".rdata", ".rds"),
                        placeholder = "CSV, XLSX, RDA, or RDS")
            ),
            column(3, br(),
              actionButton("upload_map_btn", "Map & Import",
                           icon = icon("file-import"),
                           class = "btn-primary btn-sm mt-1")
            )
          )
        )
      )
    ),

    # ── Study Details ─────────────────────────────────────────────────────────
    div(class = "section-card",
      div(class = "card-header-custom",
        h5("Study Details"),
        div(
          actionButton("add_study",      "Add study",
                       icon = icon("plus"),    class = "btn-sm btn-light"),
          actionButton("add_col_study",  "Add column",
                       icon = icon("columns"), class = "btn-sm btn-outline-light ms-2"),
          actionButton("manage_study",   "Manage columns",
                       icon = icon("sliders"), class = "btn-sm btn-outline-light ms-2"),
          actionButton("del_study",      "Remove row",
                       icon = icon("minus"),   class = "btn-sm btn-danger ms-2"),
          actionButton("clear_study",    "Clear section",
                       icon = icon("trash"),   class = "btn-sm btn-danger ms-2")
        )
      ),
      div(class = "card-body border border-top-0 rounded-bottom p-3",
        div(style = "display:flex; align-items:center; justify-content:space-between; margin-bottom:6px;",
          p(class = "text-muted small mb-0",
            "Author and Year are used as study identifiers across all sections."),
          checkboxInput("incl_multiarm", "Include multi-arm group column",
                        value = FALSE, width = "auto")
        ),
        DTOutput("study_tbl")
      )
    ),

    # ── Outcome Details ───────────────────────────────────────────────────────
    div(class = "section-card",
      div(class = "card-header-custom",
        h5("Outcome Details"),
        div(
          actionButton("add_outcome_btn",    "Add outcome",
                       icon = icon("plus"),       class = "btn-sm btn-light"),
          actionButton("manage_outcome",     "Manage columns",
                       icon = icon("sliders"),    class = "btn-sm btn-outline-light ms-2"),
          actionButton("del_outcome_row",    "Remove row",
                       icon = icon("minus"),      class = "btn-sm btn-danger ms-2"),
          actionButton("del_outcome_group",  "Remove outcome",
                       icon = icon("table-columns"), class = "btn-sm btn-warning ms-2"),
          actionButton("clear_outcome",      "Clear section",
                       icon = icon("trash"),      class = "btn-sm btn-danger ms-2")
        )
      ),
      div(class = "card-body border border-top-0 rounded-bottom p-3",
        div(class = "btn-bar",
          span(class = "label", "Show type:"),
          actionButton("otype_all",  "All",            class = "btn-sm btn-primary"),
          actionButton("otype_bin",  "Binary",         class = "btn-sm btn-outline-secondary"),
          actionButton("otype_con",  "Continuous",     class = "btn-sm btn-outline-secondary"),
          actionButton("otype_cnt",  "Count",          class = "btn-sm btn-outline-secondary"),
          actionButton("otype_tte",  "Time-to-event",  class = "btn-sm btn-outline-secondary"),
          actionButton("otype_prop", "Proportion",     class = "btn-sm btn-outline-secondary"),
          actionButton("sync_outcome", "Sync studies",
                       icon = icon("arrows-rotate"), class = "btn-sm btn-outline-secondary ms-3",
                       onclick = "$('.dataTable input:focus, .dataTable select:focus, .dataTable textarea:focus').blur();")
        ),
        p(class = "text-muted small mb-2",
          "One row per study. Use Add Outcome to define outcomes — columns are added to the right."),
        DTOutput("outcome_tbl")
      )
    ),

    # ── Risk of Bias ──────────────────────────────────────────────────────────
    div(class = "section-card",
      div(class = "card-header-custom",
        div(style = "display:flex; align-items:center; gap:6px;",
          h5("Risk of Bias"),
          actionButton("rob_help", "", icon = icon("circle-question"),
                       class = "btn-help", title = "RoB tool guide")
        ),
        div(style = "display:flex; gap:8px; align-items:center;",
          span("Tool:", style = "font-size:0.85rem;"),
          selectInput("rob_tool", NULL,
                      choices = c("RoB 2" = "rob2",
                                  "RoB 2 (Cluster RCT)" = "rob2_crt",
                                  "RoB 2 (Crossover)" = "rob2_xo",
                                  "ROBINS-I" = "robins_i",
                                  "ROBINS-II" = "robins_ii",
                                  "ROBINS-E" = "robins_e",
                                  "Newcastle-Ottawa" = "newcastle_ottawa",
                                  "QUADAS-2" = "quadas2"),
                      width = "220px")
        )
      ),
      div(class = "card-body border border-top-0 rounded-bottom p-3",
        uiOutput("rob_btn_bar"),
        uiOutput("rob_sel_info"),
        div(class = "btn-bar mb-2",
          actionButton("sync_rob", "Sync studies",
                       icon = icon("arrows-rotate"), class = "btn-sm btn-outline-secondary",
                       onclick = "$('.dataTable input:focus, .dataTable select:focus, .dataTable textarea:focus').blur();"),
          actionButton("manage_rob", "Manage columns",
                       icon = icon("sliders"), class = "btn-sm btn-outline-secondary ms-2"),
          actionButton("del_rob_row", "Remove row",
                       icon = icon("minus"), class = "btn-sm btn-danger ms-2"),
          actionButton("clear_rob", "Clear section",
                       icon = icon("trash"), class = "btn-sm btn-danger ms-2")
        ),
        DTOutput("rob_tbl")
      )
    ),

    # ── INSPECT-SR ────────────────────────────────────────────────────────────
    div(class = "section-card",
      div(class = "card-header-custom",
        div(style = "display:flex; align-items:center; gap:6px;",
          h5("INSPECT-SR"),
          actionButton("inspect_help", "", icon = icon("circle-question"),
                       class = "btn-help", title = "INSPECT-SR scoring guide")
        )
      ),
      div(class = "card-body border border-top-0 rounded-bottom p-3",
        div(class = "btn-bar",
          span(class = "label", "Insert:"),
          actionButton("btn_isp_none",    "No concerns",
                       class = "btn-sm btn-inspect-none"),
          actionButton("btn_isp_some",    "Some concerns",
                       class = "btn-sm btn-inspect-some"),
          actionButton("btn_isp_serious", "Serious concerns",
                       class = "btn-sm btn-inspect-serious"),
          actionButton("btn_isp_clear",   "Clear",
                       class = "btn-sm btn-outline-secondary"),
          actionButton("sync_inspect", "Sync studies",
                       icon = icon("arrows-rotate"), class = "btn-sm btn-outline-secondary ms-3",
                       onclick = "$('.dataTable input:focus, .dataTable select:focus, .dataTable textarea:focus').blur();"),
          actionButton("manage_inspect", "Manage columns",
                       icon = icon("sliders"), class = "btn-sm btn-outline-secondary ms-2"),
          actionButton("del_inspect_row", "Remove row",
                       icon = icon("minus"), class = "btn-sm btn-danger ms-2"),
          actionButton("clear_inspect", "Clear section",
                       icon = icon("trash"), class = "btn-sm btn-danger ms-2")
        ),
        uiOutput("inspect_sel_info"),
        p(class = "text-muted small mb-2",
          "D4.3, D4.6, D4.8, D4.9 include input columns for bayesma automated tests."),
        DTOutput("inspect_tbl")
      )
    ),

    # ── Download ──────────────────────────────────────────────────────────────
    div(class = "section-card",
      div(class = "card-header-custom", h5("Download")),
      div(class = "card-body border border-top-0 rounded-bottom p-3",
        fluidRow(
          column(6,
            h6("Entered data"),
            radioButtons("dl_format", "Format:",
              choices  = c("RDS (.rds)" = "rds", "CSV (.csv)" = "csv",
                           "Excel (.xlsx)" = "xlsx"),
              selected = "rds", inline = TRUE),
            uiOutput("dl_buttons"),
            uiOutput("dl_hint")
          ),
          column(5,
            h6("Templates"),
            p(class = "text-muted small mb-1",
              "Blank templates showing required columns."),
            radioButtons("tpl_format", "Format:",
              choices  = c("RDS (.rds)" = "rds", "CSV (.csv)" = "csv",
                           "Excel (.xlsx)" = "xlsx"),
              selected = "rds", inline = TRUE),
            div(class = "dl-row mb-2",
              downloadButton("tpl_study",   "Study Details",
                             class = "btn-sm btn-primary"),
              downloadButton("tpl_outcome", "Outcome Data",
                             class = "btn-sm btn-primary"),
              downloadButton("tpl_rob",     "Risk of Bias",
                             class = "btn-sm btn-primary"),
              downloadButton("tpl_inspect", "INSPECT-SR",
                             class = "btn-sm btn-primary")
            ),
            downloadButton("tpl_all", "Download all templates",
                           class = "btn-sm btn-outline-primary")
          )
        )
      )
    ),

    tags$br()
  )
)

# ── Server ─────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  rv <- reactiveValues(
    study        = empty_df(BASE_STUDY_COLS),
    outcome      = empty_df(c("Author", "Year")),
    outcome_defs = list(),
    rob          = empty_df(c("Author", "Year", ROB_DOMAIN_COLS[["rob2"]])),
    inspect      = empty_df(ALL_INSPECT_COLS),
    types        = list(study = list(), outcome = list(),
                        rob = list(), inspect = list()),
    sel_rob      = NULL,
    upload_rda_dfs = list(),
    sel_inspect  = NULL,
    upload_df    = NULL
  )

  outcome_type_filter  <- reactiveVal("all")
  new_outcome_type     <- reactiveVal("binary")
  imp_tab_selected     <- reactiveVal(NULL)
  imp_outcomes         <- reactiveVal(list(list(idx = 1L, type = "binary")))
  imp_outcome_next_idx <- reactiveVal(2L)
  imp_rob_tool         <- reactiveVal("rob2")

  MAX_IMP_OUTCOMES <- 8L

  for (.i in seq_len(MAX_IMP_OUTCOMES)) {
    local({
      ii <- .i
      for (.t in names(OUTCOME_TYPE_LABELS)) {
        local({
          tt <- .t
          observeEvent(input[[paste0("imp_otype_", ii, "_", tt)]], {
            cur <- imp_outcomes()
            pos <- which(vapply(cur, `[[`, integer(1), "idx") == ii)
            if (length(pos)) { cur[[pos]]$type <- tt; imp_outcomes(cur) }
          }, ignoreInit = TRUE)
        })
      }
      observeEvent(input[[paste0("imp_rm_outcome_", ii)]], {
        cur <- imp_outcomes()
        if (length(cur) > 1)
          imp_outcomes(cur[vapply(cur, `[[`, integer(1), "idx") != ii])
      }, ignoreInit = TRUE)
    })
  }

  for (.t in names(ROB_DOMAIN_COLS)) {
    local({
      tt <- .t
      observeEvent(input[[paste0("imp_rob_tool_", tt)]], {
        imp_rob_tool(tt)
      }, ignoreInit = TRUE)
    })
  }

  observeEvent(input$imp_add_outcome, {
    idx <- imp_outcome_next_idx()
    imp_outcomes(c(imp_outcomes(), list(list(idx = idx, type = "binary"))))
    imp_outcome_next_idx(idx + 1L)
  }, ignoreInit = TRUE)

  # ── Upload ────────────────────────────────────────────────────────────────

  observeEvent(input$upload_map_btn, {
    req(input$upload_file)
    result <- tryCatch(read_uploaded(input$upload_file), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 8)
      NULL
    })
    if (is.null(result)) return()

    if (result$type == "rda") {
      rv$upload_df     <- result$dfs[[1]]
      rv$upload_rda_dfs <- result$dfs
    } else {
      rv$upload_df      <- result$df
      rv$upload_rda_dfs <- list()
    }

    imp_outcomes(list(list(idx = 1L, type = "binary")))
    imp_outcome_next_idx(2L)
    imp_rob_tool(input$rob_tool %||% "rob2")
    imp_tab_selected(NULL)

    showModal(modalDialog(
      title = "Map columns to sections",
      size  = "l",
      checkboxGroupInput("import_sections", "Import into:",
        choices  = c("Study Details" = "study", "Outcome Details" = "outcome",
                     "Risk of Bias"  = "rob",   "INSPECT-SR"      = "inspect"),
        selected = c("study", "outcome", "rob", "inspect"),
        inline   = TRUE),
      tags$hr(),
      uiOutput("upload_mapping_ui"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_import", "Import", class = "btn-primary")
      )
    ))
  })

  output$upload_mapping_ui <- renderUI({
    req(rv$upload_df, input$import_sections)
    file_cols  <- c("(Skip)", names(rv$upload_df))
    rda_names  <- names(rv$upload_rda_dfs)
    outcomes   <- imp_outcomes()
    rob_tool_v <- imp_rob_tool()

    rda_selector <- if (length(rda_names) > 1)
      selectInput("import_rda_df", "Dataset to map from (.rda contains multiple):",
                  choices = rda_names, selected = rda_names[1], width = "320px")
    else NULL

    tabs <- lapply(input$import_sections, function(sec) {
      if (sec == "outcome") {
        outcome_blocks <- lapply(seq_along(outcomes), function(pos) {
          o      <- outcomes[[pos]]
          ii     <- o$idx
          ty     <- o$type
          ty_cols <- OUTCOME_TYPE_COLS[[ty]]
          cur_lbl <- isolate(input[[paste0("imp_label_", ii)]]) %||% ""

          col_rows <- lapply(ty_cols, function(col) {
            matched <- best_col_match(col, names(rv$upload_df))
            div(class = "map-row",
              div(class = "map-label", col),
              div(style = "flex:1;",
                selectInput(paste0("map_outcome_", ii, "_", col), NULL,
                            choices = file_cols, selected = matched, width = "100%")))
          })

          div(class = "imp-outcome-block",
            div(class = "imp-outcome-header",
              textInput(paste0("imp_label_", ii),
                        "Outcome label (prefix — only needed for multiple outcomes):",
                        value = cur_lbl, width = "320px",
                        placeholder = "e.g. Primary — leave blank if only one outcome"),
              if (pos > 1)
                actionButton(paste0("imp_rm_outcome_", ii), "",
                             icon = icon("trash"), class = "btn-sm btn-outline-danger mt-4")
            ),
            div(class = "btn-bar mb-1",
              span(class = "label", "Type:"),
              lapply(names(OUTCOME_TYPE_LABELS), function(t)
                actionButton(paste0("imp_otype_", ii, "_", t), OUTCOME_TYPE_LABELS[[t]],
                  class = paste("btn-sm",
                    if (t == ty) "btn-primary" else "btn-outline-secondary")))
            ),
            div(style = "max-height:220px; overflow-y:auto; margin-bottom:4px;",
                do.call(tagList, col_rows))
          )
        })

        tabPanel("Outcome Details",
          div(style = "padding-top:8px;",
            do.call(tagList, outcome_blocks),
            actionButton("imp_add_outcome", "Add outcome",
                         icon = icon("plus"), class = "btn-sm btn-outline-primary mt-1")
          )
        )
      } else if (sec == "rob") {
        rob_labels <- c(rob2 = "RoB 2", rob2_crt = "RoB 2 (CRT)", rob2_xo = "RoB 2 (XO)",
                        robins_i = "ROBINS-I", robins_ii = "ROBINS-II", robins_e = "ROBINS-E",
                        newcastle_ottawa = "Newcastle-Ottawa", quadas2 = "QUADAS-2")
        expected <- c("Author", "Year", ROB_DOMAIN_COLS[[rob_tool_v]])
        rows <- lapply(expected, function(col) {
          matched <- best_col_match(col, names(rv$upload_df))
          div(class = "map-row",
            div(class = "map-label", col),
            div(style = "flex:1;",
              selectInput(paste0("map_rob_", col), NULL,
                          choices = file_cols, selected = matched, width = "100%")))
        })
        tabPanel("Risk of Bias",
          div(style = "padding-top:8px;",
            div(class = "btn-bar mb-2",
              span(class = "label", "Tool:"),
              lapply(names(rob_labels), function(tid)
                actionButton(paste0("imp_rob_tool_", tid), rob_labels[[tid]],
                  class = paste("btn-sm",
                    if (tid == rob_tool_v) "btn-primary" else "btn-outline-secondary")))
            ),
            div(style = "max-height:340px; overflow-y:auto;", do.call(tagList, rows))
          )
        )
      } else {
        expected <- switch(sec,
          study   = c("Author", "Year", "Study_Group",
                      setdiff(BASE_STUDY_COLS, c("Author", "Year"))),
          inspect = ALL_INSPECT_COLS
        )
        rows <- lapply(expected, function(col) {
          matched <- best_col_match(col, names(rv$upload_df))
          div(class = "map-row",
            div(class = "map-label", col),
            div(style = "flex:1;",
              selectInput(paste0("map_", sec, "_", col), NULL,
                          choices = file_cols, selected = matched, width = "100%")))
        })
        tabPanel(
          title = switch(sec, study = "Study Details", inspect = "INSPECT-SR"),
          div(style = "max-height:400px; overflow-y:auto; padding-top:8px;",
              do.call(tagList, rows))
        )
      }
    })

    tagList(
      rda_selector,
      do.call(tabsetPanel, c(list(id = "mapping_tabs", selected = imp_tab_selected()), tabs))
    )
  })

  observeEvent(input$mapping_tabs,
    imp_tab_selected(input$mapping_tabs), ignoreNULL = TRUE, ignoreInit = TRUE)

  observeEvent(input$import_rda_df, {
    sel <- input$import_rda_df
    if (!is.null(sel) && sel %in% names(rv$upload_rda_dfs))
      rv$upload_df <- rv$upload_rda_dfs[[sel]]
  }, ignoreInit = TRUE)

  observeEvent(input$confirm_import, {
    req(rv$upload_df, input$import_sections)
    tryCatch({
      df <- rv$upload_df
      df[] <- lapply(df, as.character)
      n    <- nrow(df)

      for (sec in input$import_sections) {
        if (sec == "outcome") {
          for (o in imp_outcomes()) {
            ii        <- o$idx
            type_val  <- o$type
            label_val <- trimws(input[[paste0("imp_label_", ii)]] %||% "")
            type_cols <- OUTCOME_TYPE_COLS[[type_val]]
            new_cols  <- outcome_prefix_cols(label_val, type_val)

            auth_src <- best_col_match("Author", names(df))
            year_src <- best_col_match("Year",   names(df))
            imp <- data.frame(
              Author = if (auth_src != "(Skip)") df[[auth_src]] else rep("", n),
              Year   = if (year_src != "(Skip)") df[[year_src]] else rep("", n),
              stringsAsFactors = FALSE
            )
            if (nzchar(label_val)) {
              imp[[paste0(label_val, "_Name")]] <- label_val
              imp[[paste0(label_val, "_Type")]] <- type_val
            }
            for (col in type_cols) {
              col_name <- if (nzchar(label_val)) paste0(label_val, "_", col) else col
              src      <- input[[paste0("map_outcome_", ii, "_", col)]]
              imp[[col_name]] <-
                if (!is.null(src) && src != "(Skip)" && src %in% names(df)) df[[src]]
                else rep("", n)
            }

            existing <- rv$outcome
            nr_exist <- nrow(existing)
            for (nc in new_cols)
              if (!nc %in% names(existing)) existing[[nc]] <- rep("", nr_exist)
            ex_keys  <- paste(existing$Author, existing$Year, sep = "|||")
            imp_keys <- paste(imp$Author, imp$Year, sep = "|||")
            for (i in seq_len(nrow(imp))) {
              k   <- imp_keys[i]
              idx <- which(ex_keys == k)
              if (length(idx)) {
                for (nc in new_cols) existing[idx[1], nc] <- imp[i, nc]
              } else {
                nr <- blank_row(names(existing))
                nr$Author <- imp$Author[i]; nr$Year <- imp$Year[i]
                for (nc in new_cols) nr[[nc]] <- imp[i, nc]
                existing <- bind_rows(existing, nr)
                ex_keys  <- c(ex_keys, k)
              }
            }
            rv$outcome <- existing
            existing_labels <- vapply(rv$outcome_defs, `[[`, character(1), "label")
            if (!label_val %in% existing_labels)
              rv$outcome_defs <- c(rv$outcome_defs,
                list(list(label = label_val, type = type_val, cols = new_cols)))
          }
        } else {
          expected <- switch(sec,
            study   = BASE_STUDY_COLS,
            rob     = c("Author", "Year", ROB_DOMAIN_COLS[[imp_rob_tool()]]),
            inspect = ALL_INSPECT_COLS
          )
          new_df <- as.data.frame(
            matrix("", nrow = n, ncol = length(expected)),
            stringsAsFactors = FALSE
          )
          names(new_df) <- expected
          for (col in expected) {
            src <- input[[paste0("map_", sec, "_", col)]]
            if (!is.null(src) && src != "(Skip)" && src %in% names(df))
              new_df[[col]] <- df[[src]]
          }
          switch(sec,
            study = {
              rv$study <- new_df
              sg_src <- input[["map_study_Study_Group"]]
              if (!is.null(sg_src) && sg_src != "(Skip)") {
                updateCheckboxInput(session, "incl_multiarm", value = TRUE)
                if (!"Study_Group" %in% names(rv$outcome) && nrow(rv$outcome) > 0) {
                  idx_yr    <- which(names(rv$outcome) == "Year")
                  left      <- rv$outcome[, seq_len(idx_yr), drop = FALSE]
                  right     <- rv$outcome[, seq(idx_yr + 1, ncol(rv$outcome)), drop = FALSE]
                  sg_vals   <- rep("", nrow(rv$outcome))
                  out_keys  <- paste(rv$outcome$Author, rv$outcome$Year, sep = "|||")
                  for (i in seq_len(nrow(new_df))) {
                    k   <- paste(new_df$Author[i], new_df$Year[i], sep = "|||")
                    pos <- which(out_keys == k)
                    if (length(pos)) sg_vals[pos[1]] <- new_df$Study_Group[i]
                  }
                  rv$outcome <- cbind(
                    left,
                    data.frame(Study_Group = sg_vals, stringsAsFactors = FALSE),
                    right
                  )
                }
              }
            },
            rob     = { rv$rob     <- new_df },
            inspect = { rv$inspect <- new_df }
          )
        }
      }
      removeModal()
      showNotification(
        paste0("Imported: ", paste(input$import_sections, collapse = ", "), "."),
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("Import failed:", conditionMessage(e)),
                       type = "error", duration = 12)
    })
  })

  # ── Study Details ─────────────────────────────────────────────────────────

  observeEvent(input$incl_multiarm, {
    has_col <- "Study_Group" %in% names(rv$study)
    if (isTRUE(input$incl_multiarm) && !has_col) {
      n <- nrow(rv$study)
      idx <- which(names(rv$study) == "Year")
      left  <- rv$study[, seq_len(idx),          drop = FALSE]
      right <- rv$study[, seq(idx + 1, ncol(rv$study)), drop = FALSE]
      rv$study <- cbind(left,
                        data.frame(Study_Group = rep("", n), stringsAsFactors = FALSE),
                        right)
    } else if (isFALSE(input$incl_multiarm) && has_col) {
      rv$study <- rv$study[, names(rv$study) != "Study_Group", drop = FALSE]
    }
  }, ignoreInit = TRUE)

  output$study_tbl <- renderDT(dt_base(rv$study, "study_tbl"))

  observeEvent(input$study_tbl_cell_edit, {
    e <- input$study_tbl_cell_edit
    rv$study[e$row, e$col + 1] <- DT::coerceValue(e$value, rv$study[[e$col + 1]][e$row])
  }, priority = 100)

  observeEvent(input$add_study, {
    sel <- input$study_tbl_rows_selected
    if (isTRUE(input$incl_multiarm) && length(sel) > 0) {
      ref <- rv$study[sel[1], , drop = FALSE]
      pre <- list(Author = ref$Author, Year = ref$Year)
      if ("Study_Group" %in% names(ref)) pre$Study_Group <- ref$Study_Group
      rv$study <- bind_rows(rv$study, blank_row(names(rv$study), pre))
    } else {
      rv$study <- bind_rows(rv$study, blank_row(names(rv$study)))
    }
  })

  observeEvent(input$del_study, {
    sel <- input$study_tbl_rows_selected
    if (length(sel)) rv$study <- rv$study[-sel, , drop = FALSE]
  })

  observeEvent(input$add_col_study, {
    showModal(modalDialog(
      title     = "Add column to Study Details",
      textInput("new_col_name", "Column name", placeholder = "e.g. Country"),
      footer    = tagList(modalButton("Cancel"),
                          actionButton("confirm_add_col", "Add", class = "btn-primary")),
      easyClose = TRUE
    ))
  })

  observeEvent(input$confirm_add_col, {
    col_name <- trimws(input$new_col_name)
    if (nzchar(col_name) && !col_name %in% names(rv$study))
      rv$study[[col_name]] <- rep("", nrow(rv$study))
    removeModal()
  })

  # ── Outcome Details ───────────────────────────────────────────────────────

  observeEvent(input$otype_all,  outcome_type_filter("all"))
  observeEvent(input$otype_bin,  outcome_type_filter("binary"))
  observeEvent(input$otype_con,  outcome_type_filter("continuous"))
  observeEvent(input$otype_cnt,  outcome_type_filter("count"))
  observeEvent(input$otype_tte,  outcome_type_filter("tte"))
  observeEvent(input$otype_prop, outcome_type_filter("proportion"))

  outcome_display <- reactive({
    df   <- rv$outcome
    filt <- outcome_type_filter()
    if (filt == "all" || ncol(df) <= 2 || length(rv$outcome_defs) == 0) return(df)
    match_cols <- unlist(lapply(rv$outcome_defs, function(d) {
      if (d$type == filt) d$cols else character(0)
    }))
    keep <- c("Author", "Year", intersect(match_cols, names(df)))
    df[, keep, drop = FALSE]
  })

  output$outcome_tbl <- renderDT(dt_base(outcome_display(), "outcome_tbl"))

  observeEvent(input$outcome_tbl_cell_edit, {
    e        <- input$outcome_tbl_cell_edit
    col_name <- names(outcome_display())[e$col + 1]
    rv$outcome[e$row, col_name] <-
      DT::coerceValue(e$value, rv$outcome[[col_name]][e$row])
  }, priority = 100)

  observeEvent(input$del_outcome_row, {
    sel <- input$outcome_tbl_rows_selected
    if (length(sel)) rv$outcome <- rv$outcome[-sel, , drop = FALSE]
  })

  observeEvent(input$del_outcome_group, {
    if (length(rv$outcome_defs) == 0) {
      showNotification("No defined outcomes to remove.", type = "warning")
      return()
    }
    choices <- setNames(
      seq_along(rv$outcome_defs),
      sapply(rv$outcome_defs, function(d)
        paste0(d$label, " (", OUTCOME_TYPE_LABELS[[d$type]], ")"))
    )
    showModal(modalDialog(
      title = "Remove outcome",
      selectInput("del_outcome_choice", "Select outcome to remove:",
                  choices = choices, size = min(6, length(choices)),
                  selectize = FALSE),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_del_outcome", "Remove", class = "btn-danger")
      ),
      easyClose = TRUE
    ))
  })

  observeEvent(input$confirm_del_outcome, {
    idx <- as.integer(input$del_outcome_choice)
    if (!length(idx) || is.na(idx)) { removeModal(); return() }
    def      <- rv$outcome_defs[[idx]]
    drop_cols <- intersect(def$cols, names(rv$outcome))
    rv$outcome      <- rv$outcome[, setdiff(names(rv$outcome), drop_cols), drop = FALSE]
    rv$outcome_defs <- rv$outcome_defs[-idx]
    removeModal()
    showNotification(paste0("Removed outcome: ", def$label), type = "message")
  })

  observeEvent(input$sync_outcome, {
    out <- sync_outcome_rows(rv$outcome, rv$study)
    if (nrow(out) > 0) {
      for (def in rv$outcome_defs) {
        if (nzchar(def$label)) {
          name_col <- paste0(def$label, "_Name")
          type_col <- paste0(def$label, "_Type")
          if (name_col %in% names(out)) out[[name_col]] <- def$label
          if (type_col %in% names(out)) out[[type_col]] <- def$type
        }
      }
    }
    rv$outcome <- out
  })

  observeEvent(input$add_outcome_btn, {
    new_outcome_type("binary")
    showModal(modalDialog(
      title = "Add outcome",
      textInput("new_outcome_label",
                "Outcome label (column prefix — only needed for multiple outcomes):",
                value = "",
                placeholder = "e.g. Primary, HAM-D — leave blank if only one outcome"),
      div(class = "mb-2",
        tags$label("Outcome type:"),
        uiOutput("new_outcome_type_btns")
      ),
      uiOutput("new_outcome_cols_preview"),
      footer    = tagList(modalButton("Cancel"),
                          actionButton("confirm_add_outcome", "Add outcome",
                                       class = "btn-primary")),
      easyClose = TRUE
    ))
  })

  output$new_outcome_type_btns <- renderUI({
    sel <- new_outcome_type()
    div(class = "btn-bar mt-1",
      lapply(names(OUTCOME_TYPE_LABELS), function(t) {
        actionButton(paste0("new_otype_", t), OUTCOME_TYPE_LABELS[[t]],
          class = paste("btn-sm",
            if (t == sel) "btn-primary" else "btn-outline-secondary"))
      })
    )
  })

  lapply(names(OUTCOME_TYPE_LABELS), function(t) {
    observeEvent(input[[paste0("new_otype_", t)]], new_outcome_type(t), ignoreInit = TRUE)
  })

  output$new_outcome_cols_preview <- renderUI({
    cols <- OUTCOME_TYPE_COLS[[new_outcome_type()]]
    p(class = "text-muted small mt-2",
      "Columns to be added: ", tags$em(paste(cols, collapse = ", ")))
  })

  observeEvent(input$confirm_add_outcome, {
    label_val <- trimws(input$new_outcome_label %||% "")
    type_val  <- new_outcome_type()
    new_cols  <- outcome_prefix_cols(label_val, type_val)

    existing_labels <- sapply(rv$outcome_defs, function(d) d$label)
    if (label_val %in% existing_labels) {
      showNotification(
        if (nzchar(label_val))
          paste0("'", label_val, "' already exists.")
        else
          "An unlabelled outcome already exists. Use a prefix for additional outcomes.",
        type = "warning")
      return()
    }
    existing <- rv$outcome
    nr_exist <- nrow(existing)
    for (nc in new_cols)
      if (!nc %in% names(existing)) existing[[nc]] <- rep("", nr_exist)
    if (nzchar(label_val)) {
      existing[[paste0(label_val, "_Name")]] <- rep(label_val, nr_exist)
      existing[[paste0(label_val, "_Type")]] <- rep(type_val,  nr_exist)
    }
    rv$outcome      <- existing
    rv$outcome_defs <- c(rv$outcome_defs,
      list(list(label = label_val, type = type_val, cols = new_cols)))
    removeModal()
  })

  # ── Risk of Bias ─────────────────────────────────────────────────────────

  observeEvent(input$rob_tool, {
    cols  <- c("Author", "Year", ROB_DOMAIN_COLS[[input$rob_tool]])
    id_df <- rv$rob[, intersect(names(rv$rob), c("Author", "Year")), drop = FALSE]
    n     <- nrow(id_df)
    xcols <- setdiff(cols, c("Author", "Year"))
    xdf   <- if (n > 0)
      as.data.frame(setNames(replicate(length(xcols), rep("", n), simplify = FALSE), xcols),
                   stringsAsFactors = FALSE)
    else empty_df(xcols)
    rv$rob     <- if (n > 0) rbind(empty_df(cols), cbind(id_df, xdf)) else empty_df(cols)
    rv$sel_rob <- NULL
  }, ignoreInit = TRUE)

  output$rob_btn_bar <- renderUI({
    tool <- req(input$rob_tool)
    if (tool == "newcastle_ottawa") {
      div(class = "btn-bar",
        span(class = "label", "Enter numeric scores directly in the table."),
        span(class = "text-muted small", "(Selection 0-4, Comparability 0-2, Outcome 0-3)"))
    } else {
      rob2_btns <- list(
        list(id = "btn_rob_low",  label = "Low",           cls = "btn-rob-low"),
        list(id = "btn_rob_some", label = "Some concerns", cls = "btn-rob-some"),
        list(id = "btn_rob_high", label = "High",          cls = "btn-rob-high")
      )
      robins_btns <- list(
        list(id = "btn_rob_low",      label = "Low",            cls = "btn-rob-low"),
        list(id = "btn_rob_moderate", label = "Moderate",       cls = "btn-rob-moderate"),
        list(id = "btn_rob_serious",  label = "Serious",        cls = "btn-rob-serious"),
        list(id = "btn_rob_critical", label = "Critical",       cls = "btn-rob-critical"),
        list(id = "btn_rob_ni",       label = "No information", cls = "btn-rob-ni")
      )
      btns <- switch(tool,
        rob2     = rob2_btns,
        rob2_crt = rob2_btns,
        rob2_xo  = rob2_btns,
        robins_i  = robins_btns,
        robins_ii = robins_btns,
        robins_e  = robins_btns,
        quadas2 = list(
          list(id = "btn_rob_low",     label = "Low",     cls = "btn-rob-low"),
          list(id = "btn_rob_high",    label = "High",    cls = "btn-rob-high"),
          list(id = "btn_rob_unclear", label = "Unclear", cls = "btn-rob-unclear")
        )
      )
      tagList(
        div(class = "btn-bar",
          span(class = "label", "Insert:"),
          lapply(btns, function(b)
            actionButton(b$id, b$label, class = paste("btn-sm", b$cls))),
          actionButton("btn_rob_clear", "Clear", class = "btn-sm btn-outline-secondary")
        ),
        p(class = "text-muted small mb-1",
          "Click a domain cell in the table, then click a button above to fill it.")
      )
    }
  })

  output$rob_sel_info <- renderUI({
    sel <- rv$sel_rob
    if (is.null(sel)) return(NULL)
    lbl <- ROB_DOMAIN_LABELS[[input$rob_tool]][[sel$col_name]] %||% sel$col_name
    div(class = "selection-info",
        "Selected: ", tags$strong(sel$study), " — ", lbl)
  })

  output$rob_tbl <- renderDT(dt_base(rv$rob, "rob_tbl"))

  observeEvent(input$rob_tbl_cell_edit, {
    e <- input$rob_tbl_cell_edit
    rv$rob[e$row, e$col + 1] <- DT::coerceValue(e$value, rv$rob[[e$col + 1]][e$row])
  }, priority = 100)

  observeEvent(input$rob_tbl_cell_clicked, {
    ci <- input$rob_tbl_cell_clicked
    if (is.null(ci$row) || length(ci$row) == 0) return()
    col_idx <- ci$col + 1L
    if (col_idx <= 2L) return()
    rv$sel_rob <- list(
      row      = ci$row,
      col_name = names(rv$rob)[col_idx],
      study    = study_label(rv$rob$Author[ci$row], rv$rob$Year[ci$row])
    )
  })

  fill_rob <- function(val) {
    sel <- rv$sel_rob
    if (!is.null(sel) && sel$row <= nrow(rv$rob))
      rv$rob[sel$row, sel$col_name] <- val
  }

  observeEvent(input$btn_rob_low,      fill_rob("Low"))
  observeEvent(input$btn_rob_some,     fill_rob("Some concerns"))
  observeEvent(input$btn_rob_high,     fill_rob("High"))
  observeEvent(input$btn_rob_moderate, fill_rob("Moderate"))
  observeEvent(input$btn_rob_serious,  fill_rob("Serious"))
  observeEvent(input$btn_rob_critical, fill_rob("Critical"))
  observeEvent(input$btn_rob_unclear,  fill_rob("Unclear"))
  observeEvent(input$btn_rob_ni,       fill_rob("No information"))
  observeEvent(input$btn_rob_clear,    fill_rob(""))

  observeEvent(input$sync_rob, {
    rv$rob <- sync_rows(rv$rob, rv$study, c("Author", "Year"),
                        ROB_DOMAIN_COLS[[input$rob_tool]])
  })

  # ── INSPECT-SR ────────────────────────────────────────────────────────────

  output$inspect_sel_info <- renderUI({
    sel <- rv$sel_inspect
    if (is.null(sel)) return(NULL)
    lbl <- ALL_INSPECT_LABELS[[sel$col_name]] %||% sel$col_name
    div(class = "selection-info",
        "Selected: ", tags$strong(sel$study), " — ", lbl)
  })

  output$inspect_tbl <- renderDT({
    df   <- rv$inspect
    nmap <- ALL_INSPECT_LABELS[intersect(names(ALL_INSPECT_LABELS), names(df))]
    names(df)[names(df) %in% names(nmap)] <- nmap[names(df)[names(df) %in% names(nmap)]]
    dt_base(df, "inspect_tbl")
  })

  observeEvent(input$inspect_tbl_cell_edit, {
    e <- input$inspect_tbl_cell_edit
    rv$inspect[e$row, e$col + 1] <-
      DT::coerceValue(e$value, rv$inspect[[e$col + 1]][e$row])
  }, priority = 100)

  observeEvent(input$inspect_tbl_cell_clicked, {
    ci <- input$inspect_tbl_cell_clicked
    if (is.null(ci$row) || length(ci$row) == 0) return()
    col_idx <- ci$col + 1L
    if (col_idx <= 2L) return()
    rv$sel_inspect <- list(
      row      = ci$row,
      col_name = names(rv$inspect)[col_idx],
      study    = study_label(rv$inspect$Author[ci$row], rv$inspect$Year[ci$row])
    )
  })

  fill_inspect <- function(val) {
    sel <- rv$sel_inspect
    if (!is.null(sel) && sel$row <= nrow(rv$inspect))
      rv$inspect[sel$row, sel$col_name] <- val
  }

  observeEvent(input$btn_isp_none,    fill_inspect("No concerns"))
  observeEvent(input$btn_isp_some,    fill_inspect("Some concerns"))
  observeEvent(input$btn_isp_serious, fill_inspect("Serious concerns"))
  observeEvent(input$btn_isp_clear,   fill_inspect(""))

  observeEvent(input$sync_inspect, {
    rv$inspect <- sync_inspect_rows(rv$inspect, rv$study)
  })

  # ── Downloads ─────────────────────────────────────────────────────────────

  dl_ext <- function() switch(input$dl_format, rds = ".rds", csv = ".csv", xlsx = ".xlsx")

  output$dl_buttons <- renderUI({
    individual <- div(class = "dl-row mb-2",
      downloadButton("dl_study",   "Study Details",  class = "btn-sm btn-outline-primary"),
      downloadButton("dl_outcome", "Outcome Data",   class = "btn-sm btn-outline-primary"),
      downloadButton("dl_rob",     "Risk of Bias",   class = "btn-sm btn-outline-primary"),
      downloadButton("dl_inspect", "INSPECT-SR",     class = "btn-sm btn-outline-primary")
    )
    tagList(
      individual,
      downloadButton("dl_all", "Download all sections",
                     class = "btn-sm btn-primary")
    )
  })

  output$dl_hint <- renderUI({
    msg <- switch(input$dl_format,
      rds  = "One .rds per section — reload with readRDS().",
      xlsx = "One .xlsx per section.",
      csv  = "One .csv per section."
    )
    p(class = "text-muted small mt-1 mb-0", msg)
  })

  write_section <- function(df, file) {
    switch(input$dl_format,
      rds  = saveRDS(df, file),
      csv  = write.csv(df, file, row.names = FALSE),
      xlsx = writexl::write_xlsx(df, file)
    )
  }

  # ── Manage columns / Clear section / Row delete ───────────────────────────

  PROTECTED_COLS <- c("Author", "Year", "Study_Group")

  open_manage_modal <- function(section) {
    df    <- rv[[section]]
    types <- rv$types[[section]] %||% list()
    cols  <- names(df)
    if (length(cols) == 0) {
      showNotification("No columns to manage in this section.", type = "warning")
      return()
    }
    rows <- lapply(seq_along(cols), function(i) {
      col <- cols[i]
      protected <- col %in% PROTECTED_COLS
      div(class = "row g-2 align-items-center mb-1",
        div(class = "col-6",
          textInput(paste0("mc_name_", section, "_", i),
                    label = if (i == 1) "Column name" else NULL,
                    value = col, width = "100%")
        ),
        div(class = "col-5",
          selectInput(paste0("mc_type_", section, "_", i),
                      label = if (i == 1) "Type" else NULL,
                      choices = COL_TYPE_CHOICES,
                      selected = types[[col]] %||% "text",
                      width = "100%")
        ),
        div(class = "col-1 text-end pt-3",
          if (protected) tags$small("locked", class = "text-muted")
          else checkboxInput(paste0("mc_drop_", section, "_", i),
                             label = "drop", value = FALSE, width = "auto")
        )
      )
    })
    showModal(modalDialog(
      title = paste("Manage columns —", tools::toTitleCase(section)),
      size  = "l",
      easyClose = TRUE,
      tagList(
        p(class = "small text-muted",
          "Rename, set type, or drop columns. ",
          strong(paste(PROTECTED_COLS, collapse = ", ")),
          " cannot be dropped."),
        do.call(tagList, rows),
        tags$input(type = "hidden", id = paste0("mc_section"), value = section)
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("apply_manage_cols", "Apply", class = "btn-primary")
      )
    ))
    session$userData$mc_section <- section
  }

  observeEvent(input$manage_study,   open_manage_modal("study"))
  observeEvent(input$manage_outcome, open_manage_modal("outcome"))
  observeEvent(input$manage_rob,     open_manage_modal("rob"))
  observeEvent(input$manage_inspect, open_manage_modal("inspect"))

  observeEvent(input$apply_manage_cols, {
    section <- session$userData$mc_section
    if (is.null(section)) { removeModal(); return() }
    df    <- rv[[section]]
    types <- rv$types[[section]] %||% list()
    cols  <- names(df)
    new_names <- cols
    drop_idx  <- integer(0)
    new_types <- list()
    for (i in seq_along(cols)) {
      old <- cols[i]
      nm  <- trimws(input[[paste0("mc_name_", section, "_", i)]] %||% old)
      tp  <- input[[paste0("mc_type_", section, "_", i)]] %||% "text"
      drop <- isTRUE(input[[paste0("mc_drop_", section, "_", i)]])
      if (drop && !old %in% PROTECTED_COLS) {
        drop_idx <- c(drop_idx, i); next
      }
      if (!nzchar(nm) || (old %in% PROTECTED_COLS && nm != old)) nm <- old
      new_names[i] <- nm
      new_types[[nm]] <- tp
    }
    if (length(drop_idx) > 0) {
      df <- df[, -drop_idx, drop = FALSE]
      new_names <- new_names[-drop_idx]
    }
    if (anyDuplicated(new_names)) {
      showNotification("Column names must be unique. Changes not applied.",
                       type = "error")
      return()
    }
    names(df) <- new_names
    df <- apply_col_types(df, new_types)
    rv[[section]]      <- df
    rv$types[[section]] <- new_types
    removeModal()
    showNotification(paste("Updated", section, "columns."), type = "message")
  })

  # ── Clear section (with confirmation) ─────────────────────────────────────

  open_clear_modal <- function(section) {
    showModal(modalDialog(
      title = paste("Clear", tools::toTitleCase(section), "section?"),
      tagList(
        p("This removes ", strong("all rows"), " from this section. ",
          "Other sections are not affected. This cannot be undone."),
        tags$input(type = "hidden", id = "clear_section_value", value = section)
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_clear_section", "Clear", class = "btn-danger")
      ),
      easyClose = TRUE
    ))
    session$userData$clear_section <- section
  }

  observeEvent(input$clear_study,   open_clear_modal("study"))
  observeEvent(input$clear_outcome, open_clear_modal("outcome"))
  observeEvent(input$clear_rob,     open_clear_modal("rob"))
  observeEvent(input$clear_inspect, open_clear_modal("inspect"))

  observeEvent(input$confirm_clear_section, {
    section <- session$userData$clear_section
    if (is.null(section)) { removeModal(); return() }
    cols <- names(rv[[section]])
    rv[[section]] <- empty_df(cols)
    removeModal()
    showNotification(paste("Cleared", section, "section."), type = "message")
  })

  # ── Row delete for RoB / Inspect ──────────────────────────────────────────

  observeEvent(input$del_rob_row, {
    sel <- input$rob_tbl_rows_selected
    if (length(sel)) rv$rob <- rv$rob[-sel, , drop = FALSE]
  })
  observeEvent(input$del_inspect_row, {
    sel <- input$inspect_tbl_rows_selected
    if (length(sel)) rv$inspect <- rv$inspect[-sel, , drop = FALSE]
  })

  # ── Downloads (apply per-section types) ───────────────────────────────────

  output$dl_study <- downloadHandler(
    filename = function() paste0("study_details", dl_ext()),
    content  = function(file)
      write_section(apply_col_types(add_study_col(rv$study), rv$types$study), file)
  )
  output$dl_outcome <- downloadHandler(
    filename = function() paste0("outcome_data", dl_ext()),
    content  = function(file)
      write_section(apply_col_types(add_study_col(rv$outcome), rv$types$outcome), file)
  )
  output$dl_rob <- downloadHandler(
    filename = function() paste0("rob_data", dl_ext()),
    content  = function(file)
      write_section(apply_col_types(add_study_col(rv$rob), rv$types$rob), file)
  )
  output$dl_inspect <- downloadHandler(
    filename = function() paste0("inspect_sr", dl_ext()),
    content  = function(file)
      write_section(apply_col_types(add_study_col(rv$inspect), rv$types$inspect), file)
  )

  output$dl_all <- downloadHandler(
    filename = function() paste0("bayesma_data", dl_ext()),
    content  = function(file) {
      merged <- merge_all_sections(
        apply_col_types(rv$study,   rv$types$study),
        apply_col_types(rv$outcome, rv$types$outcome),
        apply_col_types(rv$rob,     rv$types$rob),
        apply_col_types(rv$inspect, rv$types$inspect)
      )
      write_section(merged, file)
    }
  )

  # ── Templates ─────────────────────────────────────────────────────────────

  tpl_ext <- function() switch(input$tpl_format, csv = ".csv", xlsx = ".xlsx", rds = ".rds")

  write_tpl <- function(df, file) {
    switch(input$tpl_format,
      csv  = write.csv(df, file, row.names = FALSE),
      xlsx = writexl::write_xlsx(df, file),
      rds  = saveRDS(df, file)
    )
  }

  output$tpl_study <- downloadHandler(
    filename = function() paste0("template_study", tpl_ext()),
    content  = function(file) write_tpl(empty_df(BASE_STUDY_COLS), file)
  )

  output$tpl_outcome <- downloadHandler(
    filename = function() paste0("template_outcome", tpl_ext()),
    content  = function(file) {
      all_cols <- unique(unlist(OUTCOME_TYPE_COLS))
      write_tpl(empty_df(c("Author", "Year", "Outcome_Name", "Outcome_Type", all_cols)), file)
    }
  )

  output$tpl_rob <- downloadHandler(
    filename = function() paste0("template_rob_", input$rob_tool, tpl_ext()),
    content  = function(file)
      write_tpl(empty_df(c("Author", "Year", ROB_DOMAIN_COLS[[input$rob_tool]])), file)
  )

  output$tpl_inspect <- downloadHandler(
    filename = function() paste0("template_inspect_sr", tpl_ext()),
    content  = function(file) write_tpl(empty_df(ALL_INSPECT_COLS), file)
  )

  output$tpl_all <- downloadHandler(
    filename = function() "templates_all.xlsx",
    content  = function(file) {
      all_outcome_cols <- unique(unlist(OUTCOME_TYPE_COLS))
      writexl::write_xlsx(
        list(
          study_details = empty_df(BASE_STUDY_COLS),
          outcome_data  = empty_df(c("Author", "Year", "Outcome_Name", "Outcome_Type",
                                     all_outcome_cols)),
          rob_data      = empty_df(c("Author", "Year",
                                     ROB_DOMAIN_COLS[[input$rob_tool %||% "rob2"]])),
          inspect_sr    = empty_df(ALL_INSPECT_COLS)
        ),
        file
      )
    }
  )

  # ── Help modals ───────────────────────────────────────────────────────────

  rob_help_tbl <- function(rows, footer = NULL) {
    tagList(
      tags$table(class = "table table-sm table-bordered",
        tags$thead(tags$tr(
          tags$th("Domain"), tags$th("Key questions"), tags$th("Ratings")
        )),
        tags$tbody(lapply(rows, function(r) {
          tags$tr(
            tags$td(tags$strong(r[[1]])),
            tags$td(r[[2]]),
            tags$td(tags$code(r[[3]]))
          )
        }))
      ),
      if (!is.null(footer)) p(class = "text-muted small mb-0", footer)
    )
  }

  observeEvent(input$rob_help, {
    tool <- input$rob_tool %||% "rob2"
    content <- switch(tool,
      rob2 = tagList(
        p(class = "text-muted small",
          "RoB 2 assesses randomised controlled trials across 5 domains."),
        rob_help_tbl(
          list(
            list("D1: Randomisation process",
                 "Was the allocation sequence random? Was allocation adequately concealed? Were there baseline imbalances suggesting a problem?",
                 "Low / Some concerns / High"),
            list("D2: Deviations from intended interventions",
                 "Were participants or carers aware of their assigned intervention? Did deviations from the intended intervention occur that could affect the outcome?",
                 "Low / Some concerns / High / NA"),
            list("D3: Missing outcome data",
                 "Were data available for all (or nearly all) participants? Could missingness in outcome data depend on the true value of the outcome?",
                 "Low / Some concerns / High"),
            list("D4: Measurement of the outcome",
                 "Was the method of measuring the outcome appropriate? Were outcome assessors blinded to the intervention assignment?",
                 "Low / Some concerns / High"),
            list("D5: Selection of reported results",
                 "Was the trial pre-registered before enrolment? Is the reported analysis consistent with the pre-registered plan?",
                 "Low / Some concerns / High"),
            list("Overall", "Overall risk of bias judgement across all domains.",
                 "Low / Some concerns / High")
          ),
          HTML("Reference: Higgins et al. (2019). A revised tool to assess risk of bias in randomised trials (RoB 2). <em>BMJ</em>, 366.")
        )
      ),
      rob2_crt = tagList(
        p(class = "text-muted small",
          "RoB 2-CRT is an extension of RoB 2 for cluster-randomised trials. It adds a domain (D2) specific to the identification and recruitment of participants within clusters — a unique source of bias when cluster assignment is known before participants are enrolled."),
        rob_help_tbl(
          list(
            list("D1: Randomisation process",
                 "Was the sequence used to randomise clusters random? Was allocation of clusters adequately concealed? Were there baseline imbalances at the cluster or individual level?",
                 "Low / Some concerns / High"),
            list("D2: Identification and recruitment of participants",
                 "Were participants identified and recruited into the trial after clusters were randomised? Could knowledge of cluster assignment have influenced who was recruited within each cluster?",
                 "Low / Some concerns / High / NA"),
            list("D3: Deviations from intended interventions",
                 "Were participants/carers/cluster-level staff aware of their assigned intervention? Did deviations from the intended intervention occur that could affect the outcome?",
                 "Low / Some concerns / High / NA"),
            list("D4: Missing outcome data",
                 "Were data available for all (or nearly all) randomised clusters and participants? Could missingness in outcome data depend on the true value of the outcome?",
                 "Low / Some concerns / High"),
            list("D5: Measurement of the outcome",
                 "Was the method of measuring the outcome appropriate? Were outcome assessors blinded to cluster assignment?",
                 "Low / Some concerns / High"),
            list("D6: Selection of reported results",
                 "Was the trial pre-registered? Is the reported analysis consistent with the pre-registered plan? Was there selection among multiple eligible outcome measurements or analyses?",
                 "Low / Some concerns / High"),
            list("Overall", "Overall risk of bias judgement across all domains.",
                 "Low / Some concerns / High")
          ),
          HTML("Reference: Higgins et al. (2021). Revised Cochrane risk-of-bias tool for cluster-randomised trials (RoB 2-CRT). <em>Methods in Medicine</em>.")
        )
      ),
      rob2_xo = tagList(
        p(class = "text-muted small",
          "RoB 2 for crossover trials uses the same 5-domain structure as standard RoB 2 but adds carry-over and period effect signalling questions within D1. Each period is treated as a separate observation; carry-over from a prior treatment period is the key additional bias concern."),
        rob_help_tbl(
          list(
            list("D1: Randomisation process (incl. carry-over & period effects)",
                 "Was the order of intervention periods randomly assigned? Was allocation adequately concealed? Were there baseline differences between periods suggesting carry-over or period effects? Was there an adequate washout period between interventions?",
                 "Low / Some concerns / High"),
            list("D2: Deviations from intended interventions",
                 "Were participants or care providers aware of their assigned intervention? Did deviations from the intended intervention occur that could affect the outcome?",
                 "Low / Some concerns / High / NA"),
            list("D3: Missing outcome data",
                 "Were outcome data available for all (or nearly all) participants across all periods? Could missingness in outcome data depend on the true value of the outcome?",
                 "Low / Some concerns / High"),
            list("D4: Measurement of the outcome",
                 "Was the outcome measured appropriately and consistently across periods? Were outcome assessors blinded to the intervention assignment?",
                 "Low / Some concerns / High"),
            list("D5: Selection of reported results",
                 "Was the trial pre-registered? Is the reported analysis consistent with the pre-registered plan? Was there selection among multiple eligible analysis approaches?",
                 "Low / Some concerns / High"),
            list("Overall", "Overall risk of bias judgement across all domains.",
                 "Low / Some concerns / High")
          ),
          HTML("Reference: Higgins et al. (2023). RoB 2: A revised tool for assessing risk of bias in randomised trials — crossover version. Cochrane Methods.")
        )
      ),
      robins_i = tagList(
        p(class = "text-muted small",
          "ROBINS-I assesses non-randomised studies of interventions across 7 domains."),
        rob_help_tbl(
          list(
            list("D1: Confounding",
                 "Were there confounding variables associated with both the intervention and the outcome that were not adequately controlled?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D2: Selection of participants",
                 "Was selection of participants into the study (or into the analysis) related to intervention and outcome?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D3: Classification of interventions",
                 "Was the intervention status correctly classified and consistently applied?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D4: Deviations from intended interventions",
                 "Did co-interventions differ between groups? Did participants adhere to their assigned intervention?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D5: Missing data",
                 "Were data missing for a substantial proportion of participants? Could missingness be related to the true value of the outcome?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D6: Measurement of outcomes",
                 "Was outcome measurement blinded to intervention status? Could assessment differ between intervention groups?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D7: Selection of reported results",
                 "Were multiple analyses or outcome measurements reported selectively based on results?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("Overall", "Overall risk of bias judgement.",
                 "Low / Moderate / Serious / Critical / NI")
          ),
          HTML("NI = No information. Reference: Sterne et al. (2016). ROBINS-I. <em>BMJ</em>, 355.")
        )
      ),
      robins_ii = tagList(
        p(class = "text-muted small",
          "ROBINS-II is the revised successor to ROBINS-I for non-randomised studies of interventions. It retains the 7-domain structure but has substantially updated signalling questions, particularly for confounding (D1), which now explicitly distinguishes baseline from time-varying confounding and requires specification of the target trial."),
        rob_help_tbl(
          list(
            list("D1: Confounding",
                 "Was the analysis appropriately adjusted for baseline confounders? Were time-varying confounders accounted for? Is there a clearly specified target trial that the study approximates?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D2: Selection of participants",
                 "Was the study population selected in a way that introduced bias relative to the target trial? Were eligible participants excluded based on factors related to the intervention and outcome?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D3: Classification of interventions",
                 "Were intervention groups correctly and consistently classified? Was measurement of the intervention influenced by knowledge of the outcome (or vice versa)?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D4: Deviations from intended interventions",
                 "Did co-interventions differ between groups? Did non-adherence to the assigned intervention occur in ways that could affect the outcome?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D5: Missing data",
                 "Were outcome or covariate data missing for a substantial proportion of participants? Could missingness be related to the true value of the outcome?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D6: Measurement of outcomes",
                 "Was outcome measurement blinded to intervention status? Could outcome measurement differ between intervention groups?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D7: Selection of reported results",
                 "Were multiple analyses, subgroups, or outcome definitions selectively reported based on results?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("Overall", "Overall risk of bias judgement.",
                 "Low / Moderate / Serious / Critical / NI")
          ),
          HTML("NI = No information. Reference: Sterne et al. (2024). ROBINS-II: revised tool for assessing risk of bias in non-randomised studies of interventions. <em>BMJ</em>.")
        )
      ),
      robins_e = tagList(
        p(class = "text-muted small",
          "ROBINS-E assesses non-randomised studies of exposures (rather than interventions). The domain structure mirrors ROBINS-I/II but uses exposure-specific language. It is appropriate for observational studies estimating the effect of an exposure on an outcome without experimental manipulation."),
        rob_help_tbl(
          list(
            list("D1: Confounding",
                 "Were baseline confounders appropriately measured and controlled? Were time-varying confounders accounted for, including those on the causal pathway?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D2: Selection of participants",
                 "Was selection into the study (or into the analysis) related to both exposure and outcome in ways that could introduce bias?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D3: Classification of exposures",
                 "Was the exposure correctly and consistently classified? Could measurement error in exposure classification have differed between groups?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D4: Departures from intended exposures",
                 "Were there changes to exposure status after baseline that could affect the outcome? Were co-exposures balanced between groups?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D5: Missing data",
                 "Were data missing for a substantial proportion of participants? Could missingness in outcome or covariate data be related to the true value of the outcome?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D6: Measurement of outcomes",
                 "Was the outcome measured appropriately and consistently? Could outcome measurement have differed by exposure group?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("D7: Selection of reported results",
                 "Were multiple analyses or outcome definitions selectively reported? Is the reported result chosen from a larger set of analyses?",
                 "Low / Moderate / Serious / Critical / NI"),
            list("Overall", "Overall risk of bias judgement.",
                 "Low / Moderate / Serious / Critical / NI")
          ),
          HTML("NI = No information. Reference: Higgins et al. (2024). ROBINS-E: a tool for assessing risk of bias in non-randomised studies of exposures. <em>BMJ</em>.")
        )
      ),
      newcastle_ottawa = tagList(
        p(class = "text-muted small",
          "The Newcastle-Ottawa Scale rates cohort and case-control studies. Scores ≥7 are generally considered low risk."),
        tags$table(class = "table table-sm table-bordered",
          tags$thead(tags$tr(
            tags$th("Domain"), tags$th("Items assessed"), tags$th("Max score")
          )),
          tags$tbody(
            tags$tr(tags$td(tags$strong("Selection")),
              tags$td("Representativeness of the exposed cohort; selection of the non-exposed cohort; ascertainment of exposure; demonstration that outcome was not present at start."),
              tags$td(tags$code("0–4"))),
            tags$tr(tags$td(tags$strong("Comparability")),
              tags$td("Comparability of cohorts on the basis of study design or analysis (e.g., controlling for age, sex, or other confounders)."),
              tags$td(tags$code("0–2"))),
            tags$tr(tags$td(tags$strong("Outcome")),
              tags$td("Assessment of outcome; sufficient follow-up length; adequacy of follow-up of cohorts."),
              tags$td(tags$code("0–3"))),
            tags$tr(tags$td(tags$strong("Total")),
              tags$td(HTML("≥7 = low risk of bias; 4–6 = moderate; ≤3 = high.")),
              tags$td(tags$code("0–9")))
          )
        ),
        p(class = "text-muted small mb-0",
          "Reference: Wells et al. (2000). The Newcastle-Ottawa Scale for assessing the quality of non-randomised studies.")
      ),
      quadas2 = tagList(
        p(class = "text-muted small",
          "QUADAS-2 assesses diagnostic test accuracy studies across 4 risk of bias domains and 3 applicability concern domains."),
        tags$table(class = "table table-sm table-bordered",
          tags$thead(tags$tr(
            tags$th("Domain"), tags$th("Key questions"), tags$th("Ratings")
          )),
          tags$tbody(
            tags$tr(tags$td(tags$strong("Risk of Bias"), colspan = 3,
                            style = "background:#f0f4f8; font-size:0.85rem;")),
            tags$tr(tags$td(tags$strong("PD1: Patient selection")),
              tags$td("Was a consecutive or random sample enrolled? Were inappropriate exclusions avoided?"),
              tags$td(tags$code("Low / High / Unclear"))),
            tags$tr(tags$td(tags$strong("PD2: Index test")),
              tags$td("Were index test results interpreted without knowledge of the reference standard?"),
              tags$td(tags$code("Low / High / Unclear"))),
            tags$tr(tags$td(tags$strong("PD3: Reference standard")),
              tags$td("Is the reference standard likely to correctly classify the target condition? Were results interpreted blinded to the index test?"),
              tags$td(tags$code("Low / High / Unclear"))),
            tags$tr(tags$td(tags$strong("PD4: Flow & timing")),
              tags$td("Was the time between index test and reference standard appropriate? Did all patients receive the same reference standard?"),
              tags$td(tags$code("Low / High / Unclear"))),
            tags$tr(tags$td(tags$strong("Applicability concerns"), colspan = 3,
                            style = "background:#f0f4f8; font-size:0.85rem;")),
            tags$tr(tags$td(tags$strong("AC1: Patient selection")),
              tags$td("Are there concerns that included patients do not match the review question?"),
              tags$td(tags$code("Low / High / Unclear"))),
            tags$tr(tags$td(tags$strong("AC2: Index test")),
              tags$td("Are there concerns that the index test, its conduct, or interpretation differ from the review question?"),
              tags$td(tags$code("Low / High / Unclear"))),
            tags$tr(tags$td(tags$strong("AC3: Reference standard")),
              tags$td("Are there concerns that the target condition as defined by the reference standard does not match the review question?"),
              tags$td(tags$code("Low / High / Unclear")))
          )
        ),
        p(class = "text-muted small mb-0",
          HTML("Reference: Whiting et al. (2011). QUADAS-2. <em>Ann Intern Med</em>, 155(8).")
        )
      )
    )
    showModal(modalDialog(
      title = tagList(icon("circle-question"), " Risk of Bias Tool Guide"),
      size  = "l",
      content,
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  observeEvent(input$inspect_help, {
    irow <- function(item, desc, score0, score1, score2) {
      tags$tr(
        tags$td(tags$strong(item)),
        tags$td(style = "font-size:0.85rem;", desc),
        tags$td(tags$code("0"), " ", score0),
        tags$td(tags$code("1"), " ", score1),
        tags$td(tags$code("2"), " ", score2)
      )
    }
    isection <- function(label) {
      tags$tr(
        tags$td(tags$strong(label), colspan = 5,
                style = "background:#f0f4f8; font-size:0.85rem;")
      )
    }
    showModal(modalDialog(
      title = tagList(icon("circle-question"), " INSPECT-SR Scoring Guide"),
      size  = "xl",
      tagList(
        p(class = "text-muted small",
          HTML("INSPECT-SR assesses study integrity across 4 domains. Score each item: <code>0</code> = No concerns, <code>1</code> = Some concerns, <code>2</code> = Serious concerns.")),
        tags$div(style = "overflow-x:auto;",
          tags$table(class = "table table-sm table-bordered",
            tags$thead(tags$tr(
              tags$th("Item"), tags$th("Description"),
              tags$th("0"), tags$th("1"), tags$th("2")
            )),
            tags$tbody(
              isection("Domain 1: Registration"),
              irow("D1.1 Prospective registration",
                   "Was the trial registered before participant enrolment began?",
                   "Registered before enrolment", "Retrospectively registered",
                   "Not registered / unclear"),
              irow("D1.2 Registration consistency",
                   "Is the published report consistent with the trial registration?",
                   "No discrepancies", "Minor discrepancies",
                   "Major discrepancies / unclear"),
              irow("D1.3 SAP consistency",
                   "Is the statistical analysis plan consistent with the registration?",
                   "Consistent SAP", "Minor inconsistencies",
                   "Major inconsistencies / no SAP"),

              isection("Domain 2: Protocol Conduct"),
              irow("D2.1 Ethical approval",
                   "Was documented ethical approval obtained?",
                   "Yes, with evidence", "Unclear",
                   "No evidence of approval"),
              irow("D2.2 Primary registry",
                   "Was the trial registered in a primary registry (e.g., ClinicalTrials.gov, ANZCTR)?",
                   "Yes", "Secondary registry only",
                   "Not registered"),
              irow("D2.3 Registered before enrolment",
                   "Was registration completed before participant enrolment began?",
                   "Yes", "After enrolment started",
                   "After completion / unclear"),
              irow("D2.4 Protocol vs. conduct",
                   "Is the conduct of the trial consistent with the registered protocol?",
                   "Consistent", "Minor inconsistencies",
                   "Major inconsistencies"),
              irow("D2.5 Other conduct concerns",
                   "Are there other concerns about how the trial was conducted?",
                   "No", "Some concerns",
                   "Serious concerns"),

              isection("Domain 3: Reporting"),
              irow("D3.1 All outcomes reported",
                   "Are all pre-specified outcomes reported in the publication?",
                   "All reported", "Minor omissions",
                   "Major omissions"),
              irow("D3.2 Selective reporting",
                   "Is there evidence that outcome reporting was selective based on results?",
                   "No evidence", "Some evidence",
                   "Clear evidence"),

              isection("Domain 4: Data Integrity"),
              irow("D4.1 Flow consistency",
                   "Is the participant flow (CONSORT diagram) consistent throughout the paper?",
                   "Consistent", "Minor inconsistencies",
                   "Major inconsistencies"),
              irow("D4.2 Analysis populations",
                   "Are analysis populations (ITT, per-protocol) appropriately defined and used?",
                   "Appropriately defined", "Some issues",
                   "Not defined / inappropriate"),
              irow("D4.3 Carlisle test",
                   "Are baseline p-values consistent with the expected uniform distribution? (Enter p-values for automated test.)",
                   "Consistent", "Borderline",
                   "Inconsistent / suspicious"),
              irow("D4.4 Data errors",
                   "Are there apparent errors in the reported data (e.g., impossible values, arithmetic errors)?",
                   "No", "Possible errors",
                   "Clear errors"),
              irow("D4.5 Data duplication",
                   "Is there evidence that data from this study appeared in another publication?",
                   "No evidence", "Possible",
                   "Clear duplication"),
              irow("D4.6 N-consistency",
                   "Are sample sizes consistent across tables, figures, and text? (Enter N values for automated check.)",
                   "Consistent", "Minor inconsistencies",
                   "Major inconsistencies"),
              irow("D4.7 Per guidelines",
                   "Were analyses conducted per relevant reporting guidelines (e.g., CONSORT)?",
                   "Yes", "Mostly",
                   "No / unclear"),
              irow("D4.8 GRIM test",
                   "Are reported means consistent with the sample size and number of decimal places? (Enter mean, n, and decimals for automated test.)",
                   "GRIM-consistent", "Borderline",
                   "GRIM-inconsistent"),
              irow("D4.9 p-value verification",
                   "Are reported p-values consistent with the test statistics and degrees of freedom? (Enter statistic, df, and p-value for automated check.)",
                   "Consistent", "Borderline",
                   "Inconsistent"),
              irow("D4.10 Retraction / correction",
                   "Has the paper been retracted, or have corrections been issued?",
                   "No", "Correction issued",
                   "Retracted"),
              irow("D4.11 Other data concerns",
                   "Are there other concerns about the integrity of the reported data?",
                   "No", "Some concerns",
                   "Serious concerns")
            )
          )
        ),
        p(class = "text-muted small mb-0",
          HTML('Reference: Riddle et al. (2025). INSPECT-SR: an instrument for assessing the credibility of randomised trial reports. <em>medRxiv</em>. <a href="https://doi.org/10.1101/2025.09.03.25334905" target="_blank">doi:10.1101/2025.09.03.25334905</a>'))
      ),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
}

shinyApp(ui, server)
