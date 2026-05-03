library(shiny)
library(bslib)
library(DT)
library(dplyr)

# ── Constants ─────────────────���────────────────────────────���─────────────────

STUDY_COLS <- c("Author", "Year", "Subgroup", "Control", "Intervention")

OUTCOME_EXTRA <- list(
  binary     = c("N_Control", "N_Intervention",
                 "Event_Control", "Event_Intervention"),
  continuous = c("N_Control", "N_Intervention",
                 "Mean_Control", "SD_Control",
                 "Mean_Intervention", "SD_Intervention"),
  count      = c("N_Control", "N_Intervention",
                 "Count_Control", "Exposure_Control",
                 "Count_Intervention", "Exposure_Intervention")
)

ROB_DOMAIN_COLS <- list(
  rob2             = c("D1", "D2", "D3", "D4", "D5", "Overall"),
  robins_i         = c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "Overall"),
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
  newcastle_ottawa = c(
    Selection = "Selection (0–4)",
    Comparability = "Comparability (0–2)",
    Outcome = "Outcome (0–3)",
    Total = "Total (0–9)"
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

INSPECT_COLS <- c(
  "d1_1", "d1_2", "d1_3",
  "d2_1", "d2_2", "d2_3", "d2_4", "d2_5",
  "d3_1", "d3_2",
  "d4_1", "d4_2", "d4_4", "d4_5", "d4_7", "d4_10", "d4_11"
)

INSPECT_LABELS <- c(
  d1_1 = "D1.1 Registration prospective",
  d1_2 = "D1.2 Registration consistent",
  d1_3 = "D1.3 SAP consistent",
  d2_1 = "D2.1 Ethical approval",
  d2_2 = "D2.2 Primary registry",
  d2_3 = "D2.3 Registered before enrolment",
  d2_4 = "D2.4 Protocol vs conduct",
  d2_5 = "D2.5 Other conduct concerns",
  d3_1 = "D3.1 All outcomes reported",
  d3_2 = "D3.2 Selective reporting",
  d4_1 = "D4.1 Flow consistent",
  d4_2 = "D4.2 Analysis populations",
  d4_4 = "D4.4 Data errors",
  d4_5 = "D4.5 Data duplication",
  d4_7 = "D4.7 Per guidelines",
  d4_10 = "D4.10 Retraction",
  d4_11 = "D4.11 Other concerns"
)

# ── Helpers ────────────────────────────────────────────��─────────────────────

empty_df <- function(cols) {
  structure(
    replicate(length(cols), character(0), simplify = FALSE),
    names = cols,
    class = "data.frame",
    row.names = integer(0)
  )
}

blank_row <- function(cols, vals = list()) {
  row <- setNames(as.list(rep("", length(cols))), cols)
  for (nm in names(vals)) row[[nm]] <- vals[[nm]]
  as.data.frame(row, stringsAsFactors = FALSE)
}

study_key <- function(df) paste(df$Author, df$Year, sep = "|||")

sync_rows <- function(target, study_df, target_id_cols, extra_cols) {
  all_cols <- c(target_id_cols, extra_cols)
  need <- empty_df(all_cols)
  if (nrow(study_df) == 0) return(target)

  existing_keys <- if (nrow(target) > 0)
    paste(target[[target_id_cols[1]]], target[[target_id_cols[2]]], sep = "|||")
  else character(0)

  for (i in seq_len(nrow(study_df))) {
    k <- paste(study_df$Author[i], study_df$Year[i], sep = "|||")
    if (!k %in% existing_keys) {
      need <- bind_rows(need, blank_row(all_cols, list(
        Author = study_df$Author[i],
        Year   = study_df$Year[i]
      )))
    }
  }
  bind_rows(target, need)
}

sync_inspect_rows <- function(target, study_df) {
  all_cols <- c("study", INSPECT_COLS)
  need <- empty_df(all_cols)
  if (nrow(study_df) == 0) return(target)

  existing <- if (nrow(target) > 0) target$study else character(0)

  for (i in seq_len(nrow(study_df))) {
    lbl <- paste0(study_df$Author[i], " (", study_df$Year[i], ")")
    if (!lbl %in% existing) {
      need <- bind_rows(need, blank_row(all_cols, list(study = lbl)))
    }
  }
  bind_rows(target, need)
}

dt_base <- function(df, id, sel_mode = "single", editable = TRUE) {
  datatable(
    df,
    elementId  = id,
    editable   = editable,
    selection  = sel_mode,
    rownames   = FALSE,
    options    = list(
      pageLength    = 50,
      scrollX       = TRUE,
      dom           = "tip",
      autoWidth     = FALSE,
      scrollCollapse = TRUE
    )
  )
}

# ── CSS ──────────────���───────────────────────���───────────────────────────────

app_css <- "
body { font-size: 0.9rem; }
.section-card { margin-bottom: 1.5rem; }
.card-header-custom {
  background: #4F6D8A; color: white;
  padding: 0.6rem 1rem;
  display: flex; align-items: center; justify-content: space-between;
}
.card-header-custom h5 { margin: 0; font-size: 1rem; }
.card-header-custom .form-select, .card-header-custom .form-check-inline {
  font-size: 0.85rem;
}
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
.btn-nos-score    { background:#3191bf; color:white; border:none; }
.btn-inspect-none    { background:#77bb41; color:white; border:none; }
.btn-inspect-none:hover { background:#5e9430; color:white; }
.btn-inspect-some    { background:#f5c518; color:#2E3A4A; border:none; }
.btn-inspect-some:hover { background:#d4a814; color:#2E3A4A; }
.btn-inspect-serious { background:#e32400; color:white; border:none; }
.btn-inspect-serious:hover { background:#c21e00; color:white; }
.selection-info {
  padding: 6px 12px; background:#eaf4fb;
  border-left: 3px solid #3191bf;
  font-size: 0.85rem; margin-bottom: 8px;
  border-radius: 0 4px 4px 0;
}
.btn-bar { display:flex; flex-wrap:wrap; gap:6px; align-items:center; margin-bottom:8px; }
.btn-bar .label { font-weight:600; font-size:0.85rem; margin-right:4px; }
.dl-row { display:flex; flex-wrap:wrap; gap:10px; align-items:flex-start; }
"

# ── UI ────────���──────────────────────────────────���───────────────────────────

ui <- fluidPage(
  theme = bs_theme(bootswatch = "cosmo", primary = "#3191bf"),
  tags$head(tags$style(HTML(app_css))),

  # Page header
  div(
    style = "background:#4F6D8A; color:white; padding:14px 20px; margin-bottom:20px;",
    h3("bayesma — Data Entry", style = "margin:0;"),
    p(style = "margin:4px 0 0; font-size:0.85rem; opacity:0.85;",
      "Enter study data below. Use the Download section to save .rda files.")
  ),

  div(
    style = "max-width:1400px; margin:0 auto; padding:0 16px;",

    # ── 1. Study Details ────��─────────────────────────────────────────────────
    div(
      class = "section-card",
      div(
        class = "card-header-custom",
        h5("1  Study Details"),
        div(
          actionButton("add_study", "Add study", icon = icon("plus"),
                       class = "btn-sm btn-light"),
          actionButton("del_study", "Remove selected", icon = icon("minus"),
                       class = "btn-sm btn-danger ms-2")
        )
      ),
      div(
        class = "card-body border border-top-0 rounded-bottom p-3",
        p(class = "text-muted small mb-2",
          "Author and Year act as study identifiers — all other sections sync from these two columns."),
        DTOutput("study_tbl")
      )
    ),

    # ── 2. Outcome Details ────────���───────────────────────────────────────���──
    div(
      class = "section-card",
      div(
        class = "card-header-custom",
        h5("2  Outcome Details"),
        div(
          style = "display:flex; gap:14px; align-items:center;",
          span("Outcome type:", style = "font-size:0.85rem;"),
          radioButtons(
            "outcome_type", NULL,
            choices  = c("Binary" = "binary",
                         "Continuous" = "continuous",
                         "Count" = "count"),
            selected = "binary",
            inline   = TRUE
          )
        )
      ),
      div(
        class = "card-body border border-top-0 rounded-bottom p-3",
        div(
          class = "btn-bar mb-2",
          actionButton("sync_outcome", "Sync studies",
                       icon = icon("arrows-rotate"),
                       class = "btn-sm btn-outline-secondary")
        ),
        DTOutput("outcome_tbl")
      )
    ),

    # ── 3. Risk of Bias ──────────��───────────────────────────────────────────
    div(
      class = "section-card",
      div(
        class = "card-header-custom",
        h5("3  Risk of Bias"),
        div(
          style = "display:flex; gap:8px; align-items:center;",
          span("Tool:", style = "font-size:0.85rem;"),
          selectInput(
            "rob_tool", NULL,
            choices = c("RoB 2"            = "rob2",
                        "ROBINS-I"         = "robins_i",
                        "Newcastle-Ottawa" = "newcastle_ottawa",
                        "QUADAS-2"         = "quadas2"),
            width = "180px"
          )
        )
      ),
      div(
        class = "card-body border border-top-0 rounded-bottom p-3",
        uiOutput("rob_btn_bar"),
        uiOutput("rob_sel_info"),
        div(
          class = "btn-bar mb-2",
          actionButton("sync_rob", "Sync studies",
                       icon = icon("arrows-rotate"),
                       class = "btn-sm btn-outline-secondary")
        ),
        DTOutput("rob_tbl")
      )
    ),

    # ── 4. INSPECT-SR ──────────────────────────────────────────────��─────────
    div(
      class = "section-card",
      div(class = "card-header-custom", h5("4  INSPECT-SR  (manual judgements)")),
      div(
        class = "card-body border border-top-0 rounded-bottom p-3",
        p(class = "text-muted small mb-2",
          "Items D4.3 (Carlisle), D4.6 (N-consistency), D4.8 (GRIM), and D4.9 (p-value) ",
          "are computed automatically by ",
          code("inspect_sr()"), " — they are not entered here."),
        div(
          class = "btn-bar",
          span(class = "label", "Insert:"),
          actionButton("btn_isp_none",    "No concerns",      class = "btn-sm btn-inspect-none"),
          actionButton("btn_isp_some",    "Some concerns",    class = "btn-sm btn-inspect-some"),
          actionButton("btn_isp_serious", "Serious concerns", class = "btn-sm btn-inspect-serious"),
          actionButton("btn_isp_clear",   "Clear",            class = "btn-sm btn-outline-secondary")
        ),
        uiOutput("inspect_sel_info"),
        div(
          class = "btn-bar mb-2",
          actionButton("sync_inspect", "Sync studies",
                       icon = icon("arrows-rotate"),
                       class = "btn-sm btn-outline-secondary")
        ),
        DTOutput("inspect_tbl")
      )
    ),

    # ── Download ─────────���───────────────────────────────��───────────────────
    div(
      class = "section-card",
      div(class = "card-header-custom", h5("Download")),
      div(
        class = "card-body border border-top-0 rounded-bottom p-3",
        fluidRow(
          column(
            6,
            h6("Individual files"),
            div(
              class = "dl-row",
              downloadButton("dl_study",   "Study Details (.rda)",  class = "btn-sm btn-primary"),
              downloadButton("dl_outcome", "Outcome Data (.rda)",   class = "btn-sm btn-primary"),
              downloadButton("dl_rob",     "Risk of Bias (.rda)",   class = "btn-sm btn-primary"),
              downloadButton("dl_inspect", "INSPECT-SR (.rda)",     class = "btn-sm btn-primary")
            )
          ),
          column(
            6,
            h6("Everything in one file"),
            p(class = "text-muted small",
              "Saves four named objects: ",
              code("study_details"), ", ", code("outcome_data"), ", ",
              code("rob_data"), ", ", code("inspect_sr"), "."),
            downloadButton("dl_all", "All datasets (.rda)", class = "btn-success")
          )
        )
      )
    ),

    tags$br()
  )
)

# ── Server ──────────────────────────────────────────���────────────────────────

server <- function(input, output, session) {

  # ── Reactive stores ──────��────────────────────────────────────���───────────
  rv <- reactiveValues(
    study   = empty_df(STUDY_COLS),
    outcome = empty_df(c("Author", "Year", OUTCOME_EXTRA[["binary"]])),
    rob     = empty_df(c("Author", "Year", ROB_DOMAIN_COLS[["rob2"]])),
    inspect = empty_df(c("study", INSPECT_COLS)),
    sel_rob     = NULL,   # list(row, col_name, label)
    sel_inspect = NULL
  )

  # ── 1. Study Details ──────────────────────────────────────────────────────

  output$study_tbl <- renderDT(dt_base(rv$study, "study_tbl"))

  observeEvent(input$study_tbl_cell_edit, {
    e <- input$study_tbl_cell_edit
    rv$study[e$row, e$col + 1] <- DT::coerceValue(e$value, rv$study[[e$col + 1]][e$row])
  })

  observeEvent(input$add_study, {
    rv$study <- bind_rows(rv$study, blank_row(STUDY_COLS))
  })

  observeEvent(input$del_study, {
    sel <- input$study_tbl_rows_selected
    if (length(sel)) rv$study <- rv$study[-sel, , drop = FALSE]
  })

  # ── 2. Outcome Details ──────────────���─────────────────────────────────────

  observeEvent(input$outcome_type, {
    cols <- c("Author", "Year", OUTCOME_EXTRA[[input$outcome_type]])
    # Rebuild, preserving Author/Year that already exist
    existing <- rv$outcome |> select(any_of(c("Author", "Year")))
    new_df   <- empty_df(cols)
    if (nrow(existing) > 0) {
      new_df <- bind_rows(new_df, bind_cols(
        existing,
        empty_df(setdiff(cols, c("Author", "Year")))[
          rep(1L, nrow(existing)), , drop = FALSE]
      ))
    }
    rv$outcome <- new_df
  }, ignoreInit = TRUE)

  output$outcome_tbl <- renderDT({
    dt_base(rv$outcome, "outcome_tbl")
  })

  observeEvent(input$outcome_tbl_cell_edit, {
    e <- input$outcome_tbl_cell_edit
    rv$outcome[e$row, e$col + 1] <- DT::coerceValue(e$value, rv$outcome[[e$col + 1]][e$row])
  })

  observeEvent(input$sync_outcome, {
    cols <- c("Author", "Year", OUTCOME_EXTRA[[input$outcome_type]])
    rv$outcome <- sync_rows(rv$outcome, rv$study, c("Author", "Year"),
                            OUTCOME_EXTRA[[input$outcome_type]])
  })

  # ── 3. Risk of Bias ────────────────────────────────────────��──────────────

  # Rebuild rob table when tool changes
  observeEvent(input$rob_tool, {
    cols <- c("Author", "Year", ROB_DOMAIN_COLS[[input$rob_tool]])
    existing <- rv$rob |> select(any_of(c("Author", "Year")))
    new_df   <- empty_df(cols)
    if (nrow(existing) > 0) {
      new_df <- bind_rows(new_df, bind_cols(
        existing,
        empty_df(setdiff(cols, c("Author", "Year")))[
          rep(1L, nrow(existing)), , drop = FALSE]
      ))
    }
    rv$rob      <- new_df
    rv$sel_rob  <- NULL
  }, ignoreInit = TRUE)

  # Button bar — changes based on tool
  output$rob_btn_bar <- renderUI({
    tool <- req(input$rob_tool)
    if (tool == "newcastle_ottawa") {
      div(class = "btn-bar",
          span(class = "label", "Enter numeric scores directly in the table."),
          span(class = "text-muted small",
               "(Selection 0–4, Comparability 0–2, Outcome 0–3)"))
    } else {
      btns <- switch(tool,
        rob2 = list(
          list(id = "btn_rob_low",      label = "Low",           cls = "btn-rob-low"),
          list(id = "btn_rob_some",     label = "Some concerns", cls = "btn-rob-some"),
          list(id = "btn_rob_high",     label = "High",          cls = "btn-rob-high")
        ),
        robins_i = list(
          list(id = "btn_rob_low",      label = "Low",      cls = "btn-rob-low"),
          list(id = "btn_rob_moderate", label = "Moderate", cls = "btn-rob-moderate"),
          list(id = "btn_rob_serious",  label = "Serious",  cls = "btn-rob-serious"),
          list(id = "btn_rob_critical", label = "Critical", cls = "btn-rob-critical")
        ),
        quadas2 = list(
          list(id = "btn_rob_low",      label = "Low",     cls = "btn-rob-low"),
          list(id = "btn_rob_high",     label = "High",    cls = "btn-rob-high"),
          list(id = "btn_rob_unclear",  label = "Unclear", cls = "btn-rob-unclear")
        )
      )
      tagList(
        div(class = "btn-bar",
            span(class = "label", "Insert:"),
            lapply(btns, function(b)
              actionButton(b$id, b$label,
                           class = paste("btn-sm", b$cls))),
            actionButton("btn_rob_clear", "Clear",
                         class = "btn-sm btn-outline-secondary")
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
  })

  observeEvent(input$rob_tbl_cell_clicked, {
    ci <- input$rob_tbl_cell_clicked
    if (is.null(ci$row) || length(ci$row) == 0) return()
    col_idx <- ci$col + 1L
    if (col_idx <= 2L) return()           # skip Author / Year
    col_name <- names(rv$rob)[col_idx]
    study_lbl <- paste(rv$rob$Author[ci$row], rv$rob$Year[ci$row])
    rv$sel_rob <- list(row = ci$row, col_name = col_name, study = study_lbl)
  })

  # Helper to fill selected RoB cell
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
  observeEvent(input$btn_rob_clear,    fill_rob(""))

  observeEvent(input$sync_rob, {
    rv$rob <- sync_rows(rv$rob, rv$study, c("Author", "Year"),
                        ROB_DOMAIN_COLS[[input$rob_tool]])
  })

  # ── 4. INSPECT-SR ─────────────��───────────────────────────────────────────

  output$inspect_sel_info <- renderUI({
    sel <- rv$sel_inspect
    if (is.null(sel)) return(NULL)
    lbl <- INSPECT_LABELS[[sel$col_name]] %||% sel$col_name
    div(class = "selection-info",
        "Selected: ", tags$strong(sel$study), " — ", lbl)
  })

  output$inspect_tbl <- renderDT({
    df <- rv$inspect
    # Short column headers in the display
    nmap <- INSPECT_LABELS[intersect(names(INSPECT_LABELS), names(df))]
    names(df)[names(df) %in% names(nmap)] <- nmap[names(df)[names(df) %in% names(nmap)]]
    dt_base(df, "inspect_tbl")
  })

  observeEvent(input$inspect_tbl_cell_edit, {
    e <- input$inspect_tbl_cell_edit
    rv$inspect[e$row, e$col + 1] <- DT::coerceValue(
      e$value, rv$inspect[[e$col + 1]][e$row]
    )
  })

  observeEvent(input$inspect_tbl_cell_clicked, {
    ci <- input$inspect_tbl_cell_clicked
    if (is.null(ci$row) || length(ci$row) == 0) return()
    col_idx <- ci$col + 1L
    if (col_idx <= 1L) return()           # skip study column
    col_name <- names(rv$inspect)[col_idx]
    rv$sel_inspect <- list(row = ci$row, col_name = col_name,
                           study = rv$inspect$study[ci$row])
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

  # ── Downloads ───────────────────��──────────────────────────���──────────────

  save_rda <- function(obj_name, data) {
    function(file) {
      assign(obj_name, data)
      save(list = obj_name, file = file, envir = environment())
    }
  }

  output$dl_study <- downloadHandler(
    filename = "study_details.rda",
    content  = function(file) {
      study_details <- rv$study
      save(study_details, file = file)
    }
  )
  output$dl_outcome <- downloadHandler(
    filename = function() paste0("outcome_data_", input$outcome_type, ".rda"),
    content  = function(file) {
      outcome_data <- rv$outcome
      save(outcome_data, file = file)
    }
  )
  output$dl_rob <- downloadHandler(
    filename = "rob_data.rda",
    content  = function(file) {
      rob_data <- rv$rob
      save(rob_data, file = file)
    }
  )
  output$dl_inspect <- downloadHandler(
    filename = "inspect_sr.rda",
    content  = function(file) {
      inspect_sr <- rv$inspect
      save(inspect_sr, file = file)
    }
  )
  output$dl_all <- downloadHandler(
    filename = "bayesma_data.rda",
    content  = function(file) {
      study_details <- rv$study
      outcome_data  <- rv$outcome
      rob_data      <- rv$rob
      inspect_sr    <- rv$inspect
      save(study_details, outcome_data, rob_data, inspect_sr, file = file)
    }
  )
}

shinyApp(ui, server)
