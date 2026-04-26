#' Generate a Bayesian meta-analysis report
#'
#' Renders a parameterised Quarto template into a publication-style statistical
#' report. Section structure follows the bayesma workflow guide; prose style
#' follows a typical systematic review (Introduction, Methods, Results, Overall
#' Interpretation). Most sections are toggleable via `run_*` arguments.
#'
#' The function writes a `.qmd` file plus a companion `.rds` (containing
#' `data`, `rob_data`, and `priors`) to the directory of `output_file`. With
#' `render = TRUE` it then renders the `.qmd` to HTML via the Quarto CLI.
#'
#' @param data A data frame with one row per study (or per arm).
#' @param studyvar Character. Column name with study identifiers.
#' @param yearvar Character. Column name with publication year.
#' @param estimand Character. Target effect: `"OR"`, `"RR"`, `"HR"`, `"IRR"`,
#'   `"MD"`, `"SMD"`, `"RD"` / `"ARR"`, `"ATE"`, `"ATT"`, or `"CATE"`. The
#'   marginal estimands (RD/ATE/ATT/CATE) are computed via [bayesma_marginal()]
#'   on the natural scale.
#' @param cate_covariate Character. Study-level covariate name. Required when
#'   `estimand = "CATE"`.
#' @param baseline_risk Numeric in (0, 1) or `"study_mean"`. Reference baseline
#'   used to back-transform two-stage binomial fits to RD/ATE/ATT.
#' @param event_ctrl,event_int,n_ctrl,n_int Character. Column names for event
#'   counts and arm sizes (binary measures).
#' @param mean_ctrl,mean_int,sd_ctrl,sd_int Character. Column names for arm
#'   means and SDs (continuous measures).
#' @param outcome_name Character. Human-readable outcome label used in titles
#'   and prose.
#' @param null_range Numeric length-2 vector. Null/ROPE range on the natural
#'   scale (e.g. `c(0.9, 1.1)` for an OR).
#' @param title,subtitle,author Character. Report metadata.
#' @param prospero_id Character. Optional PROSPERO registration identifier.
#' @param background,hypothesis Character (scalar or vector). Prose for the
#'   Introduction and Methods sections. Each element becomes a paragraph.
#' @param inclusion_criteria,exclusion_criteria Character vectors. Each element
#'   becomes a bullet point.
#' @param quality_tool Character. Risk-of-bias tool name (default Cochrane RoB 2).
#' @param software_note Character. Optional override for the Software paragraph.
#' @param rob_data Optional pre-scored risk-of-bias data frame for [rob_plot()].
#'   Used when `run_rob = TRUE`.
#' @param priors Optional named list of prior settings for the sensitivity
#'   panel. Each element is a list with `mu_prior` and `tau_prior`. Defaults
#'   to `vague` / `weak_reg` / `informative`.
#' @param run_inspect_sr,run_rob,run_egger,run_robma,run_metareg,run_subgroup,run_sensitivity_priors
#'   Logical. Section toggles.
#' @param moderators Character vector. Moderator columns for meta-regression.
#'   Required when `run_metareg = TRUE`.
#' @param subgroup_var Character vector. Subgroup column(s). Required when
#'   `run_subgroup = TRUE`.
#' @param re_min_k Optional numeric. Forwarded to [bayesma()] (and used per
#'   subgroup level): downgrade random-effects to common-effect when the
#'   number of studies in scope is below this threshold. Default `NULL`.
#' @param primary_model Character. `"auto"` (default; selected via
#'   [compare_models()]) or one of the candidate model labels (e.g.
#'   `"two_stage_re"`, `"two_stage_re_hksj"`).
#' @param output_file Character. Path for the generated `.qmd`.
#' @param render Logical. If `TRUE`, also render the `.qmd` to HTML.
#'
#' @return Invisibly, the path to the generated `.qmd` file.
#'
#' @examples
#' \dontrun{
#' bayesma_report(
#'   data         = eeg_data,
#'   studyvar     = "Author",
#'   yearvar      = "Year",
#'   estimand     = "OR",
#'   event_ctrl   = "c_event", event_int = "i_event",
#'   n_ctrl       = "c_n",     n_int     = "i_n",
#'   outcome_name = "Postoperative Delirium",
#'   null_range   = c(0.9, 1.1),
#'   title        = "EEG-guided anaesthesia and POD",
#'   author       = "Dr Benjamin Moran",
#'   output_file  = "eeg_pod_report.qmd",
#'   render       = TRUE
#' )
#' }
#'
#' @export
bayesma_report <- function(
    data,
    studyvar,
    yearvar,
    estimand = c("OR", "RR", "HR", "IRR", "MD", "SMD",
                 "RD", "ARR", "ATE", "ATT", "CATE"),
    cate_covariate = NULL,
    baseline_risk  = NULL,
    event_ctrl = NULL, event_int = NULL,
    mean_ctrl  = NULL, mean_int  = NULL,
    sd_ctrl    = NULL, sd_int    = NULL,
    n_ctrl     = NULL, n_int     = NULL,
    outcome_name = "the primary outcome",
    null_range,
    title    = "Bayesian meta-analysis report",
    subtitle = "Generated with bayesma",
    author   = "",
    prospero_id = "",
    background  = "",
    hypothesis  = "",
    inclusion_criteria = character(),
    exclusion_criteria = character(),
    quality_tool  = "Cochrane RoB 2",
    software_note = "",
    rob_data = NULL,
    priors   = NULL,
    run_inspect_sr         = FALSE,
    run_rob                = FALSE,
    run_egger              = TRUE,
    run_robma              = TRUE,
    run_metareg            = FALSE,
    run_subgroup           = FALSE,
    run_sensitivity_priors = TRUE,
    moderators   = character(),
    subgroup_var = character(),
    primary_model = "auto",
    re_min_k = NULL,
    output_file,
    render = FALSE
) {
  estimand <- rlang::arg_match(estimand)

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }
  if (!is.numeric(null_range) || length(null_range) != 2) {
    cli::cli_abort("{.arg null_range} must be a length-2 numeric vector.")
  }
  if (missing(output_file) || !nzchar(output_file)) {
    cli::cli_abort("{.arg output_file} is required.")
  }
  if (run_metareg && !length(moderators)) {
    cli::cli_abort("{.arg moderators} required when {.code run_metareg = TRUE}.")
  }
  if (run_subgroup && !length(subgroup_var)) {
    cli::cli_abort("{.arg subgroup_var} required when {.code run_subgroup = TRUE}.")
  }
  if (estimand == "CATE" && is.null(cate_covariate)) {
    cli::cli_abort(
      "{.arg cate_covariate} required when {.code estimand = \"CATE\"}."
    )
  }
  if (run_rob && is.null(rob_data) && !any(grepl("^rob_", names(data)))) {
    cli::cli_warn(
      "{.code run_rob = TRUE} but no {.arg rob_data} supplied and no \\
       {.code rob_*} columns in {.arg data}; the RoB section will be empty."
    )
  }

  is_binary <- estimand %in% c("OR", "RR", "HR", "IRR",
                                "RD", "ARR", "ATE", "ATT", "CATE") &&
               !is.null(event_ctrl)
  if (is_binary) {
    if (is.null(event_ctrl) || is.null(event_int) ||
        is.null(n_ctrl) || is.null(n_int)) {
      cli::cli_abort(
        "{.arg event_ctrl}, {.arg event_int}, {.arg n_ctrl}, {.arg n_int} \\
         required for estimand {.val {estimand}}."
      )
    }
  } else {
    if (is.null(mean_ctrl) || is.null(mean_int) ||
        is.null(sd_ctrl) || is.null(sd_int) ||
        is.null(n_ctrl) || is.null(n_int)) {
      cli::cli_abort(
        "{.arg mean_ctrl}, {.arg mean_int}, {.arg sd_ctrl}, {.arg sd_int}, \\
         {.arg n_ctrl}, {.arg n_int} required for estimand {.val {estimand}}."
      )
    }
  }

  template <- system.file("templates", "bayesma_report.qmd",
                          package = "bayesma")
  if (!nzchar(template)) {
    cli::cli_abort("Template not found. Reinstall {.pkg bayesma}.")
  }

  output_file <- normalizePath(output_file, mustWork = FALSE)
  output_dir  <- dirname(output_file)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  bundle_path <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(output_file)), "_bundle.rds")
  )
  saveRDS(
    list(data = data, rob_data = rob_data, priors = priors),
    bundle_path
  )

  null_value <- if (is_binary) 1 else 0

  params <- list(
    title              = title,
    subtitle           = subtitle,
    author             = author,
    prospero_id        = prospero_id,
    outcome_name       = outcome_name,
    estimand           = estimand,
    cate_covariate     = cate_covariate %||% "",
    baseline_risk      = if (is.null(baseline_risk)) NA_real_
                         else if (identical(baseline_risk, "study_mean")) -1
                         else baseline_risk,
    null_value         = null_value,
    null_range         = null_range,
    bundle_path        = bundle_path,
    studyvar           = studyvar,
    yearvar            = yearvar,
    event_ctrl         = event_ctrl %||% "",
    event_int          = event_int  %||% "",
    mean_ctrl          = mean_ctrl  %||% "",
    mean_int           = mean_int   %||% "",
    sd_ctrl            = sd_ctrl    %||% "",
    sd_int             = sd_int     %||% "",
    n_ctrl             = n_ctrl     %||% "",
    n_int              = n_int      %||% "",
    background         = paste(background, collapse = "\n\n"),
    hypothesis         = hypothesis,
    inclusion_criteria = inclusion_criteria,
    exclusion_criteria = exclusion_criteria,
    quality_tool       = quality_tool,
    software_note      = software_note,
    run_inspect_sr     = run_inspect_sr,
    run_rob            = run_rob,
    run_egger          = run_egger,
    run_robma          = run_robma,
    run_metareg        = run_metareg,
    run_subgroup       = run_subgroup,
    run_sensitivity_priors = run_sensitivity_priors,
    moderators         = moderators,
    subgroup_var       = subgroup_var,
    primary_model      = primary_model,
    re_min_k           = re_min_k %||% NA_real_
  )

  inject_params(template, output_file, params)

  cli::cli_alert_success("Wrote report skeleton: {.file {output_file}}")
  cli::cli_alert_info("Companion bundle: {.file {bundle_path}}")

  if (render) {
    if (!requireNamespace("quarto", quietly = TRUE)) {
      cli::cli_abort("{.pkg quarto} package required for {.code render = TRUE}.")
    }
    cli::cli_alert_info("Rendering with Quarto...")
    quarto::quarto_render(output_file, quiet = TRUE)
    cli::cli_alert_success(
      "Rendered: {.file {tools::file_path_sans_ext(output_file)}.html}"
    )
  }

  invisible(output_file)
}

#' @noRd
inject_params <- function(template, output_file, params) {
  lines    <- readLines(template)
  yaml_end <- which(lines == "---")[2]
  body     <- lines[(yaml_end + 1):length(lines)]

  yaml_block <- c(
    "---",
    sprintf("title: %s",    yaml_str(params$title)),
    sprintf("subtitle: %s", yaml_str(params$subtitle)),
    sprintf("author: %s",   yaml_str(params$author)),
    "date: today",
    "format:",
    "  html:",
    "    embed-resources: true",
    "    toc: true",
    "    toc-location: left",
    "    toc-depth: 4",
    "    number-sections: true",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "params:",
    sprintf("  title: %s",            yaml_str(params$title)),
    sprintf("  subtitle: %s",         yaml_str(params$subtitle)),
    sprintf("  author: %s",           yaml_str(params$author)),
    sprintf("  prospero_id: %s",      yaml_str(params$prospero_id)),
    sprintf("  outcome_name: %s",     yaml_str(params$outcome_name)),
    sprintf("  estimand: %s",         yaml_str(params$estimand)),
    sprintf("  cate_covariate: %s",   yaml_str(params$cate_covariate)),
    sprintf("  baseline_risk: %s",
            if (is.na(params$baseline_risk)) "~"
            else if (params$baseline_risk == -1) "\"study_mean\""
            else format(params$baseline_risk)),
    sprintf("  null_value: %s",       params$null_value),
    sprintf("  null_range: !expr c(%s)",
            paste(params$null_range, collapse = ", ")),
    sprintf("  bundle_path: %s",      yaml_str(params$bundle_path)),
    sprintf("  studyvar: %s",         yaml_str(params$studyvar)),
    sprintf("  yearvar: %s",          yaml_str(params$yearvar)),
    sprintf("  event_ctrl: %s",       yaml_str(params$event_ctrl)),
    sprintf("  event_int: %s",        yaml_str(params$event_int)),
    sprintf("  mean_ctrl: %s",        yaml_str(params$mean_ctrl)),
    sprintf("  mean_int: %s",         yaml_str(params$mean_int)),
    sprintf("  sd_ctrl: %s",          yaml_str(params$sd_ctrl)),
    sprintf("  sd_int: %s",           yaml_str(params$sd_int)),
    sprintf("  n_ctrl: %s",           yaml_str(params$n_ctrl)),
    sprintf("  n_int: %s",            yaml_str(params$n_int)),
    sprintf("  background: %s",       yaml_str(params$background)),
    sprintf("  hypothesis: %s",       yaml_str(params$hypothesis)),
    sprintf("  inclusion_criteria: !expr %s",
            yaml_char_vec(params$inclusion_criteria)),
    sprintf("  exclusion_criteria: !expr %s",
            yaml_char_vec(params$exclusion_criteria)),
    sprintf("  quality_tool: %s",     yaml_str(params$quality_tool)),
    sprintf("  software_note: %s",    yaml_str(params$software_note)),
    sprintf("  run_inspect_sr: %s",   tolower(params$run_inspect_sr)),
    sprintf("  run_rob: %s",          tolower(params$run_rob)),
    sprintf("  run_egger: %s",        tolower(params$run_egger)),
    sprintf("  run_robma: %s",        tolower(params$run_robma)),
    sprintf("  run_metareg: %s",      tolower(params$run_metareg)),
    sprintf("  run_subgroup: %s",     tolower(params$run_subgroup)),
    sprintf("  run_sensitivity_priors: %s",
            tolower(params$run_sensitivity_priors)),
    sprintf("  moderators: !expr %s",  yaml_char_vec(params$moderators)),
    sprintf("  subgroup_var: !expr %s",
            yaml_char_vec(params$subgroup_var)),
    sprintf("  primary_model: %s",    yaml_str(params$primary_model)),
    sprintf("  re_min_k: %s",
            if (is.na(params$re_min_k)) "~" else format(params$re_min_k)),
    "---"
  )

  writeLines(c(yaml_block, body), output_file)
  invisible(output_file)
}

#' @noRd
yaml_str <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"",   "\\\\\"",   x)
  paste0("\"", x, "\"")
}

#' @noRd
yaml_char_vec <- function(x) {
  if (!length(x)) return("character()")
  paste0("c(", paste(yaml_str(x), collapse = ", "), ")")
}
