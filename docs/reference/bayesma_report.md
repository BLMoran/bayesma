# Generate a Bayesian meta-analysis report

Renders a parameterised Quarto template into a publication-style
statistical report. Section structure follows the bayesma workflow
guide; prose style follows a typical systematic review (Introduction,
Methods, Results, Overall Interpretation). Most sections are toggleable
via `run_*` arguments.

## Usage

``` r
bayesma_report(
  data,
  studyvar,
  yearvar,
  estimand = c("OR", "RR", "HR", "IRR", "MD", "SMD", "RD", "ARR", "ATE", "ATT", "CATE"),
  cate_covariate = NULL,
  baseline_risk = NULL,
  event_ctrl = NULL,
  event_int = NULL,
  mean_ctrl = NULL,
  mean_int = NULL,
  sd_ctrl = NULL,
  sd_int = NULL,
  n_ctrl = NULL,
  n_int = NULL,
  outcome_name = "the primary outcome",
  null_range,
  title = "Bayesian meta-analysis report",
  subtitle = "Generated with bayesma",
  author = "",
  prospero_id = "",
  background = "",
  hypothesis = "",
  inclusion_criteria = character(),
  exclusion_criteria = character(),
  quality_tool = "Cochrane RoB 2",
  software_note = "",
  rob_data = NULL,
  priors = NULL,
  run_inspect_sr = FALSE,
  run_rob = FALSE,
  run_egger = TRUE,
  run_robma = TRUE,
  run_metareg = FALSE,
  run_subgroup = FALSE,
  run_sensitivity_priors = TRUE,
  moderators = character(),
  subgroup_var = character(),
  primary_model = "auto",
  re_min_k = NULL,
  output_file,
  render = FALSE
)
```

## Arguments

- data:

  A data frame with one row per study (or per arm).

- studyvar:

  Character. Column name with study identifiers.

- yearvar:

  Character. Column name with publication year.

- estimand:

  Character. Target effect: `"OR"`, `"RR"`, `"HR"`, `"IRR"`, `"MD"`,
  `"SMD"`, `"RD"` / `"ARR"`, `"ATE"`, `"ATT"`, or `"CATE"`. The marginal
  estimands (RD/ATE/ATT/CATE) are computed via
  [`bayesma_marginal()`](https://blmoran.github.io/bayesma/reference/bayesma_marginal.md)
  on the natural scale.

- cate_covariate:

  Character. Study-level covariate name. Required when
  `estimand = "CATE"`.

- baseline_risk:

  Numeric in (0, 1) or `"study_mean"`. Reference baseline used to
  back-transform two-stage binomial fits to RD/ATE/ATT.

- event_ctrl, event_int, n_ctrl, n_int:

  Character. Column names for event counts and arm sizes (binary
  measures).

- mean_ctrl, mean_int, sd_ctrl, sd_int:

  Character. Column names for arm means and SDs (continuous measures).

- outcome_name:

  Character. Human-readable outcome label used in titles and prose.

- null_range:

  Numeric length-2 vector. Null/ROPE range on the natural scale (e.g.
  `c(0.9, 1.1)` for an OR).

- title, subtitle, author:

  Character. Report metadata.

- prospero_id:

  Character. Optional PROSPERO registration identifier.

- background, hypothesis:

  Character (scalar or vector). Prose for the Introduction and Methods
  sections. Each element becomes a paragraph.

- inclusion_criteria, exclusion_criteria:

  Character vectors. Each element becomes a bullet point.

- quality_tool:

  Character. Risk-of-bias tool name (default Cochrane RoB 2).

- software_note:

  Character. Optional override for the Software paragraph.

- rob_data:

  Optional pre-scored risk-of-bias data frame for
  [`rob_plot()`](https://blmoran.github.io/bayesma/reference/rob_plot.md).
  Used when `run_rob = TRUE`.

- priors:

  Optional named list of prior settings for the sensitivity panel. Each
  element is a list with `mu_prior` and `tau_prior`. Defaults to `vague`
  / `weak_reg` / `informative`.

- run_inspect_sr, run_rob, run_egger, run_robma, run_metareg,
  run_subgroup, run_sensitivity_priors:

  Logical. Section toggles.

- moderators:

  Character vector. Moderator columns for meta-regression. Required when
  `run_metareg = TRUE`.

- subgroup_var:

  Character vector. Subgroup column(s). Required when
  `run_subgroup = TRUE`.

- primary_model:

  Character. `"auto"` (default; selected via
  [`compare_models()`](https://blmoran.github.io/bayesma/reference/compare_models.md))
  or one of the candidate model labels (e.g. `"two_stage_re"`,
  `"two_stage_re_hksj"`).

- re_min_k:

  Optional numeric. Forwarded to
  [`bayesma()`](https://blmoran.github.io/bayesma/reference/bayesma.md)
  (and used per subgroup level): downgrade random-effects to
  common-effect when the number of studies in scope is below this
  threshold. Default `NULL`.

- output_file:

  Character. Path for the generated `.qmd`.

- render:

  Logical. If `TRUE`, also render the `.qmd` to HTML.

## Value

Invisibly, the path to the generated `.qmd` file.

## Details

The function writes a `.qmd` file plus a companion `.rds` (containing
`data`, `rob_data`, and `priors`) to the directory of `output_file`.
With `render = TRUE` it then renders the `.qmd` to HTML via the Quarto
CLI.

## Examples

``` r
if (FALSE) { # \dontrun{
bayesma_report(
  data         = eeg_data,
  studyvar     = "Author",
  yearvar      = "Year",
  estimand     = "OR",
  event_ctrl   = "c_event", event_int = "i_event",
  n_ctrl       = "c_n",     n_int     = "i_n",
  outcome_name = "Postoperative Delirium",
  null_range   = c(0.9, 1.1),
  title        = "EEG-guided anaesthesia and POD",
  author       = "Dr Benjamin Moran",
  output_file  = "eeg_pod_report.qmd",
  render       = TRUE
)
} # }
```
