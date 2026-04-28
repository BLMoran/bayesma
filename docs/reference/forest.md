# Create Bayesian Forest Plot for Meta-Analysis

This function creates a Bayesian forest plot for meta-analysis from a
bayesma model object.

## Usage

``` r
forest(
  model,
  data,
  estimand,
  studyvar = NULL,
  year = NULL,
  c_n = NULL,
  i_n = NULL,
  c_event = NULL,
  i_event = NULL,
  c_mean = NULL,
  i_mean = NULL,
  c_sd = NULL,
  i_sd = NULL,
  c_time = NULL,
  i_time = NULL,
  sort_studies_by = "author",
  subgroup = FALSE,
  subgroup_var = NULL,
  sort_subgroup_by = "alphabetical",
  label_outcome = "Outcome",
  label_control = "Control",
  label_intervention = "Intervention",
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  xlim = NULL,
  x_breaks = NULL,
  add_rope = FALSE,
  rope_value = NULL,
  rope_color = "grey50",
  shrinkage_output = "density",
  null_value = NULL,
  color_palette = NULL,
  color_study_posterior_null_left = "deepskyblue",
  color_study_posterior_null_right = "violet",
  color_study_posterior = "dodgerblue",
  color_study_posterior_outline = "blue",
  color_overall_posterior = "blue",
  color_shrinkage_pointinterval = "purple",
  color_shrinkage_outline = "purple",
  color_shrinkage_fill = NULL,
  split_color_by_null = FALSE,
  color_favours_control = "firebrick",
  color_favours_intervention = "dodgerblue",
  add_arm_labels = TRUE,
  reverse_arm_labels = FALSE,
  add_pred = FALSE,
  add_pred_subgroup = FALSE,
  pred_output = c("density", "pointinterval"),
  color_pred_posterior = "orange",
  color_pred_outline = "darkorange2",
  color_pred_pointinterval = "orange",
  plot_width = 4,
  add_rob = FALSE,
  rob_tool = c("rob2", "robins_i", "quadas2", "robins_e"),
  add_rob_legend = FALSE,
  exclude_high_rob = FALSE,
  re_min_k = NULL,
  incl_shrinkage = TRUE,
  font = NULL
)
```

## Arguments

- model:

  A fitted bayesma object (class 'bayesma').

- data:

  A data frame containing the study data used for the meta-analysis.

- estimand:

  Character string specifying the effect measure or marginal estimand.
  Must be one of:

  - Relative-effect: `"OR"` (Odds Ratio), `"HR"` (Hazard Ratio), `"RR"`
    (Risk Ratio), `"IRR"` (Incidence Rate Ratio), `"MD"` (Mean
    Difference), `"SMD"` (Standardised Mean Difference).

  - Marginal: `"RD"` / `"ARR"` (Risk Difference), `"ATE"` (Average
    Treatment Effect), `"ATT"` (Average Treatment Effect on the
    Treated), `"CATE"` (Conditional Average Treatment Effect).

- studyvar:

  Column name containing study identifiers/authors. Default is NULL.
  Extracted from the model object if not provided.

- year:

  Column name containing publication years. Default is NULL.

- c_n:

  Column name containing control group sample sizes. Required for OR,
  RR, MD, SMD.

- i_n:

  Column name containing intervention group sample sizes. Required for
  OR, RR, MD, SMD.

- c_event:

  Column name containing control group event counts. Required for OR,
  RR, IRR.

- i_event:

  Column name containing intervention group event counts. Required for
  OR, RR, IRR.

- c_mean:

  Column name containing control group means. Required for MD, SMD.

- i_mean:

  Column name containing intervention group means. Required for MD, SMD.

- c_sd:

  Column name containing control group standard deviations. Required for
  MD, SMD.

- i_sd:

  Column name containing intervention group standard deviations.
  Required for MD, SMD.

- c_time:

  Column name containing control group time periods. Required for IRR.

- i_time:

  Column name containing intervention group time periods. Required for
  IRR.

- sort_studies_by:

  Character string specifying how to sort studies. Options: "author"
  (default), "year", or "effect".

- subgroup:

  Logical indicating whether to create subgroup analysis. Default is
  FALSE.

- subgroup_var:

  Character string. Name of the variable in data to use for subgroup
  analysis.

- sort_subgroup_by:

  Character string or vector specifying subgroup ordering. Options:
  "alphabetical" (default), "effect", or custom character vector of
  subgroup names.

- label_outcome:

  Character string for outcome label. Default is "Outcome".

- label_control:

  Character string for control group label. Default is "Control".

- label_intervention:

  Character string for intervention group label. Default is
  "Intervention".

- title:

  Character string for the plot title. Default is NULL (no title).

- subtitle:

  Character string for the plot subtitle. Default is NULL (no subtitle).

- title_align:

  Character string specifying title alignment. Options: "left"
  (default), "center"/"centre", "right".

- xlim:

  Numeric vector of length 2 specifying x-axis limits. Default is NULL
  (auto-scaled).

- x_breaks:

  Numeric vector specifying custom x-axis break points. Default is NULL
  (uses measure-specific defaults).

- add_rope:

  Logical indicating whether to add ROPE (Region of Practical
  Equivalence). Default is FALSE.

- rope_value:

  Numeric vector of length 2 specifying ROPE range, or single value for
  symmetric range around null. Default is NULL (uses Kruschke's
  recommendations: OR/HR/RR/IRR = c(0.9, 1.1), SMD = c(-0.1, 0.1), MD
  requires specification).

- rope_color:

  Color for ROPE shading. Default is transparent grey.

- shrinkage_output:

  Character string specifying shrinkage visualization. Options:
  "density" (default) or "pointinterval".

- null_value:

  Numeric value specifying x-axis value of the null line. Default is
  NULL (uses measure-specific defaults).

- color_palette:

  Character vector of colors for the plot. Default is NULL.

- color_study_posterior_null_left:

  Color for left side of null study posteriors. Default is
  "deepskyblue".

- color_study_posterior_null_right:

  Color for right side of null study posteriors. Default is "violet".

- color_study_posterior:

  Color for study posterior densities. Default is "dodgerblue".

- color_study_posterior_outline:

  Color for study posterior outlines. Default is "blue".

- color_overall_posterior:

  Color for overall posterior. Default is "blue".

- color_shrinkage_pointinterval:

  Color for shrinkage point intervals (used when
  `shrinkage_output = "pointinterval"`). Default is "purple".

- color_shrinkage_outline:

  Color for shrinkage plot outlines. Default is "purple".

- color_shrinkage_fill:

  Color for shrinkage plot fill. Default is NULL.

- split_color_by_null:

  Logical. If TRUE, posterior densities are split and coloured based on
  whether values fall above or below the null value.

- color_favours_control:

  Colour used for density regions favouring the control group when
  `split_color_by_null = TRUE`.

- color_favours_intervention:

  Colour used for density regions favouring the intervention group when
  `split_color_by_null = TRUE`.

- add_arm_labels:

  Logical indicating whether to display "Favours Control" / "Favours
  Intervention" labels above the density plot. Default is TRUE.

- reverse_arm_labels:

  Logical indicating whether to swap the positions of the "Favours
  Control" and "Favours Intervention" labels. Default is FALSE.

- add_pred:

  Logical indicating whether to add a prediction interval row beneath
  the Pooled Effect. Default is FALSE.

- add_pred_subgroup:

  Logical indicating whether to add prediction interval rows for each
  subgroup when `subgroup = TRUE`. If `FALSE` (the default), a
  prediction row is only added for the overall pooled effect. Ignored
  when `subgroup = FALSE`.

- pred_output:

  Character string specifying the visualisation for the prediction
  interval row. Options: "density" (default) or "pointinterval".

- color_pred_posterior:

  Color for the prediction interval density fill. Default is "orange".

- color_pred_outline:

  Color for the prediction interval density outline. Default is
  "darkorange2".

- color_pred_pointinterval:

  Color for the prediction interval point interval. Default is "orange".

- plot_width:

  Numeric value specifying the relative width of the plot component.
  Default is 4.

- add_rob:

  Logical indicating whether to add Risk of Bias assessment. Default is
  FALSE.

- rob_tool:

  Character string specifying RoB tool. Options: "rob2" (default).

- add_rob_legend:

  Logical indicating whether to add RoB legend. Default is FALSE.

- exclude_high_rob:

  Logical indicating whether to exclude high risk of bias studies and
  refit the model. Default is FALSE.

- re_min_k:

  Optional numeric. If the number of studies in a (sub)group is below
  this threshold, the model is downgraded to common-effect for that
  refit. Overrides any `re_min_k` stored in the original model's call
  arguments. Default is NULL (uses the value from the original fit, if
  any).

- incl_shrinkage:

  Logical. If `TRUE` (default), the shrinkage (posterior) effect column
  is shown alongside the observed effect column. Set to `FALSE` to
  display only the observed effect column.

- font:

  Character string specifying the font family to use throughout the
  plot. Default is NULL (uses system defaults).

## Value

A patchwork object containing the complete forest plot with study
information table, density plots, and effect size table.
