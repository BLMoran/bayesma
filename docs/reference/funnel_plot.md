# Create a Funnel Plot for Bayesian Meta-Analysis

Creates a publication-ready funnel plot for meta-analysis, displaying
study-level effect sizes against their precision (standard error).
Supports both brmsfit and bayesma model objects. The plot aesthetic
matches the bayesfoRest package style.

## Usage

``` r
funnel_plot(
  model,
  data,
  measure,
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
  subgroup = FALSE,
  subgroup_var = NULL,
  contour_lines = TRUE,
  contour_alpha = 0.15,
  pooled_line = c("posterior", "fixed", "none"),
  color_points = "dodgerblue",
  color_points_outline = "blue",
  color_pooled = "blue",
  color_null = "black",
  color_contour = "grey50",
  point_size = 3,
  xlim = NULL,
  x_breaks = NULL,
  null_value = NULL,
  add_null_line = FALSE,
  label_studies = FALSE,
  label_size = 3,
  title = NULL,
  subtitle = NULL,
  title_align = "left",
  font = NULL
)
```

## Arguments

- model:

  A fitted model object. Either a brmsfit object (class 'brmsfit') or a
  bayesma object (class 'bayesma').

- data:

  A data frame containing the study data used for the meta-analysis.

- measure:

  Character string specifying the effect measure. Must be one of: "OR"
  (Odds Ratio), "HR" (Hazard Ratio), "RR" (Risk Ratio), "IRR" (Incidence
  Rate Ratio), "MD" (Mean Difference), or "SMD" (Standardized Mean
  Difference).

- studyvar:

  Column name containing study identifiers/authors. Default is NULL.

- year:

  Column name containing publication years. Default is NULL.

- c_n:

  Column name containing control group sample sizes.

- i_n:

  Column name containing intervention group sample sizes.

- c_event:

  Column name containing control group event counts.

- i_event:

  Column name containing intervention group event counts.

- c_mean:

  Column name containing control group means.

- i_mean:

  Column name containing intervention group means.

- c_sd:

  Column name containing control group standard deviations.

- i_sd:

  Column name containing intervention group standard deviations.

- c_time:

  Column name containing control group time periods.

- i_time:

  Column name containing intervention group time periods.

- subgroup:

  Logical indicating whether to colour points by subgroup. Default is
  FALSE.

- subgroup_var:

  Character string. Name of the variable in data to use for subgroup
  colouring.

- contour_lines:

  Logical indicating whether to add significance contour lines (at p =
  0.10, 0.05, 0.01). Default is TRUE.

- contour_alpha:

  Numeric. Alpha transparency for contour shading. Default is 0.15.

- pooled_line:

  Character string specifying the pooled estimate line style. Options:
  "posterior" (default, uses posterior median), "fixed" (uses the
  fixed/common effect), or "none" (no pooled line).

- color_points:

  Color for study points. Default is "dodgerblue".

- color_points_outline:

  Color for study point outlines. Default is "blue".

- color_pooled:

  Color for the pooled effect line. Default is "blue".

- color_null:

  Color for the null effect line. Default is "black".

- color_contour:

  Color for significance contour fills. Default is "grey50".

- point_size:

  Numeric. Size of the study points. Default is 3.

- xlim:

  Numeric vector of length 2 specifying x-axis limits. Default is NULL
  (auto-scaled).

- x_breaks:

  Numeric vector specifying custom x-axis break points. Default is NULL.

- null_value:

  Numeric. X-axis value of the null line. Default is NULL
  (measure-specific).

- add_null_line:

  Logical. Whether to draw a vertical line at the null value. Default is
  FALSE.

- label_studies:

  Logical indicating whether to label study points. Default is FALSE.

- label_size:

  Numeric. Size of study labels. Default is 3.

- title:

  Character string for the plot title. Default is NULL.

- subtitle:

  Character string for the plot subtitle. Default is NULL.

- title_align:

  Character string specifying title alignment. Options: "left"
  (default), "center"/"centre", "right".

- font:

  Character string specifying the font family. Default is NULL.

## Value

A ggplot object containing the funnel plot.
