# Bubble Plot for Meta-Regression

Creates a bubble plot showing the relationship between a continuous
moderator and the effect size, with point sizes proportional to study
precision (inverse variance).

## Usage

``` r
bubble_plot(
  object,
  mod,
  ci = TRUE,
  ci_level = 0.95,
  size_scale = 1,
  show_studies = FALSE,
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  color_palette = c("#4292C6", "#08519C"),
  theme = ggplot2::theme_minimal()
)
```

## Arguments

- object:

  A `bayesma_reg` object.

- mod:

  Character. Name of the moderator variable to plot. Must be a
  continuous moderator.

- ci:

  Logical. Show confidence/credible band for the regression line
  (default: TRUE).

- ci_level:

  Numeric. Credible interval level (default: 0.95).

- size_scale:

  Numeric. Scaling factor for bubble sizes (default: 1).

- show_studies:

  Logical. Label study points (default: FALSE).

- xlab, ylab:

  Character. Axis labels. If NULL, uses variable names.

- title:

  Character. Plot title.

- color_palette:

  Character vector of length 2. Colors for points and regression line.

- theme:

  A ggplot2 theme (default: `theme_minimal()`).

## Value

A ggplot object.

## Details

The bubble plot is a standard visualization for meta-regression with
continuous moderators. Each study is represented by a circle:

- **Position**: x = moderator value, y = effect size

- **Size**: Proportional to study weight (1 / variance)

- **Line**: Regression line showing the moderator effect

- **Band**: Credible interval for the regression line

For centered moderators, the x-axis shows the centered values by
default. Use `centered = FALSE` in the original
[`meta_reg()`](https://blmoran.github.io/bayesma/reference/meta_reg.md)
call if you want the original scale.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- meta_reg(data, studyvar = "author", yi = "yi", vi = "vi",
                mods = ~ year + dose)

# Basic bubble plot
bubble_plot(fit, mod = "year")

# Customized plot
bubble_plot(fit, mod = "dose",
            xlab = "Dose (mg)",
            ylab = "Log Odds Ratio",
            title = "Dose-Response Relationship")
} # }
```
