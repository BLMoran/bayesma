# Create a gt Table for RoBMA Results

Produces a publication-ready `gt` table summarising the Robust Bayesian
Meta-Analysis results. Includes component summary (posterior
probabilities and Bayes factors), model-averaged estimates, and
direction probabilities.

## Usage

``` r
robma_table(
  x,
  digits = 3,
  include_components = TRUE,
  include_estimates = TRUE,
  include_direction = TRUE,
  include_models = FALSE,
  exponentiate = FALSE
)
```

## Arguments

- x:

  A `bayesma_robma` object from
  [`robma()`](https://blmoran.github.io/bayesma/reference/robma.md).

- digits:

  Integer. Number of decimal places for numeric values. Default is 3.

- include_components:

  Logical. If `TRUE` (default), includes the components summary table
  (Effect, Heterogeneity, Bias).

- include_estimates:

  Logical. If `TRUE` (default), includes the model-averaged parameter
  estimates table.

- include_direction:

  Logical. If `TRUE` (default), includes direction probabilities.

- include_models:

  Logical. If `TRUE`, includes individual model posterior probabilities
  (bridge method only). Default is `FALSE`.

- exponentiate:

  Logical. If `TRUE` and the effect is on a log

## Value

A `gt` table object.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- robma(data, studyvar = "study", ...)
robma_table(fit)

# Options
robma_table(fit, include_models = TRUE)
robma_table(fit, digits = 4)

# Save
robma_table(fit) |> gt::gtsave("robma_results.html")
robma_table(fit) |> gt::gtsave("robma_results.docx")
} # }
```
