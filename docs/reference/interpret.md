# Interpret a Bayesian meta-analysis workflow

Generates a comprehensive narrative interpretation across one or more
bayesma fits – overall effects, heterogeneity, publication bias, model
averaging, sensitivity, model comparison, and convergence diagnostics.
The function auto-detects what is provided and assembles only the
relevant sections.

## Usage

``` r
interpret(
  ...,
  null_range = c(-0.1, 0.1),
  effect_label = NULL,
  credible_level = 0.95,
  quiet = FALSE
)
```

## Arguments

- ...:

  One or more fitted bayesma objects (named or unnamed). Accepted
  classes: `bayesma`, `bayesma_mv`, `bayesma_robma`, `bayesma_egger`,
  `bayesma_metareg`, `bayesma_robma_sensitivity`, `bayesma_comparison`.

- null_range:

  Optional length-2 numeric vector for direction/ROPE probabilities.
  Defaults to `c(-0.1, 0.1)` on the natural scale.

- effect_label:

  Optional character override for the effect label used in narrative
  text (e.g. "log_or"). Inherited from fits when `NULL`.

- credible_level:

  Credible interval width used in summaries. Default `0.95`.

- quiet:

  If `TRUE`, suppresses progress messages during assembly.

## Value

An object of class `bayesma_interpretation` – a list with one element
per detected section plus a `meta` slot. The
[`print()`](https://rdrr.io/r/base/print.html) method renders the full
narrative report.

## Examples

``` r
if (FALSE) { # \dontrun{
fit  <- bayesma::bayesma(data = dat, yi = "yi", sei = "sei")
rob  <- bayesma::robma(data = dat, yi = "yi", sei = "sei")
egg  <- bayesma::egger(data = dat, yi = "yi", sei = "sei")
sens <- bayesma::robma_sensitivity(data = dat, priors = my_priors)
interpret(fit, rob, egg, sens)
} # }
```
