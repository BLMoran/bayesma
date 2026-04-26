# Compute R-squared for Meta-Regression

Calculates the proportion of between-study heterogeneity explained by
moderators. Requires fitting a model without moderators for comparison.

## Usage

``` r
r_squared(
  object,
  null_model = NULL,
  method = c("tau", "tau_sq"),
  summary = c("median", "mean")
)
```

## Arguments

- object:

  A `bayesma_reg` object (model with moderators).

- null_model:

  Optional. A `bayesma` object without moderators. If NULL,

  the function will fit one automatically using stored call arguments.

- method:

  Character. Method for computing R². Options:

  - `"tau"`: Based on reduction in tau (default). R² = 1 - (tau_reg /
    tau_null)

  - `"tau_sq"`: Based on reduction in tau². R² = 1 - (tau²_reg /
    tau²_null

- summary:

  Character. How to summarize posterior. `"median"` (default) or
  `"mean"`.

## Value

A list with:

- `R2`: Point estimate of R²

- `R2_ci`: 95% credible interval for R²

- `R2_draws`: Full posterior distribution of R²

- `tau_null`: Tau from null model

- `tau_reg`: Tau from regression model

- `method`: Method used

## Details

R² in meta-regression represents the proportion of between-study
variance (heterogeneity) explained by the moderators. It is calculated
as:

\$\$R^2 = 1 - \frac{\tau^2\_{regression}}{\tau^2\_{null}}\$\$

where \\\tau^2\_{null}\\ is the heterogeneity from a model without
moderators, and \\\tau^2\_{regression}\\ is the residual heterogeneity
after accounting for moderators.

Note that R² can be negative if the moderators do not explain any
heterogeneity (or if the model with moderators fits worse due to
overfitting with weak priors).

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit meta-regression
fit_reg <- meta_reg(data, study = "author", yi = "yi", vi = "vi",
                    mods = ~ year + quality)

# Compute R²
r2 <- r_squared(fit_reg)
print(r2)

# With pre-fitted null model
fit_null <- bayesma(data, study = "author", ...)
r2 <- r_squared(fit_reg, null_model = fit_null)
} # }
```
