# Generate Stan code for a bivariate meta-analysis specification

Generate Stan code for a bivariate meta-analysis specification

## Usage

``` r
bayesma_mv_stan_code(spec, format = TRUE)
```

## Arguments

- spec:

  A `bayesma_mv_spec` object.

- format:

  Logical. If `TRUE` (default), run the generated program through Stan's
  `stanc --auto-format` for consistent indentation and spacing. Falls
  back to the raw program if the formatter is unavailable.

## Value

An object of class `"bayesma_mv_stan_code"`.
