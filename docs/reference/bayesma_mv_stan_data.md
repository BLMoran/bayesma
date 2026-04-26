# Build the Stan data list for a bivariate meta-analysis specification

Build the Stan data list for a bivariate meta-analysis specification

## Usage

``` r
bayesma_mv_stan_data(spec)
```

## Arguments

- spec:

  A `bayesma_mv_spec` object.

## Value

A named list suitable for `cmdstanr::CmdStanModel$sample(data = ...)`.
