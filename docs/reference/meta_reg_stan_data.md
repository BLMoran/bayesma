# Build the Stan data list for a meta-regression specification

Build the Stan data list for a meta-regression specification

## Usage

``` r
meta_reg_stan_data(spec)
```

## Arguments

- spec:

  A `meta_reg_spec` object.

## Value

A named list suitable for `cmdstanr::CmdStanModel$sample(data = ...)`.
