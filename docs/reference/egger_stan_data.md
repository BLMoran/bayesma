# Build the Stan data list for an Egger's test specification

Build the Stan data list for an Egger's test specification

## Usage

``` r
egger_stan_data(spec)
```

## Arguments

- spec:

  An `egger_spec` object.

## Value

A named list suitable for `cmdstanr::CmdStanModel$sample(data = ...)`.
