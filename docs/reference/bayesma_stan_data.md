# Build the Stan data list for a bayesma specification

Constructs the list passed to
[`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).
When `spec$custom_data` is provided, those elements override or augment
the automatically built list.

## Usage

``` r
bayesma_stan_data(spec)
```

## Arguments

- spec:

  A `bayesma_spec` object.

## Value

A named list suitable for `cmdstanr::CmdStanModel$sample(data = ...)`.
