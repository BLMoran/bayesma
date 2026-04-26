# Build cmdstanr data list(s) for a RoBMA specification

Returns a named list keyed by model label (bridge) or by bias indicator
(ss). Each element is a list suitable for
`cmdstanr::CmdStanModel$sample()`. Analytic bridge components contribute
a `NULL` entry.

## Usage

``` r
robma_stan_data(spec)
```

## Arguments

- spec:

  A `bayesma_robma_spec`.

## Value

An object of class `"bayesma_robma_stan_data"`.

## Details

`spec$custom_data` overlays the auto-built data list. If `custom_data`
is a list keyed by label, per-label overrides are applied; for the ss
case a non-nested list is merged at the top level of the single program.
