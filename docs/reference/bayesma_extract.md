# Extract tidy effect components from a bayesma fit

Summaries, draws, per-study rows, and the prediction interval for a
fitted bayesma model. When `spec$custom_model` is set, extraction is
best-effort: `mu` is summarised if it exists in the fit, otherwise only
the raw summary table and draws are returned.

## Usage

``` r
bayesma_extract(fit, spec)
```

## Arguments

- fit:

  A `bayesma_fit` object.

- spec:

  A `bayesma_spec` object.

## Value

An object of class `"bayesma_effects"`.
