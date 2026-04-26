# Generate Stan code for a bayesma specification

Returns both the full Stan program and its decomposition into the
standard blocks (`functions`, `data`, `transformed_data`, `parameters`,
`transformed_parameters`, `model`, `generated_quantities`). Missing
blocks are empty strings.

## Usage

``` r
bayesma_stan_code(spec, format = TRUE)
```

## Arguments

- spec:

  A `bayesma_spec` object.

- format:

  Logical. If `TRUE` (default), run the generated program through Stan's
  `stanc --auto-format` for consistent indentation and spacing. Falls
  back to the raw program if the formatter is unavailable.

## Value

An object of class `"bayesma_stan_code"` – a list with elements
`functions`, `data`, `transformed_data`, `parameters`,
`transformed_parameters`, `model`, `generated_quantities`, and `full`.

## Details

If `spec$custom_model` is non-NULL, the user's Stan code is returned
verbatim (and parsed into blocks for inspection).
