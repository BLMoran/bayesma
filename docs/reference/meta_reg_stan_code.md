# Generate Stan code for a meta-regression specification

Generate Stan code for a meta-regression specification

## Usage

``` r
meta_reg_stan_code(spec, format = TRUE)
```

## Arguments

- spec:

  A `meta_reg_spec` object.

- format:

  Logical. If `TRUE` (default), run the generated program through Stan's
  `stanc --auto-format` for consistent indentation and spacing. Falls
  back to the raw program if the formatter is unavailable.

## Value

An object of class `"meta_reg_stan_code"`.
