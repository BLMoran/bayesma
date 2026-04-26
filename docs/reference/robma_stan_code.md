# Generate Stan programs for a RoBMA specification

For `method = "bridge"` returns a named list keyed by model label. Each
entry is either a `bayesma_stan_code`-shaped list (parsed blocks +
`full`) or, for analytic component models with no Stan program, a list
with `analytic = TRUE` and `full = NA_character_`.

## Usage

``` r
robma_stan_code(spec, format = TRUE)
```

## Arguments

- spec:

  A `bayesma_robma_spec`.

- format:

  Logical. Run each program through `stanc --auto-format`.

## Value

An object of class `"bayesma_robma_stan_code"`.

## Details

For `method = "ss"` returns a single-entry list keyed by
`bias_indicator`.

If `spec$custom_model` is set, the user's Stan program(s) replace the
generated code: a character scalar is used for the ss case; a named list
is merged per-label for the bridge case.
