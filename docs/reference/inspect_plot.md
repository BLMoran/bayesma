# Create INSPECT-SR Trustworthiness Plot

Create INSPECT-SR Trustworthiness Plot

## Usage

``` r
inspect_plot(
  data,
  sort_studies_by = "study",
  add_legend = FALSE,
  incl_checks = FALSE,
  font = NULL,
  title = NULL,
  title_align = "left",
  subtitle = NULL
)
```

## Arguments

- data:

  An `inspect_sr` object or a data frame with columns Study, D1, D2, D3,
  D4, Overall.

- sort_studies_by:

  `"study"` (default, alphabetical) or `"overall"` (by severity).

- add_legend:

  Logical. Add a legend panel (default FALSE).

- incl_checks:

  Logical. If `TRUE`, expand the table to show every INSPECT-SR item
  (1.1, 1.2, …, 4.11) grouped under domain spanners, with a per-domain
  "Overall" and a final study-level "Overall". Requires `data` to be an
  `inspect_sr` object. Default `FALSE`.

- font:

  Character. Font family (NULL = default).

- title, subtitle:

  Character. Optional.

- title_align:

  `"left"` (default), `"center"`/`"centre"`, or `"right"`.

## Value

A patchwork object.
