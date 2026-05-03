# Create a Table of Included Studies

Produces a formatted
[`gt::gt()`](https://gt.rstudio.com/reference/gt.html) table from a
study-level data frame. Handles continuous and binary outcomes, subgroup
grouping, and traffic-light risk-of-bias formatting for RoB 2, ROBINS-I,
NOS, and QUADAS-2.

## Usage

``` r
study_table(
  data,
  rob_tool = c("auto", "none", "rob2", "robins_i", "newcastle_ottawa", "quadas2"),
  outcome_type = c("auto", "continuous", "binary"),
  subgroup_col = NULL,
  font = NULL
)
```

## Arguments

- data:

  A data frame or tibble with one row per study. Expected columns depend
  on `outcome_type` and `rob_tool` — see the bundled templates
  ([`download_template()`](https://blmoran.github.io/bayesma/reference/download_template.md))
  for the exact layout.

- rob_tool:

  Risk-of-bias tool. One of `"auto"` (detect from column names),
  `"none"`, `"rob2"`, `"robins_i"`, `"newcastle_ottawa"`, or
  `"quadas2"`.

- outcome_type:

  Whether outcomes are `"continuous"`, `"binary"`, or `"auto"` (detect
  from column names).

- subgroup_col:

  Unquoted name of the subgroup column, or `NULL` (no grouping).

- font:

  Font family for the table. `NULL` uses the gt default.

## Value

A [`gt::gt()`](https://gt.rstudio.com/reference/gt.html) table object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(binary_outcome)
study_table(binary_outcome)

data(cont_outcome)
study_table(cont_outcome)
} # }
```
