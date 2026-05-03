# Launch the bayesma Data Entry App

Opens a local Shiny app for entering study details, risk-of-bias
assessments, and INSPECT-SR judgements. Data can be exported as CSV
files ready for use with
[`study_table()`](https://blmoran.github.io/bayesma/reference/study_table.md)
and
[`inspect_sr()`](https://blmoran.github.io/bayesma/reference/inspect_sr.md).

## Usage

``` r
launch_data_entry()
```

## Value

Called for its side effect; returns `NULL` invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
launch_data_entry()
} # }
```
