# Download a Data Entry Template

Writes a CSV template to disk so you know exactly which columns are
expected for each input data frame. Templates have header rows only;
fill in one row per study.

## Usage

``` r
download_template(
  type = c("continuous", "binary", "rob2", "robins_i", "newcastle_ottawa", "quadas2",
    "inspect_sr"),
  path = NULL
)
```

## Arguments

- type:

  Template type. One of `"continuous"`, `"binary"`, `"rob2"`,
  `"robins_i"`, `"newcastle_ottawa"`, `"quadas2"`, or `"inspect_sr"`.

- path:

  Directory to write the template into. Defaults to the current working
  directory.

## Value

Invisibly returns the path to the written file.

## Examples

``` r
if (FALSE) { # \dontrun{
download_template("binary")
download_template("rob2", path = "~/Desktop")
} # }
```
