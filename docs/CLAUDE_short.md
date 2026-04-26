# R coding preferences

## Defaults

- Tidyverse syntax throughout.
- Native pipe `|>`, not `%>%`.
- Assignment `<-`, not `=`.
- [`purrr::map()`](https://purrr.tidyverse.org/reference/map.html)
  family over `for` loops.
- Explicit `pkg::fn()` namespacing in every function body.
- `snake_case` names. Dots in function names are reserved for S3
  dispatch.

## Response style

- Lead with the answer. Bullets over prose.
- Yes/no question → yes/no first, context after.
- Don’t explain every detail unless asked.

## Code style

- Minimal comments. Names carry the meaning.
- No running narration (`# First, filter...`, `# Then fit...`).
- No meta/changelog comments (`# Updated`, `# NEW`,
  `# Changed from...`).
- No “why I chose this” commentary in chunks — put it in the response
  text around the chunk.

## Function design

- Primary data is the first argument.
- Required args have no defaults. Optional args always do.
- `...` between required and optional args.
- Prefix internal args with `.` when `...` forwards user data (`.by`,
  `.f`).
- String options →
  [`rlang::arg_match()`](https://rlang.r-lib.org/reference/arg_match.html)
  with enumerated default.
- Strategy variations → objects or enums, not boolean flags.
- Side-effect functions return input
  [`invisible()`](https://rdrr.io/r/base/invisible.html).
- Output type predictable from input types, not values.

## Tidyeval

- Bare column → embrace `{{ var }}`.
- String column → `.data[[var]]`.
- Character vector → `dplyr::all_of(vars)` in `across()`/`select()`.
- Bare names via `...` → forward `...` directly.
- Dynamic output names → `"{{ var }}" :=` or
  [`rlang::englue()`](https://rlang.r-lib.org/reference/englue.html).
- Don’t parse strings as code. Don’t
  [`get()`](https://rdrr.io/r/base/get.html) inside a data mask.

## Modern dplyr / tidyr

- Per-op grouping → `.by`, not `group_by()/ungroup()` sandwich.
- Joins → `dplyr::join_by(id == user_id)`, not `by = c("a" = "b")`.
- Dropping rows →
  [`dplyr::filter_out()`](https://dplyr.tidyverse.org/reference/filter.html),
  not [`filter()`](https://rdrr.io/r/stats/filter.html) with negation.
- OR in filter →
  [`dplyr::when_any()`](https://dplyr.tidyverse.org/reference/when-any-all.html),
  not nested parens.
- Value recode →
  [`dplyr::recode_values()`](https://dplyr.tidyverse.org/reference/recode-and-replace-values.html)
  / `replace_values()`, not chained `==` in `case_when()`.
- Add `.unmatched = "error"` when all cases are meant to be covered.
- Row-bind many frames → `purrr::map() |> purrr::list_rbind()`, not
  `map_dfr()`.
- Multi-row summaries →
  [`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html).
- [`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)/`pivot_wider()`,
  not `gather`/`spread`.

## Errors

``` r
validate_x <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (bad(x)) {
    cli::cli_abort("{.arg {arg}} must be ...", call = call)
  }
  invisible(x)
}
```

User-facing errors always go through
[`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
with `call = caller_env()`.

## Classes

- New structured class → S7 with a validator.
- Vector-like column class → vctrs.
- Writing S3 methods for existing generics (`print`, `format`,
  `summary`, `as_tibble`) → S3, registered via
  `S7::method(print, MyClass) <- ...`.
- Use `@` for S7 property access, not `$`.

## Progress + Bayesian output

Multi-step functions: print `"Step N: <verb> <what>"` via
[`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html),
gated by a `quiet = FALSE` argument.

Bayesian fits: suppress compilation/sampling chatter.

``` r
# cmdstanr
model$sample(data = d, refresh = 0, show_messages = FALSE, show_exceptions = FALSE)

# brms
brms::brm(y ~ x, data = d, backend = "cmdstanr", refresh = 0, silent = 2)
```

## Outputs

- Tables user-facing →
  [`gt::gt()`](https://gt.rstudio.com/reference/gt.html), not
  `kableExtra`.
- Plot composition → `patchwork` (`p1 + p2`, `p1 / p2`), not
  `gridExtra`/`cowplot`.
- Vignettes → `.qmd` with `%\VignetteEngine{quarto::html}`, not `.Rmd`.

## Anti-patterns to avoid

``` r
data %>% dplyr::filter(x > 0)                    # use |>
dplyr::inner_join(x, y, by = c("a" = "b"))       # use join_by()
purrr::map_dfr(xs, f)                            # use map() |> list_rbind()
data |> dplyr::group_by(g) |> summarise(...) |> ungroup()   # use .by
sapply(xs, f)                                    # use purrr::map_*()
dplyr::filter(!cond & !is.na(cond))              # use filter_out()
dplyr::case_when(x == 1 ~ "a", x == 2 ~ "b")     # use recode_values()
dplyr::case_match(x, ...)                        # soft-deprecated
structure(list(...), class = "range")            # use S7::new_class()
x$start                                          # use x@start on S7
brms::brm(y ~ x, data = d)                       # add refresh=0, silent=2
knitr::kable(df) |> kableExtra::kable_styling()  # use gt::gt()
options(digits = 3)                              # use withr::local_options()
```
