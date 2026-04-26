# CLAUDE.md

Guidance for writing R code, with emphasis on function design and
package development. Grounded in the tidy tools manifesto, tidy design
principles, and modern tidyverse / tidy evaluation practice.

## User defaults

- Write R using tidyverse syntax.
- Prefer
  [`purrr::map()`](https://purrr.tidyverse.org/reference/map.html)
  family over `for` loops.
- Use explicit `pkg::fn()` namespacing inside function bodies unless
  told otherwise.
- Use the native pipe `|>`, not `%>%`.
- Assignment is `<-`, not `=`.

## Response style

- Don’t explain every detail of something unless asked.
- Bullets over verbose prose.
- Lead with the answer; support it with just enough reasoning.
- If the user asks a yes/no or short question, answer that first, then
  add context only if needed.

## Code style

- Minimal comments. Only when a line needs one to make sense — not a
  running narration.
- No meta-commentary in code (`# Updated`, `# Modified`, `# New:`,
  `# Changed from previous`). Code should read as if it was written
  once, not as a changelog.
- No “why I chose this” commentary in code chunks. If it matters, put it
  in the response text around the chunk.
- Prefer self-explanatory function and variable names over explanatory
  comments.

``` r
# Bad — running narration
fit_model <- function(data) {
  # First, filter out missing values
  clean <- data |> dplyr::filter_out(is.na(y))
  # Now fit the model
  mod <- lm(y ~ x, data = clean)
  # Return the model
  mod
}

# Good — names carry the meaning
fit_model <- function(data) {
  data |>
    dplyr::filter_out(is.na(y)) |>
    lm(y ~ x, data = _)
}
```

## Progress messages in multi-step functions

For functions that run several discrete steps (fitting multiple models,
long pipelines, simulation workflows), print a short progress message at
each step using
[`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html).
Format: `"Step N: <verb> <what>"`.

``` r
fit_all_models <- function(data) {
  cli::cli_inform("Step 1: Compiling Model 1: Null Model")
  m1 <- fit_null(data)

  cli::cli_inform("Step 2: Compiling Model 2: Main Effects")
  m2 <- fit_main(data)

  cli::cli_inform("Step 3: Compiling Model 3: Interaction")
  m3 <- fit_interaction(data)

  list(null = m1, main = m2, interaction = m3)
}
```

Gate these behind a `quiet = FALSE` argument when the function might be
called inside a larger pipeline — the user can suppress them when
they’re noise.

``` r
fit_all_models <- function(data, quiet = FALSE) {
  inform <- if (quiet) \(...) invisible() else cli::cli_inform
  inform("Step 1: Compiling Model 1: Null Model")
  ...
}
```

## Bayesian model output

Suppress the default cmdstanr / brms compilation and sampling chatter:

    Compiling Stan model...
    Running MCMC with 4 parallel chains...
    Chain 3 finished in 1.1 seconds.
    ...
    Total execution time: 3.3 seconds.

Keep the user’s own step messages visible; silence everything else.

``` r
# cmdstanr
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,            # suppress iteration progress
  show_messages = FALSE,  # suppress informational messages
  show_exceptions = FALSE # suppress Stan exception chatter
)

# brms
fit <- brms::brm(
  y ~ x,
  data = data,
  backend = "cmdstanr",
  refresh = 0,
  silent = 2              # 2 = suppress compilation + sampling messages
)

# rstan
fit <- rstan::sampling(
  model,
  data = stan_data,
  refresh = 0,
  verbose = FALSE
)
```

If the progress really is useful (a long chain that might hang), make it
opt-in via a `verbose = FALSE` argument rather than on by default.

## The four manifesto principles

Every piece of code should visibly serve at least one of these:

1.  **Reuse existing data structures.** Prefer tibbles/data frames for
    rectangular data; prefer existing atomic vector types for 1-D data.
    Build custom classes (S7 by default) on top of these only when
    nothing existing fits.
2.  **Compose simple functions with the pipe.** Each function does one
    thing; the pipe is how they combine. If a function can’t be
    described in one sentence, split it.
3.  **Embrace functional programming.** Immutable inputs,
    copy-on-modify, generics over mutable state,
    [`purrr::map()`](https://purrr.tidyverse.org/reference/map.html)
    over `for`.
4.  **Design for humans.** Thinking time is the bottleneck, not compute
    time. Name things for the reader.

## Function naming

- **Verbs in imperative mood** for actions: `fit_model()`,
  `summarise()`, `parse_date()`.
- **Nouns** are acceptable for builders that return an object:
  `tibble()`, `recipe()`, `geom_point()`.
- **Prefix families, don’t suffix them.** `str_detect()`,
  `str_replace()`, `str_split()` — autocomplete works on prefixes.
  Reserve suffixes for *variations on a theme*:
  [`purrr::map_int()`](https://purrr.tidyverse.org/reference/map.html),
  [`purrr::map_chr()`](https://purrr.tidyverse.org/reference/map.html).
- **Length inversely proportional to frequency.** Short names (`n()`,
  [`c()`](https://rdrr.io/r/base/c.html)) earn their brevity by being
  used everywhere. Rare operations deserve descriptive names
  (`validate_model_specification()`).
- `snake_case` throughout. Dots (`.`) in function names are reserved for
  S3 method dispatch — avoid them in new function names regardless of
  the class system you’re using.

## Argument design

### Order

1.  Primary data first (makes the function pipe-friendly).
2.  Arguments that determine the *shape* of the output next.
3.  `...` after all required arguments.
4.  Optional arguments with defaults last — and because they come after
    `...`, users must name them.

``` r
# Signature template
my_function <- function(data, primary_arg, ..., option_a = default, option_b = NULL) {
  ...
}
```

### Required vs. optional

- **Required arguments have no defaults.** If there’s no sensible
  default, don’t invent one — let the missing-argument error fire.
- **Optional arguments always have defaults.**
- **Keep defaults short.** If a default needs computation, use `NULL`
  and resolve in the body with `%||%`:

``` r
my_function <- function(x, weights = NULL) {
  weights <- weights %||% rep(1, length(x))
  ...
}
```

### Enumerate options

For string arguments with a fixed set of valid values, use a
character-vector default and
[`rlang::arg_match()`](https://rlang.r-lib.org/reference/arg_match.html).
The first element is the default; users get a helpful error on typos.

``` r
summarise_method <- function(x, method = c("mean", "median", "trimmed")) {
  method <- rlang::arg_match(method)
  ...
}
```

### Standard argument names

Use these when they fit — consistency across packages is worth more than
creative naming.

| Purpose                  | Use       |
|--------------------------|-----------|
| Primary input data frame | `data`    |
| New data for prediction  | `newdata` |
| Missing value handling   | `na_rm`   |
| Case weights             | `weights` |
| Predictors matrix/df     | `x`       |
| Outcome vector           | `y`       |

### Dot prefix for meta-arguments

When `...` forwards user data, prefix the function’s own arguments with
`.` to avoid name collisions (the dplyr / purrr convention): `.data`,
`.by`, `.f`, `.progress`.

## Tidy evaluation

Two families of non-standard evaluation exist in the tidyverse. Know
which one a function uses before writing against it.

- **Data-masking**:
  [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html),
  [`filter()`](https://rdrr.io/r/stats/filter.html), `summarise()`,
  `arrange()`, `group_by()`.
- **Tidy selection**:
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html),
  `relocate()`, `across()`,
  [`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html).

### Forwarding user inputs

| Situation | Tool |
|----|----|
| User passes a bare column name | Embrace: `{{ var }}` |
| User passes a string column name | Pronoun: `.data[[var]]` |
| User passes a character vector of names | `dplyr::all_of(vars)` inside `across()` / `select()` |
| User passes bare names via `...` | Forward `...` directly; no special syntax |
| User passes a single-arg tidy selection via `...` | Wrap in `c(...)` |

``` r
# Embrace pattern — bare column name
summarise_var <- function(data, var) {
  data |> dplyr::summarise(mean = mean({{ var }}, na.rm = TRUE))
}

# String column name — use .data pronoun
summarise_string <- function(data, var) {
  data |> dplyr::summarise(mean = mean(.data[[var]], na.rm = TRUE))
}

# Forward dots
my_group_by <- function(data, ...) {
  data |> dplyr::group_by(...)
}

# Character vector of columns
my_select <- function(data, vars) {
  data |> dplyr::select(dplyr::all_of(vars))
}
```

### Output names with glue syntax

`"{{ var }}"` interpolates the embraced expression into the output name.
Use `:=` when the LHS is dynamic.

``` r
mean_named <- function(data, var) {
  data |> dplyr::summarise("mean_{{ var }}" := mean({{ var }}, na.rm = TRUE))
}

# Allow user override via englue()
mean_named <- function(data, var, name = rlang::englue("mean_{{ var }}")) {
  data |> dplyr::summarise("{name}" := mean({{ var }}, na.rm = TRUE))
}
```

### Injection (`!!`, `!!!`)

Use embracing (`{{ }}`) by default. Reach for `!!` / `!!!` only when:

- Building symbols from strings: `!!rlang::sym(var)` (or safer:
  `.data[[var]]`).
- Splicing a list of arguments: `!!!args`.
- Splicing a character vector of column names: `!!!rlang::syms(vars)` —
  though `dplyr::across(dplyr::all_of(vars))` is usually cleaner.

### Patterns to avoid

``` r
# Don't parse strings as code
eval(parse(text = paste("mean(", var, ")")))   # unsafe

# Don't use get() inside a data mask
dplyr::with(data, mean(get(var)))              # collision-prone

# Don't embrace a value (only embrace function arguments)
my_fn <- function(x) {
  x <- force(x)
  rlang::quo(mean({{ x }}))                    # wrong — x is a value, not an argument
}
```

## Modern dplyr / tidyr patterns

### Grouping

Prefer `.by` for per-operation grouping. It returns ungrouped data,
which is almost always what you want.

``` r
# Preferred
data |> dplyr::summarise(mean = mean(value), .by = c(group, year))

# Avoid the group_by() / ungroup() sandwich
data |>
  dplyr::group_by(group, year) |>
  dplyr::summarise(mean = mean(value)) |>
  dplyr::ungroup()
```

### Joins

Use
[`dplyr::join_by()`](https://dplyr.tidyverse.org/reference/join_by.html)
— supports equality, inequality, rolling, and overlap joins with clearer
syntax.

``` r
# Good
dplyr::inner_join(x, y, by = dplyr::join_by(id == user_id))

# Inequality / rolling
dplyr::inner_join(x, y, dplyr::join_by(id == user_id, closest(date >= start)))

# Specify matching expectations for quality control
dplyr::inner_join(x, y, by = dplyr::join_by(id), multiple = "error", unmatched = "error")
```

### Other modern idioms

- [`dplyr::pick()`](https://dplyr.tidyverse.org/reference/pick.html) —
  select columns *inside* a data-masking function.
- [`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html)
  — summaries that return multiple rows per group.
- [`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)
  / `pivot_wider()` — not `gather()` / `spread()`.
- [`tidyr::separate_wider_delim()`](https://tidyr.tidyverse.org/reference/separate_wider_delim.html)
  / `separate_wider_regex()` — not `separate()` / `extract()`.
- `purrr::map() |> purrr::list_rbind()` — not the superseded
  `map_dfr()`.

### Filtering (dplyr 1.2+)

[`filter()`](https://rdrr.io/r/stats/filter.html) is optimised for
*keeping* rows and treats `NA` as `FALSE`, which works *against* you
when you want to *drop* rows. Use
[`dplyr::filter_out()`](https://dplyr.tidyverse.org/reference/filter.html)
— same syntax, but designed around specifying which rows to drop.

``` r
# Drop deceased patients from before 2012
# filter_out() treats NAs as FALSE (don't drop), which is almost always what you want
patients |> dplyr::filter_out(deceased, date < 2012)

# Rule of thumb: if you reach for `!`, `!=`, or `& !is.na()`, use filter_out()
```

For OR-combined conditions, use
[`dplyr::when_any()`](https://dplyr.tidyverse.org/reference/when-any-all.html)
instead of nesting everything in parentheses and pipes.
[`dplyr::when_all()`](https://dplyr.tidyverse.org/reference/when-any-all.html)
is the AND counterpart (rarely needed inside
[`filter()`](https://rdrr.io/r/stats/filter.html) since comma-separated
conditions already AND).

``` r
# Rows where US/CA scored 200-300, OR where PR/RU scored 100-200
countries |>
  dplyr::filter(dplyr::when_any(
    name %in% c("US", "CA") & dplyr::between(score, 200, 300),
    name %in% c("PR", "RU") & dplyr::between(score, 100, 200)
  ))
```

`when_any()` and `when_all()` are regular vector functions — usable
anywhere, not just in [`filter()`](https://rdrr.io/r/stats/filter.html).

### Recoding and replacing (dplyr 1.2+)

The recode/replace family now has four members, organised by whether
you’re matching with conditions or values, and whether you’re building a
new column or partially updating an existing one:

|  | **Recoding (new column)** | **Replacing (same column)** |
|----|----|----|
| Match with conditions | [`dplyr::case_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html) | [`dplyr::replace_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html) |
| Match with values | [`dplyr::recode_values()`](https://dplyr.tidyverse.org/reference/recode-and-replace-values.html) | [`dplyr::replace_values()`](https://dplyr.tidyverse.org/reference/recode-and-replace-values.html) |

Reach for each based on what you’re actually doing:

``` r
# Recoding values — creating a new column from matched values
likert |>
  dplyr::mutate(
    category = dplyr::recode_values(
      score,
      1 ~ "Strongly disagree",
      2 ~ "Disagree",
      3 ~ "Neutral",
      4 ~ "Agree",
      5 ~ "Strongly agree"
    )
  )

# With a lookup table — from/to interface is cleaner than splicing
lookup <- tibble::tribble(
  ~from, ~to,
      1, "Strongly disagree",
      2, "Disagree",
      3, "Neutral",
      4, "Agree",
      5, "Strongly agree"
)

likert |>
  dplyr::mutate(category = dplyr::recode_values(score, from = lookup$from, to = lookup$to))

# Replacing values — partial update, type-stable on the input
schools |>
  dplyr::mutate(
    name = name |>
      dplyr::replace_values(
        c("UNC", "Chapel Hill")       ~ "UNC Chapel Hill",
        c("Duke", "Duke University")  ~ "Duke"
      )
  )

# replace_values() also subsumes common na_if/coalesce/replace_na patterns
dplyr::replace_values(state, NA ~ "Unknown")          # fill NAs
dplyr::replace_values(state, "Unknown" ~ NA)          # flag as missing
dplyr::replace_values(state, c(NA, "Unknown") ~ "<missing>")  # normalise both
```

**[`dplyr::case_match()`](https://dplyr.tidyverse.org/reference/case_match.html)
is soft-deprecated** in favour of `recode_values()`. Don’t reach for it
in new code.

### Safety features (dplyr 1.2+)

Use `.unmatched = "error"` on `case_when()` and `recode_values()` when
you believe you’ve covered every case. The function errors on unhandled
input rather than silently returning `NA` — exactly the behaviour you
want for defensive programming in packages.

``` r
likert |>
  dplyr::mutate(
    category = score |>
      dplyr::recode_values(
        from = lookup$from,
        to = lookup$to,
        unmatched = "error"
      )
  )
# Errors immediately if score contains a value not in lookup$from
```

`NA` values must be explicitly handled (e.g. `NA ~ NA`) when
`unmatched = "error"` is set — missing values trigger the error unless
you’ve opted in.

A few smaller 1.2 changes worth knowing:

- [`dplyr::between()`](https://dplyr.tidyverse.org/reference/between.html)
  gains a `ptype` argument for controlling output type (useful with
  ordered factors).
- [`dplyr::if_else()`](https://dplyr.tidyverse.org/reference/if_else.html),
  `case_when()`, and `coalesce()` are significantly faster thanks to a
  vctrs C rewrite.
- Returning `!= 1` row per group from `summarise()` is now defunct — use
  [`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html).
- Underscored verbs (`mutate_()`, `arrange_()`, etc.) are defunct. Use
  tidy evaluation.
- `mutate_each()` / `summarise_each()` are defunct — use
  [`dplyr::across()`](https://dplyr.tidyverse.org/reference/across.html).

## Modern ggplot2 (4.0+)

ggplot2 4.0 rewrote its internals on S7 (which dovetails nicely with the
S7-first position elsewhere in this file). The practical implications:

- **Stricter input validation.** Invalid property types error at
  construction rather than failing silently later.
  `element_text(hjust = "foo")` now errors immediately.
- **`@` is the S7-idiomatic accessor.** `plot@data` rather than
  `plot$data`. `$` still works for backwards compatibility but will be
  phased out.
- **Extensions get double dispatch via `update_ggplot()`** (successor to
  `ggplot_add()`) — relevant when writing geoms/stats/scales for a
  package.

### Theme-driven layer defaults

Non-data styling (default colours, shapes, palettes) now lives in the
theme via `theme(geom = element_geom(...))`. This centralises the look
of a plot in one place instead of scattering it across `geom_*()` calls.

``` r
ggplot2::ggplot(mpg, ggplot2::aes(displ, hwy)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
  ggplot2::theme_gray(paper = "cornsilk", ink = "navy", accent = "tomato")
```

Themes now think in `ink` (foreground), `paper` (background), and
`accent` roles rather than raw `colour`/`fill`. Access a theme property
from inside `aes()` with
[`ggplot2::from_theme()`](https://ggplot2.tidyverse.org/reference/aes_eval.html):

``` r
ggplot2::ggplot(mpg, ggplot2::aes(class, displ)) +
  ggplot2::geom_boxplot(ggplot2::aes(colour = ggplot2::from_theme(accent))) +
  ggplot2::theme(
    geom = ggplot2::element_geom(accent = "tomato", paper = "cornsilk")
  )
```

Default palettes live in the theme too, under
`palette.{aesthetic}.{type}`:

``` r
ggplot2::theme(
  palette.colour.continuous = c("chartreuse", "forestgreen"),
  palette.shape.discrete    = c("triangle", "triangle open", "triangle down open")
)
```

### Theme shortcuts

`theme()` has a lot of long argument names. The `theme_sub_*()` family
is less verbose when you’re changing several related settings:

``` r
# Instead of this
ggplot2::theme(
  axis.ticks.x        = ggplot2::element_line(colour = "red"),
  axis.ticks.length.x = grid::unit(5, "mm"),
  panel.widths        = grid::unit(5, "cm"),
  panel.spacing.x     = grid::unit(5, "mm"),
  panel.background    = ggplot2::element_rect(fill = NA)
)

# Use this
ggplot2::theme_sub_axis_x(
  ticks        = ggplot2::element_line(colour = "red"),
  ticks.length = grid::unit(5, "mm")
) +
ggplot2::theme_sub_panel(
  widths       = grid::unit(5, "cm"),
  spacing.x    = grid::unit(5, "mm"),
  background   = ggplot2::element_rect(fill = NA)
)
```

The available shortcuts are `theme_sub_axis()` (with
`_x`/`_y`/`_top`/`_bottom`/`_left`/`_right` variants),
`theme_sub_legend()`, `theme_sub_panel()`, `theme_sub_plot()`, and
`theme_sub_strip()`. Also: `ggplot2::margin_auto(1, 2)` for CSS-style
margin specification.

### Labels from column metadata

ggplot2 4 picks up the `label` attribute on columns as the default
axis/legend title. This makes `labelled`, `Hmisc`, and `gt` workflows
feel native.

``` r
attr(df$bill_dep, "label") <- "Bill depth (mm)"
ggplot2::ggplot(df, ggplot2::aes(bill_dep, bill_len)) + ggplot2::geom_point()
# x-axis automatically reads "Bill depth (mm)"
```

You can also pass a data-dictionary-style named vector via
`labs(dictionary = ...)`:

``` r
dict <- c(bill_dep = "Bill depth (mm)", bill_len = "Bill length (mm)")
ggplot2::ggplot(penguins, ggplot2::aes(bill_dep, bill_len)) +
  ggplot2::geom_point() +
  ggplot2::labs(dictionary = dict)
```

Label priority (low → high): `aes()` expression → `labs(dictionary)` →
column `label` attribute → `labs(<aes> = ...)` → `scale_*(name)` →
`guide_*(title)`. Functions are accepted anywhere a label is taken, so
`labs(y = \(x) paste0(x, " (log scale)"))` works.

### Other 4.0 features worth knowing

- **Boxplot/violin styling** — `geom_boxplot()` gains `whisker.*`,
  `box.*`, `median.*`, `staple.*` argument families. `fatten` is
  deprecated in favour of `median.linewidth`. `geom_violin()` quantiles
  are now computed from the real data, not the density.
- **Position aesthetics** — `position_nudge()` now exposes
  `nudge_x`/`nudge_y` as mappable aesthetics. `position_dodge()` gains
  an `order` aesthetic for consistent group placement.
- **`stat_manual()`** — pass a function that transforms a data frame and
  skip the whole `Stat*` class ceremony. Good for quick custom stats.
- **`coord_*(reverse = "y")`** — flip axis direction at the coord level.
  Works with `coord_sf()` where `scale_*_reverse()` doesn’t.
- **Discrete secondary axes** via
  `scale_x_discrete(sec.axis = dup_axis(...))`.
- **`facet_wrap(dir = "br")`** etc. — eight filling directions via
  two-letter codes (`tl`, `lt`, `tr`, `rt`, `bl`, `lb`, `br`, `rb`).
- **`layer(layout = ...)`** — control how data is distributed across
  facet panels (`"fixed"`, `"fixed_rows"`, integer indices).

## Outputs

### Type stability

Output type must be predictable from input *types*, never from input
*values*.

``` r
# Bad — type depends on the condition
ifelse(cond, 1L, 2.0)

# Good — always the same type
dplyr::if_else(cond, 1L, 2L)

# Use typed purrr variants when you know the return type
purrr::map_dbl(xs, f)   # not purrr::map() when a numeric vector is wanted
purrr::map_chr(xs, f)
purrr::map_lgl(xs, f)
```

### Size stability

Output size should be a predictable function of input size (same rows in
/ same rows out, or documented behaviour otherwise). Prediction
functions should return a tibble with one row per input row.

### Side-effect functions return invisibly

Functions called for side effects (writing files, printing, plotting)
must return their primary input
[`invisible()`](https://rdrr.io/r/base/invisible.html) so they can live
inside a pipeline.

``` r
write_backup <- function(data, path) {
  readr::write_csv(data, path)
  invisible(data)
}

data |>
  write_backup("raw.csv") |>
  dplyr::filter(keep) |>
  write_backup("clean.csv")
```

### Partition side effects from computation

Don’t bury [`cat()`](https://rdrr.io/r/base/cat.html),
[`options()`](https://rdrr.io/r/base/options.html), or
[`Sys.setenv()`](https://rdrr.io/r/base/Sys.setenv.html) inside a
function that also computes a result. Gate side effects behind an
explicit argument (`verbose`, `quiet`) or move them to their own
function.

## Strategies

When a function offers multiple ways of doing something, avoid boolean
flags.

``` r
# Bad — what does perl = TRUE, fixed = TRUE even mean?
grepl(pattern, x, perl = TRUE, fixed = FALSE, ignore.case = TRUE)

# Good — strategy objects carry their own arguments
stringr::str_detect(x, stringr::regex(pattern, ignore_case = TRUE))
stringr::str_detect(x, stringr::fixed(pattern))
```

Extract strategies into helper constructors (`regex()`, `fixed()`,
`coll()`) when they need different arguments. A single enum (via
[`rlang::arg_match()`](https://rlang.r-lib.org/reference/arg_match.html))
is enough when strategies share arguments.

## Functional programming over loops

Default to
[`purrr::map()`](https://purrr.tidyverse.org/reference/map.html) and
friends. Reach for `for` only for genuinely iterative, side-effecting
work where each step depends on the previous.

``` r
# Simple map
results <- purrr::map(inputs, fit_model)

# Typed output
mean_vec <- purrr::map_dbl(dfs, \(df) mean(df$value, na.rm = TRUE))

# Row-binding many data frames
combined <- purrr::map(files, readr::read_csv) |> purrr::list_rbind()

# Two inputs in lockstep
plots <- purrr::map2(datasets, titles, make_plot)

# Side effects only
purrr::walk2(datasets, paths, readr::write_csv)

# Parallel (purrr 1.1+, requires mirai daemons)
mirai::daemons(4)
results <- purrr::map(inputs, purrr::in_parallel(\(x) expensive(x)))
mirai::daemons(0)
```

## Errors

Use [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
for user-facing errors. Pass `call = caller_env()` from validators so
the error points at the user’s code, not the internal helper.

``` r
validate_positive <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.numeric(x) || any(x <= 0)) {
    cli::cli_abort(
      "{.arg {arg}} must be a positive numeric vector.",
      call = call
    )
  }
  invisible(x)
}
```

Error messages should say what’s wrong *and* what to do about it. Prefer
`cli` class styles (`{.arg}`, `{.val}`, `{.fn}`) for consistent
formatting.

## Package development

### Dependencies

Add a dependency when the functionality gain is real (string
manipulation, dates, HTTP, modelling). Don’t add one for trivial base R
replacements. The tidyverse core (`dplyr`, `purrr`, `stringr`, `tidyr`,
`rlang`, `cli`) is usually worth it for data-adjacent packages; avoid
`Imports: tidyverse` in any package.

### Namespacing in function bodies

Inside every exported or internal package function, call imports with
`pkg::fn()`. This is explicit, survives namespace changes, and makes
reviews easier. Reserve `@importFrom` for infix operators (`%||%`, `:=`)
and for generics you define methods on.

### Roxygen tags for tidyeval

``` r
#' @param var <[`data-masked`][dplyr::dplyr_data_masking]> Column to summarise.
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> Passed to [dplyr::summarise()].
#' @param cols <[`tidy-select`][dplyr::dplyr_tidy_select]> Columns to select.
```

### Classes

**Default to S7 for new classes.** S7 is designed to supersede S3 and
S4, has CRAN status, and is what the tidyverse team reaches for in new
projects. It gives you formal class definitions, automatic property
validation on construction *and* mutation, safe `@` property access,
proper multiple dispatch, and a clear migration path as S7 moves toward
base R.

``` r
Range <- S7::new_class("Range",
  properties = list(
    start = S7::class_double,
    end = S7::class_double
  ),
  validator = function(self) {
    if (length(self@start) != 1 || length(self@end) != 1) {
      "@start and @end must be length 1"
    } else if (self@end < self@start) {
      "@end must be >= @start"
    }
  }
)

# Methods
inside <- S7::new_generic("inside", "x")
S7::method(inside, Range) <- function(x, y) {
  y >= x@start & y <= x@end
}

# Multiple dispatch — real strength over S3
combine <- S7::new_generic("combine", c("x", "y"))
S7::method(combine, list(Range, Range))   <- function(x, y) ...
S7::method(combine, list(Range, S7::class_numeric)) <- function(x, y) ...
```

**Pick the class system per-class, not per-package.** A single package
can mix systems where appropriate:

- **S7** — the default for new structured objects: properties,
  validators, inheritance, multiple dispatch.
- **vctrs** — for vector-like objects that should live in a data frame
  column (custom dates, units, percentages, categoricals). Solves
  problems S7 doesn’t: coercion rules, arithmetic, `vec_c()` behaviour,
  tibble integration. Use
  [`vctrs::new_vctr()`](https://vctrs.r-lib.org/reference/new_vctr.html)
  as the constructor.
- **S3** — acceptable when writing *methods* for existing S3 generics
  (`print`, `format`, `summary`, `as.data.frame`, `as_tibble`), and when
  extending packages whose class model is S3. For a thin labelled
  wrapper with nothing to validate, S3 is also fine — but once you have
  more than 2–3 slots or any validation to do, reach for S7.
- **S4** — only when working inside Bioconductor or an existing S4
  codebase.
- **R6** — only when you genuinely need mutable reference semantics
  (rare; usually a sign the design should be reconsidered).

**Providing S3 methods for an S7 class.** S7 and S3 are designed to
coexist. Every S7 object carries an S3 class vector automatically, so S3
generics (`print`, `format`, `summary`, `plot`, `as.data.frame`,
[`tibble::as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html),
[`knitr::knit_print`](https://rdrr.io/pkg/knitr/man/knit_print.html),
[`dplyr::as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html),
etc.) find methods through ordinary S3 dispatch. This is the idiomatic
way to make S7 objects feel native in the R ecosystem.

There are two ways to register an S3 method for an S7 class. Both work;
pick based on whether the generic is S3-only or S7-aware.

``` r
# Pattern 1: register with S7::method<-() 
# Use for S3 generics that you want to feel first-class alongside your
# other S7 methods. S7 registers the S3 method for you.
S7::method(print, Range) <- function(x, ...) {
  cat("<Range: [", x@start, ", ", x@end, "]>\n", sep = "")
  invisible(x)
}

S7::method(format, Range) <- function(x, ...) {
  sprintf("[%g, %g]", x@start, x@end)
}

S7::method(summary, Range) <- function(object, ...) {
  cat("Range of width", object@end - object@start, "\n")
  invisible(object)
}
```

``` r
# Pattern 2: write a plain S3 method using the class name
# Use when you want the method to be visible as a regular S3 method
# (easier for other packages to find via methods()) or when you prefer
# the familiar S3 syntax. Both patterns dispatch identically.
print.Range <- function(x, ...) {
  cat("<Range: [", x@start, ", ", x@end, "]>\n", sep = "")
  invisible(x)
}
```

For package development, **Pattern 1** is generally preferred — it keeps
all methods for your class in one visual location and plays well with
`S7::methods()` introspection. Pattern 2 is fine for one-off methods or
when you’re gradually migrating an S3 class to S7.

Either way, remember to register the S3 method in `NAMESPACE`:

``` r
#' @export
# (roxygen adds S3method(print, Range) automatically)
```

**Common S3 generics worth implementing for user-facing S7 classes:**

| Generic | Why |
|----|----|
| `print` | Default REPL display |
| `format` | Called by `print`, used in error messages, printed inside data frames |
| `summary` | `summary(x)` at the console |
| `plot` | Base-R plotting hook |
| `as.data.frame` / [`tibble::as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html) | Interop with the rest of the tidyverse |
| [`knitr::knit_print`](https://rdrr.io/pkg/knitr/man/knit_print.html) | Controls how the object renders in R Markdown / Quarto |

**Don’t write an S3 method for a generic that already has an S7 generic
available.** If there’s an S7 version of what you need
([`S7::new_generic()`](https://rconsortium.github.io/S7/reference/new_generic.html)
or an existing one), define that instead — you get validation and
multiple dispatch. S3 methods are the bridge to the existing ecosystem,
not the default choice for new generics.

**Model objects.** Never store the training data or
[`match.call()`](https://rdrr.io/r/base/match.call.html) inside a model
object — both silently capture huge amounts of state. Store only what
[`predict()`](https://rdrr.io/r/stats/predict.html) needs (coefficients,
factor levels, preprocessing recipe).

### Testing

``` r
testthat::test_that("summarise_var returns a tibble with one row per group", {
  result <- summarise_var(mtcars, mpg, .by = cyl)
  testthat::expect_s3_class(result, "tbl_df")
  testthat::expect_equal(nrow(result), 3)
})

testthat::test_that("summarise_var supports tidyeval", {
  testthat::expect_no_error(summarise_var(mtcars, mpg * 2))
  var <- rlang::sym("mpg")
  testthat::expect_no_error(summarise_var(mtcars, !!var))
})
```

Test *behaviour*, not implementation. Test error conditions with
`testthat::expect_error(..., class = "my_error_class")`.

### Tables: gt

Use `gt` for any table output from package functions — `summary`
methods, model coefficient tables, diagnostic outputs, anything
user-facing. Don’t return raw data frames styled with
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) or
`kableExtra` for presentation tables.

``` r
summary_table <- function(fit) {
  fit |>
    broom::tidy() |>
    gt::gt() |>
    gt::fmt_number(columns = c(estimate, std.error), decimals = 3) |>
    gt::fmt_number(columns = p.value, decimals = 4) |>
    gt::cols_label(
      term      = "Term",
      estimate  = "Estimate",
      std.error = "Std. Error"
    )
}
```

`gt` tables render cleanly in HTML, LaTeX, Word, and RTF, which matters
for CRAN vignettes that need to build across formats. For tables that
belong inside a larger composed figure, `gt` integrates with `patchwork`
via
[`patchwork::wrap_table()`](https://patchwork.data-imaginist.com/reference/wrap_table.html).

### Composition: patchwork

Use `patchwork` to combine ggplots (and gt tables) into multi-panel
figures. Avoid
[`gridExtra::grid.arrange()`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html)
or `cowplot` for new code — `patchwork` has a cleaner API and active
development.

``` r
library(patchwork)

# Side by side
p1 + p2

# Stacked
p1 / p2

# Layout grammar
(p1 | p2) / p3

# Collect shared legends
(p1 | p2 | p3) + patchwork::plot_layout(guides = "collect")

# Mix a gt table in
p1 + patchwork::wrap_table(gt_table)

# Annotate
(p1 | p2) +
  patchwork::plot_annotation(
    title = "Model diagnostics",
    tag_levels = "A"
  )
```

When a package function returns a composed figure, return the
`patchwork` object directly — users can add layers or change the layout
without re-running the function.

### Vignettes: Quarto (.qmd)

Write package vignettes as `.qmd` files, not `.Rmd`. Quarto produces
nicer output, has better cross-format handling, and is the direction
both Posit and the broader R community are moving.

Minimal Quarto vignette setup in `DESCRIPTION`:

    VignetteBuilder: quarto
    Suggests:
        quarto

YAML header for a package vignette (`vignettes/intro.qmd`):

``` yaml
---
title: "Introduction to <pkgname>"
vignette: >
  %\VignetteIndexEntry{Introduction to <pkgname>}
  %\VignetteEngine{quarto::html}
  %\VignetteEncoding{UTF-8}
format:
  html:
    toc: true
---
```

The `quarto` R package (on CRAN) provides
`%\VignetteEngine{quarto::html}`. Keep one `intro.qmd` as the main
“getting started” vignette, and write task-specific ones as separate
files (`modelling.qmd`, `plots.qmd`, etc.).

Make sure `quarto` (the CLI) is available in CI — on GitHub Actions,
`quarto-dev/quarto-actions/setup@v2` handles it.

## CRAN readiness

Every package is going to CRAN, so build to CRAN standards from day one
rather than retrofitting at submission time.

### Structural requirements

- `DESCRIPTION` — `Title` in title case without “R package for…”,
  `Description` in full sentences ending in full stops, `Authors@R` with
  ORCID where available, `License` as a standard SPDX identifier
  (e.g. `MIT + file LICENSE`, `GPL (>= 3)`).
- `NAMESPACE` — generated by roxygen2, never edited by hand.
- `LICENSE` file matches the `License` field exactly. For MIT use the
  two-line `LICENSE` file pointing to `LICENSE.md`.
- `NEWS.md` — required de facto, and essential on resubmission. One
  heading per version (`# pkgname 0.1.0`), bullet points for changes,
  newest at top.
- `README.md` — generated from `README.Rmd` (or `README.qmd`) via
  [`devtools::build_readme()`](https://devtools.r-lib.org/reference/build_rmd.html).
  Include installation, minimal example, and a link to the pkgdown site.
- `cran-comments.md` — your communication with CRAN maintainers. Include
  R CMD check results and a list of reverse dependencies checked.

### Required local checks

Run all of these clean before every submission:

``` r
devtools::check()                            # local R CMD check
devtools::check(remote = TRUE, manual = TRUE) # includes URL + manual checks
urlchecker::url_check()                      # every URL in the package
devtools::spell_check()                      # DESCRIPTION, roxygen, vignettes
goodpractice::gp()                           # broader hygiene
rcmdcheck::rcmdcheck(args = "--as-cran")     # CRAN-specific checks
```

Zero errors, zero warnings, zero notes is the target. Notes that are
genuinely unavoidable (new submission, possibly-misspelled surnames) go
in `cran-comments.md` with explanation.

### Multi-platform CI

Run R CMD check on CRAN’s target platforms via GitHub Actions. Use
`usethis::use_github_action("check-standard")` to get the canonical
workflow covering:

- Linux (release)
- macOS (release)
- Windows (release)
- Linux (devel)
- Linux (oldrel-1)

Before submission, additionally test on:

``` r
devtools::check_win_devel()    # win-builder, R-devel
devtools::check_win_release()  # win-builder, R-release
devtools::check_mac_release()  # mac-builder
rhub::rhub_check()             # rhub v2, includes sanitizer/valgrind builds
```

Reverse dependency checks matter from the second submission onward:

``` r
revdepcheck::revdep_check(num_workers = 4)
```

### Common CRAN trip-ups

- **Examples must run in \< 5 seconds** on CRAN. Wrap slow examples in
  `\donttest{}` (runs on CRAN but not `R CMD check --as-cran`) or
  `\dontrun{}` (never runs). Prefer `\donttest{}` where possible — the
  examples still get checked.
- **Vignettes and tests must run in \< 10 minutes total.** Skip slow
  things on CRAN with
  [`testthat::skip_on_cran()`](https://testthat.r-lib.org/reference/skip.html).
  For Bayesian models this usually means mocking the fit or using a tiny
  toy dataset.
- **No writing to the user’s home directory, no setting
  [`options()`](https://rdrr.io/r/base/options.html) without restoring
  them, no
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
  inside the package.** Use `withr::local_*()` for any state changes.
- **No [`cat()`](https://rdrr.io/r/base/cat.html) /
  [`print()`](https://rdrr.io/r/base/print.html) at package load.** If
  you need startup messages use
  [`packageStartupMessage()`](https://rdrr.io/r/base/message.html), and
  make them suppressible.
- **URLs must resolve** and redirect properly.
  [`urlchecker::url_check()`](https://rdrr.io/pkg/urlchecker/man/url_check.html)
  catches most issues; prefer canonical URLs over shortened ones.
- **Don’t ship files starting with `.`** in `inst/` or anywhere they’ll
  be installed (CRAN rejects these).
- **`LazyData: true`** only if the package has data. Otherwise omit —
  flagged otherwise.
- **Use `\doi{}` and `\url{}`** in Rd for citations. Plain URLs trigger
  NOTEs.
- **Non-ASCII characters** in source code or data require explicit
  encoding declaration in `DESCRIPTION` (`Encoding: UTF-8`).
- **Package size** — under 5 MB installed is the soft target. Large data
  goes in a separate data package.

### Version numbering

Use semantic versioning: `MAJOR.MINOR.PATCH`. Development versions get a
fourth component (e.g. `0.1.0.9000`). Bump via
[`usethis::use_version()`](https://usethis.r-lib.org/reference/use_version.html)
and
[`usethis::use_dev_version()`](https://usethis.r-lib.org/reference/use_version.html).

### Submission workflow

1.  [`usethis::use_release_issue()`](https://usethis.r-lib.org/reference/use_release_issue.html)
    — creates a checklist issue.
2.  Work through the checklist.
3.  [`devtools::submit_cran()`](https://devtools.r-lib.org/reference/submit_cran.html)
    — handles the actual submission.
4.  Wait for the auto-check result; respond to the maintainer email.
5.  On acceptance,
    [`usethis::use_github_release()`](https://usethis.r-lib.org/reference/use_github_release.html)
    tags and releases.

## Design review checklist

Before considering a function done:

Name is a verb (or a noun for a builder) and shares a prefix with its
family.

Primary data is the first argument.

Required arguments have no defaults.

`...` sits between required and optional arguments.

Internal arguments prefixed with `.` where `...` forwards user data.

String options use
[`rlang::arg_match()`](https://rlang.r-lib.org/reference/arg_match.html)
with an enumerated default.

Output type is a function of input types, not input values.

Side-effect functions return their input
[`invisible()`](https://rdrr.io/r/base/invisible.html)ly.

No hidden reads from [`options()`](https://rdrr.io/r/base/options.html),
locale, or [`Sys.getenv()`](https://rdrr.io/r/base/Sys.getenv.html).

Strategy variations use objects or enums, not boolean flags.

User-facing errors go through
[`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)
with `call = caller_env()`.

Tidyeval arguments are documented with the appropriate roxygen tag.

Package function bodies use `pkg::fn()` namespacing.

Dropping rows uses
[`dplyr::filter_out()`](https://dplyr.tidyverse.org/reference/filter.html),
not [`filter()`](https://rdrr.io/r/stats/filter.html) with negation.

Value-based recoding uses `recode_values()` / `replace_values()`, not
chained `==` in `case_when()`.

`case_when()` / `recode_values()` use `.unmatched = "error"` when all
cases are meant to be covered.

New classes use S7 (unless writing S3 methods for existing generics, or
building a vector-like class with vctrs).

S7 classes have a `validator` covering invariants you’d otherwise check
in a constructor.

Tables use `gt`, plot composition uses `patchwork`, vignettes are
`.qmd`.

Multi-step functions print `"Step N: ..."` messages, gated by a `quiet`
argument.

No comment narration or meta-comments (`# Updated`, `# New`) in code.

Bayesian fits run with `refresh = 0` and compilation/sampling chatter
suppressed.

Slow examples wrapped in `\donttest{}`; slow tests use
[`testthat::skip_on_cran()`](https://testthat.r-lib.org/reference/skip.html).

Package passes `devtools::check(remote = TRUE, manual = TRUE)` with zero
errors/warnings/notes.

## Anti-patterns

``` r
# Magrittr pipe in new code
data %>% dplyr::filter(x > 0)                        # use |>

# Character-vector join syntax
dplyr::inner_join(x, y, by = c("a" = "b"))           # use join_by(a == b)

# Superseded purrr
purrr::map_dfr(xs, f)                                # use map() |> list_rbind()

# Persistent grouping where .by would do
data |> dplyr::group_by(g) |> dplyr::summarise(...) |> dplyr::ungroup()

# Growing a vector in a loop
out <- c(); for (i in seq_along(x)) out <- c(out, f(x[[i]]))   # use purrr::map()

# Type-unstable base R
sapply(xs, f)                                        # use purrr::map_*()

# String parsing as code
eval(parse(text = paste0("mean(", var, ")")))       # use .data[[var]] or !!sym(var)

# Boolean strategy flags
my_fn(x, use_fast = TRUE, use_parallel = FALSE)     # enum or strategy object

# Storing training data on a model
model$training_data <- df                            # keep only what predict() needs

# Filtering out with negation and NA handling
data |> dplyr::filter(!(deceased & date < 2012))     # use filter_out()
data |> dplyr::filter(!cond & !is.na(cond))          # use filter_out()

# OR conditions nested in parens
data |> dplyr::filter((a & b) | (c & d))             # use when_any()

# Chained ==s in case_when()
dplyr::case_when(x == 1 ~ "a", x == 2 ~ "b", x == 3 ~ "c")  # use recode_values()

# case_match() in new code
dplyr::case_match(x, ...)                            # soft-deprecated — use recode_values()

# Partial-update dance with case_when/if_else
dplyr::if_else(is.na(state), "Unknown", state)       # use replace_values(NA ~ "Unknown")

# case_when() without .unmatched = "error" when you're sure you've covered everything
dplyr::case_when(x == 1 ~ "a", x == 2 ~ "b")         # add .unmatched = "error" for safety

# $ access on a ggplot2 4 object
plot$data                                             # use plot@data (S7 idiom)

# Hard-coding layer defaults instead of using the theme
ggplot2::geom_point(colour = "navy", fill = "cornsilk")   # set via element_geom() in theme

# fatten in boxplots
ggplot2::geom_boxplot(fatten = 2)                    # deprecated — use median.linewidth

# Defaulting to S3 for a new structured class
structure(list(start = 1, end = 10), class = "range")  # use S7::new_class() with a validator

# $ access on an S7 object
x$start                                              # use x@start — it's validated

# Narration comments
# First, filter the data
# Then fit the model
# Return the result
# — drop them, let names carry the meaning

# Meta / changelog comments in code
# Updated: now handles NA
# NEW: added weight argument
# — these belong in NEWS.md, not in function bodies

# Loud Bayesian fits in package functions
brms::brm(y ~ x, data = d)                           # add refresh = 0, silent = 2

# kableExtra or gridExtra in new package code
knitr::kable(df) |> kableExtra::kable_styling()      # use gt::gt()
gridExtra::grid.arrange(p1, p2)                      # use patchwork: p1 + p2

# .Rmd vignettes in new packages
vignettes/intro.Rmd                                  # use intro.qmd + quarto engine

# Package writes to ~/ or changes options() without restoring
options(digits = 3)                                  # use withr::local_options()
```
