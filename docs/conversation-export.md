# Discussed and sketched a comprehensive refactoring of Stan code generation in the bayesma R package, consolidating 16 generator functions into 7 using a block-builder pattern and consolidated robma spec_component functions. Reviewed the PET-PEESE implementation’s use of sample-size-based predictors (1/sqrt(n) and 1/n), confirming it aligns with both Stanley (2017) and Pustejovsky (2017) simulation evidence showing these are preferable to SE-based predictors for avoiding mechanical correlation between d and V.

## User

Within the R/ folder are functions for a Bayesian Meta Analysis using
Stan. Every new function seems to be bigger than the ones previous. One
of the issues is the construction of the Stan model. There are numerous
different functions for each type of model that is being constructed.
I’m not sure if this is normal or not, but certain parts of the Stan
code don’t change (eg data). I was wondering if it were feasible or a
good idea to have a generate_stan_model function that takes different
arguments (eg model_type, stage, etc) that feed into the model or
whether it is better (and easier to debug) with the model construction
the way it is.

What do you think?

## Databot

The user wants my opinion on refactoring their Stan model construction
code. Let me first explore the codebase to understand the current
structure before giving advice.

Good question. Let me look at the current structure before giving you an
informed opinion.

Now I have a thorough understanding of the codebase. Let me give a
thoughtful opinion on whether to refactor into a unified
`generate_stan_model()` function.

This is a really well-informed question, and having looked through your
code, I can give you a grounded answer.

**Short answer: a unified `generate_stan_model()` is feasible and
probably worth doing, but with caveats.**

### What you’d gain

Your current setup has ~15 generator functions that all follow the same
block structure (data → parameters → transformed parameters → model →
generated quantities) and share significant overlap. A few concrete
examples of duplication:

- The `skew_normal_lpdf_custom` function block is copy-pasted
  identically between `generate_stan_code_two_stage()` and
  `generate_stan_code_one_stage()`
- The two-stage data block (`int S; vector y; vector se;`) appears in at
  least 5 places with minor conditional additions
- Generated quantities blocks (pooled, mu_new, y_rep) are largely
  standardized

A composable approach would reduce this duplication and make it easier
to add new model features without creating yet another function.

### What I’d actually recommend

Rather than one monolithic `generate_stan_model()` with many arguments,
I’d suggest a **block-builder pattern** — you’re already partway there
with helpers like `emit_prior_target()` and `emit_robust_likelihood()`.
The idea would be:

``` r
# Each block has its own builder
build_data_block(stage, likelihood, use_t_likelihood, use_multi_arm, ...)
build_parameters_block(model_type, re_dist, use_robust, ...)
build_transformed_parameters_block(...)
build_model_block(...)
build_generated_quantities_block(...)

# Then a thin orchestrator
generate_stan_model <- function(stage, model_type, likelihood, re_dist, ...) {
  parts <- c(
    build_functions_block(re_dist),
    build_data_block(stage, likelihood, ...),
    build_parameters_block(model_type, re_dist, ...),
    build_transformed_parameters_block(stage, model_type, ...),
    build_model_block(stage, likelihood, model_type, ...),
    build_generated_quantities_block(stage, model_type, ...)
  )
  paste(parts[nzchar(parts)], collapse = "\n\n")
}
```

This keeps each block independently testable and debuggable while
eliminating the duplication.

### The honest tradeoff

The risk is that a unified generator with many conditionals can become
harder to reason about than separate functions — “what Stan code does
this combination of arguments produce?” becomes a non-trivial question.
Your current approach, while repetitive, has the virtue that each
function is self-contained and readable top-to-bottom.

A few practical guidelines if you go ahead:

1.  **Write snapshot tests** — for each model configuration you support,
    store the expected Stan code as a text file and test against it.
    This is the single most important thing for debuggability.
    [`testthat::expect_snapshot()`](https://testthat.r-lib.org/reference/expect_snapshot.html)
    works well for this.
2.  **Keep the robma `spec_component_*` functions separate** for now.
    They have enough unique structure (returning
    `list(stan_code, stan_data)`) that forcing them into the same
    generator would add complexity without much benefit.
3.  **Start with the two-stage and one-stage generators**, which share
    the most overlap. The meta-regression and Egger functions can follow
    later.
4.  **Use the vector-collection pattern**
    ([`c()`](https://rdrr.io/r/base/c.html) then `paste(collapse)`)
    consistently rather than mixing it with `paste0` concatenation and
    `glue`. Your meta-regression functions already use this pattern
    well.

So: yes, it’s feasible and a good idea, but I’d approach it
incrementally — consolidate the most duplicated pieces first, keep good
test coverage, and don’t try to unify everything in one pass.

Want me to sketch out what the block builders might look like concretely
for your two-stage and one-stage generators?

## User

Yes please

## Databot

The user wants me to sketch out concrete block builders for their
two-stage and one-stage Stan generators. Let me read the actual code for
these functions to understand the exact conditionals and patterns.

Let me read the actual generator functions to build something concrete
against your real code.

Now let me also quickly check the existing helper functions for emitting
priors and robust code.

Now let me check the existing emit helpers so I know what’s available to
build on.

Good, I’ve got a thorough picture. Let me now write the sketch.

I’ve written the sketch to `R/generate_stan_model_sketch.R`. Here’s the
rationale behind the design decisions:

### Structure

The core idea is a **context object** (`ctx`) — a list that holds all
the configuration (stage, model_type, re_dist, likelihood, priors,
flags). Every builder receives this single object rather than threading
8+ arguments through every call. This makes it easy to add new flags
later without changing every function signature.

The hierarchy is:

- `generate_stan_model()` — thin orchestrator, calls block builders,
  assembles
- `build_*_block()` — one per Stan block (data, parameters, model, etc.)
- Sub-helpers like `build_re_parameters()`,
  `build_likelihood_two_stage()` — factored out where a block has
  substantial internal branching

### What’s shared vs. separate

Things that are now **shared** (were duplicated before): - The
skew_normal custom function - Data block core structure (with
stage-specific additions) - RE parameter declarations
(tau/z/nu/alpha_skew pattern) - RE prior emission logic - Generated
quantities patterns (skewness stats, mu_new, y_rep)

Things I kept **separate** deliberately: - The mixture likelihood code
for each stage (structurally distinct enough that abstracting further
would obscure more than it clarifies) - The multi-arm transformed
parameters block (complex loop logic that’s better as a readable unit) -
The robust one-stage likelihood per distribution family (each has
different approximation logic)

### What this doesn’t cover

The sketch handles `generate_stan_code_two_stage`,
`generate_stan_code_two_stage_null`, and `generate_stan_code_one_stage`.
It does **not** attempt to fold in the meta-regression generators or the
robma `spec_component_*` functions — those would be a separate phase if
you decide this pattern works.

### Next steps if you want to proceed

1.  Write **snapshot tests** for every model configuration you currently
    support, using the existing generators to produce the “expected”
    output
2.  Implement the new builders and verify they produce identical Stan
    code
3.  Swap call sites one at a time

Want me to start on any of those, or would you rather iterate on the
sketch first?

## User

So out of the 15 functions, how many will this sketch be able to handle?

## Databot

The user wants to know how many of the 15 generator functions the sketch
covers. Let me map it out based on what I learned from the exploration.

From the exploration, the 15 functions were:

1.  `generate_stan_code_two_stage()` — ✅ covered

2.  `generate_stan_code_two_stage_null()` — ✅ covered (via
    `null_model = TRUE`)

3.  `generate_stan_code_one_stage()` — ✅ covered

4.  `generate_stan_code_bias_corrected()` — ❌ not covered

5.  `generate_stan_code_mreg_two_stage()` — ❌ not covered

6.  `generate_stan_code_mreg_one_stage()` — ❌ not covered

7.  `generate_bayesian_egger_stan_binomial()` — ❌ not covered

8.  `generate_bayesian_egger_stan_generic()` — ❌ not covered 9-17.
    `spec_component_*` functions (robma) — ❌ not covered

So the sketch covers 3 out of the 15 functions, with the remaining 12
being either specialized variants for mediation regression, Egger
regression, or robustness analysis components.

**3 out of 15** — specifically:

| \# | Function | Covered? | How |
|----|----|----|----|
| 1 | `generate_stan_code_two_stage()` | Yes | `stage = "two_stage"` |
| 2 | `generate_stan_code_two_stage_null()` | Yes | `stage = "two_stage", null_model = TRUE` |
| 3 | `generate_stan_code_one_stage()` | Yes | `stage = "one_stage"` |
| 4 | `generate_stan_code_bias_corrected()` | No |  |
| 5 | `generate_stan_code_mreg_two_stage()` | No |  |
| 6 | `generate_stan_code_mreg_one_stage()` | No |  |
| 7 | `generate_bayesian_egger_stan_binomial()` | No |  |
| 8 | `generate_bayesian_egger_stan_generic()` | No |  |
| 9–17 | `spec_component_*` (9 robma functions) | No |  |

The meta-regression generators (5–6) would be the most natural next
addition since they share a lot of structure with the core
two-stage/one-stage generators — mainly adding a design matrix `X` and
coefficient vector `beta`. The Egger and robma spec functions are more
structurally distinct and would take more work to unify.

That said, even consolidating just those 3 functions is meaningful —
they’re the largest and most complex of the group, and the ones most
likely to grow as you add features.

## User

I’m not sure why the meta regression model can’t be used by this sketch?
It’s only the addition of a list of covariates, isn’t it?

## Databot

You’re right — let me check exactly what differs.

You’re basically right. The meta-regression generators are simpler than
the base ones (normal RE only, no robust, no multi-arm, no mixture) —
they just add `beta`, `X`, and per-coefficient priors. The differences
are:

1.  **Data block**: adds `int K` and `matrix[S, K] X` (or
    `matrix[N, K] X` for one-stage)
2.  **Parameters block**: adds `vector[K] beta`
3.  **Transformed parameters**: `theta = mu + X * beta + tau * z`
    instead of `mu + tau * z`
4.  **Model block**: per-coefficient beta priors (with the nice detail
    of collapsing to a single vectorised prior when all betas share the
    same prior)
5.  **Likelihood (one-stage)**: eta includes `X[i] * beta` term

These are all small, additive changes. The sketch could handle this with
a `use_meta_reg` flag (or passing a `moderators` config) that
conditionally injects `K`, `X`, `beta`, and the beta priors. Let me
update it.

Now let me update each block builder to handle the meta-regression
additions.

Now update the model block to handle beta priors:

Now add the `build_beta_priors` helper and update the likelihood to
include the `X * beta` term for one-stage meta-regression:

Now update the one-stage likelihood eta term and one-stage transformed
parameters to include `X[i] * beta`:

Now do the same for the y_rep eta in generated quantities:

Finally, for the two-stage meta-regression, the likelihood uses `theta`
which already incorporates `X * beta` via transformed parameters, so
that path needs no change. But the two-stage y_rep should use `theta[i]`
when meta-regression is enabled:

There’s also a minor detail — the two-stage meta-regression uses
`student_t_rng` for y_rep when t-likelihood is on (your current code
does this). Let me also handle the two-stage meta-reg likelihood, which
goes through `theta` and is already covered by
`build_likelihood_two_stage`. The only nuance is that the current mreg
two-stage code uses sampling syntax (`y ~ normal(theta, se)`) rather
than `target +=`, but that’s semantically equivalent — we can
standardise on `target +=`.

So with these changes, the sketch now covers **5 out of 15**:

| \# | Function | How |
|----|----|----|
| 1 | `generate_stan_code_two_stage()` | `stage = "two_stage"` |
| 2 | `generate_stan_code_two_stage_null()` | `stage = "two_stage", null_model = TRUE` |
| 3 | `generate_stan_code_one_stage()` | `stage = "one_stage"` |
| 4 | `generate_stan_code_mreg_two_stage()` | `stage = "two_stage", mreg_config = list(enabled = TRUE)` |
| 5 | `generate_stan_code_mreg_one_stage()` | `stage = "one_stage", mreg_config = list(enabled = TRUE)` |

You were right that it was a straightforward addition. The
meta-regression generators are essentially the base generators with
`beta`, `X`, and the `X * beta` term injected in a few predictable
places. The only slightly tricky bit is the per-coefficient beta prior
logic, which I factored into `build_beta_priors()`.

The updated sketch is in `R/generate_stan_model_sketch.R`. Want to look
at bringing the Egger or bias-corrected generators into the fold next,
or focus on getting tests in place for these 5?

## User

Yeah, try those 2.

## Databot

Let me read the Egger and bias-corrected generator functions to
understand their structure.

Now I have a clear picture. Let me assess these two:

1.  **`generate_stan_code_bias_corrected()`** — This is a very
    specialized model (bias correction with biased/unbiased mixture,
    inverse variance parameterization, known/unknown bias indicators).
    It has unique parameters (`B`, `p_bias`, `w_bias_raw`,
    `theta_raw_unbiased`, `theta_raw_biased`), unique transformed
    parameters (biased/unbiased theta), and a unique likelihood
    structure (mixture of biased/unbiased components). This is
    structurally very different from the core meta-analysis models.

2.  **Egger functions** — These are also fundamentally distinct, with
    their own parameter sets and likelihood structures. The binomial
    version computes transformed event rates from latent parameters,
    while the generic version uses a simpler regression approach. Both
    support multiplicative and additive heterogeneity variants, and they
    implement their own prior handling rather than relying on the
    standard `emit_prior_target` function.

These are both completely different model architectures — not variations
on the core pooling approach. Forcing them into a single
`generate_stan_model()` function would require so many conditional
branches that it would actually obscure the logic rather than clarify
it.

Having read both in detail, I’d actually push back on folding these into
the unified generator. Here’s why:

**`generate_stan_code_bias_corrected()`** is a fundamentally different
model structure. It has biased/unbiased mixture components with separate
theta vectors, inverse-variance parameterisation for tau, per-study bias
weights (`w_bias_raw`), and a known-vs-unknown-bias branching
likelihood. Almost nothing is shared with the core generators —
different data, different parameters, different transformed parameters,
different likelihood, different generated quantities.

**The Egger functions** are similarly distinct. The binomial version has
latent event rates (`p0`, `p1`), computes `sigma` from the latent rates
in transformed parameters, fits a hierarchical model for study-specific
log ORs (`delta ~ normal(d, tau)`), and then runs an Egger regression on
top. The generic version is a simple weighted regression. Both use their
own local `prior_to_stan()` helper rather than `emit_prior_target()`.
And both branch on multiplicative vs additive heterogeneity.

The core generators that the sketch unifies all share the same
fundamental structure: **pool study effects with optional random
effects, optionally adjust for moderators**. The Egger and
bias-corrected models are different statistical models entirely — not
variations on the same skeleton.

If you forced them into `generate_stan_model()`, you’d end up with
something like:

``` r
if (ctx$model_family == "egger" && ctx$egger_type == "binomial" && 
    ctx$heterogeneity == "multiplicative") {
  # completely different data block
  # completely different parameters
  # completely different everything
}
```

At that point the “unified” function is just a dispatch table, and
you’ve traded readable self-contained functions for a large function
where unrelated models are interleaved through conditionals.

**What I’d recommend instead**: leave these as standalone generators,
but extract the small bits they do share:

1.  Replace the local `prior_to_stan()` in the Egger functions with
    `emit_prior_target()` (they do the same thing)
2.  The Egger generic model is simple enough that multiplicative vs
    additive could be one function with a single conditional on the
    heterogeneity variance term, rather than two large `glue` blocks
    that differ by one line

Those are targeted cleanups that reduce duplication without forcing
unlike models into the same framework.

Refactor the Egger functions to use `emit_prior_target()` and
consolidate the multiplicative/additive variants.

Write snapshot tests for the 5 generators the sketch covers before
swapping in the new code.

Look at the robma `spec_component_*` functions to see if any of those
share enough structure to consolidate.

## User

What about the selection models and pet_peese?

## Databot

Now I have a complete picture of all the robma spec_component functions.
Let me assess which ones could be unified.

Looking at the patterns:

**Selection weight models** (`spec_component_selection_weight`,
`spec_component_selection_weight_fe`,
`spec_component_selection_weight_h0`): - These three are very similar.
The RE vs FE difference is just whether tau exists and whether sigma_i =
sqrt(tau^2 + se^2) or just se\[i\]. The H0 variant fixes mu=0. The
weight function likelihood logic (the p-value binning, omega weights,
normalisation constant) is essentially identical across all three. These
could absolutely be unified into one function with `has_re` and
`null_model` flags.

- The PET-PEESE models are also candidates for consolidation—they differ
  mainly in whether mu is included and whether the bias correction uses
  inverse square root of sample size or inverse sample size.
- The Copas and bias-corrected models follow similar patterns where
  random effects versus fixed effects comes down to whether tau is
  estimated and how sigma is computed.
- All of these are spec_component functions that feed into the robma
  model averaging framework, which is structurally distinct from the
  core meta-analysis generators.
- I’m wondering if they could be unified into their own family with a
  dedicated `generate_robma_component()` dispatcher rather than folding
  them into the main generator.
- Looking at the actual structure across PET-PEESE and selection weight
  models, there’s a consistent pattern: they all take N, y, se as core
  data, add optional mu and tau parameters, then apply a bias term to
  the normal likelihood with either se\[i\] or sigma_i as the scale.
- The selection weight models layer a weight function on top of this
  foundation, so the question becomes whether this shared structure
  could be abstracted into a more unified approach.

Actually, these models are fundamentally different enough that they
shouldn’t merge into the main generator—the selection weight models have
a completely different likelihood structure with the weight function and
normalisation constant, and the Copas model introduces yet another
selection function entirely. But within each family, there’s clear
duplication I could collapse: the PET-PEESE variants could become one
function parameterized by `null_model` and `bias_covariate`, the
selection weight variants could be unified with `has_re` and
`null_model` flags, and the Copas models could follow a similar pattern.

Looking at the consolidation potential, I could reduce the 9 robma
spec_component functions down to 4 by collapsing these redundancies
within each family—that’s a meaningful simplification worth presenting
clearly.

These are interesting — they’re structurally different from the core
generators, but they share a lot **with each other**. The robma
components all follow a similar skeleton and vary along two axes:
`has_re` (RE vs FE) and `null_model` (mu estimated vs mu=0).

Here’s what I see after reading them all:

### PET-PEESE family (3 functions → 1)

`spec_component_pet_peese`, `spec_component_pet_peese_h0`, and
`spec_component_peese_h0` are nearly identical. The differences are:

|                  | mu        | bias covariate | has_re branching |
|------------------|-----------|----------------|------------------|
| `pet_peese` (H1) | estimated | `1/sqrt(n)`    | yes              |
| `pet_peese_h0`   | 0         | `1/sqrt(n)`    | yes              |
| `peese_h0`       | 0         | `1/n`          | yes              |

These could trivially collapse into one function with `null_model` and
`bias_transform = c("inv_sqrt_n", "inv_n")` arguments.

### Selection weight family (3 functions → 1)

`spec_component_selection_weight`, `spec_component_selection_weight_fe`,
and `spec_component_selection_weight_h0` share the same weight-function
likelihood logic almost verbatim — the p-value binning, omega weights,
and normalisation constant are copy-pasted across all three. The only
differences:

- Whether `tau` exists (RE vs FE)
- Whether `mu` is estimated or fixed at 0
- The variance term in the likelihood: `sigma_i = sqrt(tau^2 + se^2)` vs
  just `se[i]`
- The normalisation constant computation adjusts slightly for tau

One function with `has_re` + `null_model` would handle all three.

### Copas and bias-corrected

`spec_component_selection_copas` already handles RE/FE via `has_re` — no
change needed. `spec_component_bias_corrected` +
`spec_component_jung_h0` could merge into one with a `null_model` flag,
though the RE bias-corrected model (which delegates to
`generate_stan_code_bias_corrected()`) is substantially more complex
than the others.

### Summary

These wouldn’t go into `generate_stan_model()` — they’re different
statistical models, not variations on “pool with optional RE.” But they
could consolidate among themselves:

| Current | Proposed | Reduction |
|----|----|----|
| `pet_peese` + `pet_peese_h0` + `peese_h0` | `spec_component_pet_peese()` | 3 → 1 |
| `selection_weight` + `selection_weight_fe` + `selection_weight_h0` | `spec_component_selection_weight()` | 3 → 1 |
| `selection_copas` | unchanged | 1 → 1 |
| `bias_corrected` + `jung_h0` | `spec_component_bias_corrected()` | 2 → 1 |

That takes the 9 robma functions from 9 → 4. Combined with the earlier 5
→ 1 consolidation of the core generators, the total generator count goes
from **15 → 5** (plus the Egger functions, which I’d leave as-is with
minor cleanup).

Want me to sketch out the consolidated PET-PEESE and selection weight
functions? Those two are the most clear-cut wins.

## User

Do all of those. Thanks

## Databot

Let me re-read the current sketch file to see where to append, and
double-check one detail on the jung_h0 function.

Now let me write all three consolidated functions and append them to the
sketch file.

Now let me read the end of the file to know where to append:

The selection weight normalisation constant is the trickiest part — it
differs subtly across the three original functions depending on whether
tau is present and whether mu is 0. Let me look more carefully at the
exact expressions.

Now I can see the four variants of the normalisation constant clearly:

1.  **H1 + RE** (`selection_weight`):
    `Phi((z_bounds[k] * se[i] - mu) / sigma[i]) - Phi((z_bounds[k+1] * se[i] - mu) / sigma[i])`

2.  **H1 + FE** (`selection_weight_fe`):
    `Phi((z_bounds[k] * se[i] - mu) / se[i]) - Phi((z_bounds[k+1] * se[i] - mu) / se[i])`
    which simplifies to
    `Phi(z_bounds[k] - mu/se[i]) - Phi(z_bounds[k+1] - mu/se[i])`

3.  **H0 + RE** (`selection_weight_h0` with random effects):
    `Phi(z_bounds[k] * se[i] / sigma_i) - Phi(z_bounds[k+1] * se[i] / sigma_i)`
    using the local sigma_i rather than the vector

4.  **H0 + FE** (`selection_weight_h0` without random effects):
    `Phi(z_bounds[k]) - Phi(z_bounds[k+1])` — the simplest form with no
    mean or variance adjustments

The pattern shows how the normalization constant changes depending on
whether we’re testing under the null or alternative hypothesis and
whether random effects are included. I’m noticing a key difference in
how sigma is handled between the two variants—the H0+RE model computes
it locally in the loop as a derived quantity, while H1+RE defines it in
the transformed parameters block, which affects how the general
expression gets instantiated in each case.

For the consolidated version, I’ll standardize on using transformed
parameters for sigma when random effects are present, so the Phi
argument can follow a consistent pattern based on whether we have random
effects and whether it’s a null model.

Building the string expressions for the Phi arguments…

Now I’m realizing the mathematical equivalence—when there’s no random
effect and no null model, the generated expression simplifies to just
the z-bounds directly, which matches what the current code does. Rather
than creating separate helper functions, I’ll build the normalization
constant lines inline with string interpolation for clarity. I also need
to handle the likelihood computation, particularly how sigma gets
calculated in the random effects case.

Good — now I can see the exact normalisation constant expressions across
the 4 variants. The general form is
`Phi((z_bounds[k] * se[i] - mu) / scale)` where scale is `sigma[i]` (RE)
or `se[i]` (FE), and mu is 0 for null models. Let me fix the selection
weight function — I used placeholder helper functions that don’t exist.
I’ll also add the bias-corrected consolidation.

Now let me find the start and end of the selection weight function to
replace it entirely:

Now I’ll replace the selection weight function with the correct
normalisation constant expression, and append the bias-corrected
consolidation:

All three are in `R/generate_stan_model_sketch.R`. Here’s a summary of
how they map to the originals:

### `spec_component_pet_peese()` — replaces 3 functions

| Original call | Consolidated call |
|----|----|
| `spec_component_pet_peese(has_re=T/F)` | `spec_component_pet_peese(has_re=T/F)` |
| `spec_component_pet_peese_h0(has_re=T/F)` | `spec_component_pet_peese(has_re=T/F, null_model=TRUE)` |
| `spec_component_peese_h0(has_re=T/F)` | `spec_component_pet_peese(has_re=T/F, null_model=TRUE, bias_transform="inv_n")` |

The only structural difference between these three was: (1) whether mu
is estimated or 0, and (2) whether the bias covariate is `1/sqrt(n)` or
`1/n`. Both are now arguments.

### `spec_component_selection_weight()` — replaces 3 functions

| Original call | Consolidated call |
|----|----|
| `spec_component_selection_weight()` | `spec_component_selection_weight(has_re=TRUE)` |
| `spec_component_selection_weight_fe()` | `spec_component_selection_weight(has_re=FALSE)` |
| `spec_component_selection_weight_h0(has_re=T/F)` | `spec_component_selection_weight(has_re=T/F, null_model=TRUE)` |

The weight-function likelihood logic (p-value binning, omega assignment,
normalisation constant) was copy-pasted across all three. Now it’s
written once. The normalisation constant uses the general form
`Phi((z_bounds[k] * se[i] - mu) / scale)` which is mathematically
equivalent to the simplified forms in the originals.

### `spec_component_bias_corrected()` — replaces 2 functions

| Original call | Consolidated call |
|----|----|
| `spec_component_bias_corrected(has_re=T/F)` | `spec_component_bias_corrected(has_re=T/F)` |
| `spec_component_jung_h0(has_re=T/F)` | `spec_component_bias_corrected(has_re=T/F, null_model=TRUE)` |

The RE + H1 case still delegates to
`generate_stan_code_bias_corrected()` because it’s substantially more
complex (inverse-variance parameterisation, per-study bias weights). The
other three variants (FE+H1, RE+H0, FE+H0) are unified.

### Overall tally

With `spec_component_selection_copas` unchanged (it already handles
RE/FE via `has_re`), plus the Egger functions left as-is, the full
picture is now:

| Group            | Before | After                     |
|------------------|--------|---------------------------|
| Core generators  | 5      | 1 (`generate_stan_model`) |
| Robma components | 9      | 4                         |
| Egger            | 2      | 2 (cleanup only)          |
| **Total**        | **16** | **7**                     |

Write snapshot tests capturing the Stan code output of each original
function before swapping in the consolidated versions.

Refactor the Egger functions to use `emit_prior_target()` and
consolidate the multiplicative/additive variants.

Review the call sites in `robma_bayesma.R` to plan how to update them to
use the new consolidated function signatures.

## User

Are the PET-PEESE models using SE & SE^2 or are they using 1/sqrt(n) &
1/n?

## Databot

They’re using `1/sqrt(n)` and `1/n`. Looking at the code:

- `spec_component_pet_peese` computes
  `inv_sqrt_n[i] = 1.0 / sqrt(n_total[i])` in transformed data, then
  regresses on that
- `spec_component_peese_h0` computes `inv_n[i] = 1.0 / n_total[i]`

This is a sample-size-based approach rather than the traditional
PET-PEESE formulation which uses SE and SE^2 as the bias covariates. The
two approaches are related but not identical — `1/sqrt(n)` is
proportional to SE only under specific assumptions (e.g., equal
variances across groups, equal group sizes).

Was that intentional, or would you prefer the standard SE-based version?
The SE-based version would be simpler since `se` is already in the data
block — no need for `n_total` or a transformed data block at all.

## User

It was changed as per Stanley (2017) as a variance stabilised predictor.
Is this wrong?

## Databot

It’s not wrong — it’s a defensible choice. Stanley & Doucouliagos (2014,
2017) do argue that `1/sqrt(n)` is a variance-stabilised predictor for
the PET test, and `1/n` for PEESE, on the grounds that SE itself is a
noisy estimate (especially with small samples) and using it as both the
weight and the predictor creates a mechanical correlation. Using sample
size avoids that endogeneity problem.

That said, the literature isn’t fully settled on this. The original
PET-PEESE formulation (Stanley & Doucouliagos, 2012) used SE and SE^2,
and many implementations still do (e.g.,
[`metafor::regtest`](https://wviechtb.github.io/metafor/reference/regtest.html),
the `RoBMA` package). The sample-size-based version assumes a fairly
direct relationship between `1/sqrt(n)` and SE, which holds cleanly for
standardised mean differences but is more approximate for log odds
ratios or log risk ratios where SE depends on event rates as well as
sample size.

Since your package supports binomial and Poisson likelihoods alongside
Gaussian, it might be worth noting this assumption in the documentation
— that the `1/sqrt(n)` predictor is most directly motivated for
continuous outcomes. But for the robma context where these are component
models being averaged, the approximation is likely fine.

Short answer: you’re on solid ground citing Stanley (2017). Just be
transparent about the choice.

## User

What about something like this:
<https://jepusto.com/posts/PET-PEESE-performance/>

## Databot

This is Pustejovsky’s 2017 blog post comparing traditional PET-PEESE
(using SE and SE²) vs sample-size variants (SPET/SPEESE using 1/√n and
1/n). The key findings support what bayesma is doing:

1.  Traditional PET/PEESE using SE/SE² have wildly incorrect Type-I
    error rates because d and V are mechanically correlated

2.  Sample-size variants (SPET/SPEESE using 1/√n and 1/n) maintain
    correct Type-I error rates

3.  SPET/SPEESE are exactly unbiased in the absence of publication bias

4.  Under publication bias, SPEESE emerges as the recommended default
    even though no single estimator consistently outperforms the others

So Pustejovsky’s analysis directly validates the sample-size approach
that bayesma is using.

This actually reinforces your choice. Pustejovsky’s simulations show
exactly why you’d want the sample-size variants:

1.  **Traditional PET/PEESE** (using SE and SE²) have badly inflated
    Type-I error rates even with no publication bias, because `d` and
    `V` are mechanically correlated — you need `d` to compute `V`
2.  **SPET/SPEESE** (using `1/sqrt(n)` and `1/n`) maintain correct
    Type-I error rates and are exactly unbiased when there’s no
    publication bias
3.  Under actual publication bias, neither variant consistently
    dominates, but Pustejovsky’s “gun to my head” recommendation is
    SPEESE (`1/n`)

So your implementation aligns with both Stanley (2017) and Pustejovsky’s
independent simulation evidence. The one thing worth noting from the
post is his observation that the `1/sqrt(n)` predictor is specifically
motivated for standardized mean differences from two-group designs —
which circles back to the caveat I mentioned about log ORs and log RRs,
where the SE-to-sample-size relationship is more complex. But within the
robma model-averaging framework, where these are component models
weighted by evidence, this is a reasonable approximation.

If anything, you might consider documenting both the Stanley (2017) and
Pustejovsky (2017) references in your PET-PEESE function documentation
to explain the choice.

## User

I thhink RoBMA does a z-transformation first to avaoid the problem. Here
is the package site:

<https://fbartos.github.io/RoBMA/>

## Databot

Your Posit AI credits have been depleted.
