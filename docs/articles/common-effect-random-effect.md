# Common-effect and random-effects models

## Introduction

The most widely used meta-analysis models are the **common-effect** and
**random-effects** models.

The common-effect model assumes all studies estimate the same true
effect, while the random-effects model allows the true effect to vary
across studies.

These models form the foundation for most extensions implemented in
**bayesma**.

## Common-effect model

Let

- (y_i) = effect estimate from study (i)
- (s_i) = standard error of the estimate

The sampling model is

\[ y_i (, s_i^2) \]

where () is the shared effect.

A prior is placed on the effect:

\[ (0, \_^2) \]

### Interpretation

All studies are assumed to estimate the same underlying effect.

## Random-effects model

The random-effects model introduces heterogeneity between studies.

\[ y_i (\_i, s_i^2) \]

\[ \_i (, ^2) \]

where

- () = overall effect
- () = between-study heterogeneity

### Assumptions

- study effects are exchangeable
- heterogeneity follows a normal distribution
- sampling errors are known
