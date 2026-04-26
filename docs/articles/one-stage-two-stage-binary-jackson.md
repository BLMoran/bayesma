# One-stage and two-stage models for binary outcomes

## Introduction

Binary outcome meta-analyses can be performed using either **two-stage**
or **one-stage** approaches.

The two-stage method first computes study-level effect sizes, then
performs meta-analysis.

The one-stage method models the original likelihood directly.

## Two-stage model

First compute the log odds ratio:

\[ y_i = (OR_i) \]

Then perform random-effects meta-analysis:

\[ y_i (\_i, s_i^2) \]

\[ \_i (, ^2) \]

## One-stage model (Jackson et al.)

For study (i) and arm (j):

\[ e\_{ij} (n\_{ij}, \_{ij}) \]

with logit link

\[ (\_{ij}) = *i + j+ z*{ij}\_i \]

where

- (\_i) = baseline risk
- () = treatment effect
- (\_i) = heterogeneity term

### Assumptions

- binomial likelihood
- logit link
- study-specific baseline risks
