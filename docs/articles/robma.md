# RoBMA

## Introduction

RoBMA (Robust Bayesian Meta-Analysis) performs model averaging across
multiple bias-adjusted meta-analysis models.

Instead of choosing a single model, RoBMA averages posterior
distributions weighted by model probabilities.

\[ p(\|y) = \_k p(\|y, M_k) p(M_k\|y) \]

This approach integrates uncertainty about publication bias and model
specification.
