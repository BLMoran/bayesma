# Bias-corrected meta-analysis model

## Introduction

The bias-corrected model accounts for systematic bias in some studies.

This approach extends the model proposed by Verde and later developed by
Jung.

## Model

The observed effect is

\[ y_i (\_i^B, s_i^2) \]

where

\[ \_i^B = (1 - I_i)\_i + I_i(\_i + \_i) \]

- (I_i) indicates whether study (i) is biased
- (\_i) represents the bias magnitude

## Hierarchical structure

\[ \_i (, ^2) \]

\[ *i (*, \_^2) \]

### Interpretation

The model separates **true effects** from **bias effects**.
