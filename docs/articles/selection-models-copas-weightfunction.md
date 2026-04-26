# Selection models for publication bias

## Introduction

Selection models adjust for **publication bias** by modelling the
probability that a study is observed.

Two common approaches are:

- Copas selection model
- Weight-function models (Vevea & Hedges)

## Copas model

Studies are observed if

\[ Z_i \> 0 \]

where

\[ Z_i = \_0 + + \_i \]

This links study precision with probability of publication.

## Weight-function models

Selection probabilities depend on p-values.

\[ w(p_i) =
``` math
\begin{cases}
\omega_1 & 0 < p_i \le a_1 \\
\omega_j & a_{j-1} < p_i \le a_j
\end{cases}
```

\]

These weights modify the likelihood.
