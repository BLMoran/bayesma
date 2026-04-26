# Alternative random-effects distributions

## Introduction

The standard random-effects model places a Gaussian distribution on the
study-level true effects:

``` math

\theta_i \sim \mathcal{N}(\mu, \tau^2)
```

This assumption is convenient but not always appropriate. When the
distribution of true effects is heavy-tailed, asymmetric, or bimodal —
for example, because a minority of studies target fundamentally
different populations or use qualitatively different implementations —
the Gaussian assumption can distort the pooled estimate and understate
predictive uncertainty.

**bayesma** supports three alternatives: Student-$`t`$, skew-normal, and
two-component mixture.

## Gaussian (default)

``` math

\theta_i \sim \mathcal{N}(\mu, \tau^2)
```

The Gaussian is the natural starting point. It implies symmetric
heterogeneity and exponentially light tails: very large deviations from
$`\mu`$ are treated as extremely unlikely. Use it as the default; switch
to an alternative only when there is substantive or data-driven reason.

## Student-$`t`$

``` math

\theta_i \sim t_\nu(\mu, \tau^2)
```

The Student-$`t`$ distribution adds a degrees-of-freedom parameter
$`\nu`$ that controls tail heaviness. As $`\nu \to \infty`$ it converges
to the Gaussian; at $`\nu = 3`$ the tails are substantially heavier.
Heavy tails mean that large study-level deviations from $`\mu`$ are more
plausible, reducing the influence of outlier studies on the pooled
estimate.

A prior on $`\nu`$ must be specified. **bayesma** uses

``` math

\nu \sim \text{Gamma}(2, 0.1)
```

by default, which assigns most mass to $`\nu \in (3, 30)`$ while
allowing both very heavy and near-Gaussian tails. This can be overridden
via `nu_prior`.

The Student-$`t`$ random-effects model is appropriate when:

- a small number of studies have effects far from the bulk
- the source of those deviations is unknown and possibly real (not
  artefact)
- the analyst wants to downweight outliers without excluding them

## Skew-normal

``` math

\theta_i \sim \text{SN}(\mu, \tau, \alpha)
```

The skew-normal distribution generalises the Gaussian with a shape
parameter $`\alpha`$ that controls the direction and degree of
asymmetry. When $`\alpha = 0`$ it reduces to
$`\mathcal{N}(\mu, \tau^2)`$. Positive $`\alpha`$ yields a right-skewed
distribution of effects; negative $`\alpha`$ a left-skewed distribution.

The skew-normal is appropriate when:

- there is a directional floor or ceiling effect (e.g., effects cannot
  be negative)
- a literature has a mixture of small and large positive effects but few
  negative effects
- theoretical reasons exist to expect asymmetric heterogeneity

A prior on the shape parameter is required. **bayesma** uses
$`\alpha \sim \mathcal{N}(0, 1)`$ by default.

## Two-component mixture

``` math

\theta_i \sim \pi \cdot \mathcal{N}(\mu_1, \tau_1^2) + (1 - \pi) \cdot \mathcal{N}(\mu_2, \tau_2^2)
```

The mixture model allows the distribution of true effects to be bimodal.
The mixing weight $`\pi`$ estimates the proportion of studies belonging
to each component. This is the most flexible option and the hardest to
fit reliably.

The mixture is appropriate when:

- substantive theory predicts two distinct populations of studies (e.g.,
  short vs long follow-up; low vs high dose)
- the forest plot shows a clear gap or bimodal distribution
- the analyst suspects a qualitative moderator that was not measured

**Caution**: the two-component mixture requires adequate $`k`$
(typically $`k \geq 20`$) and informative priors on $`\pi`$ and the
component parameters. With small $`k`$, the components are poorly
identified and the posterior is strongly prior-dependent. A prior
sensitivity analysis is essential.

Priors:

``` math

\pi \sim \text{Beta}(2, 2), \qquad
\mu_j \sim \mathcal{N}(0, 1), \qquad
\tau_j \sim \text{Half-Cauchy}(0, 0.5)
```

## Choosing a distribution

See [Assessment of RE
Distributions](https://blmoran.github.io/bayesma/articles/re-distribution-assessment.md)
for a model comparison workflow. As a heuristic:

| Situation                                                 | Distribution  |
|-----------------------------------------------------------|---------------|
| No strong reason to depart from Gaussian                  | Normal        |
| One or two studies with extreme effects                   | Student-$`t`$ |
| Effects are concentrated near zero with a long right tail | Skew-normal   |
| Forest plot shows two clusters of effects                 | Mixture       |

## Effect on prediction

The choice of RE distribution affects not only the pooled estimate but
also the predictive distribution for a new study. Heavier-tailed
distributions produce wider prediction intervals. This is the primary
reason to prefer heavier tails when outliers are present: the Gaussian
understates how variable a new study’s result might be.
