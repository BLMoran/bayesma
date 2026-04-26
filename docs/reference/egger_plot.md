# Plot method for bayesma_egger

Creates diagnostic plots for the Bayesian Egger's test, including a
funnel plot with posterior latent SEs and a posterior distribution plot
for the slope parameter.

## Usage

``` r
egger_plot(
  x,
  type = c("funnel", "beta_posterior", "both"),
  show_observed = TRUE,
  contour_lines = TRUE,
  contour_alpha = 0.15,
  color_observed = "grey50",
  color_latent = "#D55E00",
  color_pooled = "blue",
  color_contour = "grey50",
  point_size = 3,
  ...
)
```

## Arguments

- x:

  A bayesma_egger object

- type:

  Character. Type of plot: "funnel" (default), "beta_posterior", or
  "both".

- show_observed:

  Logical. Show observed SEs alongside latent SEs? Default TRUE.

- contour_lines:

  Logical. Add significance contour regions? Default TRUE.

- contour_alpha:

  Numeric. Alpha transparency for contour shading. Default 0.15.

- color_observed:

  Color for observed SE points. Default "grey50".

- color_latent:

  Color for latent SE points. Default "#D55E00".

- color_pooled:

  Color for pooled estimate line. Default "blue".

- color_contour:

  Color for contour fills. Default "grey50".

- point_size:

  Numeric. Size of study points. Default 3.

- ...:

  Additional arguments (currently unused)
