# Spectral clustering specification based on manydist dissimilarities

Spectral clustering specification based on manydist dissimilarities

## Usage

``` r
spectral_dist(num_clusters = NULL, sigma = NULL, nstart = 50)
```

## Arguments

- num_clusters:

  Number of clusters.

- sigma:

  Optional bandwidth for the Gaussian affinity. If \`NULL\`, the median
  pairwise distance is used.

- nstart:

  Number of random starts for k-means.

## Value

A \`spectral_dist_spec\` object.

## Examples

``` r


if (FALSE) { # \dontrun{
library(manydist)
library(palmerpenguins)
library(recipes)
library(generics)

data <- penguins |>
  dplyr::select(-species) |>
  tidyr::drop_na()

rec <- recipes::recipe(~ ., data = data) |>
  step_mdist(all_predictors(), preset = "gower", output = "pairwise")

spec <- spectral_dist(num_clusters = 3)

fit_obj <- generics::fit(spec, recipe = rec, data = data)

print(fit_obj)
predict(fit_obj)
predict(fit_obj, type = "embed")
} # }
```
