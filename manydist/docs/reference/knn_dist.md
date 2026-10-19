# Direct nearest-neighbour prediction

Predict outcomes from precomputed test-to-training distances, or compute
those distances from tabular training and test predictors before
predicting.

## Usage

``` r
knn_dist(
  x,
  y = NULL,
  new_data = NULL,
  k = 5,
  dist_args = list(),
  type = c("class", "prob", "numeric"),
  response = NULL
)
```

## Arguments

- x:

  Tabular training predictors when `new_data` is supplied. Otherwise, a
  numeric test-to-training distance matrix, data frame, or `MDist`
  object.

- y:

  Training outcomes, in the same order as the distance columns or
  training rows. Use a factor for classification and numeric outcomes
  for regression. Missing outcomes are not supported. May be omitted
  when `response` selects an outcome column from tabular training data.

- new_data:

  Optional tabular test predictors. When supplied,
  [`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md)
  computes their distances to the training data in `x`. Leave `NULL`
  when `x` already contains test-to-training distances.

- k:

  Integer number of neighbours, between one and the number of training
  observations.

- dist_args:

  Named list of additional arguments passed to
  [`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md),
  such as `list(preset = "gower")`. Used only with tabular input.

- type:

  Prediction type: `"class"` for class labels, `"prob"` for class
  probabilities, or `"numeric"` for regression predictions.

- response:

  Optional outcome column name (a string or unquoted column name) in
  tabular `x`. Supply either `y` or `response`, not both. The selected
  column is excluded from the training and test predictors before
  distance computation. This selects the prediction outcome, not the
  `response` argument of
  [`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md).

## Value

A factor for `type = "class"`, a data frame with one column per outcome
level for `type = "prob"`, or a numeric vector for `type = "numeric"`.
Each result has one row or element per test observation.

## Details

Each distance row corresponds to a test observation and each column to a
training observation. Distances must be finite and nonnegative.
Classification uses an unweighted majority vote; regression uses the
mean of the neighbours' outcomes. Distance ties are resolved in training
order and class-vote ties in factor-level order. Class probabilities are
neighbour class proportions.

No training-to-training matrix needs to be supplied. For tabular input,
[`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md)
uses `x` as training data and `new_data` as test data. The function
computes distances on each call; it does not return a fitted model. For
fitted tidymodels workflows, use
[`nearest_neighbor_dist()`](https://alfonsoiodicede.github.io/manydist_package/reference/nearest_neighbor_dist.md)
with
[`step_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/step_mdist.md)
instead.

## See also

[`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md),
[`nearest_neighbor_dist()`](https://alfonsoiodicede.github.io/manydist_package/reference/nearest_neighbor_dist.md),
[`step_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/step_mdist.md)

## Examples

``` r
train <- data.frame(value = c(0, 2, 8, 10))
test <- data.frame(value = c(1, 9))
labels <- factor(c("low", "low", "high", "high"))
knn_dist(train, labels, new_data = test, k = 1,
         dist_args = list(preset = "gower"))
#> [1] low  high
#> Levels: high low
train$label <- labels
knn_dist(train, response = "label", new_data = test, k = 1,
         dist_args = list(preset = "gower"))
#> [1] low  high
#> Levels: high low
d <- mdist(train["value"], new_data = test, preset = "gower")
knn_dist(d, labels, k = 1)
#> [1] low  high
#> Levels: high low
knn_dist(d, labels, k = 2, type = "prob")
#>   high low
#> 1    0   1
#> 2    1   0
knn_dist(d, c(0, 2, 8, 10), k = 2, type = "numeric")
#> [1] 1 9
```
