# Distance-based nearest-neighbour model

\`nearest_neighbor_dist()\` defines a parsnip model specification for
nearest-neighbour prediction using precomputed distance representations,
such as those produced by \[step_mdist()\]. It can be used for
classification or regression workflows in which predictors have first
been transformed into distances to the training observations.

## Usage

``` r
fit_knn_dist(x, y, k = 5, dist_fun = NULL, dist_args = list())

predict_knn_dist_class(object, new_data, type = c("class", "prob"))

predict_knn_dist_prob(object, new_data, ...)

predict_knn_dist_reg(object, new_data)

nearest_neighbor_dist(mode = "classification", neighbors = NULL)

# S3 method for class 'nearest_neighbor_dist'
tunable(x, ...)
```

## Arguments

- x:

  A \`nearest_neighbor_dist\` model specification.

- y:

  A vector of outcomes.

- k:

  Number of neighbours.

- dist_fun:

  Optional distance function taking arguments \`x\` and \`new_data\`. If
  \`NULL\`, inputs are assumed to already be distance matrices.

- dist_args:

  A named list of additional arguments passed to \`dist_fun\`.

- object:

  A fitted \`knn_dist\` object created by \`fit_knn_dist()\`.

- new_data:

  A data frame or matrix of precomputed distances, or raw predictors
  when \`dist_fun\` is supplied.

- type:

  Prediction type. For classification, \`"class"\` returns class
  predictions and \`"prob"\` returns class probabilities.

- ...:

  Additional arguments passed to lower-level methods or currently not
  used.

- mode:

  Character string specifying the model mode. Available values are
  \`"classification"\` and \`"regression"\`.

- neighbors:

  Number of neighbours. This can be an integer or a tunable parameter,
  for example \`tune::tune()\`.

## Value

A parsnip model specification of class \`"nearest_neighbor_dist"\`.

## Details

This model is intended to be used together with \[step_mdist()\] in a
\[recipes::recipe()\]. The recipe creates distance columns named
\`dist_1\`, \`dist_2\`, and so on; \`nearest_neighbor_dist()\` then
applies k-nearest neighbours to that distance representation.

The model uses a manydist-specific parsnip engine. In the usual
workflow, distances are computed by \[step_mdist()\] with \`output =
"distance_to_training"\`. The resulting distance columns are then passed
to \`nearest_neighbor_dist()\`.

Lower-level engine functions such as \`fit_knn_dist()\` and
\`predict_knn_dist\_\*()\` are exported for parsnip registration, but
users normally do not need to call them directly.

## See also

\[step_mdist()\], \[mdist()\]

## Examples

``` r
if (requireNamespace("palmerpenguins", quietly = TRUE)) {
  data("penguins", package = "palmerpenguins")

  penguins_small <- palmerpenguins::penguins |>
    dplyr::select(
      species, bill_length_mm, bill_depth_mm, flipper_length_mm,
      body_mass_g, island, sex
    ) |>
    tidyr::drop_na()

  set.seed(123)
  penguin_split <- rsample::initial_split(
    penguins_small,
    prop = 0.75,
    strata = species
  )

  penguin_train <- rsample::training(penguin_split)
  penguin_test  <- rsample::testing(penguin_split)

  rec <- recipes::recipe(species ~ ., data = penguin_train) |>
    step_mdist(
      recipes::all_predictors(),
      preset = "gower",
      output = "distance_to_training"
    )

  spec <- nearest_neighbor_dist(
    mode = "classification",
    neighbors = 5
  ) |>
    parsnip::set_engine("manydist")

  wf <- workflows::workflow() |>
    workflows::add_recipe(rec) |>
    workflows::add_model(spec)

  fit <- workflows::fit(wf, data = penguin_train)

  predict(fit, new_data = penguin_test) |>
    dplyr::slice_head(n = 5)
}
#> # A tibble: 5 × 1
#>   .pred_class
#>   <fct>      
#> 1 Adelie     
#> 2 Adelie     
#> 3 Adelie     
#> 4 Adelie     
#> 5 Adelie     
```
