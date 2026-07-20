# Add manydist dissimilarities to a recipe

\`step_mdist()\` is a \[recipes::recipe()\] step that replaces selected
predictors by a distance-based representation computed with \[mdist()\].
It is designed for distance-based learning workflows, especially
nearest-neighbour prediction and clustering models that operate on
dissimilarity matrices.

## Usage

``` r
step_mdist(
  recipe,
  ...,
  role = "predictor",
  trained = FALSE,
  output = "distance_to_training",
  preset = "custom",
  method_cat = "tot_var_dist",
  method_num = NULL,
  commensurable = FALSE,
  ncomp = NULL,
  threshold = NULL,
  columns = NULL,
  train_predictors = NULL,
  preprocessor = NULL,
  skip = FALSE,
  id = recipes::rand_id("mdist")
)
```

## Arguments

- recipe:

  A recipe object.

- ...:

  Selector(s) for the predictor columns to be used in \[mdist()\]. These
  are passed to \[recipes::recipes_eval_select()\] during preparation.

- role:

  Role for the new distance columns. The default is \`"predictor"\`.

- trained:

  Logical for recipes internals. Do not set manually.

- output:

  Character string specifying the type of distance output.
  \`"distance_to_training"\` returns distances from the baked data to
  the training observations and is the usual choice for prediction
  workflows. \`"pairwise"\` returns the within-training pairwise
  dissimilarity matrix and is intended for training-only distance-based
  clustering workflows.

- preset:

  Character string specifying the distance preset passed to \[mdist()\].
  Available values include \`"custom"\`, \`"gower"\`,
  \`"unbiased_dependent"\`, \`"u_dep"\`, \`"u_indep"\`, \`"u_mix"\`,
  \`"hl"\`, \`"gudmm"\`, \`"dkss"\`, \`"mod_gower"\`, and
  \`"euclidean"\`. Preset parameters are normally fixed by the selected
  preset. The exception is \`method_num\` for \`preset = "euclidean"\`
  when all predictors selected by the step are numeric.

- method_cat:

  Character string specifying the categorical-variable dissimilarity
  passed to \[mdist()\] when \`preset = "custom"\`. Common values
  include \`"matching"\` and \`"tvd"\`. Use \[all_dist_method_specs()\]
  to inspect available methods.

- method_num:

  Character string or \`NULL\` specifying numerical-variable
  preprocessing. Supported values are \`"none"\` for no preprocessing,
  \`"std"\` for standard-deviation scaling, \`"range"\` for range
  scaling, \`"robust"\` for inter-quartile-range-based scaling, and
  \`"pc_scores"\` for principal-component score scaling. With \`preset =
  "custom"\`, \`NULL\` defaults to \`"none"\`. With \`preset =
  "euclidean"\`, the default is \`"std"\`; when all predictors selected
  by the step are numeric, an explicitly supplied value can override
  that default.

- commensurable:

  Logical. If \`TRUE\`, dissimilarities are scaled so that the average
  contribution of each variable to the overall distance is equal to 1,
  when supported by the selected distance specification.

- ncomp:

  Integer or \`NULL\`. Number of principal components to retain when
  \`method_num = "pc_scores"\`. If \`NULL\`, all available components
  are used unless \`threshold\` is supplied and supported by the
  underlying method.

- threshold:

  Numeric or \`NULL\`. Optional cumulative variance threshold used when
  \`method_num = "pc_scores"\`.

- columns:

  Names of columns selected at prep time. Used internally by recipes.

- train_predictors:

  Training predictors stored at prep time. Used internally by recipes to
  compute distances from new observations to the training observations.

- preprocessor:

  Internal fitted manydist preprocessor.

- skip:

  Logical. Standard recipes argument indicating whether the step should
  be skipped when baking new data.

- id:

  Character string. Unique step identifier.

## Value

An updated recipe with a manydist step.

## Details

The step can produce either distances from new observations to the
training observations, or the within-training pairwise dissimilarity
matrix. The former is the usual choice for supervised prediction
workflows; the latter is useful for distance-based clustering workflows
fitted on the training data.

During \[recipes::prep()\], \`step_mdist()\` stores the selected
training predictors and fits the internal manydist preprocessor. During
\[recipes::bake()\], the selected predictors are removed and replaced by
distance columns named \`dist_1\`, \`dist_2\`, and so on.

With \`output = "distance_to_training"\`, baking the training data
returns the training pairwise distances, while baking new data returns
distances from each new observation to each training observation. This
rectangular representation is suitable for nearest-neighbour prediction
models.

With \`output = "pairwise"\`, the step returns the within-training
pairwise dissimilarity matrix. Baking genuinely new data is not
supported in this mode, because the output is intended for training-only
clustering workflows such as \[pam_dist()\] or \[spectral_dist()\].

## See also

\[mdist()\], \[nearest_neighbor_dist()\], \[pam_dist()\],
\[spectral_dist()\], \[all_dist_method_specs()\]

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

  # Distance-to-training representation for prediction workflows
  rec <- recipes::recipe(species ~ ., data = penguins_small) |>
    step_mdist(
      recipes::all_predictors(),
      preset = "gower",
      output = "distance_to_training"
    )

  rec_prep <- recipes::prep(rec, training = penguins_small)
  baked <- recipes::bake(rec_prep, new_data = penguins_small)

  baked |> dplyr::slice_head(n=5)

  # Pairwise representation for clustering workflows
  rec_pairwise <- recipes::recipe(~ ., data = penguins_small) |>
    step_mdist(
      recipes::all_predictors(),
      preset = "gower",
      output = "pairwise"
    )

  rec_pairwise_prep <- recipes::prep(rec_pairwise, training = penguins_small)
  pairwise_dist <- recipes::bake(rec_pairwise_prep, new_data = penguins_small)

  pairwise_dist
}
#> Warning: When `preset` is not 'custom', distance-related arguments are ignored: `commensurable`, `method_num`. Set `preset = 'custom'` to specify them manually.
#> Warning: When `preset` is not 'custom', distance-related arguments are ignored: `commensurable`, `method_num`. Set `preset = 'custom'` to specify them manually.
#> # A tibble: 333 × 333
#>    dist_1 dist_2 dist_3 dist_4 dist_5 dist_6 dist_7 dist_8 dist_9 dist_10
#>     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
#>  1 0      0.181  0.215  0.206  0.0591 0.164  0.0864 0.196  0.0713  0.131 
#>  2 0.181  0      0.0580 0.0777 0.214  0.0290 0.238  0.0452 0.224   0.284 
#>  3 0.215  0.0580 0      0.0536 0.220  0.0595 0.232  0.0444 0.238   0.278 
#>  4 0.206  0.0777 0.0536 0      0.194  0.0729 0.214  0.0883 0.204   0.234 
#>  5 0.0591 0.214  0.220  0.194  0      0.215  0.0703 0.240  0.0222  0.0821
#>  6 0.164  0.0290 0.0595 0.0729 0.215  0      0.251  0.0341 0.233   0.293 
#>  7 0.0864 0.238  0.232  0.214  0.0703 0.251  0      0.277  0.0747  0.0676
#>  8 0.196  0.0452 0.0444 0.0883 0.240  0.0341 0.277  0      0.263   0.323 
#>  9 0.0713 0.224  0.238  0.204  0.0222 0.233  0.0747 0.263  0       0.0632
#> 10 0.131  0.284  0.278  0.234  0.0821 0.293  0.0676 0.323  0.0632  0     
#> # ℹ 323 more rows
#> # ℹ 323 more variables: dist_11 <dbl>, dist_12 <dbl>, dist_13 <dbl>,
#> #   dist_14 <dbl>, dist_15 <dbl>, dist_16 <dbl>, dist_17 <dbl>, dist_18 <dbl>,
#> #   dist_19 <dbl>, dist_20 <dbl>, dist_21 <dbl>, dist_22 <dbl>, dist_23 <dbl>,
#> #   dist_24 <dbl>, dist_25 <dbl>, dist_26 <dbl>, dist_27 <dbl>, dist_28 <dbl>,
#> #   dist_29 <dbl>, dist_30 <dbl>, dist_31 <dbl>, dist_32 <dbl>, dist_33 <dbl>,
#> #   dist_34 <dbl>, dist_35 <dbl>, dist_36 <dbl>, dist_37 <dbl>, …
```
