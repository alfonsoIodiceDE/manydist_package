# Mixed-type dissimilarities for distance-based learning

Computes a dissimilarity object for numerical, categorical, or
mixed-type data. The function combines continuous and categorical
components according to either a predefined \`preset\` or a user-defined
custom specification.

## Usage

``` r
mdist(
  x,
  new_data = NULL,
  response = NULL,
  method_cat = "tvd",
  method_num = "std",
  commensurable = TRUE,
  ncomp = NULL,
  threshold = NULL,
  preset = "custom",
  interaction = FALSE,
  prop_nn = 0.1,
  score = "ba",
  decision = "prior_corrected",
  gower_average = TRUE
)
```

## Arguments

- x:

  A data frame or matrix containing the training observations. Columns
  can be numeric, factors, or a mixture of both.

- new_data:

  Optional data frame or matrix containing new observations. If
  supplied, distances are computed from rows of \`new_data\` to rows of
  \`x\`, producing a rectangular test-to-training dissimilarity matrix.

- response:

  Optional response variable used for response-aware categorical
  dissimilarities. It can be supplied as an unquoted column name or as a
  character string. The response column is removed from the predictors
  before computing distances.

- method_cat:

  Character string specifying the categorical-variable dissimilarity
  used when \`preset = "custom"\`. Common values include \`"matching"\`
  and \`"tvd"\`. Use \[all_dist_method_specs()\] to inspect available
  methods.

- method_num:

  Character string specifying numerical-variable preprocessing.
  Supported values are \`"none"\` for no preprocessing, \`"std"\` for
  standard-deviation scaling, \`"range"\` for range scaling,
  \`"robust"\` for inter-quartile-range-based scaling, and
  \`"pc_scores"\` for principal-component score scaling. This argument
  is used directly when \`preset = "custom"\`. When \`preset =
  "euclidean"\` and all predictors are numeric, it can also override the
  default \`"std"\` preprocessing; in particular, \`"none"\` gives
  ordinary Euclidean distance on the original variables.

- commensurable:

  Logical. If \`TRUE\`, dissimilarities are scaled so that the average
  contribution of each variable to the overall distance is equal to 1.

- ncomp:

  Integer or \`NULL\`. Number of principal components to retain when
  \`method_num = "pc_scores"\`. If \`NULL\`, all available components
  are used unless \`threshold\` is supplied and supported by the
  underlying method.

- threshold:

  Numeric or \`NULL\`. Optional cumulative variance threshold used when
  \`method_num = "pc_scores"\`.

- preset:

  Character string specifying a predefined distance specification.
  Available values include \`"custom"\`, \`"gower"\`,
  \`"unbiased_dependent"\`, \`"u_dep"\`, \`"u_indep"\`, \`"u_mix"\`,
  \`"hl"\`, \`"gudmm"\`, \`"dkss"\`, \`"mod_gower"\`, and
  \`"euclidean"\`. When \`preset\` is not \`"custom"\`, arguments such
  as \`method_cat\`, \`method_num\`, \`commensurable\`, and
  \`interaction\` are normally handled by the preset and user-supplied
  values are ignored. The exception is \`method_num\` for \`preset =
  "euclidean"\` when all predictors are numeric.

- interaction:

  Logical. If \`TRUE\`, adds an interaction-aware continuous-categorical
  component based on local predictive separability.

- prop_nn:

  Numeric. Proportion of nearest neighbours used when \`interaction =
  TRUE\`.

- score:

  Character string specifying the score used when \`interaction =
  TRUE\`. Available values include \`"ba"\` for balanced accuracy and
  \`"logloss"\`.

- decision:

  Character string specifying the decision rule used when \`score =
  "ba"\`. The default is \`"prior_corrected"\`.

- gower_average:

  Logical; only used when \`preset = "gower"\`. If \`TRUE\`, returns the
  standard Gower dissimilarity averaged over variables, matching the
  scale of \[cluster::daisy()\] with \`metric = "gower"\`. If \`FALSE\`,
  returns the sum of per-variable Gower contributions, equivalent to
  multiplying the averaged Gower dissimilarity by the number of active
  variables.

## Value

An object of class \`"MDist"\`. The object contains the computed
dissimilarity in its \`\$distance\` field, the selected \`preset\`, the
training data, and a list of parameters describing the fitted distance
specification. Square train-train dissimilarities are stored as
\`"dissimilarity"\`/\`"dist"\` objects; rectangular test-to-training
dissimilarities are stored as \`"dissimilarity"\`/\`"matrix"\` objects.

## Details

\`mdist()\` is the main distance-construction function in \`manydist\`.
It can return ordinary train-train dissimilarities or rectangular
test-to-training dissimilarities when \`new_data\` is supplied. The
resulting object stores both the dissimilarity matrix and metadata about
the distance specification that was used.

With \`preset = "custom"\`, users manually choose the numerical
preprocessing, categorical dissimilarity, commensurability, and optional
interaction term.

The \`"gower"\` preset follows the usual Gower construction based on
range scaling for continuous variables and matching dissimilarities for
categorical variables. The \`gower_average\` argument controls whether
the result is averaged over variables or returned as a sum of
variable-wise contributions.

\#' The \`"u_dep"\`, \`"unbiased_dependent"\`, \`"u_indep"\`, and
\`"u_mix"\` presets are convenience specifications for unbiased or
commensurable mixed-variable dissimilarities. The \`"euclidean"\` preset
computes Euclidean distance after standardizing numerical variables and
one-hot encoding and standardizing categorical variables. For
numerical-only inputs, \`"std"\` remains the default, but \`method_num\`
can be overridden. Setting \`method_num = "none"\` produces the same
numerical distances as \[stats::dist()\] applied directly to the
original variables. The \`"gudmm"\`, \`"dkss"\`, and \`"mod_gower"\`
presets provide additional mixed-type distance constructions. Some
presets currently support only train-train distances and will stop if
\`new_data\` is supplied.

Use \[all_dist_method_specs()\] to inspect the available distance
components and method specifications.

## See also

\[step_mdist()\], \[all_dist_method_specs()\]

## Examples

``` r
if (requireNamespace("palmerpenguins", quietly = TRUE)) {
  data("penguins", package = "palmerpenguins")

  penguins_small <- palmerpenguins::penguins |>
    dplyr::select(
      bill_length_mm, bill_depth_mm, flipper_length_mm,
      body_mass_g, species, island, sex
    ) |>
    tidyr::drop_na()

  # Gower distance on mixed-type data
  d_gower <- mdist(penguins_small, preset = "gower")
  d_gower

  # Custom mixed-type specification
  d_custom <- mdist(
    penguins_small,
    preset = "custom",
    method_cat = "matching",
    method_num = "std",
    commensurable = TRUE
  )

  d_custom

  # Train-to-new-data distances
  penguin_split <- rsample::initial_split(penguins_small, prop = 0.75)
  penguin_train <- rsample::training(penguin_split)
  penguin_test  <- rsample::testing(penguin_split)

  d_new <- mdist(
    penguin_train,
    new_data = penguin_test,
    preset = "gower"
  )

  d_new
}
#> MDist object
#>   preset : gower 
#>   number of training observations : 249 
#>   number of test observations     : 84 
#>   number of continuous variables   : 4 
#>   number of categorical variables   : 3 
#>   parameters:
#>     - commensurability adjustment: FALSE
```
