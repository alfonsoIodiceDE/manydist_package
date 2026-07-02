# Distance specifications

## 1 Overview

`manydist` provides tools for constructing dissimilarities/distances for
numerical, categorical, and mixed-type data.

The main function is
[`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md):
the user can either specify the distance via a `preset` or customize it.
The custom specification requires:

- `method_cat` controls the categorical dissimilarity;
- `method_num` controls the preprocessing of numerical variables;
- `commensurable` controls whether variable-wise contributions should
  be, on average, the same;

## 2 Setup

``` r

library(manydist)
library(dplyr)
library(tidyr)
library(recipes)
```

We use the `palmerpenguins` data for illustration.

``` r

if (requireNamespace("palmerpenguins", quietly = TRUE)) {
  data("penguins", package = "palmerpenguins")

  penguins_small <- palmerpenguins::penguins |>
    dplyr::select(-year)|>
    tidyr::drop_na()
}
```

The data contain numerical variables, such as bill length and body mass,
and categorical variables, such as island and sex.

``` r

skim_cols <- c(
  "skim_type",
  "skim_variable",
  "complete_rate",
  "numeric.mean",
  "numeric.sd",
  "factor.n_unique"
)

skimr::skim(penguins_small) |>
  dplyr::as_tibble() |>
  dplyr::select(dplyr::any_of(skim_cols)) |>
  kableExtra::kbl(format = "html", digits = 2) |>
  kableExtra::kable_styling(full_width = FALSE)
```

| skim_type | skim_variable | complete_rate | numeric.mean | numeric.sd | factor.n_unique |
|:---|:---|---:|---:|---:|---:|
| factor | species | 1 | NA | NA | 3 |
| factor | island | 1 | NA | NA | 3 |
| factor | sex | 1 | NA | NA | 2 |
| numeric | bill_length_mm | 1 | 43.99 | 5.47 | NA |
| numeric | bill_depth_mm | 1 | 17.16 | 1.97 | NA |
| numeric | flipper_length_mm | 1 | 200.97 | 14.02 | NA |
| numeric | body_mass_g | 1 | 4207.06 | 805.22 | NA |

## 3 Computing a mixed-type dissimilarity

The simplest way to compute a mixed-type dissimilarity is to use a
`preset`.

For example, the `"gower"` preset computes a Gower-type dissimilarity
for mixed numerical and categorical data.

``` r

d_gow <- penguins_small |> dplyr::select(-species) |> mdist(preset = "gower")

d_gow
```

    MDist object
      preset : gower
      number of observations : 333
      number of continuous variables   : 4
      number of categorical variables   : 2
      parameters:
        - commensurability adjustment: FALSE

The result is an `MDist` object. Internally, `MDist` is implemented as
an R6 object, which stores the computed dissimilarity together with
metadata and methods for extracting or reusing the result. The
dissimilarity can be extracted as a standard `dist` object using
`to_dist()`. For display, we convert it to a matrix and show only the
first few rows and columns.

``` r

d_gow_dist <- d_gow$to_dist()

d_gow_mat <- as.matrix(d_gow_dist)

d_gow_mat[1:5, 1:5] |>
  round(digits = 2) 
```

         1    2    3    4    5
    1 0.00 0.21 0.25 0.24 0.07
    2 0.21 0.00 0.07 0.09 0.25
    3 0.25 0.07 0.00 0.06 0.26
    4 0.24 0.09 0.06 0.00 0.23
    5 0.07 0.25 0.26 0.23 0.00

The `"gower"` preset gives the same dissimilarities as
[`cluster::daisy()`](https://rdrr.io/pkg/cluster/man/daisy.html) with
`metric = "gower"`.

``` r

d_daisy <- cluster::daisy(
  penguins_small |> dplyr::select(-species),
  metric = "gower"
)

all.equal(
  as.matrix(d_gow$to_dist()),
  as.matrix(d_daisy),
  tolerance = 1e-10
)
```

    [1] TRUE

The `preset` argument provides shortcuts for common distance
specifications. Some presets are classical mixed-type distances, while
others are designed to construct unbiased or association-aware
dissimilarities.

| preset | description |
|:---|:---|
| gower | Gower-type mixed-data dissimilarity. |
| euclidean | Euclidean distance after one-hot encoding categorical variables. |
| u_indep | Unbiased distance with independent variable-wise contributions. |
| u_dep | Unbiased association-aware distance. |
| u_mix | Mixed unbiased association-aware distance. |
| hl | Heterogeneous Euclidean-type distance. |
| gudmm | Generalized distance for mixed data. |
| dkss | Distance for mixed data based on distributional comparisons. |
| mod_gower | Modified Gower-type dissimilarity. |

For the complete list of supported specifications and their argument
combinations, see:

``` r

?mdist
```

## 4 Custom specification

Presets are convenient, but users can also define custom specifications.

In a custom specification, `method_cat` controls the categorical
dissimilarity and `method_num` controls the preprocessing applied to
numerical variables.

``` r

d_custom <- mdist(
  penguins_small |> dplyr::select(-species),
  preset = "custom",
  method_cat = "matching",
  method_num = "std",
  commensurable = TRUE
)

d_custom
```

    MDist object
      preset : custom
      number of observations : 333
      number of continuous variables   : 4
      number of categorical variables   : 2
      parameters:
        - categorical method: matching
        - numerical preprocessing: std
        - commensurability adjustment: TRUE
        - number of principal components: 6
        - inertia threshold: NULL

Here, categorical variables are compared using matching dissimilarities,
numerical variables are standardized, and variable-wise contributions
are made commensurable. Note that pairwise distances are always computed
via Manhattan distance, with the exection of `euclidean` and `hl`
presets.

## 5 Response-aware specifications

Some distance specifications can use a response variable when
constructing the dissimilarity. In these cases, the response is supplied
through the `response` argument of
[`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md).

The response is not treated as an ordinary predictor. Instead, it is
removed from the predictor set and used only by the components that
support response-aware construction.

Response-aware construction is available in two ways:

1.  through response-aware presets;
2.  through response-aware categorical methods in custom specifications.

The response-aware presets can be inspected with
`all_dist_method_specs(mode = "response_aware_only")`.

``` r

response_aware_presets <- all_dist_method_specs(mode = "response_aware_only") |>
  dplyr::filter(spec_type == "preset") |>
  dplyr::distinct(preset) |>
  dplyr::arrange(preset)

response_aware_presets |>
  kableExtra::kbl(
    format = "html",
    col.names = "Response-aware preset"
  ) |>
  kableExtra::kable_styling(full_width = FALSE)
```

| Response-aware preset |
|:----------------------|
| u_dep                 |
| u_mix                 |

For example, the `"u_dep"` preset can use the response when constructing
the dissimilarity.

``` r

d_response <- mdist(
  penguins_small,
  response = species,
  preset = "u_dep"
)

d_response
```

    MDist object
      preset : u_dep
      number of observations : 333
      number of continuous variables   : 4
      number of categorical variables   : 2
      parameters:
        - categorical method: tvd
        - numerical preprocessing: pc_scores
        - commensurability adjustment: TRUE

The response can also be supplied as a character string.

``` r

d_response_chr <- mdist(
  penguins_small,
  response = "species",
  preset = "u_dep"
)

d_response_chr
```

    MDist object
      preset : u_dep
      number of observations : 333
      number of continuous variables   : 4
      number of categorical variables   : 2
      parameters:
        - categorical method: tvd
        - numerical preprocessing: pc_scores
        - commensurability adjustment: TRUE

In both calls, `species` is used as the response variable and is not
included among the predictors.

Response-aware construction can also be used in custom specifications.
In this case, `preset = "custom"` is combined with a response-aware
categorical method.

The following table lists the response-aware categorical methods
implemented directly in `manydist`.

``` r

dist_methods_tbl() |>
  dplyr::filter(
    argument == "method_cat",
    response_aware,
    engine == "manydist") |> 
  dplyr::select(method) |> 
  kableExtra::kbl(full_width = FALSE)
```

| method    |
|:----------|
| gifi_chi2 |
| le_and_ho |
| tvd       |

In addition to the methods implemented directly in `manydist`,
categorical dissimilarities available through `philentropy` can also be
used in response-aware custom specifications.
