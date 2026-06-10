# Compare LOVO diagnostics across multiple distance specifications

This function compares \*\*leave-one-variable-out (LOVO)\*\* diagnostics
across multiple distance definitions supported by `manydist`.

## Usage

``` r
compare_lovo_mdist(
  x,
  methods,
  dims = 2,
  keep_dist = FALSE,
  .progress = FALSE,
  ...
)
```

## Arguments

- x:

  A data frame or tibble containing the predictors.

- methods:

  A \*\*named list\*\* describing the distance specifications to
  compare. Each element must be a list of arguments passed to
  [`lovo_mdist`](https://alfonsoiodicede.github.io/manydist_package/reference/lovo_mdist.md).

  For example:


      methods = list(
        gower = list(preset = "gower"),
        u_dep = list(preset = "unbiased_dependent"),
        custom = list(
          preset = "custom",
          method_cat = "matching",
          method_num = "std",
          commensurable = TRUE
        )
      )

- dims:

  Number of dimensions used for the MDS configuration when computing
  congruence-based diagnostics.

- keep_dist:

  Logical; if `TRUE`, distance matrices from the LOVO computations are
  retained. This increases memory usage.

- .progress:

  Logical; if `TRUE`, progress messages are printed while computing LOVO
  diagnostics for each method.

- ...:

  Additional arguments passed to
  [`lovo_mdist`](https://alfonsoiodicede.github.io/manydist_package/reference/lovo_mdist.md)
  unless overridden in the method-specific argument list.

  These may include optional clustering diagnostics, for example
  `cluster_k`, `cluster_methods`, and `hclust_method`.

## Value

An object of class `MDistLOVOCompare` containing:

- results:

  A tibble with one row per method-variable combination.

- methods:

  The list of distance specifications used.

- dims:

  Number of MDS dimensions used.

- n_obs:

  Number of observations in the dataset.

## Details

For each distance specification, the function:

1.  Computes the full mixed-type distance using
    [`mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md).

2.  Recomputes the distance repeatedly leaving out one variable at a
    time.

3.  Measures the impact of each variable using metrics such as mean
    absolute deviation (MAD), congruence-based diagnostics, and, when
    requested, clustering-based agreement measures.

The results are combined across methods and returned as an
`MDistLOVOCompare` object, which supports
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html).

## See also

[`lovo_mdist`](https://alfonsoiodicede.github.io/manydist_package/reference/lovo_mdist.md),
[`mdist`](https://alfonsoiodicede.github.io/manydist_package/reference/mdist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(manydist)
library(palmerpenguins)

data <- penguins |>
  dplyr::select(-species) |>
  tidyr::drop_na()

cmp <- compare_lovo_mdist(
  x = data,
  methods = list(
    gower = list(preset = "gower"),
    u_dep = list(preset = "unbiased_dependent")
  ),
  cluster_k = 3
)

summary(cmp)
autoplot(cmp, metric = "mad_importance")
autoplot(cmp, metric = "pam_importance")
} # }
```
