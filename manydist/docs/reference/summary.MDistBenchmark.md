# Summarize a distance benchmark

Displays the pairwise diagnostics already computed by
[`benchmark_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/benchmark_mdist.md),
without refitting distances or recomputing diagnostics.

## Usage

``` r
# S3 method for class 'MDistBenchmark'
summary(object, metric = NULL, digits = 3L, ..., n = Inf, width = NULL)
```

## Arguments

- object:

  An object returned by
  [`benchmark_mdist()`](https://alfonsoiodicede.github.io/manydist_package/reference/benchmark_mdist.md).

- metric:

  Character vector selecting diagnostics to display. `NULL` (the
  default) displays all available diagnostics. Choices include `"mad"`,
  `"relative_distance"`, `"mds_congruence"`, `"alienation"`, and
  available ARI columns such as `"ari_pam"`. `"ari"` selects all
  available ARI columns.

- digits:

  Number of decimal places in the displayed tables.

- ...:

  Additional arguments passed to the matrix print method.

- n:

  Maximum number of methods to display in each table. The default,
  `Inf`, displays every successful method. The returned tibble is not
  limited.

- width:

  Width of the printed output; `NULL` uses the default width.

## Value

Invisibly returns the complete tibble with one row per unique pair of
successful specifications, regardless of `metric`, `digits`, or `n`.
Columns include the method labels, MAD, symmetric relative distance, MDS
congruence, alienation, and any requested clustering ARIs. If fewer than
two specifications succeed, the tibble has zero rows.

## Examples

``` r
if (FALSE) { # \dontrun{
pairwise_results <- summary(benchmark_result)
summary(benchmark_result, metric = "relative_distance")
summary(benchmark_result, metric = c("alienation", "ari"))
} # }
```
