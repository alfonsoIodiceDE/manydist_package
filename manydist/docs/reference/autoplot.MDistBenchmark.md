# Plot pairwise distance benchmark comparisons

Draws an annotated triangular heatmap of the pairwise diagnostics from
\[benchmark_mdist()\].

## Usage

``` r
# S3 method for class 'MDistBenchmark'
autoplot(
  object,
  metric = c("relative_distance", "mad", "alienation", "mds_congruence", "ari"),
  cluster_method = NULL,
  digits = 2,
  ...
)
```

## Arguments

- object:

  An object returned by \[benchmark_mdist()\].

- metric:

  Character string selecting \`"mad"\`, \`"relative_distance"\`,
  \`"mds_congruence"\`, \`"alienation"\`, or \`"ari"\`. A specific ARI
  column, such as \`"ari_pam"\`, can also be supplied.

- cluster_method:

  Optional clustering method used when \`metric = "ari"\`. If \`NULL\`,
  all available ARI metrics are shown, with one facet per clustering
  method.

- digits:

  Number of decimal places used for cell labels.

- ...:

  Currently unused.

## Value

A \`ggplot\` object.
