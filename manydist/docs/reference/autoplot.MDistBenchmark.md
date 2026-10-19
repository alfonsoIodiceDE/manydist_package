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

## Details

The default fill gradient runs from coral (\`#E76F51\`) for lower values
to blue (\`#008CFF\`) for higher values, with regular-weight white cell
labels. To use another palette, add a fill scale to the returned plot.
Self-comparisons on the main diagonal are left blank. ARI plots use
fixed colour limits of 0 and 1, so the same ARI has the same colour
across plots. Negative ARIs retain their numerical labels but use the
low-end colour.
