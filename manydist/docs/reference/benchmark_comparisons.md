# Extract pairwise distance benchmark comparisons

Returns the pairwise diagnostics computed by \[benchmark_mdist()\].

## Usage

``` r
benchmark_comparisons(x)
```

## Arguments

- x:

  An object returned by \[benchmark_mdist()\].

## Value

A tibble with one row per unique pair of successful distance
specifications. It contains pairwise MAD, symmetric relative distance,
MDS congruence, alienation, and—when requested—one adjusted Rand index
column per clustering method.

## Examples

``` r
if (FALSE) { # \dontrun{
benchmark_comparisons(benchmark_result)
} # }
```
