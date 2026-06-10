# List available \`manydist\` methods

Returns a tibble describing the categorical methods, numerical
preprocessing options, and preset specifications currently available in
\`manydist\`.

## Usage

``` r
dist_methods_tbl()
```

## Value

A tibble with one row per available method. The columns are:

- method:

  Name of the method or preset.

- argument:

  The \`mdist()\` argument to which the method belongs:
  \`"method_cat"\`, \`"method_num"\`, or \`"preset"\`.

- data_type:

  Type of variables targeted by the method.

- distance_basis:

  Broad methodological family.

- response_aware:

  Logical; whether the method can use a response variable.

- engine:

  Implementation source.

## Details

The table is used by helpers such as \[all_dist_method_specs()\] and
\[benchmark_mdist()\] to construct valid \`mdist()\` specifications.

## Examples

``` r
dist_methods_tbl()
#> # A tibble: 79 × 6
#>    method        argument   data_type   distance_basis response_aware engine    
#>    <chr>         <chr>      <chr>       <chr>          <lgl>          <chr>     
#>  1 gifi_chi2     method_cat categorical association    TRUE           manydist  
#>  2 le_and_ho     method_cat categorical association    TRUE           manydist  
#>  3 tvd           method_cat categorical association    TRUE           manydist  
#>  4 additive_symm method_cat categorical association    TRUE           philentro…
#>  5 avg           method_cat categorical association    TRUE           philentro…
#>  6 bhattacharyya method_cat categorical association    TRUE           philentro…
#>  7 canberra      method_cat categorical association    TRUE           philentro…
#>  8 chebyshev     method_cat categorical association    TRUE           philentro…
#>  9 clark         method_cat categorical association    TRUE           philentro…
#> 10 cosine        method_cat categorical association    TRUE           philentro…
#> # ℹ 69 more rows

dist_methods_tbl() |>
  dplyr::count(argument)
#> # A tibble: 3 × 2
#>   argument       n
#>   <chr>      <int>
#> 1 method_cat    64
#> 2 method_num     5
#> 3 preset        10
```
