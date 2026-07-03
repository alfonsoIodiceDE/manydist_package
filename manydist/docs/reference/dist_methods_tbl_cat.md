# List available categorical dissimilarities

Convenience wrapper around \[dist_methods_tbl()\] returning only methods
that can be used as categorical dissimilarities through \`method_cat\`.

## Usage

``` r
dist_methods_tbl_cat()
```

## Value

A tibble containing the rows of \[dist_methods_tbl()\] with \`argument
== "method_cat"\`.

## Examples

``` r
dist_methods_tbl_cat()
#> # A tibble: 64 × 6
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
#> # ℹ 54 more rows
```
