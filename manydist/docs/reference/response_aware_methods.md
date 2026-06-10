# List response-aware methods

Returns the methods in \[dist_methods_tbl()\] that can use a response
variable. The output can optionally be restricted by data type or by
\`mdist()\` argument.

## Usage

``` r
response_aware_methods(data_type = NULL, argument = NULL)
```

## Arguments

- data_type:

  Optional character vector used to restrict the returned methods by
  data type, for example \`"categorical"\` or \`"mixed"\`.

- argument:

  Optional character vector used to restrict the returned methods by
  argument, for example \`"method_cat"\` or \`"preset"\`.

## Value

A character vector of response-aware method names.

## Examples

``` r
response_aware_methods()
#>  [1] "gifi_chi2"         "le_and_ho"         "tvd"              
#>  [4] "additive_symm"     "avg"               "bhattacharyya"    
#>  [7] "canberra"          "chebyshev"         "clark"            
#> [10] "cosine"            "czekanowski"       "dice"             
#> [13] "divergence"        "euclidean"         "fidelity"         
#> [16] "gower"             "harmonic_mean"     "hassebrook"       
#> [19] "hellinger"         "inner_product"     "intersection"     
#> [22] "jaccard"           "jeffreys"          "jensen-shannon"   
#> [25] "jensen_difference" "k_divergence"      "kulczynski_d"     
#> [28] "kulczynski_s"      "kullback-leibler"  "kumar-johnson"    
#> [31] "lorentzian"        "manhattan"         "matusita"         
#> [34] "minkowski"         "motyka"            "neyman"           
#> [37] "non-intersection"  "pearson"           "prob_symm"        
#> [40] "ruzicka"           "soergel"           "sorensen"         
#> [43] "squared_chi"       "squared_chord"     "squared_euclidean"
#> [46] "taneja"            "tanimoto"          "topsoe"           
#> [49] "wavehedges"        "u_dep"             "u_mix"            
response_aware_methods(argument = "method_cat")
#>  [1] "gifi_chi2"         "le_and_ho"         "tvd"              
#>  [4] "additive_symm"     "avg"               "bhattacharyya"    
#>  [7] "canberra"          "chebyshev"         "clark"            
#> [10] "cosine"            "czekanowski"       "dice"             
#> [13] "divergence"        "euclidean"         "fidelity"         
#> [16] "gower"             "harmonic_mean"     "hassebrook"       
#> [19] "hellinger"         "inner_product"     "intersection"     
#> [22] "jaccard"           "jeffreys"          "jensen-shannon"   
#> [25] "jensen_difference" "k_divergence"      "kulczynski_d"     
#> [28] "kulczynski_s"      "kullback-leibler"  "kumar-johnson"    
#> [31] "lorentzian"        "manhattan"         "matusita"         
#> [34] "minkowski"         "motyka"            "neyman"           
#> [37] "non-intersection"  "pearson"           "prob_symm"        
#> [40] "ruzicka"           "soergel"           "sorensen"         
#> [43] "squared_chi"       "squared_chord"     "squared_euclidean"
#> [46] "taneja"            "tanimoto"          "topsoe"           
#> [49] "wavehedges"        "u_dep"             "u_mix"            
response_aware_methods(argument = "preset")
#>  [1] "gifi_chi2"         "le_and_ho"         "tvd"              
#>  [4] "additive_symm"     "avg"               "bhattacharyya"    
#>  [7] "canberra"          "chebyshev"         "clark"            
#> [10] "cosine"            "czekanowski"       "dice"             
#> [13] "divergence"        "euclidean"         "fidelity"         
#> [16] "gower"             "harmonic_mean"     "hassebrook"       
#> [19] "hellinger"         "inner_product"     "intersection"     
#> [22] "jaccard"           "jeffreys"          "jensen-shannon"   
#> [25] "jensen_difference" "k_divergence"      "kulczynski_d"     
#> [28] "kulczynski_s"      "kullback-leibler"  "kumar-johnson"    
#> [31] "lorentzian"        "manhattan"         "matusita"         
#> [34] "minkowski"         "motyka"            "neyman"           
#> [37] "non-intersection"  "pearson"           "prob_symm"        
#> [40] "ruzicka"           "soergel"           "sorensen"         
#> [43] "squared_chi"       "squared_chord"     "squared_euclidean"
#> [46] "taneja"            "tanimoto"          "topsoe"           
#> [49] "wavehedges"        "u_dep"             "u_mix"            
```
