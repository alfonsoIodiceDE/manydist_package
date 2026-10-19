# Create a grid of \`mdist()\` method specifications

Construct a tibble of valid method specifications that can be passed to
\[benchmark_mdist()\] or iterated over manually to benchmark \[mdist()\]
across multiple configurations.

## Usage

``` r
all_dist_method_specs(
  method_cat = NULL,
  method_num = NULL,
  preset = NULL,
  commensurable = NULL
)
```

## Arguments

- method_cat:

  Optional character vector restricting the categorical methods used in
  custom specifications.

- method_num:

  Optional character vector restricting the numerical preprocessing
  methods used in custom specifications.

- preset:

  \`NULL\` (the default) returns custom combinations only. A character
  vector selects named presets; \`"all"\` selects all named presets.
  Include \`"custom"\` explicitly to also return custom combinations,
  for example \`c("all", "custom")\`. Component filters apply only to
  custom rows.

- commensurable:

  Optional logical vector restricting the commensurable values used in
  custom specifications.

## Value

A tibble where each row represents one valid \`mdist()\` specification.
The tibble contains \`preset\`, \`method_cat\`, \`method_num\`, and
\`commensurable\`.

## Details

Select named presets or construct custom specifications from the
categorical methods and numerical preprocessing options listed in
\[dist_methods_tbl()\].

Rows with \`preset = "custom"\` use the explicit categorical, numerical,
and commensurability options. Named presets determine those settings
internally.

With \`commensurable = TRUE\`, custom grids retain only \`"std"\` and
\`"pc_scores"\` numerical preprocessing. Variable-wise affine scaling
cancels under commensurability, so \`"robust"\` is omitted as redundant
with \`"std"\`. The existing exclusions of \`"none"\` and \`"range"\`
are retained. Range scaling clips new observations outside the training
range, so its equivalence on training data does not extend to all new
observations. These restrictions affect grid construction, not the
options accepted by \[mdist()\].

Named preset rows retain \`NA\` in their component columns: their
settings are determined by the preset. Supplying component filters
without requesting custom rows gives an error rather than silently
ignoring those filters.

## Examples

``` r
all_dist_method_specs()
#> # A tibble: 413 × 4
#>    preset method_cat method_num commensurable
#>    <chr>  <chr>      <chr>      <lgl>        
#>  1 custom HL         none       FALSE        
#>  2 custom HL         pc_scores  FALSE        
#>  3 custom HL         pc_scores  TRUE         
#>  4 custom HL         range      FALSE        
#>  5 custom HL         robust     FALSE        
#>  6 custom HL         std        FALSE        
#>  7 custom HL         std        TRUE         
#>  8 custom HLeucl     none       FALSE        
#>  9 custom HLeucl     pc_scores  FALSE        
#> 10 custom HLeucl     pc_scores  TRUE         
#> # ℹ 403 more rows
all_dist_method_specs(preset = c("gower", "u_indep", "u_dep"))
#> # A tibble: 3 × 4
#>   preset  method_cat method_num commensurable
#>   <chr>   <chr>      <chr>      <lgl>        
#> 1 gower   NA         NA         NA           
#> 2 u_dep   NA         NA         NA           
#> 3 u_indep NA         NA         NA           
all_dist_method_specs(preset = "all")
#> # A tibble: 10 × 4
#>    preset    method_cat method_num commensurable
#>    <chr>     <chr>      <chr>      <lgl>        
#>  1 euclidean NA         NA         NA           
#>  2 gower     NA         NA         NA           
#>  3 hl        NA         NA         NA           
#>  4 u_dep     NA         NA         NA           
#>  5 u_dep_bw  NA         NA         NA           
#>  6 u_indep   NA         NA         NA           
#>  7 u_mix     NA         NA         NA           
#>  8 dkss      NA         NA         NA           
#>  9 gudmm     NA         NA         NA           
#> 10 mod_gower NA         NA         NA           
all_dist_method_specs(preset = c("all", "custom"))
#> # A tibble: 423 × 4
#>    preset    method_cat method_num commensurable
#>    <chr>     <chr>      <chr>      <lgl>        
#>  1 euclidean NA         NA         NA           
#>  2 gower     NA         NA         NA           
#>  3 hl        NA         NA         NA           
#>  4 u_dep     NA         NA         NA           
#>  5 u_dep_bw  NA         NA         NA           
#>  6 u_indep   NA         NA         NA           
#>  7 u_mix     NA         NA         NA           
#>  8 dkss      NA         NA         NA           
#>  9 gudmm     NA         NA         NA           
#> 10 mod_gower NA         NA         NA           
#> # ℹ 413 more rows
all_dist_method_specs(method_cat = c("tvd", "le_and_ho"))
#> # A tibble: 14 × 4
#>    preset method_cat method_num commensurable
#>    <chr>  <chr>      <chr>      <lgl>        
#>  1 custom le_and_ho  none       FALSE        
#>  2 custom le_and_ho  pc_scores  FALSE        
#>  3 custom le_and_ho  pc_scores  TRUE         
#>  4 custom le_and_ho  range      FALSE        
#>  5 custom le_and_ho  robust     FALSE        
#>  6 custom le_and_ho  std        FALSE        
#>  7 custom le_and_ho  std        TRUE         
#>  8 custom tvd        none       FALSE        
#>  9 custom tvd        pc_scores  FALSE        
#> 10 custom tvd        pc_scores  TRUE         
#> 11 custom tvd        range      FALSE        
#> 12 custom tvd        robust     FALSE        
#> 13 custom tvd        std        FALSE        
#> 14 custom tvd        std        TRUE         
```
