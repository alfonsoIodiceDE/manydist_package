# Create a grid of \`mdist()\` method specifications

Construct a tibble of valid method specifications that can be passed to
\[benchmark_mdist()\] or iterated over manually to benchmark \[mdist()\]
across multiple configurations.

## Usage

``` r
all_dist_method_specs(
  mode = c("full", "presets_only", "response_aware_only"),
  method_cat = NULL,
  method_num = NULL,
  preset = NULL,
  commensurable = NULL
)
```

## Arguments

- mode:

  Character string controlling the initial pool of specifications.
  Supported values are \`"full"\`, \`"presets_only"\`, and
  \`"response_aware_only"\`.

- method_cat:

  Optional character vector restricting the categorical methods used in
  component-based specifications.

- method_num:

  Optional character vector restricting the numerical preprocessing
  methods used in component-based specifications.

- preset:

  Optional character vector restricting preset-based specifications.

- commensurable:

  Optional logical vector restricting the commensurable values used in
  component-based specifications.

## Value

A tibble where each row represents one valid \`mdist()\` specification.
The tibble contains \`spec_type\`, \`preset\`, \`method_cat\`,
\`method_num\`, and \`commensurable\`.

## Details

The resulting tibble combines preset-based specifications and custom
component-based specifications built from the currently available
categorical methods and numerical preprocessing options listed in
\[dist_methods_tbl()\].

\`mode\` defines the initial candidate pool. Any explicit argument
filters supplied by the user are applied afterwards and therefore
restrict the selected pool further.

With \`mode = "response_aware_only"\`, preset specifications and
categorical methods are restricted to response-aware methods.

## Examples

``` r
all_dist_method_specs()
#> # A tibble: 482 × 5
#>    spec_type preset    method_cat method_num commensurable
#>    <chr>     <chr>     <chr>      <chr>      <lgl>        
#>  1 preset    euclidean NA         NA         NA           
#>  2 preset    gower     NA         NA         NA           
#>  3 preset    hl        NA         NA         NA           
#>  4 preset    u_dep     NA         NA         NA           
#>  5 preset    u_indep   NA         NA         NA           
#>  6 preset    u_mix     NA         NA         NA           
#>  7 preset    dkss      NA         NA         NA           
#>  8 preset    gudmm     NA         NA         NA           
#>  9 preset    mod_gower NA         NA         NA           
#> 10 preset    custom    NA         NA         NA           
#> # ℹ 472 more rows
all_dist_method_specs(mode = "presets_only")
#> # A tibble: 10 × 5
#>    spec_type preset    method_cat method_num commensurable
#>    <chr>     <chr>     <chr>      <chr>      <lgl>        
#>  1 preset    euclidean NA         NA         NA           
#>  2 preset    gower     NA         NA         NA           
#>  3 preset    hl        NA         NA         NA           
#>  4 preset    u_dep     NA         NA         NA           
#>  5 preset    u_indep   NA         NA         NA           
#>  6 preset    u_mix     NA         NA         NA           
#>  7 preset    dkss      NA         NA         NA           
#>  8 preset    gudmm     NA         NA         NA           
#>  9 preset    mod_gower NA         NA         NA           
#> 10 preset    custom    NA         NA         NA           
all_dist_method_specs(mode = "response_aware_only")
#> # A tibble: 354 × 5
#>    spec_type preset method_cat    method_num commensurable
#>    <chr>     <chr>  <chr>         <chr>      <lgl>        
#>  1 preset    u_dep  NA            NA         NA           
#>  2 preset    u_mix  NA            NA         NA           
#>  3 component custom additive_symm none       FALSE        
#>  4 component custom additive_symm pc_scores  FALSE        
#>  5 component custom additive_symm pc_scores  TRUE         
#>  6 component custom additive_symm range      FALSE        
#>  7 component custom additive_symm robust     FALSE        
#>  8 component custom additive_symm robust     TRUE         
#>  9 component custom additive_symm std        FALSE        
#> 10 component custom additive_symm std        TRUE         
#> # ℹ 344 more rows
all_dist_method_specs(mode = "full", method_cat = c("tvd", "le_and_ho"))
#> # A tibble: 26 × 5
#>    spec_type preset    method_cat method_num commensurable
#>    <chr>     <chr>     <chr>      <chr>      <lgl>        
#>  1 preset    euclidean NA         NA         NA           
#>  2 preset    gower     NA         NA         NA           
#>  3 preset    hl        NA         NA         NA           
#>  4 preset    u_dep     NA         NA         NA           
#>  5 preset    u_indep   NA         NA         NA           
#>  6 preset    u_mix     NA         NA         NA           
#>  7 preset    dkss      NA         NA         NA           
#>  8 preset    gudmm     NA         NA         NA           
#>  9 preset    mod_gower NA         NA         NA           
#> 10 preset    custom    NA         NA         NA           
#> # ℹ 16 more rows
```
