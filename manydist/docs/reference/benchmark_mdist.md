# Benchmark \`mdist()\` over multiple method specifications

Applies \[mdist()\] repeatedly over a tibble of distance-method
specifications, typically generated with \[all_dist_method_specs()\].

## Usage

``` r
benchmark_mdist(x, response = NULL, specs = all_dist_method_specs())
```

## Arguments

- x:

  A data frame or tibble of predictors, optionally including the
  response column.

- response:

  Optional response column inside \`x\`, supplied either unquoted or as
  a character string.

- specs:

  A tibble of method specifications. By default, this is generated with
  \[all_dist_method_specs()\]. It must contain the columns
  \`spec_type\`, \`preset\`, \`method_cat\`, \`method_num\`, and
  \`commensurable\`.

## Value

A tibble containing the supplied specifications together with:

- result:

  The corresponding output of \[mdist()\], or an error object if the
  specification failed.

- ok:

  Logical indicator; \`TRUE\` if the run completed successfully,
  \`FALSE\` otherwise.

- error:

  Error message for failed runs, \`NA\` otherwise.

## Details

Each row of \`specs\` is interpreted as one valid \`mdist()\`
configuration. Preset-based and custom component-based specifications
are both supported. Failed specifications are caught and returned in the
output rather than stopping the full benchmark.

Preset specifications use the \`preset\` column and ignore
\`method_cat\`, \`method_num\`, and \`commensurable\`. Component
specifications are evaluated as \`preset = "custom"\` and use
\`method_cat\`, \`method_num\`, and \`commensurable\`.

This function is intended for benchmarking, validation, and sensitivity
analyses across multiple distance specifications.

## Examples

``` r
if (requireNamespace("palmerpenguins", quietly = TRUE)) {
  data("penguins", package = "palmerpenguins")

  penguins_small <- palmerpenguins::penguins |>
    dplyr::select(
      species, bill_length_mm, bill_depth_mm, flipper_length_mm,
      body_mass_g, island, sex
    ) |>
    tidyr::drop_na()

  specs <- all_dist_method_specs(mode = "presets_only")

  res <- benchmark_mdist(
    penguins_small,
    response = species,
    specs = specs
  )

  res |>
    dplyr::select(spec_type, preset, ok, error)
}
#> Warning: For method(s) 'HLeucl', category dissimilarities do not depend on conditional profiles; `response` was therefore ignored.
#> Warning: For method(s) 'matching', category dissimilarities do not depend on conditional profiles; `response` was therefore ignored.
#> # A tibble: 10 × 4
#>    spec_type preset    ok    error
#>    <chr>     <chr>     <lgl> <chr>
#>  1 preset    euclidean TRUE  NA   
#>  2 preset    gower     TRUE  NA   
#>  3 preset    hl        TRUE  NA   
#>  4 preset    u_dep     TRUE  NA   
#>  5 preset    u_indep   TRUE  NA   
#>  6 preset    u_mix     TRUE  NA   
#>  7 preset    dkss      TRUE  NA   
#>  8 preset    gudmm     TRUE  NA   
#>  9 preset    mod_gower TRUE  NA   
#> 10 preset    custom    TRUE  NA   
```
