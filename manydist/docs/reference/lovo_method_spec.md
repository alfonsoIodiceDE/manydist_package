# Create a LOVO method specification for \`compare_lovo_mdist()\`

Helper to build one method specification for \[compare_lovo_mdist()\].
This allows tidy-style specification of \`response\`, e.g. \`response =
Name\`, while storing the method definition as a regular list.

## Usage

``` r
lovo_method_spec(response = NULL, ...)
```

## Arguments

- response:

  Optional response column, supplied either unquoted (e.g. \`Name\`) or
  quoted (e.g. \`"Name"\`).

- ...:

  Additional arguments passed on to \[lovo_mdist()\] through
  \[compare_lovo_mdist()\].

## Value

A named list of arguments suitable for one element of the \`methods\`
argument in \[compare_lovo_mdist()\].

## Examples

``` r
if (FALSE) { # \dontrun{
methods <- list(
  tvd_sup = lovo_method_spec(
    response = species,
    preset = "custom",
    method_cat = "tvd",
    method_num = "std",
    commensurable = TRUE,
    response_used = TRUE
  ),
  tvd_unsup = lovo_method_spec(
    response = species,
    preset = "custom",
    method_cat = "tvd",
    method_num = "std",
    commensurable = TRUE,
    response_used = FALSE
  )
)
} # }
```
