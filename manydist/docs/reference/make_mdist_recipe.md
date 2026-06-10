# Build a recipe that computes mdist-based distance features

This helper creates a tidymodels recipe using step_mdist(). It supports
both mdist presets and custom param sets.

## Usage

``` r
make_mdist_recipe(df, mdist_type, mdist_preset, param_set, outcome)
```

## Arguments

- df:

  A data frame.

- mdist_type:

  "preset" or "custom".

- mdist_preset:

  Name of the preset (if mdist_type == "preset").

- param_set:

  A list of custom mdist arguments (if mdist_type == "custom").

- outcome:

  Name of the outcome variable.

## Value

A \`recipes::recipe()\` object.
