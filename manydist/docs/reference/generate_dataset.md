# generate toy mixed datasets for the supplementary material: figures 1 to 3

generate toy mixed datasets for the supplementary material: figures 1 to
3

## Usage

``` r
generate_dataset(
  n,
  porig,
  pn,
  pnnoise,
  pcnoise,
  sigma,
  qoptions,
  seed = NULL,
  mode = "per_variable"
)
```

## Arguments

- n:

  number of observations

- porig:

  number of original informative continuous variables

- pn:

  number of continuous variables (total, before adding noise variables)

- pnnoise:

  number of extra numeric noise variables

- pcnoise:

  number of extra categorical noise variables

- sigma:

  sd for noise added to informative variables

- qoptions:

  number of bins for categorization (vector if per_variable, scalar if
  shared)

- seed:

  optional seed

- mode:

  either "per_variable" or "shared"
