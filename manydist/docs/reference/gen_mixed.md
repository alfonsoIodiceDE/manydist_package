# Generate mixed-type clustering data with signal and noise

Generates synthetic mixed datasets with controllable numerical and
categorical signal/noise structure and balanced cluster sizes.

## Usage

``` r
gen_mixed(
  k_true,
  clustSizeEq = 50,
  numsignal = 2,
  numnoise = 2,
  catsignal = 2,
  catnoise = 2,
  q = 5,
  q_err = 9,
  numsep = 0.1,
  catsep = 0.5,
  seed = NULL,
  error_type = c("normal", "chisq"),
  error_df = 2,
  error_scale = 1
)
```

## Arguments

- k_true:

  Number of clusters

- clustSizeEq:

  Observations per cluster

- numsignal:

  Number of numerical signal variables

- numnoise:

  Number of numerical noise variables

- catsignal:

  Number of categorical signal variables

- catnoise:

  Number of categorical noise variables

- q:

  Number of categories for signal categorical variables

- q_err:

  Number of categories for categorical noise

- numsep:

  Separation for numerical signal

- catsep:

  Separation for categorical signal

- seed:

  Optional seed

- error_type:

  Error distribution ("normal" or "chisq")

- error_df:

  Degrees of freedom for chi-square noise

- error_scale:

  Scale of noise

## Value

A list with elements `df`, `X_num`, `X_cat`, `y`, `num_cols`, `cat_cols`
