
<img src="man/figures/manydist_logo.png" align="left" width="220" alt="manydist logo"/>

# manydist

**Distance-Based Learning for Mixed-Type Data**

`manydist` provides tools for constructing, computing, and using
distance measures for numerical, categorical, and mixed-type data.

<br clear="left"/>

## Installation

``` r
# install.packages("palmerpenguins")
install.packages("manydist")
```

You can also install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("alfonsoIodiceDE/manydist_package", subdir = "manydist")
```

### Initialization

``` r
library(manydist)
library(palmerpenguins)
library(dplyr)
library(tidyr)
```

## Getting started

The main function in `manydist` is `mdist()`. It computes
dissimilarities for numerical, categorical, and mixed-type data.

``` r
penguins_small <- penguins |>
  select(-year,-sex) |> drop_na()
```

A mixed-type dissimilarity can be computed with a preset. For example,
`preset = "gower"` computes a Gower-type dissimilarity.

``` r
D <- penguins_small |> select(-species) |> 
  mdist(preset = "gower")

as.matrix(D$distance)[1:5,1:5] |> round(digits=2)
```

    ##      1    2    3    4    5
    ## 1 0.00 0.05 0.10 0.09 0.08
    ## 2 0.05 0.00 0.08 0.11 0.10
    ## 3 0.10 0.08 0.00 0.08 0.11
    ## 4 0.09 0.11 0.08 0.00 0.07
    ## 5 0.08 0.10 0.11 0.07 0.00

Some presets can use a response variable. When a response is supplied,
it is used for response-aware distance construction and is not treated
as an ordinary predictor.

``` r
D_resp <- penguins_small |> mdist(
  response = species,
  preset = "u_dep"
)
```

    ## Warning: 'tvd' requires >1 categorical variable. Switching to 'matching'.

    ## Warning: For method(s) 'matching', category dissimilarities do not depend on
    ## conditional profiles; `response` was therefore ignored.

``` r
D_resp
```

    ## MDist object
    ##   preset : u_dep 
    ##   number of observations : 342 
    ##   number of continuous variables   : 4 
    ##   number of categorical variables   : 1 
    ##   parameters:
    ##     - categorical method: matching
    ##     - numerical preprocessing: pc_scores
    ##     - commensurability adjustment: TRUE

Custom specifications can be defined by combining a categorical method,
a numerical preprocessing method, and a commensurability rule.

``` r
D_custom <- mdist(penguins_small |> select(-species),
                  preset = "custom",
                  method_cat = "matching",
                  method_num = "std",
                  commensurable = TRUE
)

D_custom
```

    ## MDist object
    ##   preset : custom 
    ##   number of observations : 342 
    ##   number of continuous variables   : 4 
    ##   number of categorical variables   : 1 
    ##   parameters:
    ##     - categorical method: matching
    ##     - numerical preprocessing: std
    ##     - commensurability adjustment: TRUE
    ##     - number of principal components: 5
    ##     - inertia threshold: NULL
