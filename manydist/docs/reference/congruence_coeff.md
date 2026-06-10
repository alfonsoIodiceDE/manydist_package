# Congruence coefficient between two configurations

Computes the congruence coefficient between two data configurations
using the Frobenius inner product of their pairwise distance matrices.

## Usage

``` r
congruence_coeff(L1, L2)
```

## Arguments

- L1:

  A numeric matrix or data frame (rows = observations)

- L2:

  A numeric matrix or data frame with the same number of rows as `L1`

## Value

A scalar in \\\[-1,1\]\\ measuring similarity between the two
configurations.

## Details

The congruence coefficient is defined as \$\$ \frac{\langle D_1, D_2
\rangle_F} {\sqrt{\langle D_1, D_1 \rangle_F \langle D_2, D_2
\rangle_F}} \$\$ where \\D_1\\ and \\D_2\\ are the pairwise distance
matrices derived from `L1` and `L2`.
