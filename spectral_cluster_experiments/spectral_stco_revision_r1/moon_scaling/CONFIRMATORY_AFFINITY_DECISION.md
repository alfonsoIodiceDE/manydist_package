# Confirmatory affinity decision — 2026-09-21

The revised experiment will use the global Ng--Jordan--Weiss Gaussian
affinity. For each distance matrix, define the single dataset-and-method-level
scale

```text
sigma_0 = median of all positive upper-triangular pairwise dissimilarities.
```

The primary analysis uses `sigma = sigma_0`. A pre-specified sensitivity
analysis uses `sigma = 0.5 * sigma_0` and `sigma = 2 * sigma_0`, applying the
same rule to every method. These are global bandwidths; no observation-specific
scale is used, and no best bandwidth is selected by method, dataset, or
replicate.

The self-tuned seventh-neighbour affinity evaluated in the exploratory smoke
is retained only as an internal robustness diagnostic. It is not part of the
confirmatory Simulation 3 specification or the primary revised-paper claim.

