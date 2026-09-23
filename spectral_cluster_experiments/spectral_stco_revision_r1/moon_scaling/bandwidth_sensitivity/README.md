# Gaussian-bandwidth sensitivity sidecar

This sidecar preserves the original Simulation 3 smoke and its artifacts. It
reconstructs each signal dataset and all five distance methods under the same
source, seed, gate, and clustering contracts, then varies only the global
Gaussian bandwidth:

```text
sigma(c) = c * median of the positive pairwise dissimilarities
c in {0.5, 1, 2}
```

The value `c = 1` remains the pre-specified primary analysis. The other values
are a symmetric, modest sensitivity analysis; they are not candidates from
which a dataset-specific best result may be selected. All methods receive the
same rule and multipliers, evaluated on their own dissimilarity scale.

The runner requires the parent smoke artifact, verifies that `c = 1`
reproduces its Gaussian ARI for every method, hashes all computational sources
including the loaded `manydist/R` files, and rejects a dirty package source
tree. Existing artifacts are never overwritten.

From this directory:

```sh
Rscript tests/run_tests.R
Rscript run_task.R --configuration reference --replicate 1
Rscript summarize_smoke.R
```

Run the task once for each of `reference`, `numeric_medium`, `numeric_high`,
`categorical_high`, and `joint_high`. This remains a one-replicate exploratory
check. The multiplier grid is now suitable for freezing prospectively in the
confirmatory experiment.

