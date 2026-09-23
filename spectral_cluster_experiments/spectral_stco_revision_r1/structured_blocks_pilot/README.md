# Structured-block pilot for revised Simulation 3

This pilot is separate from the historical Simulation 1 audit and from the
full-run manifest.  It evaluates a proposed replacement with correlated
numerical variables and associated multilevel categorical variables.

The 24 categorical variables comprise eight binary, eight three-level, and
eight four-level factors.  Within each arity, noisy shared nominal labels are
arranged in different four-variable blocks in each association regime.
The group-specific block arrangements produce distinct association
structures, while modest group-specific multinomial profiles make those
structures identifiable from individual observations.  Numerical and
categorical random variables are still generated independently conditional
on their respective regime labels; there is no direct numerical-to-
categorical generating equation.

The paired arms contain exactly the same observations within each block:

- `aligned` matches numerical group `G` to categorical association regime
  `H = G`;
- `decoupled` assigns every numerical group equal numbers from all three
  categorical regimes.

Thus, only cross-block alignment changes.  The pilot uses the agreed package
baseline, gated active-average construction, global Gaussian affinity, and
package implementations of Gower, modified Gower, and Euclidean one-hot.

Run from the repository root:

```sh
Rscript spectral_cluster_experiments/revision_r1/structured_blocks_pilot/tests/run_tests.R
Rscript spectral_cluster_experiments/revision_r1/structured_blocks_pilot/run_pilot.R
```

The five-pair result and the generator-selection audit are recorded in
`PILOT_RESULTS_2026-09-21.md`.  The proposed confirmatory design is frozen in
`FULL_EXPERIMENT_PLAN.md`; the active full-run manifest has not yet been
expanded, so no confirmatory tasks have been launched.
