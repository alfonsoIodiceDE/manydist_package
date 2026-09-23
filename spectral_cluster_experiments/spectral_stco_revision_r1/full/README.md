# Enriched R1 experiment runner

This directory contains the launch scaffold for the corrected and enriched
simulation study. It is separate from every historical experiment directory
and never writes to `data/` or `data_old/`.

The canonical manifest is `design.yml`. It expands to 1,650 independent task
artifacts: 50 independently generated datasets per condition. Multiple
K-means starts are used to select one fit within each dataset and are never
counted as independent replicates.

## Frozen paper-local method

The no-interaction distance is delegated to a pure package call:
`manydist::mdist(preset = "custom", method_num = "pc_scores", method_cat =
"tvd", commensurable = FALSE, interaction = FALSE)`. Only the gated extension
is paper-local:

```text
A       = {j : the multiplicity-adjusted gate selects variable j}
D_int,A = |A|^{-1} sum_{j in A} D_int,j, or the zero matrix if A is empty
s_num   = mean of the positive upper-triangular entries of D_num
D_base  = D_num + D_cat
D_R1    = D_base + gamma s_num D_int,A, with gamma = 1 primary
```

`g_j` is a pre-specified permutation gate. Complete categorical rows are
permuted together, and each observed variable-level maximum is compared with
the permutation distribution of the maximum over all variables and category
pairs. The primary configuration uses 99 permutations and family-wise level
0.01. Consequently, when no interaction is detected, `D_R1` and `D_base` are
exactly the same matrix. Averaging only the selected components prevents
irrelevant gated-out categorical variables from diluting a detected signal;
scaling by `s_num` puts the added term on the numerical block's empirical
scale.

## Submitted-study continuity reruns

The final five scenario families rerun the data-generating schemes stated in
the submitted manuscript: the four variable-balance/noise configurations of
Simulation 1 and the cross-type construction of Simulation 2. Together their
two sample-size cells and 50 replicates add 500 tasks. They are appended to
the manifest, so the original task IDs 1--1,150 remain unchanged.

These generators follow the manuscript specification rather than silently
adopting discrepancies in the exploratory historical scripts. In Simulation
1, the historical code formed categorical variables from unshifted latent
numerics, used different predictor indexing, and used sample sizes 500/1,000
rather than the manuscript's 1,000/2,000. Literal use of the shifted values in
the manuscript's Bernoulli model can make the categorical signal nearly
deterministic, so this must be examined when interpreting the rerun. In
Simulation 2, the scripts used 0.5 rather than 0.1, an intercept of 8 rather
than 3, and four rather than three categorical variables. A separate post-hoc
historical-code audit remains appropriate.

That post-hoc audit is now supported by two separately named generator
functions in `R/scenarios.R`: `full_make_legacy_simulation_1()` reconstructs
the mechanisms in `SetupSim*.R`, and `full_make_legacy_simulation_2()`
reconstructs `simulation xBIG.R`. They preserve executable-code details rather
than correcting them silently, including the unshifted numerical predictors
and copied-zero categorical predictor columns in Simulation 1, and the fourth
categorical variable, coefficient 0.5, offset 8, and sorted sampled-index
assignment in Simulation 2. They generate new seeded observations from the
historical mechanisms; the exact old observations cannot be recovered because
the historical random seeds and row-level datasets were not saved. These
generators are tested but intentionally not yet added to `design.yml`, so the
frozen 1,650-task manifest and all existing task IDs remain unchanged.

Simulation 1 also requested equicorrelation -0.5 for the third cluster. That
matrix is not positive semidefinite at the declared dimensions. The continuity
generator records this fact and clips its negative eigenvalues to zero, which
reproduces the effective eigen-based legacy draw without presenting the
submitted covariance as valid.

## Preflight

From the repository root:

```sh
Rscript spectral_cluster_experiments/revision_r1/tests/run_tests.R
Rscript spectral_cluster_experiments/revision_r1/full/tests/run_tests.R
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --count-tasks
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --list-tasks
```

A smoke run uses reduced sample sizes and a configuration-hashed output
directory:

```sh
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 1 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 251 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 301 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 851 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 1151 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 1251 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 1351 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 1451 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 1551 --smoke
Rscript spectral_cluster_experiments/revision_r1/full/merge_results.R --include-smoke
```

## Full execution

Each invocation runs exactly one task and writes it atomically. Existing task
artifacts are left unchanged, so interrupted or parallel runs can be resumed.

```sh
Rscript spectral_cluster_experiments/revision_r1/full/run_task.R --task-id 1
```

On a machine with enough memory, four local workers can be launched with:

```sh
seq 1 1650 | xargs -P 4 -I {} Rscript \
  spectral_cluster_experiments/revision_r1/full/run_task.R --task-id {}
```

Start conservatively with one or two workers if memory is limited. The new
`n = 2,000` continuity tasks and their dense distance/affinity matrices are
the peak-memory cases.

After the task run:

```sh
Rscript spectral_cluster_experiments/revision_r1/full/merge_results.R
```

The merger records missing tasks and method-level failures; it does not hide
or impute them.

## Comparator boundary

The core manifest currently enables the gated active-average method, its pure
package no-interaction ablation, Gower, standardized Euclidean one-hot,
modified Gower, DKSS, and k-prototypes. The superseded package interaction
path is explicitly disabled rather than included as a revised method.
Ahmad--Dey, KAMILA, and a mixed model-based comparator remain disabled with
explicit blockers until their implementations and restart/model-selection
contracts are verified. Thus the core run is launchable, but comparator gate
G4 is not yet complete.
