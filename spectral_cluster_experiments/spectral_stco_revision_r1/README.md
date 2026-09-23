# Statistics & Computing revision: gated active-average distance

This directory freezes the paper-specific composition used for the first
revision. It does **not** add a public `manydist` preset or change the package
vignette. Reusable distance primitives remain in `manydist`; the code here
only orchestrates them for this manuscript and returns the intermediate
objects needed to answer the reviewers.

## Mathematical specification

Let `p` be the number of numerical variables and `q` the number of nominal
variables.

The no-interaction baseline is the unmodified package result

```r
manydist::mdist(
  x = data,
  preset = "custom",
  method_num = "pc_scores",
  method_cat = "tvd",
  commensurable = FALSE,
  interaction = FALSE,
  ncomp = p
)
```

where `p` is the number of numerical variables. Separate package calls on the
numerical and categorical blocks give `D_num` and `D_cat`, and the code asserts
that their sum is exactly the full package baseline `D_base`.

For nominal variable `j`, the existing `manydist` KNN routine supplies the
prior-corrected balanced-accuracy interaction level dissimilarity derived from
`D_num`; mapping it to observations gives `D_int,j`. A gate indicator `g_j` is
obtained with a joint-label permutation test:
the rows of the complete categorical block are permuted together, preserving
all categorical--categorical associations, and the null statistic is the
maximum KNN separability over every categorical variable and category pair.
The primary design uses 99 permutations and family-wise level 0.01.

Let `A = {j : g_j = 1}`. The interaction contribution is averaged over the
selected variables only and placed on the empirical numerical-block scale:

```text
D_int,A = |A|^{-1} sum_{j in A} D_int,j
s_num   = mean of the positive upper-triangular entries of D_num
D_R1    = D_base + gamma s_num D_int,A
```

The primary setting is `gamma = 1`. If `A` is empty, `D_int,A` is the exact
zero matrix and `D_R1` is bit-for-bit identical to `D_base`. Averaging only
over active variables prevents gated-out nuisance factors from diluting a
detected component. This is the construction used by Simulation 3 and the
aligned full runner. The earlier `u_dep_bw + rho/q` implementation is retained
only as historical audit code and is not an enabled revised-paper method.

## Usage

From the repository root:

```r
source("spectral_cluster_experiments/revision_r1/R/manydist_bridge.R")
source("spectral_cluster_experiments/revision_r1/R/distance_r1.R")
source("spectral_cluster_experiments/revision_r1/moon_scaling/R/gated_distance.R")

fit <- scmix_gated_active_distance(
  data = analysis_data[c("x1", "x2", "x3", "c1", "c2")],
  gamma = 1,
  prop_nn = 0.10
)

fit$distance
fit$distance_no_interaction
fit$components$interaction_level_delta$c1
fit$gate
fit$diagnostics
```

`distance_no_interaction` is the pure package call and is the paired ablation.
When every permutation gate is zero, `distance` is bit-for-bit the same matrix.

Pass the analysis columns explicitly whenever the data also contain truth,
identifier, or grouping columns. Missing data and ordered factors are rejected
in this pilot implementation rather than handled silently.

## Provenance and internal API boundary

The configuration pins `manydist` version 0.5.2 and source commit
`ab27241f1a0846bb13087f804e23fd52746d1578`. The bridge records the installed
package path/version, detected source commit, and source-tree status. Calls to
the non-exported `cat_delta()` and `delta_int_knn()` functions are isolated in
`R/manydist_bridge.R`; this is the only file that should change if those
internal APIs change.

Run the lightweight contract tests with:

```sh
Rscript spectral_cluster_experiments/revision_r1/tests/run_tests.R
```

## Experiment suites

`moon_scaling/` contains Simulation 3 and the mechanism checks that motivated
the active-average construction. `full/` is the aligned multi-scenario runner;
its manifest, task artifacts, merger, and tests all use the same construction.

`pilot/` and `pilot_code_reconstruction/` are retained as dated audit trails of
earlier candidate formulas. Their readouts and results must not be described
as the final revised method or combined with artifacts from the aligned
runner. Configuration hashes keep these generations separate mechanically.
