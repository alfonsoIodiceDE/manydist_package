# Simulation 3: RSS26 moon-scaling mechanism study

This is a paper-local, higher-dimensional extension of the RSS26 interaction
moons demonstration. It is separate from the frozen R1 full-study manifest
and does not modify the package implementation.

The study asks a deliberately narrow question: when the six clusters are
defined by the cross-product of moon geometry and an arc-band category, does
the gated interaction term improve spectral clustering, remain selective as
the numerical and categorical dimensions grow, and reduce exactly to the
package no-interaction baseline in a matched cross-type null?

## Frozen method comparison

The no-interaction ablation is a direct package call:

```r
manydist::mdist(
  x,
  preset = "custom",
  method_num = "pc_scores",
  method_cat = "tvd",
  commensurable = FALSE,
  interaction = FALSE,
  ncomp = p_numeric
)
```

The gated method begins with that exact matrix. It calls the package's
internal interaction-delta routine through the existing revision bridge,
applies the pre-specified maximum permutation gate, averages only the active
observation-level interaction components, and adds that average with weight
`gamma * mean(D_num)`. If no gate opens, it returns the package baseline
without numerical modification. The suite never calls
`mdist(interaction = TRUE)`.

The smoke also contains the three distance competitors used in the submitted
Simulation 1 and 2 tables:

- Gower (`preset = "gower"`);
- modified Gower (`preset = "mod_gower"`); and
- standardized Euclidean distance after one-hot encoding
  (`preset = "euclidean"`).

The last method is labelled explicitly because the historical experiment code
called it "naive", whereas the manuscript prose describes a different
additive Euclidean-plus-matching construction. DKSS is not part of this smoke:
it was omitted from the submitted simulation tables because its runtime was
prohibitive. K-prototypes and other real-data comparators can be evaluated in
the confirmatory study after their contracts are frozen.

## Generator and configurations

The first two numerical columns are the noisy two-moon coordinates. Additional
columns are noisy smooth projections of the same latent coordinates, so the
dimension sequence is nested without copying the response labels. The first
categorical variable is the three-level arc band. Up to four three-level
nuisance variables are generated with an OA(9, 4, 3, 2) within every
moon-by-band cell, giving exact pairwise categorical independence.

The signal and null datasets are paired. The null applies one common row
permutation to the complete categorical block; numerical data, truth,
categorical margins, and the joint distribution inside the categorical block
are unchanged. Five nested configurations cover `(p, q) = (2, 2), (6, 2),
(10, 2), (2, 5), (10, 5)`.

## Reproduction

From this directory, run the contract tests:

```sh
Rscript tests/run_tests.R
```

Run one immutable smoke artifact at a time:

```sh
Rscript run_task.R --configuration reference --state signal --replicate 1
Rscript run_task.R --configuration reference --state null --replicate 1
```

Valid configuration names are `reference`, `numeric_medium`, `numeric_high`,
`categorical_high`, and `joint_high`. Existing artifacts are never
overwritten. After all ten signal/null tasks exist, validate and aggregate
them with:

```sh
Rscript summarize_smoke.R
```

The complete frozen settings are in `design.yml`. The single-replicate smoke
is a software and mechanism check, not an estimate of expected performance.

The separately versioned `bandwidth_sensitivity/` sidecar leaves these
artifacts unchanged and evaluates the global Gaussian scale multipliers
`c = 0.5, 1, 2`, with `c = 1` retained as the primary analysis.

`CONFIRMATORY_AFFINITY_DECISION.md` records the subsequent decision to use
only the global NJW Gaussian construction in the confirmatory experiment. The
self-tuned smoke is retained as an exploratory diagnostic, not as a
confirmatory analysis.
