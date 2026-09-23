# Historical-code generator smoke — 2026-09-21

The separately named `legacy_simulation_1` and `legacy_simulation_2`
generators reconstruct the executable data-generating mechanisms in
`spectral_clustering_cristina_generated_data`. They do not replace the
manuscript-specification generators and are not yet enabled in `design.yml`.

## Generator checks

Simulation 1 generated valid datasets for all six historical cells:
`n = 500, 1000` crossed with `(p_numeric, q_categorical) = (15,15),
(10,20), (20,10)`. Every dataset had 30 variables and exact cluster sizes
20%, 30%, and 50%.

Simulation 2 generated valid datasets at `n = 500` and `n = 1000`, with six
numerical and four categorical variables. Its categorical margins were exact:

- `C1`: 40%, 20%, 40%;
- `C2`: 80%, 20%;
- `C3`: 60%, 40%;
- `C4`: 40%, 30%, 30%.

Unit tests verify fixed-seed reproducibility, factor structure, the recursive
numerical equations, truth assignment, and invalid-dimension rejection.

## One-seed Simulation 2 mechanism check

The original-size check used `n = 500`, data seed `26092101`, the global
Gaussian affinity, and the five method roles used for the Simulation 3 smoke.
The proposed method was the construction agreed after the original full-runner
scaffold: a pure `manydist::mdist(preset = "custom", method_num =
"pc_scores", method_cat = "tvd", commensurable = FALSE, interaction = FALSE)`
baseline plus the multiplicity-gated active-component average, scaled by the
mean numerical dissimilarity. The gate used 99 permutations, family-wise
level 0.01, `prop_nn = 0.10`, prior-corrected balanced accuracy, and `gamma =
1`.

| Method | ARI |
| --- | ---: |
| Gated active-average interaction | 0.8131 |
| Pure package no interaction | 0.2458 |
| Gower | 0.2486 |
| Modified Gower | 0.6633 |
| Standardized Euclidean one-hot | 0.2701 |

The gate selected `C1` and `C2`, which define the coefficient matrix and true
groups, and rejected nuisance variables `C3` and `C4`. The added interaction
term had mean distinct-pair dissimilarity 1.363, versus 5.083 for the numerical
block.

This is encouraging but is only one pre-confirmatory seed. At `n = 100`, the
same gate stayed closed and returned the package baseline exactly, showing why
the historical sample size—not a miniature substitute—should be used for the
replicated pilot.

## Runner alignment

The `full/` runner has now been aligned with the active-average construction
used above and in Simulation 3. Its manifest and result schema use `gamma`, the
superseded package interaction method is disabled, and the self-tuned affinity
is excluded from the confirmatory specification. Unit tests verify the
active-average equation and exact recovery of the package baseline when no
gate is active. A complete end-to-end task also passed with pinned package and
source-commit provenance.

Earlier smoke artifacts are retained under their old design hashes as
historical diagnostics; the runner never merges artifacts across hashes.
