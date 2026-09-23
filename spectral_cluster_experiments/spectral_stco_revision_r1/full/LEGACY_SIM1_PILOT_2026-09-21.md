# Historical-code Simulation 1 pilot — 2026-09-21

This pilot used the code-faithful balanced Simulation 1 generator with
`n = 500`, 15 numerical variables, 15 binary categorical variables, and three
clusters in proportions 20%, 30%, and 50%. It used three fixed data seeds:
`26092201`, `26092301`, and `26092401`.

The comparison was frozen to the five roles used in the Simulation 3 smoke:
the gated active-average method, its pure-package no-interaction baseline,
Gower, modified Gower, and standardized Euclidean one-hot. All methods used
the same global Gaussian affinity with the median positive pairwise distance
as bandwidth. The gate used 99 permutations, family-wise level 0.01,
`prop_nn = 0.10`, prior-corrected balanced accuracy, and `gamma = 1`.

## ARI results

| Method | Seed 26092201 | Seed 26092301 | Seed 26092401 | Mean | SD |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gated active average | 0.727 | 0.733 | 0.733 | 0.731 | 0.004 |
| Pure package no interaction | 0.766 | 0.843 | 0.632 | 0.747 | 0.106 |
| Gower | 0.715 | 0.699 | 0.685 | 0.700 | 0.015 |
| Modified Gower | 0.632 | 0.583 | 0.618 | 0.611 | 0.025 |
| Standardized Euclidean one-hot | 0.609 | 0.580 | 0.347 | 0.512 | 0.143 |

The paired gated-minus-baseline differences were `-0.039`, `-0.109`, and
`+0.101`, for a three-seed mean of `-0.016`. The gated result was much more
stable in these three draws, but three replicates are far too few to make a
variance claim.

## Gate interpretation

All 15 categorical variables were selected in all three datasets, with the
smallest attainable adjusted Monte Carlo p-value (`0.01`). This is not a
false-positive failure: the historical generator constructs every categorical
variable directly from the numerical block. Thus the study commonly called
the "no-interaction" Simulation 1 contains strong real cross-type dependence;
it is not a valid null experiment for the gate.

The historical saved results reinforce the intended role of this simulation:
the `Udep Int.` and `Udep` ARIs are exactly equal in every one of the 50 paired
replicates for all six variable-balance/sample-size cells. Simulation 1 should
therefore be presented as a continuity and broad clustering-performance study,
not as evidence that interaction modelling improves ARI. The matched-null arm
of Simulation 3 is the appropriate check that the revised construction
collapses exactly to its no-interaction baseline when the gate does not detect
cross-type structure.

## Pilot conclusion

The aligned gated method remains competitive and has the highest three-seed
mean among the external distance competitors, but it does not systematically
improve on its no-interaction ablation in this Simulation 1 cell. That is a
scientifically coherent result: the detected dependence is real but largely
redundant with already strong marginal cluster information. A larger paired
run is needed to estimate the small mean contrast and its variability; the
method or weight should not be retuned in response to these three pilot seeds.

The separate estimated-partition contribution analysis is recorded in
`LEGACY_SIM1_POSTHOC_2026-09-21.md`.  It refines the interpretation above:
the interaction geometry is partly redundant with the baseline but is
strongly and reproducibly aligned with the gated partition; its improvement
in recovery of the declared simulation labels is what is not systematic.
