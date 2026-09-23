# Five-replicate spectral-only pilot

## Frozen pilot

- Five Monte Carlo replicates per retained condition.
- Simulation 1: reference signal and matched null.
- Simulation 2: interaction condition at `n = 500`.
- Simulation 3: `n_per_group = 120`, profile strength `0.35`, aligned and
  exactly decoupled arms.
- The gate used 199 permutations and family-wise alpha 0.01.
- Every method used the same global-Gaussian NJW spectral clustering pipeline
  at the primary bandwidth.
- Methods: gated active average, pure-package no-interaction baseline, Gower,
  modified Gower, Euclidean one-hot, and the Ahmad--Dey pairwise
  dissimilarity. DKSS and direct k-prototypes were not included.
- Two task workers completed all 25 tasks in 682.8 seconds.

## Gate and paired checks

| Experiment | State | Active interactions | Exact fallback | Mean gated minus baseline ARI | Positive fraction |
| --- | --- | ---: | :---: | ---: | ---: |
| Simulation 1 | signal | 1 in every replicate | no | +0.400 | 1.00 |
| Simulation 1 | null | 0 in every replicate | yes, 5/5 | 0.000 | 0.00 |
| Simulation 2 | interaction | 2 in every replicate | no | +0.502 | 1.00 |
| Simulation 3 | aligned | 2, 3, 9, 6, 7 | no | +0.032 | 0.60 |
| Simulation 3 | decoupled | 0 in every replicate | yes, 5/5 | 0.000 | 0.00 |

The Simulation 3 aligned paired differences were `-0.116`, `-0.135`,
`+0.024`, `+0.274`, and `+0.113`. The effect is therefore variable rather
than uniformly favorable, but its pilot mean is positive. Every matched null
or decoupled task returned the no-interaction distance and partition exactly.

## Mean ARI

| Method | Sim. 1 signal | Sim. 1 null | Sim. 2 interaction | Sim. 3 aligned | Sim. 3 decoupled |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gated active average | 0.883 | 0.483 | 0.754 | 0.389 | 0.307 |
| Package no interaction | 0.483 | 0.483 | 0.252 | 0.357 | 0.307 |
| Gower | 0.564 | -0.001 | 0.258 | 0.385 | -0.002 |
| Modified Gower | 0.625 | 0.507 | 0.707 | 0.368 | -0.001 |
| Euclidean one-hot | 0.597 | 0.001 | 0.284 | 0.459 | -0.002 |
| Ahmad--Dey distance | 0.801 | 0.365 | 0.377 | 0.377 | -0.001 |

Ahmad--Dey did not win a condition on average. Its apparent win in the
one-replicate null smoke was not stable across the five pilot replicates. The
gated method had the highest mean in both interaction-focused Simulations 1
and 2. Euclidean one-hot had the highest aligned Simulation 3 mean, but it
collapsed when the identical categorical block was decoupled from the truth.

## Decision

All predeclared pilot criteria passed:

1. all null and decoupled gates closed and recovered the package baseline;
2. all signal and aligned gates activated;
3. the gated-minus-baseline mean was positive in Simulations 1 and 2, with a
   gain in every replicate;
4. the Simulation 3 aligned mean contribution was positive, with gains in
   three of five replicates;
5. all 150 method-level evaluations produced finite ARIs and no errors.

The full paired study is justified. Its interpretation must emphasize
conditional improvement and exact fallback, not universal dominance. The
large Simulation 3 variance is a reason to run the pre-specified full Monte
Carlo study rather than a reason to select favorable seeds.
