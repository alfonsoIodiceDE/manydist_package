# Permutation sensitivity: 199 versus 499

## Design

The same 25 pilot datasets, data seeds, clustering seeds, methods, and primary
global bandwidth were rerun. The only change was increasing the gate from 199
to 499 permutations. Two workers completed the run in 1,500.4 seconds.

## Stable parts of the study

- Simulation 1 had identical active sets and gated ARIs at 199 and 499
  permutations in all ten signal/null tasks.
- Simulation 2 had identical active sets and gated ARIs in all five tasks.
- Every Simulation 1 null and Simulation 3 decoupled task kept an empty gate
  and recovered the package baseline exactly at both permutation counts.
- All non-gated competitor results were necessarily unchanged because the
  datasets and clustering seeds were identical.

Thus 199 permutations were sufficient for the Simulation 1 and Simulation 2
pilot decisions.

## Simulation 3 aligned sensitivity

| Replicate | Active at 199 | Active at 499 | Active-set Jaccard | Gated ARI at 199 | Gated ARI at 499 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 2 | 4 | 0.500 | 0.242 | 0.392 |
| 2 | 3 | 5 | 0.600 | 0.231 | 0.381 |
| 3 | 9 | 9 | 1.000 | 0.485 | 0.485 |
| 4 | 6 | 3 | 0.500 | 0.498 | 0.461 |
| 5 | 7 | 6 | 0.857 | 0.491 | 0.506 |

Eight variable-level gate decisions changed across four of the five aligned
tasks. Decisions at 199 permutations were separated only by the adjacent
Monte Carlo values 0.010 and 0.015. At 499 permutations, active-set adjusted
p-values were at most 0.010 and inactive minima were between 0.012 and 0.020.

The higher-resolution gate materially stabilized the clustering result:

| Quantity | 199 permutations | 499 permutations |
| --- | ---: | ---: |
| Mean gated ARI | 0.389 | 0.445 |
| SD gated ARI | 0.140 | 0.056 |
| Mean no-interaction ARI | 0.357 | 0.357 |
| Mean gated minus baseline | +0.032 | +0.088 |
| SD gated minus baseline | 0.170 | 0.095 |
| Positive paired fraction | 0.60 | 1.00 |

## Decision

A single permutation count is computationally wasteful. The evidence supports
a pre-specified hybrid final design:

- 199 permutations for Simulations 1 and 2, where active sets and results were
  identical at 199 and 499;
- 499 permutations for Simulation 3, where the 199-permutation grid was too
  coarse near alpha 0.01.

The final manuscript must report this allocation transparently. It uses the
same family-wise alpha in every experiment; only the Monte Carlo resolution
differs. The full study should continue to emphasize the paired estimand and
exact fallback rather than universal superiority over all competitors.
