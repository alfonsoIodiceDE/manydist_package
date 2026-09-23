# Full experiment plan following the structured-block pilot

## Manuscript order

Internal generator names should remain stable; only manuscript display labels
change.

1. Revised Simulation 1: the controlled interaction/matched-null mechanism
   study currently stored under `moon_scaling`.
2. Revised Simulation 2: the code-faithful historical Simulation 2 continuity
   study.
3. Revised Simulation 3: the paired structured multilevel-block study defined
   here.

The historical Simulation 1 reconstruction remains an audit/supplementary
analysis and is not silently relabelled as the new generator.

## Revised Simulation 3 confirmatory cells

Freeze the pilot generator before launching.  The primary cell is:

- three balanced groups;
- 12 numerical variables;
- 24 categorical variables with levels 2, 3, and 4;
- numerical loading 0.72 and mean scale 0.65;
- categorical copy probability 0.85;
- categorical profile strength 0.35;
- aligned and exactly decoupled paired arms;
- 50 common-random-number replicates per arm.

Add two prespecified robustness blocks:

1. sample size: 120 and 240 observations per group at profile strength 0.35;
2. categorical profile strength: 0.20, 0.35, and 0.50 at 120 observations per
   group.

After removing the duplicated primary cell, this is eight conditions and 400
task-level replicates.  Every aligned/decoupled pair must reuse the same block
observations and differ only in their alignment.

## Frozen analysis

- Pure package baseline: `preset = "custom"`, `method_num = "pc_scores"`,
  `method_cat = "tvd"`, `commensurable = FALSE`, `interaction = FALSE`.
- Proposed method: multiplicity-gated active interaction average, `gamma = 1`
  and `prop_nn = 0.10`.
- Confirmatory gate: 999 permutations, family-wise alpha 0.01, complete-row
  categorical permutations, prior-corrected balanced accuracy.
- Spectral clustering: global Gaussian affinity with the median positive
  off-diagonal distance, known `k = 3`.
- A modest bandwidth sensitivity at multipliers 0.5, 1, and 2 should be run on
  the primary cell, not multiplied over every robustness condition.

Core competitors remain Gower, modified Gower, standardized Euclidean
one-hot, DKSS, and k-prototypes, using package implementations wherever
available.  Pending external implementations should remain disabled unless
their provenance and parameter contract are resolved before launch.

## Required outputs

For every replicate retain:

- ARI for every method;
- active gate count and adjusted Monte Carlo p-values;
- exact-baseline nesting indicator;
- gated-minus-baseline ARI and estimated-partition agreement;
- interaction alignment with the identified gated partition;
- empirical categorical marginal differences and association-matrix
  separation;
- runtime and method failures.

Primary summaries are paired aligned-versus-decoupled contrasts, activation
frequency, exact-nesting frequency, and paired gated-minus-baseline ARI with
Monte Carlo uncertainty.  The post-hoc alignment measure is descriptive and
must not be used to retune the gate or `gamma`.

## Execution organization

Each arm and replicate should be one resumable task with a configuration hash.
Run paired seeds on the same worker batch where practical.  The 24-variable,
999-permutation gate is the dominant cost, so the full runner should parallelize
at task level and merge only artifacts with the same frozen configuration
hash.  Run one end-to-end smoke task from each arm before dispatching the full
array.

