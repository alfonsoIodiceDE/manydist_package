# R1 pilot readout

Run date: 2026-09-18

This is a design and implementation check, not the empirical evidence for the
revision. The run used eight independently generated datasets per scenario,
180 observations in each moons scenario, paired spectral-clustering
initialization streams, and the pre-specified neighbour-proportion grid
`0.05, 0.10, 0.20`.

## Frozen candidate checked by the pilot

The numerical block is the `manydist` `u_dep_bw` construction applied to the
numerical variables: standardized and whitened PC scores, Manhattan distance,
then block scaling to mean distinct-pair contribution `p`.

For categorical variable `j`, the final pilot uses the bounded raw components

```text
D_cat,j = D_tvd,j + rho g_j D_int,j,    rho = 1/q.
```

Neither categorical component is mean-normalized. The binary `g_j` is obtained
from 99 joint-label permutations using the maximum over every variable and
category pair, at family-wise level `alpha = 0.01`. Complete categorical rows
are permuted together, preserving categorical--categorical associations. The
paired no-interaction ablation is `D_num,bw + sum_j D_tvd,j`; when every gate is
zero, the two distance matrices are exactly identical.

## Key mean ARIs

Ranges below are over `prop_nn = 0.05, 0.10, 0.20`.

| Scenario | Affinity | AB-BW R1 interaction | AB-BW R1 no interaction | Reading |
|---|---:|---:|---:|---|
| Interaction moons | Gaussian | 0.950-0.952 | 0.431 | Large interaction gain |
| Interaction moons | Self-tuning | 1.000 | 0.767 | Large interaction gain |
| Null interaction | Gaussian | 0.659 | 0.659 | Exactly unchanged |
| Null interaction | Self-tuning | 0.845 | 0.845 | Exactly unchanged |
| Marginal-signal control | Gaussian | 1.000 | 0.988 | Signal retained |
| Marginal-signal control | Self-tuning | 1.000 | 1.000 | Signal retained |

The interaction result is stable across the three pilot neighbour proportions.
The permutation gate selected the signal band in all 24
replicate-by-neighbour settings, selected no nuisance variable, and selected no
variable in any of the 24 matched-null settings. Consequently all 24 null
interaction distance matrices were bit-for-bit equal to their no-interaction
counterparts. This is a design check over only eight datasets per setting, not
a guarantee of perfect finite-sample selection in the full study.

## Diagnostic lesson from the submitted normalization

The diagnostic reconstruction that normalizes every nonzero interaction
matrix can inflate chance-level KNN variation. In the null scenario its mean
self-tuning ARI is `0.584-0.713`, compared with `0.809-0.832` for the raw R1
blend before gating and `0.845` for both the gated R1 method and the
no-interaction ablation. This is why the enriched experiments use the raw,
nested, permutation-gated composition rather than the separately normalized
submitted formula.

## Recorded failure

The full-data package call `mdist(..., preset = "u_dep_bw")` fails in both
moons scenarios because the exactly null categorical TVD block produces
non-finite values after package-level categorical commensuration. All 96
failures are isolated to that diagnostic comparator; the paper-specific R1
method calls `u_dep_bw` on the numerical block only and completed every fit.
The exact messages are stored in `results/pilot_failures.csv`.

## Decision

Proceed to the enriched experiments, but keep the claim conditional:
cross-type interaction improves clustering when the data contain a strong,
predictively separable cross-type structure; it should not be described as a
universal improvement. The gate gives the desired nested behavior in this
pilot, while its 1% level makes clear that exact null recovery remains
probabilistic in future samples. The full study must retain the
null-interaction and irrelevant-variable controls, the neighbour and `rho`
sensitivity grids, both affinities, and the raw interaction and gate
diagnostics.

## Artifacts

- `run_pilot.R`: reproducible runner.
- `results/pilot_results.csv`: replicate-level results and captured errors.
- `results/pilot_summary.csv`: grouped means, medians, standard deviations,
  and successful-replicate counts.
- `results/pilot_failures.csv`: comparator failures.
- `results/pilot_diagnostics.rds`: fitted component diagnostics.
- `results/pilot_gate_diagnostics.csv`: observed maximum statistics,
  adjusted permutation values, and gate decisions.
- `results/session_info.txt`: R and package session information.
