# Simulation 3 moon-scaling smoke readout — 2026-09-21

All ten pre-specified signal/null tasks completed for design tag `71a554b19e`
at `n = 360`. There were no spectral-clustering `ifault` values. This is one
paired replicate per configuration: it is a software, gate, and mechanism
check, not an estimate of expected ARI.

The comparison contains the same five distance roles as the submitted
Simulation 1 and 2 tables: the interaction method, its no-interaction
ablation, Gower, modified Gower, and the historical standardized-Euclidean
one-hot comparator. The last method was called "naive" by the historical
code, but that label is not used here because the manuscript prose describes
a different additive Euclidean-plus-matching construction.

## Gate and nesting checks

In every signal task, the multiplicity-adjusted gate selected the designed
`band` variable and rejected every nuisance variable. In every matched null,
it selected nothing. The no-active-gate branch returned the package baseline
matrix exactly; the equal null ARIs below are therefore not merely similar
fits from two independently constructed distances.

All five signal detections have adjusted Monte Carlo p-value 0.01, which is
the smallest attainable value with 99 permutations. That is adequate for this
software smoke but too discrete for a stable confirmatory assessment; the
confirmatory permutation count must be frozen at a larger value before launch.

| Configuration | p numeric | q categorical | Signal selection | Null selection | Mean added term / mean numerical block, signal |
| --- | ---: | ---: | --- | --- | ---: |
| reference | 2 | 2 | band only | none | 0.556 |
| numeric medium | 6 | 2 | band only | none | 0.558 |
| numeric high | 10 | 2 | band only | none | 0.565 |
| categorical high | 2 | 5 | band only | none | 0.556 |
| joint high | 10 | 5 | band only | none | 0.565 |

The four added nuisance variables did not dilute the single active component:
the gated results for `(p, q) = (2, 5)` are identical to `(2, 2)`, and those
for `(10, 5)` are identical to `(10, 2)`. This is the intended consequence of
averaging over active gated components rather than all available categorical
variables.

## Signal ARI: Gaussian affinity

| Configuration | Gated | No interaction | Gower | Modified Gower | Standardized Euclidean one-hot |
| --- | ---: | ---: | ---: | ---: | ---: |
| reference | 0.905 | 0.490 | 0.199 | 0.629 | 0.510 |
| numeric medium | 0.904 | 0.475 | 0.552 | 0.534 | 0.598 |
| numeric high | 0.816 | 0.689 | 0.522 | 0.515 | 0.574 |
| categorical high | 0.905 | 0.490 | 0.174 | 0.629 | 0.178 |
| joint high | 0.816 | 0.689 | 0.626 | 0.515 | 0.504 |

The gated method has the highest Gaussian-affinity ARI in every signal cell.
Its paired gain over the package no-interaction baseline ranges from 0.127 to
0.429.

## Signal ARI: self-tuned affinity

| Configuration | Gated | No interaction | Gower | Modified Gower | Standardized Euclidean one-hot |
| --- | ---: | ---: | ---: | ---: | ---: |
| reference | 1.000 | 0.453 | 0.181 | 1.000 | 0.255 |
| numeric medium | 1.000 | 0.456 | 0.558 | 1.000 | 0.368 |
| numeric high | 1.000 | 0.486 | 0.913 | 0.684 | 0.616 |
| categorical high | 1.000 | 0.453 | -0.009 | 1.000 | -0.010 |
| joint high | 1.000 | 0.486 | 0.432 | 0.684 | 0.412 |

The gated method reaches ARI 1.000 in every cell. Modified Gower ties it in
three of the five configurations, which is important to report rather than
presenting the positive control as universal superiority. In the two
high-numerical-dimensional cells, the gated method remains ahead of every
competitor.

## Matched-null check

| Configuration | Gaussian: gated = baseline | Self-tuned: gated = baseline | Active gates |
| --- | ---: | ---: | ---: |
| reference | 0.490 | 0.453 | 0 |
| numeric medium | 0.475 | 0.456 | 0 |
| numeric high | 0.689 | 0.486 | 0 |
| categorical high | 0.490 | 0.453 | 0 |
| joint high | 0.689 | 0.486 | 0 |

Thus the observed interaction gain is specific to the paired cross-type
signal. It is not obtained by paying an unavoidable penalty or bonus in the
matched null. Null ARI is not expected to be near zero: the numerical moons
remain intact and spectral clustering is still asked for six groups. The
negative-control criteria are zero gate activations and exact recovery of the
baseline distance, not poor clustering.

## Interpretation and launch decision

The smoke supports the intended narrow narrative:

1. the package-only no-interaction construction is a genuine nested baseline;
2. the gate detects the designed cross-type signal without selecting the
   balanced nuisance factors;
3. the interaction term improves ARI when the six groups require both moon
   geometry and arc band; and
4. the improvement persists when numerical and categorical dimensions grow.

It does not establish expected performance or general superiority. The data
generator is intentionally a mechanism positive control, the result contains
one random seed, and affinity choice matters for several competitors. The
orthogonal-array construction also makes the categorical TVD block exactly
zero, so this study isolates interaction recovery; it does not test a setting
where marginal categorical and interaction information must both be blended.
A later mixed-signal sensitivity can address that distinct question.

The current 99-permutation implementation took approximately 7–18 seconds per
signal dataset, compared with less than two seconds for the package baseline
at this sample size. These are engineering diagnostics, not fair runtime
benchmarks: the gated wrapper performs separate numerical and categorical
package calls to assert the baseline decomposition before running the gate.
That contract should be optimized and frozen before any manuscript runtime
claim.

**Decision:** the technical and scientific smoke gate passes. A confirmatory
run is warranted with a replicate count and seed schedule frozen before any
additional results are inspected. No method or weight should be tuned in
response to this single replicate. DKSS can remain a separately timed
sensitivity analysis because it was not in the submitted Simulation 1 and 2
tables and was historically omitted for runtime. Before enabling the
confirmatory runner, extend the provenance fingerprint to the loaded package
source files and use a tag derived from every computational source, not only
the design file. The current ten artifacts have one common source-hash set,
match the current suite files, and were produced from a clean `manydist/R`
tree, so this hardening does not invalidate the smoke.

Finally, the active-component average is the candidate revised construction,
not the equation currently written in the manuscript. If it is adopted, the
method section, theoretical rationale, and all reported revised experiments
must use this same frozen definition.

## Subsequent Gaussian-bandwidth check

A separately versioned sidecar subsequently evaluated
`sigma = c * median(D)` for `c = 0.5, 1, 2`, applying the same rule to every
method and retaining `c = 1` as primary. The gated method had the highest ARI
in all 15 configuration-by-bandwidth signal cells. See
`bandwidth_sensitivity/SMOKE_READOUT_2026-09-21.md`. The parent design and the
ten artifacts summarized above were not changed.
