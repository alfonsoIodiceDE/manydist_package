# Gaussian-bandwidth sensitivity smoke — 2026-09-21

The sidecar reproduced every parent Gaussian ARI exactly at the primary
multiplier `c = 1`. It then applied `c = 0.5` and `c = 2` to all five methods
without changing their distances, data, gates, clustering seeds, or any other
setting. All five artifacts completed without a spectral-clustering `ifault`.

## Signal ARI

| Configuration | c | Gated | No interaction | Gower | Modified Gower | Standardized Euclidean one-hot |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| reference | 0.5 | 0.904 | 0.501 | 0.656 | 0.594 | 0.656 |
| reference | 1 | 0.905 | 0.490 | 0.199 | 0.629 | 0.510 |
| reference | 2 | 0.905 | 0.356 | 0.199 | 0.632 | 0.510 |
| numeric medium | 0.5 | 0.904 | 0.479 | 0.696 | 0.531 | 0.801 |
| numeric medium | 1 | 0.904 | 0.475 | 0.552 | 0.534 | 0.598 |
| numeric medium | 2 | 0.904 | 0.396 | 0.426 | 0.490 | 0.493 |
| numeric high | 0.5 | 0.904 | 0.556 | 0.696 | 0.504 | 0.868 |
| numeric high | 1 | 0.816 | 0.689 | 0.522 | 0.515 | 0.574 |
| numeric high | 2 | 0.894 | 0.766 | 0.427 | 0.495 | 0.598 |
| categorical high | 0.5 | 0.904 | 0.501 | -0.012 | 0.595 | -0.009 |
| categorical high | 1 | 0.905 | 0.490 | 0.174 | 0.629 | 0.178 |
| categorical high | 2 | 0.905 | 0.356 | 0.314 | 0.632 | 0.155 |
| joint high | 0.5 | 0.904 | 0.556 | 0.631 | 0.504 | 0.578 |
| joint high | 1 | 0.816 | 0.689 | 0.626 | 0.515 | 0.504 |
| joint high | 2 | 0.894 | 0.766 | 0.628 | 0.492 | 0.548 |

The gated method has the highest ARI in all 15 configuration-by-bandwidth
cells. Its advantage over the no-interaction baseline is positive throughout,
ranging from 0.127 to 0.549. Its smallest advantage over the strongest external
competitor is 0.036, in the `numeric_high`, `c = 0.5` cell, where standardized
Euclidean one-hot obtains 0.868 and the gated method obtains 0.904.

The gated result is essentially unchanged over the bandwidth grid for the
reference, numeric-medium, and categorical-high configurations. In the two
`p = 10` configurations it is non-monotone—0.904, 0.816, and 0.894—but both
sensitivity values improve on the primary value. Therefore the primary
`c = 1` result is not a favourable bandwidth selected for the proposed method.

Competitor behaviour is more bandwidth-dependent. For example, narrowing the
kernel raises Gower from 0.199 to 0.656 in the reference configuration and
raises standardized Euclidean one-hot from 0.574 to 0.868 in the
numeric-high configuration. This confirms that bandwidth must be reported and
that comparisons should use an identical pre-specified rule for every method.

## Interpretation

This exploratory check supports retaining the NJW-style global Gaussian
affinity as the primary analysis with `c = 1`, while reporting `c = 0.5` and
`c = 2` as supplementary sensitivity values. It does not justify selecting a
different multiplier by method, dataset, or replicate. Because this is one
replicate from a deliberately favourable mechanism study, stability of mean
ARI and method rankings must still be evaluated in the confirmatory run.

