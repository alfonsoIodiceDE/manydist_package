# Structured multilevel-block pilot — 2026-09-21

## Frozen pilot construction

The final pilot uses five paired seeds and three balanced groups of 120
observations.  Each dataset has 12 correlated numerical variables and 24
nominal variables: eight binary, eight three-level, and eight four-level.

The numerical groups have different positive-definite loading structures and
modest location differences.  Within each categorical arity, variables share
noisy nominal labels in different group-specific four-variable blocks.  A
moderate group-profile strength of 0.35 makes the association regime
identifiable from an individual row without making the categorical problem
trivial.  Numerical and categorical draws are independent conditional on
their regime labels; no categorical variable is generated directly from a
numerical observation.

The paired arms reuse every numerical row and every categorical row exactly
once:

- `aligned`: numerical group `G` equals categorical association regime `H`;
- `decoupled`: every `G` contains exactly 40 rows from every `H`.

Only the cross-block alignment differs.  The gate used 199 permutations and
family-wise alpha 0.01.  Clustering used the agreed global Gaussian affinity.

## Generator-selection audit

The initial equal-margin Gaussian-copula construction produced categorical-
only ARI near 0.03 with 12 variables.  Increasing to 24--36 variables and
stronger correlations raised it only to approximately 0.09--0.13.  An
equal-margin nominal association-block construction was also unrecoverable
by the additive observation distances.

The final generator therefore retains the group-specific nominal association
blocks and adds modest group-specific multinomial profiles.  Profile strengths
0.20, 0.35, and 0.50 were checked using categorical-only TVD before rerunning
the gate or competitors.  Their two-seed mean aligned ARIs were 0.082, 0.454,
and 0.719, respectively; 0.35 was frozen as the moderate setting.  In every
case the decoupled categorical ARI was approximately zero.

## Gate and nesting results

| Seed | Active variables: aligned | Active variables: decoupled | Exact fallback: decoupled |
| --- | ---: | ---: | :---: |
| 26092501 | 13 | 0 | yes |
| 26092601 | 13 | 0 | yes |
| 26092701 | 14 | 0 | yes |
| 26092801 | 5 | 0 | yes |
| 26092901 | 12 | 0 | yes |

Thus, every aligned replicate activated the interaction construction, whereas
every decoupled replicate returned the pure-package baseline bit for bit.  In
the aligned arm, the active interaction alignment ratio with the identified
gated partition ranged from 1.275 to 1.372.

## Clustering results

| Method | Aligned mean ARI | Aligned SD | Decoupled mean ARI | Decoupled SD |
| --- | ---: | ---: | ---: | ---: |
| Gated active average | 0.487 | 0.069 | 0.376 | 0.035 |
| Pure package no interaction | 0.432 | 0.041 | 0.376 | 0.035 |
| Numerical block only | 0.388 | 0.022 | 0.388 | 0.022 |
| Categorical TVD only | 0.425 | 0.060 | -0.004 | 0.001 |
| Gower | 0.454 | 0.017 | -0.003 | 0.002 |
| Modified Gower | 0.425 | 0.023 | -0.003 | 0.001 |
| Standardized Euclidean one-hot | 0.528 | 0.045 | -0.003 | 0.002 |

The aligned gated-minus-baseline ARI differences were `+0.124`, `-0.018`,
`+0.046`, `+0.077`, and `+0.046`: four gains in five paired datasets, with
mean difference `+0.055` and SD `0.052`.  In the decoupled arm every paired
difference was exactly zero.

Euclidean one-hot has the highest aligned mean in this small pilot, but it
collapses when the categorical association regime is decoupled from the
numerical truth.  The gated construction instead improves its own baseline
when alignment is detected and preserves the baseline exactly when it is not.
This paired robustness result is the primary reason to promote the generator;
the experiment should not be described as showing universal dominance over
every competitor.

This behaviour is explained by the Euclidean preset's geometry.  The 24
categorical variables expand to 72 dummy columns, compared with 12 numerical
columns, and the package standardizes every dummy column.  Across the five
datasets, the mean categorical-block Euclidean distance was 11.91 and the
mean numerical-block distance was 4.75, a ratio of 2.51.  The categorical
block consequently supplied 85.7% of the average squared Euclidean distance.
That emphasis is advantageous when categorical regimes are aligned with the
truth and destructive when the identical categorical observations are
decoupled.  It is therefore a substantive robustness result rather than an
unexpected inconsistency.

## Decision

The pilot supports promotion to the planned full experiment.  It establishes:

1. genuinely multilevel categorical data with nontrivial within-block
   association;
2. moderate recoverability of both blocks in the aligned arm;
3. complete gate closure and exact nesting in the matched decoupled arm;
4. a positive mean interaction contribution in the aligned arm without an
   artificially perfect result.

Five replicates are not enough for inferential performance claims.  The full
paired experiment must determine the frequency and uncertainty of the gated
improvement and retain all competitor results.
