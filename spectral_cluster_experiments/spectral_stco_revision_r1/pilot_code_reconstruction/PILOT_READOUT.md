# Reconstructed Simulation 2 paired pilot — 2026-09-19; corrected 2026-09-20

**Correction (2026-09-20):** The first historical-comparator sidecars
(`historical_feb2026_sidecars/`) normalized categorical profiles by column
rather than by row. Their ARIs and the earlier readout are invalid. This
report uses only `historical_feb2026_sidecars_v2/`, after an imbalanced-level
regression test and an n = 500 exact comparison against historical
`mdist()` (maximum absolute distance difference 3.6e-15). The original
pilot R1 gated/ungated results did not change.

Ten paired n = 500 replicates were completed for each of two explicitly
labelled reconstructions of the available historical Simulation 2 script.
Each generated dataset was evaluated with the revised gated AB-BW distance,
its exact no-interaction ablation, and a current-package legacy-style
diagnostic. A separate sidecar evaluated the distance actually computed by
the available February 2026 source for the historical Simulation 2 call:
Euclidean distance on standardized full PC scores plus categorical TVD.
Both Gaussian and self-tuned affinities used the same pilot clustering
engine and seeds. Generator and historical-distance tests passed; all 20
n = 500 tasks and their sidecars completed without errors. A separate
n = 100 run was only a software sanity check and is excluded below.

| Index rule | Affinity | R1 gated mean ARI | R1 no-interaction mean ARI | February-source distance mean ARI | Paired gated minus no-interaction |
| --- | --- | ---: | ---: | ---: | ---: |
| literal script | Gaussian | 0.2752 | 0.2751 | 0.2712 | +0.0001 |
| literal script | self-tuned | 0.3032 | 0.3196 | 0.2746 | -0.0164 |
| corrected proportions | Gaussian | 0.2522 | 0.2521 | 0.2491 | +0.0001 |
| corrected proportions | self-tuned | 0.2876 | 0.3065 | 0.2531 | -0.0188 |

The gate detected C1 and C2 in 10/10 replicates of each variant and did not
detect nuisance C3 or C4. Thus failure to improve is not because the gate
remained off. At the default rho = 1/4, the mean added gated term was about
0.148 (literal) or 0.145 (corrected), compared with the numerical block's
mean of 6. The revised interaction signal is present but relatively small
in the final distance.

This pilot **does not support claiming that the new method wins** on the
reconstructed Simulation 2 design. In the literal variant, the paired
Gaussian gain was positive for two datasets, exactly zero for seven, and
negative for one; the self-tuned difference was positive for three, zero
for one, and negative for six. The corrected variant was similar. The R1
gated distance exceeds the February-source distance by only 0.0031–0.0040
mean ARI with Gaussian affinity and by 0.029–0.035 with self-tuning, but
the R1 *ungated* distance exceeds both with self-tuning. Thus the apparent
self-tuned improvement cannot be attributed to the gated interaction.

## Commensurability-only diagnostic (2026-09-20)

The surviving Simulation 2 script sets `commensurable = FALSE` in “Udep
Int.” but `TRUE` in “Udep.” A separate exploratory script uses `FALSE` for
both, so it cannot establish how the archived 50-replicate results were
created. To isolate this setting, we directly called the February source
`mdist()` on the **same** ten datasets per index rule, holding geometry,
affinity, truth, and clustering seeds fixed. We did not use the older
`distance_by_method` wrapper.

| Generator | Euclidean/PC, comm. off | Euclidean/PC, comm. on | Manhattan/std, comm. off | Manhattan/std, comm. on |
| --- | ---: | ---: | ---: | ---: |
| literal script, Gaussian mean ARI | 0.2712 | 0.2993 | 0.2760 | 0.3043 |
| corrected proportions, Gaussian mean ARI | 0.2491 | 0.2446 | 0.2545 | 0.2670 |

For Euclidean/PC with self-tuned affinity, commensurability on also raised
mean ARI from 0.2746 to 0.2970 (literal) and from 0.2531 to 0.2694
(corrected). Thus switching it **off is not a general explanation for higher
ARI in this pilot**; it helped only slightly in the corrected/Gaussian
Euclidean/PC cell. It changes block balance substantially: in the first
literal dataset, the categorical block's mean pairwise contribution was
0.158 with commensurability off versus 4.008 on, while the numerical
means were 3.149 versus 3.018. This is a post-hoc diagnostic, not a
replacement for a pre-specified revision study.
In particular, reproducing the *two different settings written in the
surviving script* gives 0.2712 for its “Udep Int.” distance versus 0.3043
for its “Udep” distance on the literal/Gaussian pilot—opposite to the
archived ranking. That confirms these pilot datasets do not recreate the
archived experiment, rather than identifying a single cause of its gap.

There is a more consequential historical comparison problem. The available
Simulation 2 script labelled its first arm “Udep Int.”, but called the
February package's custom Euclidean/PC-score branch, where the
`interaction = TRUE` argument was ignored. A direct test against source
commit `79563a7` found exact equality between `interaction = TRUE` and
`FALSE` for this call. The labelled “Udep” arm also changed numerical
geometry, scaling, and commensurability, so the archived contrast does not
isolate an interaction effect. The saved 50-replicate driver and its exact
package environment have not been recovered; therefore this code-path
finding must not be represented as proof of what produced the archived
0.55/0.52 table. The later saved script additionally passes `prop_nn` and
`score`, which the February `mdist()` signature did not accept.

Interpretation is limited in three ways:

1. The exact archived 50-replicate generator is missing. The surviving
   script has an indexing error; the two variants bracket that uncertainty
   but do not prove which generated the submitted results.
2. Base-R sampling matches the statistical design but not the unavailable
   SimDesign/sampling random-number streams seed-for-seed.
3. The February-source distance is an exact reconstruction of the available
   February package *distance formula* on these pilot data, verified against
   the old numerical/categorical component functions. It is **not** an exact
   reproduction of the archived experiment, because its 50-replicate driver,
   package environment, and random-number streams remain unverified. Its
   pilot ARIs around 0.25–0.28 must not be compared directly with the
   submitted table's 0.55/0.52. The current-package legacy-style arm is
   still retained in the original pilot artifacts as a diagnostic but is
   omitted from the table because it is not the submitted distance.

The isolated first n = 500 replicate took about 12 seconds; median elapsed
time across the ten literal and ten corrected tasks was about 20 and 21
seconds, respectively, with three concurrent R workers. This does not
benchmark the largest n = 2,000 study tasks.

**Next decision:** establish what produced the archived table (especially
the 50-replicate driver and package version) before treating it as evidence
for interaction. Independently, the small relative weight of the gated
term warrants a theory-led, pre-specified sensitivity analysis or more
modest interaction claims; do not tune the primary method against these
pilot ARIs. The full 1,650-task revision study should not launch on the
assumption that the interaction method is already validated.
