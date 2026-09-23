# R1 full-manifest smoke review — 2026-09-19

The nine pre-specified representative tasks were run with reduced sample
sizes (n = 90 for the submitted-study families and Simulation 2; n = 36
for the moon scenarios). All nine completed, with no warnings, runner errors,
method-level failures, or spectral-clustering ifault values. The manifest
contains 1,650 tasks; the other 1,641 were **not** run.

Results are in results/smoke/59e3500bd1/; merged CSVs are in
results/merged/59e3500bd1/smoke/. These are one-replicate software and
design diagnostics, not estimates of expected clustering performance.

| Task | Scenario | Gates detected | R1 interaction vs no-interaction ARI, Gaussian | R1 interaction vs no-interaction ARI, self-tuned |
| --- | --- | ---: | ---: | ---: |
| 1 | corrected marginal overlap | 0/4 | 0.343 vs 0.343 | 0.248 vs 0.248 |
| 251 | nonlinear interaction | 1/2 | 0.929 vs 0.469 | 0.929 vs 0.794 |
| 301 | null interaction | 0/2 | 0.784 vs 0.784 | 0.784 vs 0.784 |
| 851 | nearest-neighbour sensitivity, signal | 1/2 for each of four settings | 1.000 vs 0.398 for every setting | 1.000 vs 1.000 for every setting |
| 1151 | submitted Simulation 1, balanced | 15/15 | 0.794 vs 0.692 | 0.967 vs 0.967 |
| 1251 | submitted Simulation 1, categorical dominant | 20/20 | 1.000 vs 1.000 | 1.000 vs 1.000 |
| 1351 | submitted Simulation 1, numeric dominant | 0/10 | 0.666 vs 0.666 | 1.000 vs 1.000 |
| 1451 | submitted Simulation 1, numeric noise | 7/15 | 0.509 vs 0.475 | 0.710 vs 0.710 |
| 1551 | submitted Simulation 2 | 2/3 | 0.281 vs 0.277 | 0.258 vs 0.258 |

The controlled moon signal/null pair supports the intended *mechanism*: the
gate activates with designed cross-type signal and stays off in the matched
null, where the two distances and clustering results coincide. It does not
establish general superiority: modified Gower reached ARI 1.000 on the
nonlinear-signal smoke task.

The submitted-study continuity scenarios are **not** faithful reproductions
of the historical code. They deliberately implement the submitted manuscript
specification. In all four Simulation 1 smoke datasets, every cluster-2
observation had level 1 and every cluster-3 observation level 0 in every
categorical variable (27 and 45 rows, respectively). The manuscript's
requested cluster-3 equicorrelation of -0.5 is indefinite for 10–20
variables. Negative-eigenvalue clipping produces a singular draw whose
continuous-variable sum is fixed, contributing to categorical collapse.
The historical code generated categories from unshifted continuous draws
before adding the displayed cluster mean shifts; its sample sizes and
predictor indexing also differ from the manuscript. Consequently, high
Simulation 1 ARIs here cannot be taken as independent confirmation of the
method or as a reproduction of the submitted tables.

Simulation 2 is the more serious efficacy warning. The two relevant category
variables passed the gate, but the gated term's mean contribution was about
0.189, compared with the numerical block's target mean of 6. The ARI gain was
negligible at this smoke size. On the same dataset DKSS with Gaussian affinity
reached 0.954 and modified Gower with self-tuned affinity reached 0.819.
The historical Simulation 2 scripts also differ from the manuscript
(coefficient 0.5 rather than 0.1, intercept 8 rather than 3, and four rather
than three categorical variables).

## Decision before full launch

Technical smoke check: **pass**. Scientific/design launch gate: **hold**.
Before a 1,650-task run, decide whether the submitted-study reruns should
follow the manuscript specification, the historical executable code, or both
as explicitly labelled sensitivity analyses. Then run a small original-size
pilot for Simulation 2 and one timing/memory pilot at the largest Simulation 1
size. Any interaction-weight or block-weight sensitivity should be declared
before interpreting the resulting ARIs, not tuned retrospectively to this
single smoke replicate. The pending comparator implementations documented in
the README are a separate launch gate.
