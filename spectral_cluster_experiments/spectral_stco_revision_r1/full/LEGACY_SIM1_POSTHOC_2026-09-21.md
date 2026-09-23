# Historical-code Simulation 1: post-hoc interaction contribution

This analysis re-created the three fixed datasets in
`LEGACY_SIM1_PILOT_2026-09-21.md` and used the same frozen settings.  Its
diagnostic quantities use only estimated partitions.  The known simulation
labels are retained solely for the external ARI evaluation and do not enter
the interaction-alignment calculation.

For an estimated partition `c`, interaction alignment is summarized by

```text
R_int(c) = mean(D_int,A between estimated clusters) /
           mean(D_int,A within estimated clusters).
```

Values above one indicate that the accepted interaction component supports
the estimated partition.  Because evaluating `D_int,A` against the gated
partition is partly circular, the primary sensitivity diagnostic removes each
active interaction component in turn, reconstructs and reclusters the
distance without it, and measures the held-out component against that
leave-one-interaction-out partition.

## Replicate results

| Seed | ARI to truth: baseline | ARI to truth: gated | Gated - baseline | ARI between estimated partitions | Interaction alignment: baseline partition | Interaction alignment: gated partition | Mean held-out alignment | Mean leave-one-out ARI to gated partition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 26092201 | 0.766 | 0.727 | -0.039 | 0.843 | 1.385 | 2.060 | 2.064 | 1.000 |
| 26092301 | 0.843 | 0.733 | -0.109 | 0.769 | 1.196 | 2.001 | 2.003 | 0.997 |
| 26092401 | 0.632 | 0.733 | +0.101 | 0.572 | 1.705 | 2.161 | 2.168 | 1.000 |

All 15 interaction components were active in every replicate.  Every one of
the 45 held-out component checks had an alignment ratio above one; the lowest
replicate-level minimum was 1.806.  Removing one component at a time produced
partitions almost identical to the full gated partition: the mean ARI between
the leave-one-out and gated partitions was 0.999 across the three datasets.

The active-average interaction and baseline distances had Spearman
correlations of 0.617, 0.606, and 0.650.  Thus, the interaction geometry is
partly redundant with the baseline but is not simply a rescaled copy of it.
It contributes enough additional geometry to change the fitted partition:
the baseline--gated partition ARIs ranged from 0.572 to 0.843.

## Interpretation

The diagnostic does **not** support describing Simulation 1 as a case in
which the observed interactions fail to contribute to the identified
partition.  The interactions are internally coherent, are strongly aligned
with the gated solution, and collectively define a solution that is highly
robust to removing any single interaction component.

What the three-replicate pilot shows instead is that this coherent interaction
solution does not consistently improve recovery of the generator's declared
labels: it reduces external ARI in two datasets and improves it in one.  This
is the distinction that matters for real unsupervised applications.  A
post-hoc diagnostic can measure contribution, coherence, and robustness with
respect to an identified partition, but without external labels it cannot
establish that the partition is closer to an unknown truth.

The result remains useful for the revised-paper narrative, but the precise
claim should be that Simulation 1 contains genuine, solution-aligned
interaction structure whose incremental value for recovering the declared
generating partition is not systematic.  The three datasets are too few for
a population-level performance or variance claim; the diagnostic should be
summarized over the full replicated experiment before manuscript conclusions
are finalized.

## Reproducible artifacts

The analysis is run by `run_legacy_sim1_posthoc.R`.  It writes:

- `results/posthoc/legacy_sim1_three_seed/replicate_diagnostics.csv`;
- `results/posthoc/legacy_sim1_three_seed/component_leave_one_out.csv`;
- `results/posthoc/legacy_sim1_three_seed/summary.csv`.

