# Simulation 2 code-reconstruction pilot

This pilot is deliberately separate from the frozen 1,650-task full manifest.
It reconstructs the *available* historical Simulation 2 script at n = 500,
but does not claim to reproduce the missing 50-replicate driver or its random
number streams exactly. The historical SimDesign and sampling packages are
not installed; the generator uses equivalent base-R Gaussian and fixed-subset
sampling designs.

The primary literal variant preserves the surviving script's unparenthesized
categorical index expressions. At n = 500, this yields C1 counts 200/59/241
and C4 counts 200/89/211. The corrected variant changes only those indices
to implement the intended proportions 200/100/200 and 200/150/150. Both use
six numerical variables, four categorical variables, the script's
coefficient matrix (0.8, 0.5, -0.4), and intercept 8 for coefficient 0.5.

The original three paired distances are:

1. the paper-local AB-BW R1 gated distance;
2. the exact no-interaction ablation of that R1 distance;
3. a **current-package legacy-style diagnostic** using PC scores,
   categorical TVD, interaction, prop_nn = 0.05 and log-loss.

The third arm is **not** the exact submitted method: the historical
implementation requested Euclidean numerical geometry while the current
manydist custom path uses Manhattan and applies interaction. The February
2026 source's custom Euclidean/PC-score branch actually ignored its
`interaction` argument. A separate, source-verified historical-distance
sidecar therefore evaluates Euclidean distance on normalized full PC
scores plus categorical TVD, using the same pilot data, spectral engine,
affinities, and clustering seeds. This is a faithful reconstruction of
that February *distance formula*, not proof that the unrecovered archived
50-replicate driver used exactly that package environment.

Generator test:

    Rscript spectral_cluster_experiments/revision_r1/pilot_code_reconstruction/tests/test_generator.R

One original-size replicate:

    Rscript spectral_cluster_experiments/revision_r1/pilot_code_reconstruction/run_task.R --replicate 1 --index-rule literal --n 500

Task outputs are immutable RDS artifacts under results/<design-hash>/n0500/.
Historical-distance outputs are separate immutable RDS sidecars under
historical_feb2026_sidecars_v2/<design-hash>/n0500/. The earlier
historical_feb2026_sidecars/ are **invalid**: they used column rather than
row normalization of categorical profiles and must not be analyzed. They
were preserved to make the correction auditable. To verify the comparator
against the February source and score one saved dataset:

    Rscript spectral_cluster_experiments/revision_r1/pilot_code_reconstruction/tests/test_historical_feb2026_distance.R
    Rscript spectral_cluster_experiments/revision_r1/pilot_code_reconstruction/score_historical_feb2026_sidecar.R spectral_cluster_experiments/revision_r1/pilot_code_reconstruction/results/59e3500bd1/n0500/literal/replicate-01.rds

No original simulation files, archived results, package presets or full
manifest were changed. See `PILOT_READOUT.md` before interpreting the
archived manuscript table as an interaction effect.
