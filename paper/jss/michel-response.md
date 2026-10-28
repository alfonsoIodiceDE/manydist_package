# Response to Michel’s comments

This is the concise, colleague-facing response log. The complete annotations are preserved in [michel-comments-0708.md](michel-comments-0708.md), and working status is maintained in [michel-comments-tracker.md](michel-comments-tracker.md).

Entries are added here after a comment has been assessed. Unresolved author decisions are explicitly marked `needs discussion`; completed changes are verified against the manuscript. Stable IDs `M-001` through `M-140` match the tracker.

## Progress

- Assessed and documented: 140 of 140
- Resolved or superseded: 98
- Partially resolved: 5
- In progress: 1
- Needs discussion: 36
- Not yet assessed: 0

## Responses

### Abstract — M-001 to M-005

| ID | Response | Status |
|---|---|---|
| M-001 | Replaced the technical property list with the more direct statement that `manydist` constructs mixed-variable distances that account for the identified biases. | resolved |
| M-002 | Recast the contribution list and removed the marked conjunction. | resolved |
| M-003 | Moved pipeline integration into a separate, copyedited sentence describing unsupervised and supervised pipelines. | resolved |
| M-004 | The Palmer penguins example is now explicitly connected to distance construction, diagnostics, and benchmarking. | resolved |
| M-005 | The World Development Indicators example now clearly names unsupervised and supervised pipelines, clustering and classification, and resample-specific preprocessing. | resolved |

**Location:** Abstract.

### Introduction: importance bias — M-007 to M-018

| ID | Response | Status |
|---|---|---|
| M-007 | Deleted the “mirror image” characterization. | resolved |
| M-008 | Removed the marked comparison to a matching categorical variable. | resolved |
| M-009 | The phrase about favouring categorical contributions remains, but the annotation has no explanatory text. We should confirm whether the suggested dominance wording is preferred. | needs discussion |
| M-010 | The containing sentence was replaced; the requested preposition change is no longer applicable. | superseded |
| M-011 | Recast the Gower explanation around a mean numerical distance below one. We should discuss whether the stronger wording “much smaller than one” and “categorical variables tend to dominate” is needed. | partially resolved |
| M-012 | Reworded the sentence to say that skewed numerical distributions amplify the effect and corrected the spelling of “distributions.” | resolved |
| M-013 | Replaced “heavy-tailed” with “long-tailed.” | resolved |
| M-014 | Deleted the marked clause about other scaling choices. | resolved |
| M-015–M-016 | Replaced the containing passage with a new statement about undesirable influences on the overall distance. | superseded |
| M-017 | Adopted the proposed framing and connected it to the general mixed-variable-distance formulation. | resolved |
| M-018 | Presented multivariate additivity and commensurability as two essential properties of the proposed unbiased distances. | resolved |

**Location:** Introduction, importance-bias paragraph and opening of the additivity/commensurability paragraph.

### Introduction: additivity and commensurability — M-019 to M-028

| ID | Response | Status |
|---|---|---|
| M-019 | “Requires” remains; consider “entails” or “specifies.” | needs discussion |
| M-020 | The marked clause explaining that each variable enters additively remains. | needs discussion |
| M-021 | Deleted the anticipatory explanation of block-dependent variable contributions. | resolved |
| M-022 | The wording still begins “the variable-specific distances”; the intended “individual variable” wording needs clarification. | needs discussion |
| M-023 | Deleted the reference to being drawn from a more dispersed distribution. | resolved |
| M-024 | Changed the definition to “are defined to be commensurable.” | resolved |
| M-025 | The sentence retains “so that”; consider splitting it and beginning the consequence with “Hence.” | needs discussion |
| M-026 | Changed “For the general additive formulation” to “In the general additive formulation.” | resolved |
| M-027 | “This is achieved” remains; consider making commensurability the grammatical subject. | needs discussion |
| M-028 | Deleted the qualification about deliberate subject-matter weighting. | resolved |

**Location:** Introduction, additivity/commensurability paragraph.

### Overall organization — M-006

| ID | Response | Status |
|---|---|---|
| M-006 | Reorganized the article around the introduction, a dedicated review of other software, the unified framework, distance construction, diagnostics, learning pipelines, and the conclusion. | resolved |

**Location:** Article-level organization.

### Introduction: package contribution and organization — M-029 to M-038

| ID | Response | Status |
|---|---|---|
| M-029 | Clarified that `manydist` implements unbiased mixed-variable distances as well as distances that do not necessarily belong to that framework. | resolved |
| M-030 | Replaced the unclear reference to “discipline” with a copyedited learning-pipeline contribution and an explicit training-only fit-and-apply explanation. | resolved |
| M-031 | The `mdist()` list now covers distance construction, commensurability, aggregation, and reuse of fitted preprocessing. The placement of association- and response-aware construction remains under discussion. | partially resolved |
| M-032 | Retained a dedicated association-aware subsection. The placement of the first introductory explanation remains unsettled. | partially resolved |
| M-033 | Marked the introductory awareness paragraph for possible removal or relocation. | needs discussion |
| M-034 | “Scale- and type-aware” remains in the paragraph marked for possible removal or relocation. | needs discussion |
| M-035 | Replaced “pairwise benchmarking” with the clearer “pairwise comparisons of candidate distances” and “comparisons between candidate distances.” | resolved |
| M-036 | Created a dedicated `Other software` section. | resolved |
| M-037 | The explicit list of clustering, ordination, and nearest-neighbour methods remains; consider replacing it with “distance-based data analysis methods.” | needs discussion |
| M-038 | Removed the sentence claiming that `manydist` builds on rather than replaces the listed tools. | resolved |

**Location:** End of the Introduction and opening of `Other software`.

### Other software — M-039 to M-060

| ID | Response | Status |
|---|---|---|
| M-039 | Changed “numerical matrices” to “numerical data.” | resolved |
| M-040 | “Extensible catalogue” remains and may need a plainer explanation. | needs discussion |
| M-041 | The `nomclust` description remains; its scope and the Boriah-family distances should be clarified. | needs discussion |
| M-042 | The broader numerical-distance capabilities of `cluster::daisy()` are not yet stated. | needs discussion |
| M-043 | The phrase “same underlying construction” remains without explaining the implementation differences. | needs discussion |
| M-044 | The importance-bias property of Gower remains in the package-comparison paragraph. | needs discussion |
| M-045 | Consider replacing “packages for mixed data” with the more precise “packages implementing mixed-variable distances.” | needs discussion |
| M-046 | The description of the explicitly weighted family remains unchanged. | needs discussion |
| M-047 | References for the Podani, Wishart, Harikumar, and Ahmad formulations still need to be added or the list shortened. | needs discussion |
| M-048 | “Learns an adaptive mixed-data distance” remains; consider Michel’s more direct data-estimation wording. | needs discussion |
| M-049 | “Automatically” remains in the `kproto()` description. | needs discussion |
| M-050 | “Related to this family” remains. | needs discussion |
| M-051 | Singular “a weighting” remains. | needs discussion |
| M-052 | The description does not yet state that `kamila` and `cluspcamix` do not yield a separate distance matrix. | needs discussion |
| M-053 | Removed the containing comparison rather than revising it. | superseded |
| M-054 | Removed “likewise” with the containing sentence. | resolved |
| M-055 | Replaced the generic refitting statement with an explicit `manydist` fit-and-apply description: data-dependent quantities are estimated from the training data only and then held fixed for test-to-training distances. The revision also states why fitting once on the complete data set would cause leakage. | resolved |
| M-056 | Deleted the broad and cryptic novelty claim and replaced it with a concrete resampling explanation. | resolved |
| M-057 | The table description of `daisy()` has not yet been broadened. | needs discussion |
| M-058 | The meaning of “variants” in the software table remains unclear. | needs discussion |
| M-059 | Removed the containing “share a consequence” claim. | superseded |
| M-060 | The three families described in the prose are not yet represented explicitly in the table. | needs discussion |

**Location:** `Other software` section and its two tables.

### Software-review conclusion and roadmap — M-061 to M-069

| ID | Response | Status |
|---|---|---|
| M-061 | Moved the article roadmap before the new `Other software` section. | resolved |
| M-062 | Deleted the marked sentence about the shortcomings of the reviewed tools. | resolved |
| M-063 | Recast the transition as “complements this ecosystem.” | resolved |
| M-064 | Retained the association-aware and response-aware terms and explained them elsewhere; the location of their first explanation is still under discussion. | partially resolved |
| M-065 | The roadmap still uses “develops” rather than the suggested direct section-by-section wording. | needs discussion |
| M-066 | “Fitted preprocessing” remains in the roadmap. | needs discussion |
| M-067 | Replaced “pairwise benchmarking” in the roadmap with the broader “benchmarking of candidate distances.” | resolved |
| M-068 | Recast the roadmap sentence as “shows how distances can be embedded in `tidymodels` workflows.” | resolved |
| M-069 | The notation still uses `I`, `Q_n`, and `Q_c`; whether to align it with the earlier paper remains to be decided. | needs discussion |

**Location:** Conclusion of `Other software`, Introduction roadmap, and opening of the framework.

### Framework organization — M-070 to M-077

| ID | Response | Status |
|---|---|---|
| M-070 | “Contribution” remains rather than “distance contribution.” | needs discussion |
| M-071 | Retained the statement about sums of metrics but added a reference to the proof in van de Velden et al. | resolved |
| M-072 | Replaced “Scale- and type-aware distances” with “Multivariate additivity and commensurability.” | resolved |
| M-073 | Deleted the constrained-weights paragraph. | resolved |
| M-074 | Removed “beyond commensurability” from the association-aware heading. | resolved |
| M-075 | The current purple trial whitens the retained PCA scores and applies a single factor to scale their complete Manhattan distance to mean $Q_n$. The final construction and corresponding package update remain to be agreed. | in progress |
| M-076 | Removed “by itself.” | resolved |
| M-077 | Replaced the earlier opening with a concise distinction between commensurability and within-block association awareness. | resolved |

**Location:** Unified framework, additivity/commensurability and association-aware subsections.

### Whitening and Figure 1 — M-078 to M-081

| ID | Response | Status |
|---|---|---|
| M-078 | Defines whitening as PCA rotation followed by division by $\sqrt{\lambda_h}$, giving uncorrelated retained coordinates with unit variance; the current trial uses it before block scaling. | resolved |
| M-079 | Added an explicit reference to Figure 1 and revised its third panel to display block-scaled whitened coordinates, retaining the concise three-panel layout. | resolved |
| M-080 | Distinguishes PCA decorrelation from whitening and explains that the common block factor does not undo whitening's direction-specific rescaling. | resolved |
| M-081 | Distinguishes PCA rotation, whitening, block-level Manhattan calibration, component-wise empirical scaling, and Mahalanobis distance. | resolved |

**Location:** Framework, association-aware figure and `Numerical variables` subsection.

### Numerical and categorical framework — M-082 to M-085

| ID | Response | Status |
|---|---|---|
| M-082 | Moved the numerical construction and transformation list before multivariate additivity, commensurability, and the association-aware discussion. | resolved |
| M-083 | The indicator-based paragraph remains. We should either remove it or retain only the sentence needed to explain package options. | needs discussion |
| M-084 | Shortened the categorical association passage, corrected the conditional-distribution notation to $\mathbf{R}^{k'\mid k}$, and now states explicitly that identically zero dissimilarities cannot be reciprocal-mean scaled. The package still needs an agreed zero-contribution or fallback rule. | partially resolved |
| M-085 | The response-aware passage still needs reorganization and its bold revision markup must be removed before submission. | needs discussion |

**Location:** Framework, `Numerical variables` and `Categorical variables` subsections.

### Distance construction — M-086 to M-094

| ID | Response | Status |
|---|---|---|
| M-086 | The marked punctuation disappeared when the section opening was rewritten. | superseded |
| M-087 | Introduced the data, numerical preprocessing, categorical dissimilarity, and commensurability arguments—including defaults and main options—before introducing presets. The `MDist` object is described afterwards. | resolved |
| M-088 | `u_dep` remains the principal worked preset, but the annotation contains no written explanation. We should confirm whether a different example was intended. | needs discussion |
| M-089 | Replaced “carries” with “contains” in the rewritten `MDist` description. | resolved |
| M-090 | Removed “consumed by downstream” and now states simply that `to_dist()` returns the standard R `dist` representation. | resolved |
| M-091 | Removed the arbitrary `hclust()` and `cmdscale()` examples and stopped at the returned `dist` representation. | resolved |
| M-092 | Added a dedicated preset subsection and table, with an explanation of the `u_` family and comparison presets. | resolved |
| M-093 | Moved direct component specification before presets and clarified that `preset = "custom"` is already the default. | resolved |
| M-094 | Explained why the presets are included and why some require dedicated implementations outside the custom three-component interface. | resolved |

**Location:** `Distance construction` section.

### LOVO diagnostics — M-095 to M-106

| ID | Response | Status |
|---|---|---|
| M-095 | Replaced the strong recommendation with an objective motivation and organized the diagnostics into direct distance measures and downstream MDS or clustering measures. | resolved |
| M-096 | Replaced “predictor” with “variable” in the LOVO description. | resolved |
| M-097 | Connected the reported diagnostics directly to the changes produced by leaving out each variable and recomputing the distance. | resolved |
| M-098 | Added an itemized explanation of MAD, normalized relative distance, MDS diagnostics, and optional clustering diagnostics. | resolved |
| M-099 | Introduced classical MDS, the selected dimensionality, and the full and reduced configurations before defining congruence and alienation. | resolved |
| M-100 | The malformed notation disappeared with the rewritten passage. | superseded |
| M-101 | Removed the unsupported “variance left unexplained” interpretation and defines alienation directly from the congruence coefficient. | resolved |
| M-102 | Separated direct from downstream diagnostics and explains that they describe different objects; the MDS measures are explicitly conditional on interest in an MDS representation. | resolved |
| M-103 | The annotation on “benchmarking” contains no written explanation. The word remains in the section title because the exported function is `benchmark_mdist()`, while the prose generally uses “comparison.” | needs discussion |
| M-104 | Added an explicit reference to the LOVO figure and explains what it shows. | resolved |
| M-105 | Simplified the plotting call while retaining explicit `metric` and `reorder` arguments to reproduce the displayed quantity and ordering. | resolved |
| M-106 | Removed the claim that the two metrics “disagree” and instead explains that they describe different objects. | superseded |

**Location:** `Distance diagnostics and benchmarking`, opening and LOVO subsection.

### Comparing distance specifications — M-107 to M-132

| ID | Response | Status |
|---|---|---|
| M-107 | Split the material into comparisons of LOVO diagnostics and comparisons of complete candidate distance specifications. | resolved |
| M-108 | Added an explicit reference to the comparative LOVO figure and explains its PAM-based quantities. | resolved |
| M-109 | Removed the “disagree” wording and adopted the direct-versus-downstream distinction. | resolved |
| M-110 | Reorganized candidate comparison to begin with distance-magnitude measures and then separately introduce MDS and clustering comparisons. | resolved |
| M-111 | Avoided the broad magnitude-versus-geometry claim and distinguishes original dissimilarities from specified downstream analyses. | resolved |
| M-112 | Replaced “successful distances” with “distances that could be computed.” | resolved |
| M-113 | Specifies classical MDS, the default dimensionality, and the representations being compared. | resolved |
| M-114 | Moved optional clustering comparisons into a separate paragraph explaining when they are run and how to interpret ARI. | resolved |
| M-115 | Most prose now uses “comparison,” but “benchmarking” remains in the section title and exported function name. Preferred article terminology remains to be agreed. | needs discussion |
| M-116 | Explains directly how candidate specifications are supplied or generated with `all_dist_method_specs()`. | resolved |
| M-117 | Updated the example to use `summary(distance_benchmark)`, which displays and returns the complete pairwise-results tibble. Removed the unnecessary separate extractor from the package. The annotation has no written explanation, so we should confirm whether further change was intended. | needs discussion |
| M-118 | Simplified the code by removing the `dplyr::mutate()` label manipulation. | resolved |
| M-119 | Defines MAD and relative distance explicitly for two candidate distances. | resolved |
| M-120 | Replaced “low-dimensional geometry” with the more precise “chosen classical MDS representations.” | resolved |
| M-121 | Uses “measure agreement” rather than “determines.” | resolved |
| M-122 | Explains that lower ARI means increasingly different assignments and therefore quantifies how much the partitions differ. | resolved |
| M-123 | States that the diagnostics describe different objects rather than alternative estimates or rankings of one quantity. | resolved |
| M-124 | Added a concrete `autoplot()` example, figure reference, and explanation of the displayed pairwise comparison. | resolved |
| M-125 | Removed the implementation-error discussion from the article. | superseded |
| M-126 | Rephrased cautiously that the diagnostics do not, by themselves, select a preferred distance. | resolved |
| M-127 | Removed the ambiguous shorthand with the containing results narrative. | superseded |
| M-128 | Removed the unsupported ratio interpretation of alienation. | superseded |
| M-129 | Defines ARI as agreement between partitions and explicitly states that it does not establish preference. | resolved |
| M-130 | Removed the vague pronoun with the containing results narrative. | superseded |
| M-131 | Removed the unclear comparison with the containing results narrative. | superseded |
| M-132 | Replaced the strong synthesis with a cautious, use-dependent interpretation of the diagnostics. | resolved |

**Location:** Comparative LOVO and candidate-distance subsections.

### Learning pipelines — M-133 to M-140

| ID | Response | Status |
|---|---|---|
| M-133 | Added an introduction explaining why data-adaptive distances must be fitted on fitting data, names the unsupervised and supervised settings, and distinguishes the outer training/test split from the analysis/assessment splits used in resampling. | resolved |
| M-134 | Replaced the general introduction with a worked penguin example: direct `mdist()` distances and `nearest_neighbor_dist()`, then the equivalent `step_mdist()` workflow. Retained the WDI clustering and supervised classification applications and removed the data-leakage experiment at the author's request. | resolved |
| M-135 | Defines a `tidymodels` recipe before introducing `step_mdist()`. | resolved |
| M-136 | Condensed the repeated fit-and-apply explanation and removed the unclear pipeline-level analogy. | resolved |
| M-137 | Names `step_mdist()` explicitly when introducing its `output` argument. | resolved |
| M-138 | The worked example distinguishes square training distances from rectangular test-to-training distances and explains their row and column alignment. The clustering application retains the separate `output = "pairwise"` representation. | resolved |
| M-139 | Introduces supervised nearest-neighbour classification before explaining its required representation and later gives it a dedicated subsection. | resolved |
| M-140 | “Fixed snapshot” remains, and the annotation contains no written explanation. We should confirm whether Michel wanted the phrase removed or clarified. | needs discussion |

**Location:** `Distance-based learning pipelines` section.

<!--
Use this structure for each addressed comment:

### M-000 — PDF page.item

**Comment:** Concise paraphrase of Michel’s concern.

**Response:** What was changed, or why no change was made.

**Location:** Current section/paragraph/line anchor in article.qmd.

**Status:** needs discussion | resolved | partially resolved | already addressed | superseded | declined
-->

### 2026-10-08 — Benchmark API follow-up

The benchmark example now uses the standard `summary()` interface for pairwise results; `benchmark_mdist()` still retains the fitted distances, specifications, and run status. The separate `benchmark_comparisons()` function has been removed, and package examples and documentation have been aligned with this interface.

All package tests passed, and the manuscript was rebuilt and the revised benchmark pages visually checked. M-117 remains open for discussion because Michel's original annotation contains no written explanation.

### 2026-10-08 — Direct kNN example in Section 6.1

Simplified the direct example to calculate only test-to-training Gower distances, with numerical ranges estimated from the training observations. It now uses `fit_knn_dist()` and `predict_knn_dist_class()`; `nearest_neighbor_dist()` is introduced with the complete recipe/workflow example. The unnecessary training-to-training matrix and manual distance-column naming have been removed. The two approaches give identical predictions for all 84 test observations (accuracy 1). This is an author-requested refinement; no Michel-comment statuses changed.

The manuscript was rebuilt successfully and Section 6.1 visually checked on pages 26–28, including the code-block page breaks and the printed `TRUE` prediction comparison.

### 2026-10-08 — Unified direct prediction interface

Section 6.1 now uses the public `knn_dist()` interface rather than separate fit/predict engine calls. It demonstrates both precomputed test-to-training distances and tabular data with `response = "species"`, which extracts the training labels and excludes the response from the training and test predictors. The subsequent `nearest_neighbor_dist()` / `step_mdist()` workflow is retained. Direct precomputed, direct tabular, and workflow predictions agree for all 84 test observations; all five package test files pass. No Michel-comment statuses changed.

The updated package was installed locally, its help examples were executed successfully, and the 40-page manuscript was rebuilt. The revised Section 6.1 was visually checked on pages 26–28, including both printed `TRUE` comparisons.

### 2026-10-08 — LOVO plot ordering refinement

At the author's request, LOVO comparison plots now reorder variables by their mean selected metric across methods within each variable type, preserving the categorical/numerical groups and their blue/red background shading. The single-method LOVO plot follows the same rule. Global `top_n` selection is unchanged. Package help and regression tests were updated; all five test files pass and the reordered comparison was visually checked. No Michel-comment statuses changed.

### 2026-10-08 — Specification-table simplification and article review

Removed the redundant `spec_type` field: named presets select predefined settings, and `preset = "custom"` selects the explicit numerical, categorical, and commensurability settings. Updated the benchmark paragraph, source vignettes, help, and tests. The complete grid contains 482 candidates after removing the redundant generic-custom row. All five test files pass and the locally installed package and 34-page rendered manuscript are updated. The author's other article edits were preserved; only two equation-reference colour-markup defects were additionally repaired.

The review in `article-review-2026-10-08.md` records both remaining manuscript points and reproduced implementation mismatches involving test-time weights, Gower normalization, and the averaging convention. These are open follow-ups, not resolved changes, and require revalidation of the fit-and-apply claims before submission. No Michel-comment statuses changed.

### 2026-10-08 — Training-only commensurability correction

Corrected the test-dependent commensurability weights identified in the review.
Each numerical or categorical contribution is now divided by its own
training-to-training mean before aggregation; the test batch does not affect
that mean. This applies to plain `mdist()` and `step_mdist()`, including
indicator-based categorical specifications. Robust preprocessing was also
corrected to use training medians and IQRs. Regression checks include adding
unrelated test observations, reordering the batch, and predicting one
observation at a time.

The existing training averaging convention is preserved. The distinct-pair
formula versus full-square averaging question and the mixed-data Gower
normalization discrepancy remain separate open issues. On the paper's penguin
split, training distances are unchanged, Gower is unaffected, and one of 84
predictions changes under `u_dep`; the other checked presets retain their
predictions. No manuscript wording or Michel-comment statuses were changed.

All six package test files pass, and the corrected package was installed
locally. The manuscript PDF was not rebuilt in this implementation task.

### 2026-10-08 — Gower and HLeucl application fixes

Corrected Gower's test-to-training normalization to average over original
predictors, rather than expanded dummy columns. Replaced the Euclidean
cross-distance backend used by `HLeucl` and `hl`, which returned zeros on
larger examples in this installation. The fix covers commensurable and
non-commensurable `HLeucl`; training transformations and training distances
remain unchanged. Removed the unused Rfast dependency.

All seven package test files pass. Larger-data regression tests verify
independent reference calculations, training/new-data agreement, single-row
application, batch independence, and recipe consistency. The broader audit
passes for all 71 specifications it can compute; one categorical measure
rejects sparse profiles with a non-finite-value error. On the paper's penguin
split, Gower predictions are unchanged and three of 84 `hl` predictions change
(accuracy 81/84 to 84/84). The self-pair averaging convention remains to be
discussed. No manuscript wording or Michel-comment statuses changed; the PDF
was not rebuilt.

Installed the corrected package and confirmed Gower and `hl` consistency in a
fresh R session. Independent numerical operator checks, including reduced PCA,
are retained as regression tests alongside the new-data checks.

### 2026-10-16 — Downstream LOVO diagnostic figure

Improved the new two-panel LOVO illustration without changing its diagnostics.
The panels now use horizontal layouts, a shared readable variable axis, and
short titles. The caption distinguishes the two-dimensional classical MDS
alienation coefficient from the three-cluster complete-linkage `1 - ARI`
diagnostic, explains that larger values indicate greater changes, and identifies
the variable-type shading and different horizontal scales. The simple
`autoplot()` calls remain visible; formatting code is omitted from the printed
article. Rebuilt and visually checked the manuscript. No Michel-comment
statuses changed.

Subsequent author-requested refinement: restored the conventional orientation
with variables on the horizontal axes. Shortened variable names and used
two-line labels at a smaller font size without rotation, retaining the
side-by-side panels. The caption now refers to different vertical scales.
Rebuilt and visually checked the 35-page PDF, including page 16 and its
following page. No diagnostic values or Michel-comment statuses changed.

The side-by-side assembly command, `patchwork::wrap_plots(a, b, nrow = 1)`, is
now also printed after the two `autoplot()` calls. Only presentation adjustments
remain hidden. Rebuilt and visually verified page 16; figure appearance and
diagnostic values are unchanged.

Removed redundant interpretation guidance from the across-specification LOVO
subsection and shortened its captions, relying on the definitions already given
for single-specification LOVO. Retained the comparison-specific qualifications:
full-versus-reduced comparisons are within each specification, species labels
are not used, and association-aware constructions may change after omission.
The author's latest edits elsewhere were preserved. This wording-only update
has not been rendered; no Michel-comment statuses changed.

### 2026-10-16 — Simplified candidate-specification interface

Removed the `mode` argument from `all_dist_method_specs()`. The default
constructs custom combinations filtered by the component arguments; preset
names select only those presets, and `preset = "all"` selects all named
presets. Adding `"custom"` to a selection combines these with filtered custom
rows. Commensurable custom combinations retain `std` and `pc_scores` only,
removing the redundant `robust` rows without changing the distances computed
by `mdist()`. There are 413 custom combinations and 10 named presets.

Updated the manuscript's candidate-grid explanation and calls, package help,
release notes, source vignettes, and replication-script calls. The paper now
qualifies affine equivalence as a training-distance property because range
clipping can affect out-of-range test observations. All eight package test
files pass; the revised package was installed locally. Other author edits
were preserved. No Michel-comment statuses changed; the PDF was not rebuilt.

### 2026-10-16 — Direct kNN prediction and Section 6.1

Simplified public `knn_dist()` by removing `dist_fun`. Tabular input now calls
`mdist()` automatically, with distance settings supplied in `dist_args`;
precomputed test-to-training distances remain supported. The internal
tidymodels engine interfaces are unchanged.

Section 6.1 now starts with tabular prediction, then demonstrates externally
computed distances, and finally gives the complete recipe/workflow. The
training-only construction, response-column exclusion, matrix orientation, and
training-outcome ordering are explained in that sequence. Other author edits
were preserved. All eight package test files pass, the package is installed
locally, and the exact manuscript chunks and installed help examples execute
successfully. All three interfaces give identical predictions for the 84 test
penguins (accuracy 1). No Michel-comment statuses changed.

Rebuilt the 35-page manuscript and visually checked pages 21–24. Both printed
prediction comparisons return `TRUE`; no unresolved references or `dist_fun`
mentions remain in the PDF. The tabular and precomputed examples precede the
workflow as requested.

### 2026-10-16 — Workflow boxes moved to Section 6.3

At the author's clarification, Section 6.1 retains its original purpose: the
fixed-five-neighbour penguin workflow demonstrates equivalence with the two
direct interfaces. All 84 test predictions agree; no penguin tuning remains.

The existing WDI pipeline in Section 6.3 now has restrained, phase-labelled
boxes: Preprocess, Specify, Assemble, Tune; then Select, Finalize, Evaluate.
The `last_fit()` idiom is preserved, combining final fitting, prediction, and
test evaluation. The distance and neighbour selection, cross-validation design,
performance summaries, and figures are unchanged. Rebuilt the 35-page PDF and
visually checked pages 23, 28, and 30. The WDI selection remains Unbiased mixed
with nine neighbours, with test balanced accuracy 0.800 and accuracy 0.717.
No package behavior or Michel-comment statuses changed.

### 2026-10-16 — Right-aligned pipeline phase labels

Refined both Section 6.3 workflow boxes by replacing the phase-name code
comments with a dedicated, right-aligned column of muted grey labels.
Each label stays with its code; the code retains the native JSS formatting.
Rebuilt the 35-page PDF and visually verified both boxes on pages 28 and 30.
The analysis and results are unchanged; no Michel-comment statuses changed.

### 2026-10-16 — Upper-triangular benchmark summaries

Replaced the printed pairwise-results tibble with one upper-triangular table
per diagnostic in `summary.MDistBenchmark()`. The `metric` argument selects
measures (all by default; `"ari"` selects available clustering ARIs), and
`digits` controls display rounding. The complete, unrounded tibble remains
the invisible return value, preserving existing extraction and analysis code.
Updated package documentation and the manuscript description; rebuilt the
36-page PDF with all five benchmark tables together on page 20. All eight
test files pass, and the updated package is installed locally. No changes
to computed diagnostics or Michel-comment statuses.

### 2026-10-26 — Carlo R2-004: clustrd function name

Corrected the mixed-data related-software table entry from
`clustrd::cluspca()` to `clustrd::cluspcamix()` in both table data files,
consistent with the existing manuscript prose. Source correction complete;
PDF rebuild pending.

### 2026-10-26 — Explicit block_comm option (Carlo R2-005: manuscript pending)

Added public `block_comm = FALSE` to `mdist()` and `step_mdist()`. For a custom
distance with `method_num = "pc_scores"` and `commensurable = TRUE`, setting
`block_comm = TRUE` enables the whitened-PC Manhattan numerical-block scaling
previously available only through `u_dep_bw`. The preset and explicit settings
produce identical distances; categorical contributions retain their existing
variable-wise treatment. All eight test files pass without warnings, and the
updated package is installed locally. The author will revise the manuscript
wording separately; R2-005 is not yet marked fully addressed. PDF, website,
and source archive remain unchanged by this implementation step.

### 2026-10-26 — Carlo R2-005: manuscript identification added

The presets table now includes `u_dep_bw`, with numerical weighting explicitly
identified as blockwise. The following paragraph explains whitening, Manhattan
distance, scaling of the complete numerical block, and the continued
variable-wise commensurability of the categorical TVD contributions. It also
names the equivalent custom settings, including `block_comm = TRUE`, in the
interface section rather than the theoretical section. R2-005 is addressed.
The manuscript PDF has been rebuilt; the earlier R2-004 function-name
correction is now also reflected in the rendered related-software table.

At the author's subsequent wording preference, both PC-based presets are
labelled "principal-component scaling" in Table 4, with "blockwise" retained
for `u_dep_bw` in the unbiased column. The accompanying paragraph retains the
explicit whitening explanation. This table-label adjustment is saved in the
source and will appear in the next PDF render.

### 2026-10-26 — Commensurable numerical default changed at author's request

The package now uses no preliminary numerical scaling by default when
commensurability is enabled (`method_num = NULL` resolves to `"none"`);
without commensurability, the custom `mdist()` default remains `"std"`.
Explicit numerical methods remain available. The `u_indep` and `u_mix`
presets and the commensurable grid representative now use `"none"`, while
PCA presets and non-commensurable preset definitions are preserved. All eight
test files pass without warnings, and the package is installed locally.
The manuscript's defaults and preset/grid descriptions are pending the
author's revisions; no further comment statuses were changed by this step.

### 2026-10-26 — Review progress recorded in annotated PDF

Reviewed the saved manuscript revisions and marked 18 clearly addressed
second-round annotations green in a separate progress copy. These cover the
section roadmap, mixed-data function name, identification of blockwise
commensurability, Gower/preset exposition, removal of the flagged comparison
chunk, small wording corrections, and use of the default comparative LOVO
plotting metric. Exact annotation IDs and supporting evidence are recorded
in `coauthor-review-progress-2026-10-26.md` and `coauthor-review-status.json`.
Partially addressed comments retain their original colours, as do all other
pending comments. The original annotated PDF is preserved; all comment text,
authors, and reply relationships remain intact in the progress copy.

The author clarified that removal from the visible manuscript is sufficient
for code-deletion requests. R2-045 is therefore also addressed and green:
the lengthy plotting-formatting code is hidden from the paper. R2-044's
broader suggestion about package defaults remains open. The progress copy
now contains 19 green annotations, with all 99 records preserved.

### 2026-10-27 — Abstract comment R2-001

The author confirms that the abstract comment is addressed. Recorded this
confirmation and queued R2-001 to turn green at the next progress-copy update;
the PDF was not regenerated now. R2-074 remains a separate, unchanged item.

### 2026-10-27 — Reply-only comments: light-blue convention

Light blue will identify comments answered without manuscript changes.
The author confirms this status for R2-006, Michel's prepared-data suggestion
at the start of Data preparation (annotated PDF p. 9). Queued its light-blue
colour for the next progress-copy update. No reply text has been supplied or
transcribed here; this records the author's status confirmation only.

### 2026-10-27 — Michel R2-014: returned object and print output

Added `print(d_gower)` before "Other distance specifications" and explained
the printed overview: the custom specification, observation and variable
counts, numerical preprocessing, categorical method, and commensurability.
The prose distinguishes the returned `MDist` object from its
`d_gower$distance` component and names `as.matrix()` for full matrix access.
The exact example runs successfully. R2-014 is addressed in the source and
queued for green; PDF rendering remains pending the author's nearby
unfinished footnote.

### 2026-10-27 — Cleaner package print output

Updated the package's `MDist` print method so that principal-component counts
and inertia thresholds are shown only for `method_num = "pc_scores"`, and
NULL thresholds are never printed. This removes irrelevant PCA entries from
the Gower example's normal `print(d_gower)` output. All nine test files pass;
the updated package is installed locally. Distance calculations and stored
parameters are unchanged.

### 2026-10-27 — Michel R2-015 and updated R2-014 example

Replaced the Gower construction example with the current default, as requested
by the author. The manuscript now shows explicit `method_num = "none"`,
`method_cat = "tvd"`, and `commensurable = TRUE`, followed by the equivalent
short call and `print(d_default)`. The object explanation now refers to
`d_default$distance`. Verified that both calls and the `u_mix` preset produce
identical distances. Updated the associated default/preset descriptions and
tables and rebuilt the 36-page PDF successfully. R2-015 is addressed and
queued for green; R2-014 remains addressed with the revised example.

### 2026-10-27 — Retain the short Gower implementation distinction

Kept the default example first and restored a short, separate Gower
construction call without printed output. The text explains that this call
sums contributions, whereas conventional Gower implementations average them;
for the complete, equally weighted example, division by the predictor count
recovers the conventional scale. Verified agreement with `cluster::daisy()`
and the `gower` preset and rebuilt the 36-page PDF. The longer comparison
chunk and its equality-check output remain omitted.

### 2026-10-27 — Additional summary illustration

Following the author's request, added `summary(d_default)` immediately after
the print explanation. The text distinguishes descriptive distance statistics
from quick preview-block checks, including the sampled triangle-inequality
check, rather than presenting these checks as an exhaustive verification.
The example uses a fixed hidden seed and restores the preceding RNG state;
no package behavior or subsequent stochastic analysis was changed.

Rebuilt the 37-page manuscript and checked the summary layout on page 12.
The complete output stays together, with the following presets subsection
starting cleanly on page 13. The progress annotation PDF remains unchanged.

### 2026-10-27 — Latest colour/status update

The progress PDF now records 26 green annotations and one light-blue
reply-only annotation; the remaining 72 retain their original colours.
Applied the previously queued abstract, default-object, and default-example
updates and recorded the additional completed wording/code changes. R2-022
and R2-025 remain partial for the reasons in the detailed progress log.
All existing comment text, including the author's manually appended answer,
is preserved. A dated backup of the preceding progress copy is retained.

The author's R2-006 answer is now available in the PDF and is transcribed
verbatim here (no manuscript change for this comment):

> I would keep this, I would not define a dataset that is already availble just to drop one variable and some missings. Also, I believe it serves as an example of minimal data preparation

Separately, corrected `penguin_train[, 1]` to `penguin_train[, -1]` in the
`u_mix` example at the author's explicit instruction. No other manuscript
changes or article PDF rendering were performed in this update.

### 2026-10-27 — Consistent non-commensurable training/new-data example

Both calls now use robust numerical preprocessing, TVD categorical
dissimilarities, and `commensurable = FALSE`. Updated the explanation to
describe training-fitted category dissimilarities and numerical preprocessing,
and explicitly noted that no commensurability weights are estimated in this
example. The former mismatch is resolved; the remaining R2-025 requests are
tracked separately. Source-only wording update; PDFs were not regenerated.
