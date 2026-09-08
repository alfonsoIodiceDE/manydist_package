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
| M-117 | Simplified the example to print `benchmark_comparisons(distance_benchmark)` directly. The annotation has no written explanation, so we should confirm whether further change was intended. | needs discussion |
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
| M-134 | Added a general workflow subsection before the more specific resampling demonstration and renamed the latter around its motivation. | resolved |
| M-135 | Defines a `tidymodels` recipe before introducing `step_mdist()`. | resolved |
| M-136 | Condensed the repeated fit-and-apply explanation and removed the unclear pipeline-level analogy. | resolved |
| M-137 | Names `step_mdist()` explicitly when introducing its `output` argument. | resolved |
| M-138 | Separately explains pairwise clustering distances, new-to-training prediction distances, and MDS as a non-workflow downstream use. | resolved |
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
