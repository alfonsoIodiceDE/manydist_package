# Response to Michel’s comments

This is the concise, colleague-facing response log. The complete annotations are preserved in [michel-comments-0708.md](michel-comments-0708.md), and working status is maintained in [michel-comments-tracker.md](michel-comments-tracker.md).

Entries are added here after a comment has been assessed. Unresolved author decisions are explicitly marked `needs discussion`; completed changes are verified against the manuscript. Stable IDs `M-001` through `M-140` match the tracker.

## Progress

- Assessed and documented: 85 of 140
- Resolved or superseded: 42
- Partially resolved: 4
- In progress: 1
- Needs discussion: 38
- Not yet assessed: 55

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
| M-073 | Marked the constrained-weights paragraph `TO BE REMOVED`; the deletion has not yet been completed. | in progress |
| M-074 | Removed “beyond commensurability” from the association-aware heading. | resolved |
| M-075 | Marked the association-aware subsection for revision. The implementation uses ordinary PCA scores after normalization and then commensurates the component-wise distances; this must be distinguished from explicit PCA whitening by $\bm{\Lambda}^{-1/2}$. | needs discussion |
| M-076 | “By itself” remains. | needs discussion |
| M-077 | The proposed reorganization and revised opening of the association-aware discussion have not yet been implemented. | needs discussion |

**Location:** Unified framework, additivity/commensurability and association-aware subsections.

### Whitening and Figure 1 — M-078 to M-081

| ID | Response | Status |
|---|---|---|
| M-078 | Added a footnote defining whitening through its equivalence with Mahalanobis distance. | resolved |
| M-079 | Confirmed that Figure 1 still has no textual cross-reference. We should add an explicit reference where the Palmer penguins illustration is introduced. | needs discussion |
| M-080 | Defined whitening mathematically and stated the full-rank and retained-subspace qualifications. | resolved |
| M-081 | The definition is now present, but the text must distinguish exact whitening/Mahalanobis equivalence from the implemented PCA-score rotation followed by empirical component-wise commensurability. This remains linked to the planned association-aware rewrite. | needs discussion |

**Location:** Framework, association-aware figure and `Numerical variables` subsection.

### Numerical and categorical framework — M-082 to M-085

| ID | Response | Status |
|---|---|---|
| M-082 | The numerical-transformation list still follows the general association-aware discussion. We should decide whether to introduce numerical preprocessing first, as suggested. | needs discussion |
| M-083 | The indicator-based paragraph remains. We should either remove it or retain only the sentence needed to explain package options. | needs discussion |
| M-084 | The categorical association passage should be shortened and recast around the implementation rather than closely repeating the methodological paper. In addition, exact categorical independence currently yields zero component means and therefore `NaN` after commensurability; the package and text need an explicit rule for this case. | needs discussion |
| M-085 | The response-aware passage still needs reorganization and its bold revision markup must be removed before submission. | needs discussion |

**Location:** Framework, `Numerical variables` and `Categorical variables` subsections.

<!--
Use this structure for each addressed comment:

### M-000 — PDF page.item

**Comment:** Concise paraphrase of Michel’s concern.

**Response:** What was changed, or why no change was made.

**Location:** Current section/paragraph/line anchor in article.qmd.

**Status:** needs discussion | resolved | partially resolved | already addressed | superseded | declined
-->
