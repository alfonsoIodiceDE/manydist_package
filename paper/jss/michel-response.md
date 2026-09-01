# Response to Michel’s comments

This is the concise, colleague-facing response log. The complete annotations are preserved in [michel-comments-0708.md](michel-comments-0708.md), and working status is maintained in [michel-comments-tracker.md](michel-comments-tracker.md).

Entries are added here after a comment has been assessed. Unresolved author decisions are explicitly marked `needs discussion`; completed changes are verified against the manuscript. Stable IDs `M-001` through `M-140` match the tracker.

## Progress

- Assessed and documented: 27 of 140
- Resolved or superseded: 18
- Partially resolved: 3
- Needs discussion: 6
- Not yet assessed: 113

## Responses

### Abstract — M-001 to M-005

| ID | Response | Status |
|---|---|---|
| M-001 | Replaced the technical property list with the more direct statement that `manydist` constructs mixed-variable distances that account for the identified biases. | resolved |
| M-002 | Recast the contribution list and removed the marked conjunction. | resolved |
| M-003 | Moved pipeline integration into a separate sentence. The sentence still needs copyediting before this point is closed. | partially resolved |
| M-004 | The Palmer penguins example is now explicitly connected to distance construction, diagnostics, and benchmarking. | resolved |
| M-005 | The World Development Indicators example now explicitly names unsupervised and supervised pipelines, clustering and classification, and resample-specific preprocessing. The sentence still needs grammatical correction. | partially resolved |

**Location:** Abstract.

### Introduction: importance bias — M-007 to M-018

| ID | Response | Status |
|---|---|---|
| M-007 | Deleted the “mirror image” characterization. | resolved |
| M-008 | Removed the marked comparison to a matching categorical variable. | resolved |
| M-009 | The phrase about favouring categorical contributions remains, but the annotation has no explanatory text. We should confirm whether the suggested dominance wording is preferred. | needs discussion |
| M-010 | The containing sentence was replaced; the requested preposition change is no longer applicable. | superseded |
| M-011 | Recast the Gower explanation around a mean numerical distance below one. We should discuss whether the stronger wording “much smaller than one” and “categorical variables tend to dominate” is needed. | partially resolved |
| M-012 | Reworded the sentence to say that skewed numerical distributions amplify the effect. A spelling typo remains to be corrected. | resolved |
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

<!--
Use this structure for each addressed comment:

### M-000 — PDF page.item

**Comment:** Concise paraphrase of Michel’s concern.

**Response:** What was changed, or why no change was made.

**Location:** Current section/paragraph/line anchor in article.qmd.

**Status:** needs discussion | resolved | partially resolved | already addressed | superseded | declined
-->
