# Michel comments tracker

This is the working log for Michel’s annotations. The extracted comments remain unchanged in [michel-comments-0708.md](michel-comments-0708.md).

## Baseline

- Started: 2026-08-31
- Manuscript: [article.qmd](article.qmd)
- Git commit: `6bc38990bf1b2647be31ce46cbb46d82839086f9`
- Baseline `article.qmd` hash: `c6e2599f6bd41ed0749b2f96c29d58aaea4c85e2`
- Source: Michel’s annotated `jss paper 0708.pdf`, represented by the immutable extraction above
- Annotations: 140 across 20 PDF pages

The author edits `article.qmd`. After each reported batch, the tracker and [michel-response.md](michel-response.md) are reconciled against the manuscript diff and the rendered paper is checked.

## Status conventions

- `open`: not yet assessed or changed
- `needs discussion`: a substantive author decision is required
- `in progress`: an edit has begun but is not yet complete or verified
- `resolved`: the comment is fully addressed and verified
- `partially resolved`: some, but not all, of the concern is addressed
- `already addressed`: the current manuscript already addresses the old-PDF comment
- `superseded`: later restructuring or rewriting made the original comment inapplicable
- `declined`: no change was made; the reason is recorded

## Progress

| Status | Count |
|---|---:|
| Open | 55 |
| Needs discussion | 38 |
| In progress | 1 |
| Partially resolved | 4 |
| Resolved | 37 |
| Superseded | 5 |
| All other statuses | 0 |
| **Total** | **140** |

## Comments

| ID | PDF item | Type | Marked-text anchor | Status | Current manuscript location | Decision / change |
|---|---:|---|---|---|---|---|
| M-001 | 1.1 | Highlight | through multivariate addi tivity and commensurability. | resolved | Abstract | Replaced the technical property list with a direct statement that the distances account for the biases. |
| M-002 | 1.2 | StrikeOut | and | resolved | Abstract | Recast the contribution list and removed the marked conjunction. |
| M-003 | 1.3 | Highlight | and the | resolved | Abstract | Moved pipeline integration into a separate, copyedited sentence describing unsupervised and supervised pipelines. |
| M-004 | 1.4 | Highlight | The Palmer penguins data | resolved | Abstract | Now explicitly links the Palmer penguins example to construction, diagnostics, and benchmarking. |
| M-005 | 1.5 | Highlight | A World Development Indicators snapshot then places the distance inside clustering and… | resolved | Abstract | Now identifies the World Development Indicators snapshot, unsupervised and supervised pipelines, clustering and classification, and resample-specific refitting in a grammatically complete sentence. |
| M-006 | 1.6 | Text | [general note] | resolved | Overall organization | Reorganized the article around the introduction, a dedicated software review, the unified framework, construction, diagnostics, pipelines, and conclusion. |
| M-007 | 2.1 | Highlight | Numerical variables face the mirror image of this problem through scaling. | resolved | Introduction — importance bias | Deleted the “mirror image” sentence. |
| M-008 | 2.2 | Highlight | of a matching categorical | resolved | Introduction — importance bias | Removed the marked comparison in the rewritten Gower discussion. |
| M-009 | 2.3 | Highlight | to favor categorical contributions | needs discussion | Introduction — importance bias | The phrase remains, and the annotation contains no written explanation; confirm whether Michel wants the alternative “categorical variables tend to dominate.” |
| M-010 | 2.4 | Highlight | to correct one variable at a time. | superseded | Introduction — importance bias | The containing sentence was replaced, so the requested preposition change is no longer applicable. |
| M-011 | 2.5 | Highlight | the resulting mean pairwise distance is typically small relative to that of a matching… | partially resolved | Introduction — importance bias | Rewritten around a mean numerical distance below one; discuss whether “much smaller than one” and “categorical variables tend to dominate” are needed for the intended argument. |
| M-012 | 2.6 | Highlight | Skewness compounds the effect: | resolved | Introduction — importance bias | Adopted “Skewed distributions for the numerical variables amplify this effect” and corrected the spelling of “distributions.” |
| M-013 | 2.7 | Highlight | heavy-tailed | resolved | Introduction — importance bias | Replaced “heavy-tailed” with “long-tailed.” |
| M-014 | 2.8 | StrikeOut | while other scaling choices can amplify or attenuate the same variable unpredictably. | resolved | Introduction — importance bias | Deleted the marked clause. |
| M-015 | 2.9 | Highlight | distortions | superseded | Introduction — importance bias | Replaced the containing passage with a new statement about undesirable factors in the overall distance. |
| M-016 | 2.10 | Highlight | Because these distortions follow from measurement type and scale rather than from… | superseded | Introduction — importance bias | Replaced the containing sentence with the reformulated problem statement and citation. |
| M-017 | 2.11 | Highlight | are hard to anticipate and to correct one variable at a time. | resolved | Introduction — importance bias | Replaced the passage with Michel’s proposed framing of the undesirable factors and the general construction. |
| M-018 | 2.12 | Highlight | van de Velden et al. (2026) formalize unbiased mixed-variable distances through two… | resolved | Introduction — additivity and commensurability | Reframed these as two essential properties of the proposed unbiased distances. |
| M-019 | 2.13 | Highlight | requires | needs discussion | Introduction — additivity and commensurability | “Requires” remains; decide whether “entails” or “specifies” is preferable. |
| M-020 | 2.14 | StrikeOut | so that each variable enters the total additively. | needs discussion | Introduction — additivity and commensurability | The marked explanatory clause remains. |
| M-021 | 2.15 | Highlight | Importantly, a variable-specific contribution need not depend on that variable alone;… | resolved | Introduction — additivity and commensurability | Deleted the potentially confusing anticipatory explanation. |
| M-022 | 2.16 | Highlight | the | needs discussion | Introduction — additivity and commensurability | Current wording still begins “the variable-specific distances”; clarify the intended “individual variable” revision. |
| M-023 | 2.17 | StrikeOut | drawn from a more dispersed dis tribution, | resolved | Introduction — additivity and commensurability | Deleted the marked phrase. |
| M-024 | 2.18 | Highlight | are commensurable | resolved | Introduction — additivity and commensurability | Changed to “are defined to be commensurable.” |
| M-025 | 2.19 | Highlight | so that | needs discussion | Introduction — additivity and commensurability | The sentence still uses “so that”; decide whether to split it and begin “Hence, no variable…”. |
| M-026 | 2.20 | Highlight | For | resolved | Introduction — additivity and commensurability | Changed “For the general additive formulation” to “In the general additive formulation.” |
| M-027 | 2.21 | Highlight | this is | needs discussion | Introduction — additivity and commensurability | “This is achieved” remains; decide whether to make commensurability the grammatical subject. |
| M-028 | 2.22 | StrikeOut | - while still allowing a variable to be up- or down-weighted deliberately when… | resolved | Introduction — additivity and commensurability | Deleted the marked qualification. |
| M-029 | 2.23 | Highlight | this framework. | resolved | Introduction — package contribution | Now states that `manydist` implements the unbiased framework as well as distances outside it. |
| M-030 | 2.24 | Highlight | and the discipline required to use an estimated distance inside a resampled workflow. | resolved | Introduction and Other software | Replaced “discipline” with a copyedited learning-pipeline contribution and an explicit training-only fit-and-apply explanation. |
| M-031 | 2.25 | Highlight | a | partially resolved | Introduction — `mdist()` capabilities | The list covers construction, commensurability, aggregation, and reuse; where to introduce association- and response-aware construction remains unsettled. |
| M-032 | 3.1 | Highlight | Association-aware distances | partially resolved | Introduction and framework | A dedicated association-aware subsection is present, but the first introduction and placement of the overview paragraph remain under discussion. |
| M-033 | 3.2 | Text | [general note] | needs discussion | Introduction — blue review paragraph | The awareness paragraph is explicitly marked for removal or relocation. |
| M-034 | 3.3 | Highlight | Scale- and type-aware | needs discussion | Introduction — blue review paragraph | The terminology remains inside the paragraph marked for removal or relocation. |
| M-035 | 3.4 | Highlight | pairwise benchmarking of candidate distances, | resolved | Introduction and roadmap | Replaced the opaque phrase with “pairwise comparisons of candidate distances” and “comparisons between candidate distances.” |
| M-036 | 3.5 | Highlight | Several R packages provide distances for homogeneous numerical, categorical, or binary… | resolved | Other software | Created a dedicated `Other software` section and moved the review there. |
| M-037 | 3.6 | Highlight | clustering, ordination, and nearest-neighbour met hods. | needs discussion | Other software | The specific method list remains; decide whether to replace it with “distance-based data analysis methods.” |
| M-038 | 3.7 | StrikeOut | manydist builds on rather than replaces these tools, computing comparable per-block… | resolved | Other software | Removed the sentence from rendered output by commenting it out. |
| M-039 | 4.1 | Highlight | matrices, | resolved | Other software — homogeneous packages | Changed “numerical matrices” to “numerical data” and corrected the punctuation. |
| M-040 | 4.2 | Highlight | extensible | needs discussion | Other software — homogeneous packages | “Extensible catalogue” remains and may need a plainer explanation. |
| M-041 | 4.3 | Highlight | provides frequency-weighted similarities. | needs discussion | Other software — homogeneous packages | The description remains; clarify `nomclust` and the Boriah-family coverage. |
| M-042 | 4.4 | Highlight | cluster::daisy(), | needs discussion | Other software — mixed packages | The text still presents `daisy()` only with the Gower implementations; its broader numerical capabilities are not yet noted. |
| M-043 | 4.5 | Highlight | of the same underlying construction | needs discussion | Other software — mixed packages | The phrase remains without explaining precisely how the implementations differ. |
| M-044 | 4.6 | Highlight | and its tendency to give relatively greater influence to categorical variables when… | needs discussion | Other software — mixed packages | The Gower property remains embedded in the package comparison. |
| M-045 | 4.7 | Highlight | for mixed data | needs discussion | Other software — mixed packages | Consider Michel’s clearer formulation “packages implementing mixed-variable distances.” |
| M-046 | 4.8 | Highlight | controls the relative influence of numerical and categorical information through… | needs discussion | Other software — mixed packages | The family description remains unchanged. |
| M-047 | 4.9 | Highlight | Gower, Podani, Wishart, Huang, Harikumar, and Ahmad | needs discussion | Other software — mixed packages | References for the named formulations have not yet been added. |
| M-048 | 4.10 | Highlight | learns an adaptive | needs discussion | Other software — mixed packages | The “adaptive distance” wording remains; Michel’s data-estimation wording has not been adopted. |
| M-049 | 4.11 | StrikeOut | automatically, | needs discussion | Other software — mixed packages | The marked word remains in the `kproto()` description. |
| M-050 | 4.12 | Highlight | Related to | needs discussion | Other software — mixed packages | “Related to this family” remains. |
| M-051 | 4.13 | Highlight | a weighting | needs discussion | Other software — mixed packages | Singular “a weighting” remains. |
| M-052 | 4.14 | Highlight | Both address the balance between variable types as part of a clustering procedure, | needs discussion | Other software — mixed packages | The text does not yet explain that the procedures do not yield a separate distance matrix. |
| M-053 | 4.15 | Highlight | it differs from the commensurable distances developed here, which explicitly equalize… | superseded | Other software — mixed packages | Removed the containing comparison rather than revising it. |
| M-054 | 4.16 | StrikeOut | likewise | resolved | Other software — mixed packages | Removed the marked word with the containing sentence. |
| M-055 | 4.17 | Highlight | the estimated quantities must be refit within a resampling scheme, exactly as recipes… | resolved | Other software — resampling paragraph | Replaced the generic refitting language with an explicit fit-and-apply description: `manydist` estimates data-dependent quantities from the training data only, holds them fixed for test-to-training distances, and explains the leakage that would result from fitting on the complete data set. |
| M-056 | 4.18 | Highlight | To our knowledge no other mixed-data distance implementation exposes a fitted… | resolved | Other software — resampling paragraph | Deleted the broad novelty claim and replaced it with an explicit description of the resampling behavior. |
| M-057 | 4.19 | Highlight | Gower dissimilarity for numerical, nominal, ordinal and binary variables. | needs discussion | Other software table | The table description of `daisy()` has not yet been broadened. |
| M-058 | 4.20 | Highlight | variants | needs discussion | Other software table | The meaning of “variants” remains unexplained. |
| M-059 | 4.21 | Highlight | share a consequence | superseded | Other software | Removed the containing claim. |
| M-060 | 4.22 | Text | [general note] | needs discussion | Other software table | The three families described in the prose are not yet made visible in the table. |
| M-061 | 5.1 | Highlight | The remainder | resolved | Introduction — roadmap | Moved the roadmap before the dedicated `Other software` section. |
| M-062 | 5.2 | StrikeOut | These tools are useful, but their variable contributions can depend on scale,… | resolved | Other software — concluding paragraph | Deleted the marked sentence. |
| M-063 | 5.3 | Highlight | them | resolved | Other software — concluding paragraph | Recast the transition as “complements this ecosystem.” |
| M-064 | 5.4 | Highlight | association-aware and response-aware | partially resolved | Introduction, Other software, and framework | The terms are retained and explained elsewhere, but the introductory overview containing their first explanation is still marked for possible removal. |
| M-065 | 5.5 | Highlight | develops | needs discussion | Introduction — roadmap | The roadmap still uses “develops” and section labels rather than Michel’s proposed direct section-by-section wording. |
| M-066 | 5.6 | Highlight | fitted preprocessing. | needs discussion | Introduction — roadmap | “Fitted preprocessing” remains. |
| M-067 | 5.7 | Highlight | pairwise benchmarking | resolved | Introduction — roadmap | Changed the roadmap wording to the broader “benchmarking of candidate distances.” |
| M-068 | 5.8 | Highlight | embeds distances | resolved | Introduction — roadmap | Recast the roadmap sentence as “shows how distances can be embedded in `tidymodels` workflows.” |
| M-069 | 5.9 | Highlight | I | needs discussion | Framework — notation | The `I`, `Q_n`, and `Q_c` notation remains unchanged. |
| M-070 | 6.1 | Highlight | the contribution | needs discussion | Framework — general setup | The text still uses “contribution” rather than “distance contribution.” |
| M-071 | 6.2 | Highlight | Each component may depend on all elements of the two observations; the index… | resolved | Framework — additivity | Retained the metric statement but added a citation to the proof in van de Velden et al. |
| M-072 | 6.3 | Highlight | Scale- and type-aware distances | resolved | Framework — additivity and commensurability | Replaced the heading with “Multivariate additivity and commensurability.” |
| M-073 | 7.1 | Highlight | The formulation also permits constrained weights, such as common weights within… | in progress | Framework — commensurability | Marked the paragraph `TO BE REMOVED`, but it remains in the source pending deletion. |
| M-074 | 7.2 | StrikeOut | beyond commensurability | resolved | Framework — association-aware distances | Removed “beyond commensurability” from the heading. |
| M-075 | 7.3 | Highlight | Association-aware distances: | needs discussion | Framework — association-aware distances | The subsection is explicitly marked for revision. The implementation rotates normalized predictors to ordinary PCA scores and subsequently commensurates each component; it does not apply the explicit $\bm{\Lambda}^{-1/2}$ whitening transformation currently shown in the manuscript. |
| M-076 | 7.4 | StrikeOut | by itself | needs discussion | Framework — association-aware distances | The marked words remain. |
| M-077 | 7.5 | Highlight | When several variables encode the same underlying source of variation, that variation… | needs discussion | Framework — association-aware distances | The proposed reorganization and revised opening have not yet been implemented. |
| M-078 | 9.1 | Highlight | whitening, | resolved | Framework — numerical variables | Added a footnote defining whitening and its Mahalanobis-distance equivalence. |
| M-079 | 9.2 | Highlight | Figure 1: | needs discussion | Framework — association-aware figure | Confirmed that the figure still has no textual cross-reference. Add a sentence such as “@fig-independence-association illustrates the change in geometry…” at the point where the example is introduced. |
| M-080 | 9.3 | Highlight | whitening | resolved | Framework — numerical variables | Added a mathematical footnote defining whitening and its equivalence to Mahalanobis distance. |
| M-081 | 9.4 | Highlight | to remove linear redundancy. | needs discussion | Framework — numerical variables | The definition is now supplied, but the text must distinguish exact whitening/Mahalanobis equivalence from the implemented PCA-score rotation followed by empirical component-wise commensurability; this remains linked to `M-075`. |
| M-082 | 10.1 | Highlight | Common transformations include: • standard-deviation scaling: (x - x)/sx, which uses… | needs discussion | Framework — numerical variables | The transformation list still follows the general association-aware discussion. Decide whether to introduce numerical preprocessing before association-aware constructions, as Michel suggests. |
| M-083 | 11.1 | Highlight | Indicator-based dissimilarities treat a binary indicator representation as numerical.… | needs discussion | Framework — categorical variables | The paragraph remains. Decide whether to delete it or reduce it to a single sentence needed to explain the available categorical presets. |
| M-084 | 11.2 | Highlight | Association-based dissimilarities incorporate relationships among categorical… | needs discussion | Framework — categorical variables | The passage remains close to the methodological paper and should be shortened around what the package implements. A technical decision is also required for exact independence: all dissimilarities and their means become zero, so the current commensurability code returns `NaN`. |
| M-085 | 11.3 | Highlight | The same construction accommodates a response. When an outcome is supplied, the… | needs discussion | Framework — categorical variables | The response-aware paragraph remains substantially unchanged and is still bold as revision markup. Reorganize it around the `response` argument, what is fitted, the single-predictor behavior, and the absence of a response-aware numerical counterpart. |
| M-086 | 12.1 | StrikeOut | [general note] | open | — | — |
| M-087 | 12.2 | Highlight | A preset identifies a complete distance specification. | open | — | — |
| M-088 | 12.3 | Highlight | u_dep | open | — | — |
| M-089 | 12.4 | Highlight | carries | open | — | — |
| M-090 | 12.5 | Highlight | consumed by downstream | open | — | — |
| M-091 | 12.6 | Highlight | hclust () and cmdscale (). | open | — | — |
| M-092 | 13.1 | Highlight | in the gower, euclidean, and hl presets. | open | — | — |
| M-093 | 14.1 | Highlight | The custom preset exposes the three choices directly through the method_num,… | open | — | — |
| M-094 | 14.2 | Highlight | Table 3 summarizes the presets | open | — | — |
| M-095 | 15.1 | Highlight | A distance is an intermediate data representation, so its sensitivity should be… | open | — | — |
| M-096 | 15.2 | Highlight | one predictor | open | — | — |
| M-097 | 15.3 | Highlight | The diagnostics | open | — | — |
| M-098 | 15.4 | Highlight | relative_distance | open | — | — |
| M-099 | 15.5 | Highlight | the alienation between the full and reduced configurations. | open | — | — |
| M-100 | 15.6 | Text | in the [O, | open | — | — |
| M-101 | 15.7 | Highlight | measuring the configuration variance left unexplained | open | — | — |
| M-102 | 15.8 | Highlight | Reporting both is deliberate: a variable can shift many pairwise distances while… | open | — | — |
| M-103 | 15.9 | Highlight | benchmarking | open | — | — |
| M-104 | 16.1 | Highlight | Figure 2: Relative | open | — | — |
| M-105 | 16.2 | Highlight | lovo_gower$autoplot( = metric "relative_distance", + reorder= + TRUE +) + = ggplot2::… | open | — | — |
| M-106 | 16.3 | Highlight | The two metrics disagree | open | — | — |
| M-107 | 17.1 | Highlight | Comparing distance specifications | open | — | — |
| M-108 | 18.1 | Highlight | Figure 3: Stability | open | — | — |
| M-109 | 18.2 | Highlight | and can disagree, | open | — | — |
| M-110 | 18.3 | Highlight | The pairwise diagnostics | open | — | — |
| M-111 | 18.4 | Highlight | changes in magnitude from changes in geometry. | open | — | — |
| M-112 | 18.5 | Highlight | every pair of successful distances. | open | — | — |
| M-113 | 18.6 | Highlight | The MDS congruence coefficient and its corresponding alienation coefficient compare… | open | — | — |
| M-114 | 18.7 | Highlight | If is specified, 0 cluster_ k benchmark_mdist also applies each requested clustering… | open | — | — |
| M-115 | 18.8 | Highlight | Benchmarking | open | — | — |
| M-116 | 18.9 | Highlight | evaluates an explicit table of distance spec ifications. | open | — | — |
| M-117 | 19.1 | Highlight | benchmark_pairs <- benchmark_comparisons(distance_benchmark) | open | — | — |
| M-118 | 19.2 | Highlight | candidate_specs <- all_dist_method_specs( mode = "presets_only", + preset = c("gower",… | open | — | — |
| M-119 | 19.3 | Highlight | under pairs of distances; | open | — | — |
| M-120 | 19.4 | Highlight | low-dimensional geometry, | open | — | — |
| M-121 | 19.5 | Highlight | determines | open | — | — |
| M-122 | 19.6 | Highlight | whether those changes alter the partition. | open | — | — |
| M-123 | 19.7 | Highlight | are therefore complementary rather than alternative rankings. | open | — | — |
| M-124 | 19.8 | Highlight | renders any pairwise diagnostic autoplot () as an annotated triangular heatmap. | open | — | — |
| M-125 | 19.9 | Highlight | An error in one specification is captured | open | — | — |
| M-126 | 19.10 | Highlight | a benchmark remains a sensitivity analysis rather than a model-selection rule. | open | — | — |
| M-127 | 19.11 | Highlight | Gower | open | — | — |
| M-128 | 19.12 | Highlight | has more than twice the alienation | open | — | — |
| M-129 | 19.13 | Highlight | a much lower PAM ARI. | open | — | — |
| M-130 | 19.14 | Highlight | This | open | — | — |
| M-131 | 20.1 | Highlight | distance, and so the same gap from the Gower | open | — | — |
| M-132 | 20.2 | Highlight | Taken together, the diagnostics indicate how much of an analysis rests on the choice… | open | — | — |
| M-133 | 20.3 | Highlight | Distance-based learning pipelines | open | — | — |
| M-134 | 20.4 | Highlight | Refitting the distance within resamples | open | — | — |
| M-135 | 20.5 | Highlight | recipe. | open | — | — |
| M-136 | 20.6 | Highlight | The fitted step stores the training data and preprocessing parameters, so the distance… | open | — | — |
| M-137 | 20.7 | Highlight | Its | open | — | — |
| M-138 | 20.8 | Highlight | required by the downstream task: for clustering and "pai rwise" for prediction.… | open | — | — |
| M-139 | 20.9 | Highlight | whereas nearest-neighbour prediction | open | — | — |
| M-140 | 21.1 | Highlight | a fixed snapshot | open | — | — |

## Verification notes

### 2026-08-31 — Abstract and partial introduction

- Reconciled `M-001`–`M-005` and `M-007`–`M-028`; `M-006` remains open because it concerns the paper’s overall organization.
- Follow-up: completed the abstract and introduction copyedits, including “which constructs,” “distributions,” “integrate … into,” the clustering/classification sentence, doubled spaces, and “mixed-variable distance.”
- `M-009`, `M-019`, `M-020`, `M-022`, `M-025`, and `M-027` are explicitly parked for discussion.
- Rendered `article.qmd` successfully to a 33-page PDF. Visual inspection of pages 1–3 found no clipping, overlap, broken references, or other layout defects.

### 2026-09-01 — Introduction, software review, and framework revisions

- Reconciled `M-006`, `M-029`–`M-078`, `M-080`, and `M-081`; `M-079` remains open because the figure-reference issue was not changed.
- The manuscript now has a dedicated `Other software` section, a revised package-summary transition, explicit training-to-test resampling language, consolidated additivity/commensurability headings, and a whitening/Mahalanobis footnote.
- The resampling paragraph after `kdml::dkss()` should explicitly name `manydist` or be moved to the package-summary paragraph; as written, it appears to describe DKSS.
- Follow-up: corrected “seamless,” the package-implementation sentence, “distance computations,” “the package’s main function,” the punctuation after “For numerical data,” and the comma splice in the blue awareness paragraph.
- `git diff --check` still reports trailing whitespace in the paragraph marked `TO BE REMOVED` and after the whitening footnote; these lines are part of unresolved review material.
- Review markers remain intentionally visible for the awareness paragraph, constrained weights, and the association-aware subsection.
- Follow-up: `M-055` is closed after the revised paragraph explicitly names `manydist`, separates fitting from application, and states the leakage consequence. The shorter resampling statement after the mixed-software table now repeats the same point and should be consolidated with the detailed paragraph.
- Technical follow-up for `M-075`/`M-081`: `manydist/R/ndist.R` uses `step_normalize()` followed by `step_pca()` and then divides each component-wise distance by its empirical mean. The manuscript currently writes the distinct whitening transformation $\mathbf{X}\mathbf{V}\bm{\Lambda}^{-1/2}$. Preserve the Mahalanobis equivalence as a definition of whitening, but do not imply that the implemented commensurable Manhattan construction is itself Mahalanobis distance.
- `M-079` is now assessed: the figure label exists, but no `@fig-independence-association` reference appears in the prose.
- `M-082`–`M-085` are now assessed and parked for discussion. A balanced two-factor test confirmed that exact categorical independence produces zero total-variation category dissimilarities and `NaN` commensurable distances because the component means are zero; the implementation and corresponding manuscript statement need an explicit zero-contribution or fallback rule.
