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
| Open | 0 |
| Needs discussion | 42 |
| In progress | 1 |
| Partially resolved | 4 |
| Resolved | 80 |
| Superseded | 13 |
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
| M-073 | 7.1 | Highlight | The formulation also permits constrained weights, such as common weights within… | resolved | Framework — commensurability | Deleted the constrained-weights paragraph. |
| M-074 | 7.2 | StrikeOut | beyond commensurability | resolved | Framework — association-aware distances | Removed “beyond commensurability” from the heading. |
| M-075 | 7.3 | Highlight | Association-aware distances: | in progress | Framework — association-aware numerical distances | The current purple trial whitens the retained PCA scores and then applies one factor to scale the complete numerical Manhattan distance to mean $Q_n$. This is distinct from the package's current component-wise empirical scaling, and the final choice remains to be agreed. |
| M-076 | 7.4 | StrikeOut | by itself | resolved | Framework — association-aware distances | Removed “by itself.” |
| M-077 | 7.5 | Highlight | When several variables encode the same underlying source of variation, that variation… | resolved | Framework — association-aware distances | Replaced the earlier opening with a concise distinction between commensurability and within-block association awareness. |
| M-078 | 9.1 | Highlight | whitening, | resolved | Framework — numerical variables | Defines whitening as PCA rotation followed by division by $\sqrt{\lambda_h}$, so that the retained coordinates are uncorrelated with unit variance. The current trial uses this transformation before block scaling. |
| M-079 | 9.2 | Highlight | Figure 1: | resolved | Framework — association-aware figure | Added an explicit textual reference to Figure 1 and revised its third panel to show the block-scaled whitened coordinates. The figure retains its concise three-panel layout. |
| M-080 | 9.3 | Highlight | whitening | resolved | Framework — numerical variables | Distinguishes decorrelation by PCA rotation from whitening, which additionally standardizes the PC variances, and explains that the subsequent common block factor does not undo this relative rescaling. |
| M-081 | 9.4 | Highlight | to remove linear redundancy. | resolved | Framework — numerical variables | The revision distinguishes PCA rotation, whitening, block-level Manhattan calibration, component-wise empirical scaling, and Mahalanobis distance; Euclidean distance on all whitened PCs is Mahalanobis, while Manhattan distance is not. |
| M-082 | 10.1 | Highlight | Common transformations include: • standard-deviation scaling: (x - x)/sx, which uses… | resolved | Framework — general setup | Moved the numerical construction and transformation list before multivariate additivity, commensurability, and the association-aware discussion. |
| M-083 | 11.1 | Highlight | Indicator-based dissimilarities treat a binary indicator representation as numerical.… | needs discussion | Framework — categorical variables | The paragraph remains. Decide whether to delete it or reduce it to a single sentence needed to explain the available categorical presets. |
| M-084 | 11.2 | Highlight | Association-based dissimilarities incorporate relationships among categorical… | partially resolved | Framework — categorical variables | The passage is shorter, corrects the conditional-distribution notation to $\mathbf{R}^{k'\mid k}$, and states the zero-mean problem explicitly: identically zero dissimilarities cannot be reciprocal-mean scaled. The package still needs an agreed zero-contribution or fallback rule. |
| M-085 | 11.3 | Highlight | The same construction accommodates a response. When an outcome is supplied, the… | needs discussion | Framework — categorical variables | The response-aware paragraph remains substantially unchanged and is still bold as revision markup. Reorganize it around the `response` argument, what is fitted, the single-predictor behavior, and the absence of a response-aware numerical counterpart. |
| M-086 | 12.1 | StrikeOut | [general note] | superseded | Distance construction — opening | The punctuation marked in the earlier passage disappeared when the section opening was rewritten. |
| M-087 | 12.2 | Highlight | A preset identifies a complete distance specification. | resolved | Distance construction — basic interface | Introduced `x`, `method_num`, `method_cat`, and `commensurable`, including defaults and main options, before discussing presets; the `MDist` object follows the preset table. |
| M-088 | 12.3 | Highlight | u_dep | needs discussion | Distance construction — preset example | `u_dep` remains the principal worked preset, but the annotation contains no written explanation. Confirm whether Michel intended a different example or only marked the term. |
| M-089 | 12.4 | Highlight | carries | resolved | Distance construction — `MDist` object | Replaced “carries” with “contains” in the rewritten object description. |
| M-090 | 12.5 | Highlight | consumed by downstream | resolved | Distance construction — `MDist` object | Removed the marked wording; the text now states simply that `to_dist()` returns the standard R `dist` representation. |
| M-091 | 12.6 | Highlight | hclust () and cmdscale (). | resolved | Distance construction — `MDist` object | Removed the arbitrary list of downstream functions and stopped at the returned R `dist` representation. |
| M-092 | 13.1 | Highlight | in the gower, euclidean, and hl presets. | resolved | Distance construction — presets | Added a dedicated preset subsection and table, together with an explanation of the `u_` family and the comparison presets. |
| M-093 | 14.1 | Highlight | The custom preset exposes the three choices directly through the method_num,… | resolved | Distance construction — basic interface | Moved component specification before presets and clarified that `preset = "custom"` is the default, so it need not be requested explicitly. |
| M-094 | 14.2 | Highlight | Table 3 summarizes the presets | resolved | Distance construction — presets | Explained why the presets are included and why some require dedicated implementations outside the three-component custom interface. |
| M-095 | 15.1 | Highlight | A distance is an intermediate data representation, so its sensitivity should be… | resolved | Diagnostics — opening | Replaced the prescriptive claim with an objective motivation and explicitly separated direct distance diagnostics from downstream MDS and clustering diagnostics. |
| M-096 | 15.2 | Highlight | one predictor | resolved | Diagnostics — LOVO | Replaced “predictor” with “variable” throughout the LOVO introduction. |
| M-097 | 15.3 | Highlight | The diagnostics | resolved | Diagnostics — LOVO | Connected the diagnostic quantities directly to the change produced by recomputing the distance after each variable is omitted. |
| M-098 | 15.4 | Highlight | relative_distance | resolved | Diagnostics — LOVO | Added an itemized explanation of `mad_importance`, normalized `relative_distance`, MDS diagnostics, and optional clustering diagnostics. |
| M-099 | 15.5 | Highlight | the alienation between the full and reduced configurations. | resolved | Diagnostics — MDS-based LOVO | Introduced classical MDS, the default dimensionality, the compared configurations, and the congruence/alienation relationship before interpretation. |
| M-100 | 15.6 | Text | in the [O, | superseded | Diagnostics — MDS-based LOVO | The malformed notation disappeared with the rewritten passage. |
| M-101 | 15.7 | Highlight | measuring the configuration variance left unexplained | resolved | Diagnostics — MDS-based LOVO | Removed the unsupported variance-explained interpretation and now defines alienation as $\sqrt{1-c_j^2}$, interpreted as change in the selected MDS representation. |
| M-102 | 15.8 | Highlight | Reporting both is deliberate: a variable can shift many pairwise distances while… | resolved | Diagnostics — LOVO | Separated direct and downstream diagnostics and states that they describe different objects; MDS measures are explicitly conditional on interest in an MDS representation. |
| M-103 | 15.9 | Highlight | benchmarking | needs discussion | Diagnostics — section title | The annotation has no written explanation. The term remains in the section title because the exported function is `benchmark_mdist()`, while the prose now generally uses “comparison.” |
| M-104 | 16.1 | Highlight | Figure 2: Relative | resolved | Diagnostics — LOVO figure | Added an explicit textual reference to the figure and explained what it shows and why it is useful. |
| M-105 | 16.2 | Highlight | lovo_gower$autoplot( = metric "relative_distance", + reorder= + TRUE +) + = ggplot2::… | resolved | Diagnostics — LOVO figure | Simplified the call to `ggplot2::autoplot(lovo_gower, ...)`; retained explicit `metric` and `reorder` arguments so the displayed diagnostic and ordering are reproducible. |
| M-106 | 16.3 | Highlight | The two metrics disagree | superseded | Diagnostics — LOVO interpretation | Removed the “disagree” framing; the revision states that the measures describe different objects. |
| M-107 | 17.1 | Highlight | Comparing distance specifications | resolved | Diagnostics — subsection structure | Split the material into “Comparing LOVO diagnostics across distance specifications” and “Comparing candidate distance specifications,” making the two operations explicit. |
| M-108 | 18.1 | Highlight | Figure 3: Stability | resolved | Diagnostics — comparative LOVO figure | Added an explicit reference and a detailed explanation of the plotted within-distance PAM LOVO quantities. |
| M-109 | 18.2 | Highlight | and can disagree, | resolved | Diagnostics — organization | Removed the marked wording and replaced it with a direct-versus-downstream distinction. |
| M-110 | 18.3 | Highlight | The pairwise diagnostics | resolved | Diagnostics — candidate comparison | Reorganized candidate comparison to begin with direct magnitude measures and then separately introduce MDS and clustering comparisons. |
| M-111 | 18.4 | Highlight | changes in magnitude from changes in geometry. | resolved | Diagnostics — candidate comparison | Avoided the geometry claim and now distinguishes original dissimilarity matrices from the results of specified MDS or clustering procedures. |
| M-112 | 18.5 | Highlight | every pair of successful distances. | resolved | Diagnostics — candidate comparison | Replaced “successful distances” with “every pair of distances that could be computed.” |
| M-113 | 18.6 | Highlight | The MDS congruence coefficient and its corresponding alienation coefficient compare… | resolved | Diagnostics — candidate comparison | Now specifies classical MDS, the default two-dimensional representation, and that congruence and alienation compare those representations. |
| M-114 | 18.7 | Highlight | If is specified, 0 cluster_ k benchmark_mdist also applies each requested clustering… | resolved | Diagnostics — candidate comparison | Moved clustering into its own paragraph and explains separately when it is run and how partitions are compared. |
| M-115 | 18.8 | Highlight | Benchmarking | needs discussion | Diagnostics — terminology | Most prose now uses “compare” or “comparison,” but “benchmarking” remains in the section title and exported function name; confirm the preferred article terminology. |
| M-116 | 18.9 | Highlight | evaluates an explicit table of distance spec ifications. | resolved | Diagnostics — candidate specifications | Replaced the opaque description with a direct explanation that specifications may be supplied or generated with `all_dist_method_specs()`. |
| M-117 | 19.1 | Highlight | benchmark_pairs <- benchmark_comparisons(distance_benchmark) | needs discussion | Diagnostics — benchmark output | Simplified the example to print `benchmark_comparisons(distance_benchmark)` directly, but the annotation contains no written explanation; confirm whether further change was intended. |
| M-118 | 19.2 | Highlight | candidate_specs <- all_dist_method_specs( mode = "presets_only", + preset = c("gower",… | resolved | Diagnostics — benchmark example | Simplified the example by removing the `dplyr::mutate()` label manipulation and retaining only preset selection and the benchmark call. |
| M-119 | 19.3 | Highlight | under pairs of distances; | resolved | Diagnostics — direct comparisons | Recast the passage explicitly for two candidate distances and defines both MAD and the symmetric relative-distance measure. |
| M-120 | 19.4 | Highlight | low-dimensional geometry, | resolved | Diagnostics — MDS comparisons | Replaced the broad geometry wording with the more precise “chosen classical MDS representations.” |
| M-121 | 19.5 | Highlight | determines | resolved | Diagnostics — clustering comparisons | Uses “measure agreement” rather than “determines.” |
| M-122 | 19.6 | Highlight | whether those changes alter the partition. | resolved | Diagnostics — clustering comparisons | Now explains that decreasing ARI indicates increasingly different assignments and measures how much partitions differ. |
| M-123 | 19.7 | Highlight | are therefore complementary rather than alternative rankings. | resolved | Diagnostics — synthesis | States more strongly that the diagnostics describe different objects and should not be interpreted as alternative estimates of one quantity. |
| M-124 | 19.8 | Highlight | renders any pairwise diagnostic autoplot () as an annotated triangular heatmap. | resolved | Diagnostics — benchmark figure | Added a concrete `autoplot()` call, an explicit figure reference, and an explanation of the displayed relative-distance comparison. |
| M-125 | 19.9 | Highlight | An error in one specification is captured | superseded | Diagnostics — candidate comparison | Removed the implementation-error discussion from the article. |
| M-126 | 19.10 | Highlight | a benchmark remains a sensitivity analysis rather than a model-selection rule. | resolved | Diagnostics — synthesis | Rephrased cautiously: the diagnostics describe different aspects of candidate distances but do not by themselves select a preferred distance. |
| M-127 | 19.11 | Highlight | Gower | superseded | Diagnostics — removed results narrative | Removed the containing numerical comparison rather than retaining the ambiguous shorthand. |
| M-128 | 19.12 | Highlight | has more than twice the alienation | superseded | Diagnostics — removed results narrative | Removed the unsupported ratio interpretation. |
| M-129 | 19.13 | Highlight | a much lower PAM ARI. | resolved | Diagnostics — clustering comparison | The new text defines ARI as agreement between partitions and explicitly states that it does not establish preference. |
| M-130 | 19.14 | Highlight | This | superseded | Diagnostics — removed results narrative | Removed the vague pronoun with the containing narrative. |
| M-131 | 20.1 | Highlight | distance, and so the same gap from the Gower | superseded | Diagnostics — removed results narrative | Removed the unclear comparison with the containing narrative. |
| M-132 | 20.2 | Highlight | Taken together, the diagnostics indicate how much of an analysis rests on the choice… | resolved | Diagnostics — synthesis | Replaced the strong conclusion with a cautious statement that usefulness depends on intended use and that the diagnostics do not alone choose a distance. |
| M-133 | 20.3 | Highlight | Distance-based learning pipelines | resolved | Pipelines — opening | Added a substantive introduction explaining why data-adaptive distances must be fitted on fitting data, named the unsupervised and supervised settings, and distinguished the outer training/test split from the analysis/assessment splits used in resampling. |
| M-134 | 20.4 | Highlight | Refitting the distance within resamples | resolved | Pipelines — organization | Added a general workflow introduction before the more specific resampling subsection and renamed the latter around its motivation. |
| M-135 | 20.5 | Highlight | recipe. | resolved | Pipelines — workflow introduction | Defines a `tidymodels` recipe before introducing `step_mdist()`. |
| M-136 | 20.6 | Highlight | The fitted step stores the training data and preprocessing parameters, so the distance… | resolved | Pipelines — workflow introduction | Condensed the repeated explanation and removed the unclear “pipeline-level counterpart” sentence. |
| M-137 | 20.7 | Highlight | Its | resolved | Pipelines — workflow representations | Names `step_mdist()` explicitly when introducing its `output` argument. |
| M-138 | 20.8 | Highlight | required by the downstream task: for clustering and "pai rwise" for prediction.… | resolved | Pipelines — workflow representations | Separately explains square pairwise distances for clustering and rectangular new-to-training distances for nearest-neighbour prediction, and mentions MDS as a non-workflow downstream use. |
| M-139 | 20.9 | Highlight | whereas nearest-neighbour prediction | resolved | Pipelines — workflow representations | Introduces supervised nearest-neighbour classification before explaining its required representation and later provides a dedicated application subsection. |
| M-140 | 21.1 | Highlight | a fixed snapshot | needs discussion | Pipelines — WDI data | The phrase “fixed snapshot” remains and the annotation contains no written explanation; confirm whether Michel intended it to be removed or clarified. |

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

### 2026-09-03 — Framework reordering and render repair

- `M-082` is resolved: the numerical construction and its common transformations now precede the association-aware discussion.
- Replaced the block-level `\new{...}` wrapper around the moved equations and list with a scoped LaTeX colour group. This preserves the blue review markup while allowing Quarto to process equation identifiers and cross-references.
- Rendered the complete 34-page manuscript successfully and visually checked pages 5–10; the moved blue block, equations, references, page breaks, and Figure 1 render correctly.

### 2026-09-08 — Distance construction through learning pipelines

- Reconciled `M-086`–`M-140`; all 140 annotations have now been assessed and documented in the response log.
- In this batch, 42 comments are resolved and 8 are superseded by the rewritten material. Five remain for discussion: `M-088`, `M-103`, `M-115`, `M-117`, and `M-140`.
- Rendered the complete 35-page manuscript successfully and visually checked pages 11–32. The revised sections, code blocks, tables, plots, and section transitions have no clipping, overlap, or missing cross-references.
- Final copy/layout fixes identified in the rendered PDF: change the Table 3 caption from “clustering function” to `mdist()`; remove the manually written “Figure” before the references to Figures 2, 3, 4, and 6; replace “subsetted” with “subset”; and change the `dkss` description from kernel product to kernel summation similarity.
- Table 4 currently breaks across pages 13–14 with only the `u_mix` row on the second page. Consider tightening or repositioning it during the final layout pass.
- Terminology follow-up: the pipeline section now reserves *training/test* for the initial split and final evaluation, and uses *analysis/assessment* for the fitting and held-out portions of each resample. “Subsetted” was corrected to “subset” in the same pass. The resulting 35-page PDF compiled successfully, and pages 21–31 were visually checked without finding new layout defects.

### 2026-09-08 — Section 3 formula audit and block-level proposal

- Replaced the numerical part of the association-aware subsection with the proposed block-level construction and marked the complete passage in purple. The text explicitly states that this is a proposal and does not yet describe the current package implementation.
- Corrected the sample-correlation normalization to $1/(I-1)$; defined the unscaled principal-coordinate Manhattan distance, its mean over distinct ordered pairs, and the factor that makes the complete numerical block have empirical mean $Q_n$.
- Added the normal-theory component means, which are proportional to $\sqrt{\lambda_h}$; clarified that the construction is block- rather than component-level commensurability; and documented the consequences of retaining $\mathcal{D}<Q_n$ components, tied eigenvalues, and omitting whitening.
- Updated Figure 1 and its interpretation to show the proposed block-scaled principal coordinates, explain why the point clouds differ principally by rotation and common scaling, and display the unequal PC contributions explicitly in a fourth panel.
- Corrected the conditional-distribution notation in the categorical subsection to $\mathbf{R}^{k'\mid k}$ and $\mathbf{R}^{y\mid k}$, stated the observed-category condition, and retained the unresolved zero-mean/fallback issue explicitly.
- Implementation follow-up remains: the package currently commensurates principal components separately, uses means that include diagonal zeros, and in some new-to-training paths estimates scale factors from the rectangular distance matrix. The implementation and affected results must be updated only after the proposal is agreed.
- Rendered the complete 35-page manuscript successfully and visually checked the revised Section 3 pages. The coloured text, equations, Figure 1, and transition into Section 4 have no clipping, overlap, or broken layout.

### 2026-09-08 — Whitened block-scaling trial

- Replaced the variance-preserving PCA trial with a whitened alternative at the author's request. PCA rotation decorrelates the numerical coordinates; division by $\sqrt{\lambda_h}$ additionally gives them unit variance; a single subsequent factor scales their complete Manhattan distance to mean $Q_n$.
- Under joint normality, each whitened component has expected absolute pairwise difference $2/\sqrt{\pi}$ and, after block scaling, expected contribution $Q_n/\mathcal{D}$. Outside normality, their empirical mean absolute contributions need not be exactly equal.
- Replaced $K$ by $\mathcal{D}$ for the number of retained components, reserving $K$ for the number of clusters in later sections.
- Revised Figure 1 makes the direction-specific rescaling visible through the approximately circular whitened configuration. After reviewing a four-panel trial, the contribution bar chart was removed and the concise three-panel layout restored.
- The methodological decision remains open. Whitening deliberately amplifies low-variance directions; the common block factor controls only the total numerical weight and does not undo that amplification.
- Condensed the surrounding explanation after the trial, retaining the definition, block balance, Mahalanobis relationship, and dimensionality choice.
- The author shortened the passage further by removing the rendered normal-theory and singularity commentary; the defining equations and principal interpretation remain intact.
- Rendered the complete manuscript successfully and visually checked the revised Section 3 pages. The equations, three-panel figure, caption, and page transitions have no clipping or overlap.
