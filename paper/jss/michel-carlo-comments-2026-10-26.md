# Michel and Carlo: JSS paper comments, October 2026

Source: [Annotated PDF](article_cc_mvdv.pdf). Extracted on 2026-10-26.

## Scope and extraction notes

- 99 annotation records across 23 of the PDF's 36 pages: Michel 72; Carlo (PDF author `carlocavicchia`) 27.
- Types: 86 highlights, 11 strikeouts, one underline, and one sticky-note reply.
- 90 annotations contain written comments; nine are markings without written comments. The reply is retained as a separate record and linked to its parent.
- Page numbers refer to this annotated PDF, not necessarily the current manuscript.
- Marked text is reconstructed from annotation coordinates; spacing and mathematical notation may differ from the original. Empty or ambiguous markings are explicitly flagged.
- Reviewer wording is preserved, including proposed code, typos, and claims that require verification. These are suggestions to evaluate, not instructions automatically applied.
- No manuscript, package, or existing response/tracker statuses have been changed for this review round.

## Suggested order of work (editorial assessment, not reviewer quotations)

1. **Benchmark interpretation and scaling** (R2-066, R2-069, R2-072, R2-075; pp. 20-22). Clarify what each measure answers. Symmetric relative distance is unchanged by a common rescaling of both matrices, but it does not remove a scale difference between them. Decide whether to normalize candidates for structural comparisons or explicitly distinguish scale differences from downstream MDS/clustering changes. Michel prefers alienation and ARI for the displayed illustration.
2. **Distance construction and response handling** (R2-005, R2-012 to R2-017, R2-022 to R2-026; pp. 8, 11-13). Explain the returned MDist object before using its fields; distinguish components from presets; name the block-scaling option; review response-aware examples, warnings, and train/test wording against actual behavior.
3. **LOVO exposition and examples** (R2-028, R2-037, R2-041 to R2-055; pp. 14-19). Separate direct-distance, MDS, and clustering diagnostics; state defaults and reordering; shorten plotting code; explain redundant diagnostic columns and equivalent specifications.
4. **Simpler learning examples** (R2-006, R2-078 to R2-095; pp. 9, 23-27). Clarify predictor-data versus precomputed-distance input, returned object versus its distance matrix, and mapping/recipe operations. Discuss prepared example datasets before making any package/data additions.
5. **Author decisions, not automatic changes** (R2-001, R2-002, R2-074, R2-089, R2-098). Decide whether to restore a short leakage demonstration, add a timing comparison, alter WDI predictors, and assign responsibility for clean-session replication materials. The earlier decision to remove the leakage experiment is not reversed by these comments.
6. **Copyediting and submission consistency** (including R2-003, R2-004, R2-018, R2-051, R2-096 to R2-099). Correct names/typos, align conclusions with shown results, remove promotional phrasing, and reconcile computational-details and bibliography versions. Statements in the comments about JSS requirements and review times have not been independently verified here.

## Technical caution before applying proposed wording

In R2-037, Michel proposes describing `ac_importance` as “1 minus the alienation coefficient.” The current implementation instead stores the alienation coefficient itself:
`sqrt(pmax(0, 1 - cc^2))` in `manydist/R/lovo_mdist.R`.
Larger alienation means greater MDS change, whereas larger congruence means greater similarity. Do not adopt the proposed formula literally without resolving that distinction. This is separate from clustering importance, which uses `1 - ARI`.

## Complete annotation inventory


## Page 1

### R2-001 — Carlo — Highlight

Marked passage:

> out-of-sample evaluations.

Reviewer comment (verbatim):

> We promise proper out-of-sample evaluation here, and assert the leakage argument again on p. 4, p. 13 and in the Conclusion, but with the old Section 6.2 removed there is no longer anywhere in the paper where the two procedures are contrasted.

### R2-002 — Carlo — Highlight

Marked passage:

> Journal of Statistical Software

Reviewer comment (verbatim):

> Not a manuscript edit, but for the submission package: JSS requires replication materials alongside the PDF and the package source, preferably as a single standalone script, and reviewers are instructed to run it. Worth agreeing who assembles it and confirming it reproduces every table and figure from a clean session, including the set.seed() calls. Also not an edit: at 36 pages we are past the 30-page mark at which JSS warns that review times lengthen considerably. No reason to cut anything on that basis, but worth knowing when we set expectations on timing.

## Page 3

### R2-003 — Carlo — Highlight

Marked passage:

> Section 6 shows how
> distances can be embedded in tidymodels workflows.

Reviewer comment (verbatim):

> Section 6.1 now opens with knn_dist(), which is deliberately outside tidymodels, before moving to recipes and workflows. The sentence could say Section 6 moves from a direct interface to tidymodels workflows.

## Page 5

### R2-004 — Carlo — Highlight

Marked passage:

> clustrd::cluspca()
> (Markos et al. 2019)

Reviewer comment (verbatim):

> This should be cluspcamix(). cluspca() is the continuous-only function; our own text on p. 5 already calls it cluspcamix() a few lines above.

## Page 8

### R2-005 — Carlo — Highlight

Marked passage:

> manydist
> also provides a groupwise commensurability option,

Reviewer comment (verbatim):

> The option is described but never named, so a reader cannot use it. Could we give the argument and its value, or say which preset reaches it?

## Page 9

### R2-006 — Michel — Highlight

Marked passage:

> The examples use the Palmer penguins data from the palmerpenguins package. We retain
> four numerical variables, two categorical variables, and the categorical response species,
> which will later be used as a response variable. Observations with missing values are removed.
> The categorical variables are stored as factors, as required for their categorical treatment by
> mdist().

Reviewer comment (verbatim):

> But why not remove all of this and store the prepared data in the package? Then: data(penguin) (or something like that if we want to change the data name) would suffice. 
> Actually: 
> penguins<-na.omit(penguins)[,-8] 
> also gives what you wrote

## Page 11

### R2-007 — Carlo — Highlight

Marked passage:

> specifiy

Reviewer comment (verbatim):

> specify

### R2-008 — Michel — StrikeOut

Marked passage:

> unnormalized sum can therefore be obtained as

Reviewer comment (verbatim):

> distance, based on all variables except the "response" variable "species" in column 1, can be obtained as

### R2-009 — Michel — Highlight

Marked passage:

> R> library("manydist")

Reviewer comment (verbatim):

> would move this to before (with prepared penguin_train and penguin_test to skip all that)

### R2-010 — Michel — StrikeOut

Marked passage: not reliably recoverable from the annotation geometry; consult the PDF.

No written comment. Strikeout indicates a proposed deletion; confirm its exact extent in the PDF.

### R2-011 — Michel — Highlight

Marked passage:

> penguin_train |> dplyr::select(-species),

Reviewer comment (verbatim):

> or just: penguin_train[,-1]

### R2-012 — Michel — Highlight

Marked passage:

> This differs from the usual Gower distance only by the normalization by the number of
> variables. The same distance on the conventional Gower scale is available directly through the
> preset:
> "gower"

Reviewer comment (verbatim):

> Note that this differs from common implementations (e.g., those in daisy and gower) as those divide the gower distance by the number of variables. In the mdist package, this "average" gower distance can be obtained by dividing by the number of variables or by selecting the gower preset. We discuss the preset in Section 4.3.

### R2-013 — Michel — StrikeOut

Marked passage:

> R> d_gower <- mdist(penguin_train |> dplyr::select(-species),preset = "gower")
> R> p <- penguin_train |> dplyr::select(-species) |> ncol()
> R> all.equal(
> +
> as.matrix(d_gower_sum$distance) / p,
> +
> as.matrix(d_gower$distance)
> +)
> [1] TRUE
> Thus, the preset provides the normalized Gower distance, corresponding to standard imple-
> mentations such as those in cluster and gower.

Reviewer comment (verbatim):

> I think this is not needed. Also: it is in a sense confusing as we use parts that have not been explained. For example: $distance gives the distance object. Etc.
> In fact, I suppose that is what we have to provide/explain here. What do we get when we do this?

### R2-014 — Michel — Highlight

Marked passage:

> Other distance specifications

Reviewer comment (verbatim):

> So before this, we need to say what the mdist call yields. For example by doing print(d_gower). Explaining some of that output (which I like) and perhaps observing that the distance object is d_gower$distance
> (I notice we do print later; should be here I think)

### R2-015 — Michel — Highlight

Marked passage:

> For example,

Reviewer comment (verbatim):

> Why not give an example: 
> d_default<-mdist(penguins_train[,-1], method_num="std", method_cat="tvd", commensurable="TRUE")
> This is in fact the default setting and could also be obtained by simplyb calling d<-mdist(penguins_train[,-1]).
> Note that when commensurable =TRUE, etc.

### R2-016 — Michel — Highlight

Marked passage:

> The Gower distance introduced above

Reviewer comment (verbatim):

> The average Gower distance implemented in daisy, cluster, gower etc., as mentioned above, is also on of the presets: preset="gower". Table 4 ...

### R2-017 — Michel — Highlight

Marked passage:

> Presets

Reviewer comment (verbatim):

> In addition to constructing a mixed distance by specifying all components (method_cat, method_num and commensurable), several mixed distance variants are available in mdist as presets. This includes variants from practice or recent literature that are not additive, and hence cannot be obtained by specifying the individual components. Table 4  provides a complete list. Note that the Gower etc.

## Page 12

### R2-018 — Carlo — Highlight

Marked passage:

> gudmm

Reviewer comment (verbatim):

> We expand gudmm as “Generalized Unified Distance Metric”, but the cited Mousavi and Sehhati (2023) title on p. 35 reads “generalized multi-aspect distance metric” — one of the two needs correcting.

### R2-019 — Carlo — Highlight

Marked passage:

> Distance

Reviewer comment (verbatim):

> Distances

### R2-020 — Michel — Highlight

Marked passage:

> d_u_mix |> print()

Reviewer comment (verbatim):

> so this already has been done

## Page 13

### R2-021 — Carlo — Highlight

Marked passage:

> setting,

Reviewer comment (verbatim):

> settings,

### R2-022 — Michel — Highlight

Marked passage:

> d_response <- mdist(penguin_train,response = species,preset = "u_mix")

Reviewer comment (verbatim):

> But I find this a bit weird. I would prefer it that when specifying a response, the user has to indicate a method_cat. That way, the user is forced to think about it (and select an association based one). I do not mind it defaulting to tvd. But, I would not want to do this using preset.
> I also noticed that I can do preset=gower. But, what happens then? I do not get a warning. However, when I do preset=hl I do get a warning...
> I think it should be: response=species, method_cat="tvd", commensurable=TRUE)

### R2-023 — Michel — StrikeOut

Marked passage:

> rectangular

No written comment. Strikeout indicates a proposed deletion; confirm its exact extent in the PDF.

### R2-024 — Michel — Highlight

Marked passage:

> the weights and the scaling quantities

Reviewer comment (verbatim):

> commensurability weights and scaling quantities? Or: using the settings also used for calculating the distances for the training observations

### R2-025 — Michel — Highlight

Marked passage:

> d_response_new <- mdist(
> +
> penguin_train,
> +
> response = species,
> +
> new_data = penguin_test |> dplyr::select(-species),
> +
> preset = "u_dep"
> +)

Reviewer comment (verbatim):

> How about doing two: distance based on gower (non-supervised). distance using tvd. Than first you do not have the response field. Second, you make clear that for the training, the response needst to be supplied whereas for the test this needs to be surpressed. Also: let's not do u_dep.
> Also, probably helpful to explain what we get: distance is a n_test x n_train matrix of distance from the test points to the training points. (You explained above when/how this can be useful. Perhaps, first the distance calculation, then say what we get and how it is useful (so move that part here)

### R2-026 — Michel — Highlight

Marked passage:

> is not supplied for
> 
> species
> penguin_test:

Reviewer comment (verbatim):

> should not/cannot? What happens if I try? I get an error (good: the error may could be made more informative. Not sure if possible (it now says that response is matched to two fields. That is wrong of course. It would be nice to get something like: response is only possible for the training data. Not important though

### R2-027 — Michel — Highlight

Marked passage:

> MDS

Reviewer comment (verbatim):

> multidimensional scaling (ref?)

## Page 14

### R2-028 — Carlo — Highlight

Marked passage:

> mad_normalized

Reviewer comment (verbatim):

> mad_normalized is printed but never defined — the bullets above introduce only mad_importance and relative_distance — and its values are identical to relative_distance in all six rows. Either drop one column from the output, or state how they differ.

### R2-029 — Carlo — Highlight

Marked passage:

> (e.g,

Reviewer comment (verbatim):

> (e.g.,

### R2-030 — Carlo — Highlight

Marked passage:

> cluster).

Reviewer comment (verbatim):

> clusters).

### R2-031 — Michel — Highlight

Marked passage:

> Diagnostics

Reviewer comment (verbatim):

> Additional diagnostics

### R2-032 — Michel — StrikeOut

Marked passage:

> is

Reviewer comment (verbatim):

> :

### R2-033 — Michel — StrikeOut

Marked passage:

> is

Reviewer comment (verbatim):

> : a normalized

### R2-034 — Michel — StrikeOut

Marked passage: not reliably recoverable from the annotation geometry; consult the PDF.

Reviewer comment (verbatim):

> That is,

### R2-035 — Michel — StrikeOut

Marked passage:

> the

No written comment. Strikeout indicates a proposed deletion; confirm its exact extent in the PDF.

### R2-036 — Michel — Highlight

Marked passage: not reliably recoverable from the annotation geometry; consult the PDF.

Reviewer comment (verbatim):

> s

### R2-037 — Michel — Underline

Marked passage: not reliably recoverable from the annotation geometry; consult the PDF.

Reviewer comment (verbatim):

> But the bullet list below is not so nice. We could perhaps do:
> For the impact on mds, we assess the similarity of the full mds configuration with that of those with a variable left-out. As measure for this similarity we consider "alienation importance" (ac_importance) which is defined as 1 minus the alienation coefficient (reference): Higher values indicate a bigger impact of a variable on the mds solutions whereas smaller values etc.  For completeness, concruence importance,  which is etc. etc. , is also calculated and stored as cc_importance. 
> For the impact on cluster analysis, we consider the adjusted rand index (ARI) of a cluster analysis with all variables to that where a variable is left out. Depending on the chosen cluster algorithm, we have pam_importance, hclust_importance and spectral_importance, which are all defined as one minus the ARI so that higher values indicate more influence on the clustering solution.

### R2-038 — Michel — Highlight

Marked passage:

> lovo_gower$results

Reviewer comment (verbatim):

> This after the explanation of the plot. Or maybe even remove?

### R2-039 — Michel — Highlight

Marked passage:

> the measure

Reviewer comment (verbatim):

> to be plotted

## Page 15

### R2-040 — Carlo — Highlight

Marked passage:

> with some defaults?).

No written comment. Highlight retained without inferring a requested change.

### R2-041 — Carlo — Highlight

Marked passage:

> These MDS diagnostics should be distinguished direct distance diagnostics.

Reviewer comment (verbatim):

> These MDS diagnostics should be distinguished from direct distance diagnostics.

### R2-042 — Michel — Highlight

Marked passage:

> ggplot2::autoplot(
> +
> lovo_gower,
> +
> metric = "relative_distance",
> +
> reorder = TRUE
> +)

Reviewer comment (verbatim):

> But this is the default. So why not state that. Also, as we just explained about downstream measures, perhaps no harm in to be clear that here we do not have this and hence we only have two choices: mad_importance or relative_distance, with relative_distance the default when plotting.

### R2-043 — Michel — Text

Reply to R2-041.

Marked passage: not reliably recoverable from the annotation geometry; consult the PDF.

Reviewer comment (verbatim):

> But also: What do we mean by these? This doesn not correspond to the text. Some connection/introduction is still needed

## Page 16

### R2-044 — Michel — Highlight

Marked passage:

> R> a <- ggplot2::autoplot(
> +
> lovo_gower,
> +
> metric = "ac_importance"
> +)
> R> b <- ggplot2::autoplot(
> +
> lovo_gower,
> +
> metric = "hclust_importance"
> +)

Reviewer comment (verbatim):

> Let's only do this. Not all the code to get the plot as is. Just:
> ggplot2::autoplot(lovo_gower, ac_importance)
> and
> ggplot2::autoplot(lovo_gower, hclust_importance)
> The resulting plots represented next to each other in Figure 2.
> 
> All the formatting you do should either be standard (for example: why change the axes labels? If we do not like what we give, we should change that). (But I think that the plots are fine without all this work)

### R2-045 — Michel — StrikeOut

Marked passage:

> R> variable_labels <- c(
> +
> island = "Island", sex = "Sex",
> +
> bill_length_mm = "Bill\nlength",
> +
> bill_depth_mm = "Bill\ndepth",
> +
> flipper_length_mm = "Flipper\nlength",
> +
> body_mass_g = "Body\nmass"
> +)
> R> a <- a + ggplot2::labs(
> +
> title = "(a) MDS configuration", y = "Alienation coefficient", x = NULL
> +)
> R> b <- b + ggplot2::labs(
> +
> title = "(b) Hierarchical clustering", y = "1 - ARI", x = NULL
> +)
> R> lovo_plot_style <- list(
> +
> ggplot2::scale_x_discrete(labels = variable_labels),
> +
> ggplot2::scale_y_continuous(limits = c(0, NA)),
> +
> ggplot2::theme(
> +
> legend.position = "none",
> +
> plot.title = ggplot2::element_text(size = 10, face = "plain"),
> +
> axis.title = ggplot2::element_text(size = 9),
> +
> axis.text = ggplot2::element_text(size = 8),
> +
> axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
> +
> panel.grid.minor = ggplot2::element_blank()
> +
> )
> +)
> R> a <- a + lovo_plot_style
> R> b <- b + lovo_plot_style
> R> patchwork::wrap_plots(a, b, nrow = 1)

No written comment. Strikeout indicates a proposed deletion; confirm its exact extent in the PDF.

## Page 17

### R2-046 — Carlo — Highlight

Marked passage:

> `Unbiased Gower`

Reviewer comment (verbatim):

> By our own statement on p. 11 — that affine rescalings cancel under commensurability, so “std”, “range” and “robust” give identical training distances — this specification is u_indep under another name. Both diagnostics in Figures 3 and 4 are scale-invariant, so this series must be numerically identical to u_indep, which then appears separately in Section 5.2 and Section 6 with nothing connecting the two names. Either we say in one line that it coincides with u_indep, which makes the commensurability point nicely, or we swap it for u_dep or u_mix so the figure shows a genuinely distinct fourth specification. I would prefer the second.

### R2-047 — Carlo — Highlight

Marked passage:

> Compare LOVO diagnostics across distance specifications

Reviewer comment (verbatim):

> Imperative where every other heading is nominal, and close enough to Section 5.2’s title two pages later to be confusing. “LOVO diagnostics across specifications” would separate them.

### R2-048 — Michel — Highlight

Marked passage: not reliably recoverable from the annotation geometry; consult the PDF.

Reviewer comment (verbatim):

> and.....(explain what happens then in the compare_function.

### R2-049 — Michel — Highlight

Marked passage:

> For example,

Reviewer comment (verbatim):

> To see how this works in practice, consider the following example where we compare the leave-one-variable-out results for different distances. First we define a list with the distance variants. Here: Gower, Euclidean, HL and unbiased Gower:
> variants etc.

### R2-050 — Michel — Highlight

Marked passage:

> ggplot2::autoplot(
> +
> lovo_comparison,
> +
> metric = "relative_distance",

Reviewer comment (verbatim):

> But as this is the defualt we can leave the metric. And/or, mention that this is the default. Then move onto the pam_importance

## Page 19

### R2-051 — Carlo — Highlight

Marked passage:

> If we want to compare the Gower, unbiased, and unbiased preset variants we can do:

Reviewer comment (verbatim):

> The three presets requested are gower, u_indep and u_dep, so this should read “the Gower, unbiased independent and unbiased dependent presets”.

### R2-052 — Carlo — Highlight

Marked passage:

> distant-based

Reviewer comment (verbatim):

> distance

### R2-053 — Michel — Highlight

Marked passage:

> Other diagnostics introduced above can be displayed by changing the
> argument.
> metric

Reviewer comment (verbatim):

> We also need to explain the reorder field. The default is FALSE. At least, when I do not use it, the PAM importance plot changes (the order of variables, sex&island and also others). For the first plot when I just do default, I get exactly the same....

### R2-054 — Michel — Highlight

Marked passage:

> Note that, for association-based specifications, leaving out a variable can change not only the
> available information but also the construction of the remaining contributions. Under
> u_dep,
> 
> for example, removing either categorical variable leaves only one categorical variable. Because
> the association-based categorical construction requires at least two categorical variables, the
> reduced specification falls back to matching.

Reviewer comment (verbatim):

> well actually, it already does this for u_mix, not? (The u_dep also influence the continous variables). Perhaps better mention u_mix here. 
> Also: we could add that a warning is given when lovo is done with an association-based distance and 2 or less variables. (Actually: When we do association-based with only 1 cat variable we also get the warning. Which is good. Not sure where/how to mention this but here could be a good place

### R2-055 — Michel — Highlight

Marked passage:

> mind.

Reviewer comment (verbatim):

> Maybe we can also add a summary remark about the why of lovo comparison. Something like:
> 
> By comparing the importances of variables in distance calculations accross distance definitions, biases or other undesirable properties can be exposed. Van de Velden et al., for example, use it to show how by construction certain mixed distance specification either over or under emphasize the contributions of categorical variables.

### R2-056 — Michel — Highlight

Marked passage:

> Candidate distance

Reviewer comment (verbatim):

> Recall that distance specifications in mdist are described (defined?) by four fields:

### R2-057 — Michel — Highlight

Marked passage:

> Named presets

Reviewer comment (verbatim):

> The named presets

### R2-058 — Michel — StrikeOut

Marked passage:

> rows with
> use the other three fields explicitly.
> preset = "custom"

Reviewer comment (verbatim):

> rows of what?

### R2-059 — Michel — Highlight

Marked passage:

> the function
> benchmark_mdist()

Reviewer comment (verbatim):

> But below we start with all_dist_method_specs? What is the organization here? I think we need to me more structured/clear about this: the benchmark_mdist can be used to compare. The variants to be compared need to be specified by supplying the four fields. Now: it seems this can only be done through the all_dist_method_specs(). Is that true? Can't I do something similar to what was done for lovo comparison? (I tried, I cannot. For the benchmark the variants need to be in a different format (tibble?). That is ok, but than perhaps we need to be explicit and state this. I suppose it is possible to construct candidate_specs in the correct format also without filtering from all_dist_method_specs. Do we need to explain? Is that trivial? Not important? Actually, for me it is not so easy to do this. For example, if I want to use the same variants as in the previous section: gower, euclidean, hl, unbiased gower, how would I do that? (Well: I managed to do that running two filters and rbind.data.frame...I guess one needs to build a data frame with the (named?) columns. But: Can we also name the variants? (Now: it uses presets,  and for customs the settings. I guess this is fine, but obviously becomes a bit lengthy)

### R2-060 — Michel — Highlight

Marked passage:

> of custom

Reviewer comment (verbatim):

> all? or no. But, what 413 combinations? otherwise seems rather random. I think we do need to explain that combinations of these fields, where we removed duplicates. That is, as commensurable=TRUE, voids all indepedendent numerical standardization methods, these candidates are not provided). (Actually I see we do say this later, so just change order/exposition a bit)

### R2-061 — Michel — Highlight

Marked passage:

> "std"

Reviewer comment (verbatim):

> better none

### R2-062 — Michel — Highlight

Marked passage:

> Alternatively, the specifications table can be constructed manually.

Reviewer comment (verbatim):

> What does that mean?

### R2-063 — Michel — Highlight

Marked passage:

> successful candidates,

Reviewer comment (verbatim):

> successful in what sense?

### R2-064 — Michel — Highlight

Marked passage:

> "ari_pam").
> metric =

Reviewer comment (verbatim):

> This is fine, but earlier we had pam_importance. It makes sense not to do that here (we do not consider lovo). But perhaps worth mentioning somewhere

### R2-065 — Michel — Highlight

Marked passage:

> x =

Reviewer comment (verbatim):

> I think previously we didn't do this. Why start?

## Page 20

### R2-066 — Carlo — Highlight

Marked passage:

> 5.663
> 5.663

Reviewer comment (verbatim):

> Both gower comparisons print identically — mad 5.663 against u_dep and against u_indep, relative_distance 1.773 for both — and this is forced rather than a copy error. Solving backwards from the three printed pairs gives mean|u_dep| = mean|u_indep| = 6.025 and mean|gower| = 0.363, so mean(u) − mean(gower) = 5.662 against the printed 5.663. Commensurability drives every unbiased preset’s mean distance to about p = 6, so any commensurable preset sits the same distance from Gower and mad reports the scale gap, not structure. The awkward part is that the paragraph below introduces relative_distance as the remedy for differing scales, but r tends to 2 as the scales diverge and ceilings at 1.886 here. The three measures that do discriminate — mds_congruence 0.942 vs 0.992, alienation 0.336 vs 0.130, ari_pam 0.481 vs 0.393 — are the ones we pass over. Either put the candidates on a common scale first, or state that mad and relative_distance are interpretable only within a shared scale and that cross-commensurability comparisons should be read from the MDS and ARI columns.

### R2-067 — Michel — Highlight

Marked passage:

> summary(distance_benchmark)

Reviewer comment (verbatim):

> I would do: metric="ari_pam" then we only get one. Write that if metric is not given, all pairwise tables are given. 
> Also: Note that the last row can be removed (it is empty by construction)

### R2-068 — Michel — Highlight

Marked passage:

> u_indep

No written comment. Highlight retained without inferring a requested change.

### R2-069 — Michel — Highlight

Marked passage:

> For two candidate distances
> and
> ,
> is the mean absolute difference between their
> mad
> d a
> d b
> pairwise dissimilarities. When comparing distance specifications, the overall scale of the
> distance values can differ considerably, a direct comparison of the absolute difference may
> therefore not be appropriate. Hence,
> is also provided, defined as:
> relative distance

Reviewer comment (verbatim):

> Also in light of Carlo's comment, I think we need to be clear beforehand what is interesting/relevant and how to read the comparisons. Personally, ARI is pretty clear I think: it compares the two cluster solutions. From below, I now also understand relative distance: It basically compares two distance matrices, assesses the absolute differences and normalizes (divides by total)...But: Why/when is this interesting? 
> Practically speaking, I would say that the impact on an analysis is relevant here. So: ARIS and congruences (although, that is already harder to understand). Even more relevant would be the impact on external information (which I suppose we will get to: classification)

## Page 21

### R2-070 — Carlo — StrikeOut

Marked passage:

> the

No written comment. Strikeout indicates a proposed deletion; confirm its exact extent in the PDF.

### R2-071 — Carlo — Highlight

Marked passage:

> algorithn(s)

Reviewer comment (verbatim):

> algorithm

### R2-072 — Michel — Highlight

Marked passage:

> by calculating the the congruence coefficient
> and its corresponding alienation coefficient.

Reviewer comment (verbatim):

> I think we only do alienation? And: we didn't explain this/show metric="alienation". And it probably needs some explanation. Also: I think we also allow mds_congruence (but: not mds_alienation...). I think it is all fine, doesn't matter much but some consistency is nice. (Also: alienation means lower is better, this must be mentioned somewhere)

### R2-073 — Michel — Highlight

Marked passage:

> ggplot2::autoplot(
> +
> distance_benchmark,
> +
> metric = "relative_distance"
> +)

Reviewer comment (verbatim):

> OK. But, if you do custom (I did) we have very long labels (the ingredients of the distance)

## Page 22

### R2-074 — Carlo — Highlight

Marked passage:

> 6. Distance-based learning pipelines

Reviewer comment (verbatim):

> This is where the removed Section 6.2 sat — the single-matrix versus per-fold refit comparison, with the optimism table and the rank flips. The removal is clean, so this is a decision rather than an oversight, but the leakage argument is still asserted in the abstract, on p. 4, on p. 13 and in the Conclusion with nothing shown. In fairness to the cut, the effect was small against the paired standard errors in two of three rows and rested on five folds. My own view is that a short version is worth restoring, framed as a demonstration of the mechanism rather than a magnitude claim, since it is the one result none of the nine packages in Table 2 can produce. If we keep it out, the abstract’s promise should be softened to match.

### R2-075 — Michel — Highlight

Marked passage:

> The same display can be requested for raw MAD, MDS congruence or alienation, and, when
> clustering has been requested, ARI.

Reviewer comment (verbatim):

> Actually: I would prefer alienation and ari

### R2-076 — Michel — Highlight

Marked passage:

> this

Reviewer comment (verbatim):

> ?

### R2-077 — Michel — Highlight

Marked passage:

> knn_dist(),

Reviewer comment (verbatim):

> a function from our package?

## Page 23

### R2-078 — Michel — Highlight

Marked passage: not reliably recoverable from the annotation geometry; consult the PDF.

Reviewer comment (verbatim):

> Why this? Why not knn_preds ? The result is not a table. Just the predictions...

### R2-079 — Michel — Highlight

Marked passage:

> knn_predictions <- tibble::tibble(
> +
> .pred_class = tabular_predictions
> +)

Reviewer comment (verbatim):

> This just makes a data frame (or tibble, seems the same here: I would do knn_preds_df<-as.data.frame(knn_preds)
> seems to give the same...(Only: different column name, your "tibble" has this strange: .pred_class name. I guess that is needed somewhere...

### R2-080 — Michel — Highlight

Marked passage:

> dplyr::bind_cols(penguin_test["species"], knn_predictions) |>
> +
> yardstick::accuracy(truth = species, estimate = .pred_class)
> # A tibble: 1 x 3
> .metric
> .estimator .estimate
> <chr>
> <chr>
> <dbl>

Reviewer comment (verbatim):

> I do not see the point of doing this. In fact, after getting the predicitions using knn_dist. I would simply do:
> 
> table(knn_preds,penguin_test$species)
>  
> This gives the confusion table. Clearly accuracy =1 .

### R2-081 — Michel — Highlight

Marked passage:

> No training-to-training
> matrix needs to be supplied:

Reviewer comment (verbatim):

> as opposed to the previous case? I do not understand...
> What we will get here is a matrix with distances between all test points (84) to all training points (249). I would state that here, than the code.

### R2-082 — Michel — Highlight

Marked passage:

> The resulting matrix has 84 rows and 249 columns.

Reviewer comment (verbatim):

> But actually, the output is not a matrix but a mdist object. With same functionality (print and summary as before). It also has $distance. That is hte 84 x 249 matrix

### R2-083 — Michel — Highlight

Marked passage:

> must be in the same order as the training outcomes.

Reviewer comment (verbatim):

> How could it be different?

### R2-084 — Michel — Highlight

Marked passage:

> distance_predictions <- knn_dist(
> +
> x = d_knn_test$distance,
> +
> y = penguin_train$species,
> +
> k = 5,
> +
> type = "class"
> +)

Reviewer comment (verbatim):

> OK. So: instead of giving the xtrain and xtest, we give the distance between xtrain and xtest. That is the summary, correct?

## Page 24

### R2-085 — Michel — Highlight

Marked passage:

> tabular

Reviewer comment (verbatim):

> I do not like this notation. Both have tabular input...something along data vs distance input...

### R2-086 — Michel — Highlight

Marked passage:

> parsnip

Reviewer comment (verbatim):

> no reference/explanation needed?

### R2-087 — Michel — Highlight

Marked passage:

> consumed

No written comment. Highlight retained without inferring a requested change.

### R2-088 — Michel — Highlight

Marked passage:

> R> penguin_knn_spec <- nearest_neighbor_dist(

No written comment. Highlight retained without inferring a requested change.

## Page 25

### R2-089 — Carlo — Highlight

Marked passage:

> R> wdi_clustering <- wdi_2022 |>
> +
> dplyr::select(-country)

Reviewer comment (verbatim):

> This keeps all three classifications, so the clustering runs on both income and lending, whereas Section 6.3 on p. 27 drops lending as administratively related to income. Our own Introduction flags exactly this redundancy for independence-based distances. Either drop lending here too, or add a line on why it is acceptable in a descriptive pass. Income is also both a clustering variable here and the target in 6.3;  defensible, but currently unremarked.

### R2-090 — Michel — Highlight

Marked passage:

> clustering_recipes <- purrr::map(
> +
> distance_presets,
> +
> ~ recipes::recipe(~ ., data = wdi_clustering) |>
> +
> step_mdist(
> +
> recipes::all_predictors(),
> +
> preset = .x,
> +
> output = "pairwise"
> +
> )
> +)

Reviewer comment (verbatim):

> This is where it happens. To me a lot of things happenning. Not just recipes. Also map....

### R2-091 — Michel — Highlight

Marked passage:

> is combined with partitioning around medoids (PAM) and spectral clustering.

Reviewer comment (verbatim):

> But that is a choice. Or not? This needs to be done later, by applying to the distance matrix.

### R2-092 — Michel — Highlight

Marked passage:

> data("wdi_2022", package = "manydist")

Reviewer comment (verbatim):

> why not data(wdi_2022) ?

## Page 26

### R2-093 — Michel — Highlight

Marked passage:

> clustering_recipes$u_mix,

Reviewer comment (verbatim):

> is this the distance matrix corresponding to u_mix?

### R2-094 — Michel — Highlight

Marked passage:

> The ARI in Table 5 instead

Reviewer comment (verbatim):

> How?

## Page 27

### R2-095 — Michel — Highlight

Marked passage:

> We now predict income group, comparing the five distance specifications under a common
> resampling design. The recipes use
> to provide the rep-
> output = "distance_to_training"
> resentation required by
> Countries labelled
> nearest_neighbor_dist().
> "Not classified"
> are removed. Country is an identifier, and lending category is excluded because it is admin-
> istratively related to income group. The two-country North America level is pooled with
> Latin America and the Caribbean to reduce the risk that a categorical level is absent from an
> analysis set.

Reviewer comment (verbatim):

> I have the impression that here we again are first going to do some preparations that are irrelevant for the package explanation. Better to have all done outside and give the widi_classification data in the package

## Page 32

### R2-096 — Carlo — Highlight

Marked passage:

> The Palmer penguins data showed how scale
> and association shape the distance geometry,

Reviewer comment (verbatim):

> With the earlier geometry illustration gone, this is carried only by the all.equal() demonstrations in Section 4.2, with no figure. Either narrow the sentence or bring the illustration back.

## Page 33

### R2-097 — Carlo — Highlight

Marked passage:

> The tidymodels integration eliminates the need to write and maintain
> manual loops by automating repetitive model-fitting and evaluation tasks within
> a consistent framework. This makes the code cleaner and easier to reproduce,
> while simplifying preprocessing, parameter tuning, resampling, and performance
> comparison across multiple models.

Reviewer comment (verbatim):

> These two sentences read as package promotion and sit oddly in a paragraph whose successor opens “The framework has deliberate limits”. They are also still bold. I would cut them, or keep only the substantive half in our own register.

### R2-098 — Carlo — Highlight

Marked passage:

> through the compiled and parallel back-ends already available for
> homogeneous data. By

Reviewer comment (verbatim):

> This concedes we do not have those back-ends, and there is no timing comparison anywhere. The JSS guide asks specifically for advantages and disadvantages against existing solutions, and gower, Rfast and parallelDist are in Table 1 precisely because the computation is O(n²). A small timing table against those three would answer it before a referee asks.

### R2-099 — Carlo — Highlight

Marked passage:

> The results in this article were obtained using R 4.5.2 with manydist 0.5.2, tidymodels 1.4.1,
> dplyr 1.2.0 and palmerpenguins 0.1.1. R itself and all packages used are available from CRAN
> at https://CRAN.R-project.org/. The development version of manydist and its documentation
> are hosted at https://alfonsoiodicede.github.io/manydist_package/.

Reviewer comment (verbatim):

> The code also calls purrr, forcats, tidyr and patchwork, which tidymodels does not cover. Extending the same sentence with their versions is enough. Also: this says manydist 0.5.2 while the package entry on p. 34 says 0.5.1 — the two should agree.

