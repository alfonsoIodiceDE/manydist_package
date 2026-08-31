# Michel's comments on the JSS paper

Source PDF: `jss paper 0708.pdf`

Extracted annotations: **140** across **20 pages**.

Annotation types: Highlight (124), StrikeOut (12), Text (4).

The quoted passages below were reconstructed from the annotation coordinates. Minor spacing differences may occur around equations or inline code.

## Page 1

### 1. Highlight

> **Marked text:** through multivariate addi tivity and commensurability.

**Michel's comment:**

Not sure we need this here. 
We could perhaps replace the last part by: constructing mixed-variable distances that correct for such biases.

### 2. StrikeOut

> **Marked text:** and

**Action:** Delete the marked text.

### 3. Highlight

> **Marked text:** and the

**Michel's comment:**

lot of commas and ands. Perhaps this last part separate? (In addition, furthermore we show/illustrate how its functionality can be integrated in the tidymodels....etc.?

### 4. Highlight

> **Marked text:** The Palmer penguins data

**Michel's comment:**

We illustrate the proposed diagnostic tools etc. using the Palmer penguins data. Furthermore,... (But then we need to more clear describe what we mean by: placing insider clustering and classification pipelines...

### 5. Highlight

> **Marked text:** A World Development Indicators snapshot then places the distance inside clustering and classification pipelines, where the preprocessing is refit within each resample.

**Michel's comment:**

Must make more clear what we mean

### 6. Text

**Michel's comment:**

I think we need to think about setup structure to keep better track of all functionality. 
1) Intro
2) Related/other software
3) The unified framework
4) Distance construction: mdist: method_cat, method_num + commensurability option. This allows for a huge list of (additive) distances; existing and non-existing. The commensurability option makes them unbiased. We also offer presets for a selection of additive ones and some non-additive ones from the literature. We should also probably say something about association-based and not.
5) Distance diagnostics: LOVO and comparison tools. But: all unsupervised. (But: association-based can/should be there).
6) Supervised settings: First review basics (we have y). Then: Make clear how our association-based distance offer a new set of distances. Then explain importance of careful (correct) model training and assessment. Explain the importance and peculiarities of data-splitting, and how this should be done to avoid leakage. Explain what methods we include. (KNN, Supervised clustering?) Show how it works. 
7) Conclusion

## Page 2

### 1. Highlight

> **Marked text:** Numerical variables face the mirror image of this problem through scaling.

**Michel's comment:**

why is that mirror image? Maybe we can delete without it affecting the rest

### 2. Highlight

> **Marked text:** of a matching categorical

**Action:** Marked without a written comment.

### 3. Highlight

> **Marked text:** to favor categorical contributions

**Action:** Marked without a written comment.

### 4. Highlight

> **Marked text:** to correct one variable at a time.

**Michel's comment:**

correct for.

### 5. Highlight

> **Marked text:** the resulting mean pairwise distance is typically small relative to that of a matching categorical variable, which is why Gower tends to favor categorical contributions (Hennig and Liao 2013; Foss et al. 2016).

**Michel's comment:**

most pairwise distances tend to be much smaller than one, which is why categorical variables tend to dominate in Gower's distance (HL, Foss)

### 6. Highlight

> **Marked text:** Skewness compounds the effect:

**Michel's comment:**

Skewed distributions for the numerical variables amplify this effect:

### 7. Highlight

> **Marked text:** heavy-tailed

**Michel's comment:**

Is it heavy tailed or long tailed? (I would think that heavy means many observations close to the tail. This may be less problematic (or have reverse effects) than skewness with a long but narrow tail

### 8. StrikeOut

> **Marked text:** while other scaling choices can amplify or attenuate the same variable unpredictably.

**Action:** Delete the marked text.

### 9. Highlight

> **Marked text:** distortions

**Michel's comment:**

distortions to the overall distance

### 10. Highlight

> **Marked text:** Because these distortions follow from measurement type and scale rather than from information content, they are hard to anticipate and to correct one variable at a time.

**Action:** Marked without a written comment.

### 11. Highlight

> **Marked text:** are hard to anticipate and to correct one variable at a time.

**Michel's comment:**

they are undesirable factors in the calculation of an overall distance measure. van de Velden et al describe these problems and formulate a general mixed variable distance that can be used to construct mixed-variable distances that do not suffer from these problems.

### 12. Highlight

> **Marked text:** van de Velden et al. (2026) formalize unbiased mixed-variable distances through two properties: multivariate additivity and commensurability.

**Michel's comment:**

Two properties (or factors?) are essential for the unbiased mixed-variable distances proposed in van de Velden et al: multivar additivity and commensurability....

### 13. Highlight

> **Marked text:** requires

**Michel's comment:**

entails that? specifies that?

### 14. StrikeOut

> **Marked text:** so that each variable enters the total additively.

**Action:** Delete the marked text.

### 15. Highlight

> **Marked text:** Importantly, a variable-specific contribution need not depend on that variable alone; it may be a function of the full same-type block, which is what later permits association-based constructions within the additive form.

**Michel's comment:**

I am not sure if we need this here. This is important but to emphasize it here, without more context, is probably just confusing.

### 16. Highlight

> **Marked text:** the

**Michel's comment:**

individual variable

### 17. StrikeOut

> **Marked text:** drawn from a more dispersed dis tribution,

**Action:** Delete the marked text.

### 18. Highlight

> **Marked text:** are commensurable

**Michel's comment:**

are considered to be. Or: are defined to be commensurable

### 19. Highlight

> **Marked text:** so that

**Michel's comment:**

. Hence, no variable contributes more (or less) to the distance based on its

### 20. Highlight

> **Marked text:** For

**Michel's comment:**

In

### 21. Highlight

> **Marked text:** this is

**Michel's comment:**

commensurability can be

### 22. StrikeOut

> **Marked text:** - while still allowing a variable to be up- or down-weighted deliberately when subject-matter knowledge justifies it.

**Action:** Delete the marked text.

### 23. Highlight

> **Marked text:** this framework.

**Michel's comment:**

unbiased mixed-variable distance framework. 
(should we also add: as well as other mixed-variable distances that do not necessarily fall into it.)

### 24. Highlight

> **Marked text:** and the discipline required to use an estimated distance inside a resampled workflow.

**Michel's comment:**

this is not clear. The dicipline? I guess that this concerns the "tidyverse" stuff. But, we have to write also for those not familiar with the jargon. I guess this is related to the last part in the abstract. We need to formulate this contribution clearly: 
As well as how distances can be implemented in a workflow allowing for a structured comparison of in- and out-of sample performances.

### 25. Highlight

> **Marked text:** a

**Michel's comment:**

that:
1. Can be used to define and calculate additive mixed variable distances comprising of numerical distances (allowing for various preprocessing options), and categorical variable distances (allowing for different category dissimilarity definitions).. 
2. Can be used to construct unbiased distances by enforcing commensurability through variable-wise weighting
3. Allows for the construction of distances incorporating association between variables as well as association with a response variable. 
4. Retains the preprocessing steps (such as preprocessing steps, category dissimilarities and variable specific weights) so that the same distance can be applied to "new observations" without requiring re-estimation. 

I think I combined the old 1 and 3. The fact that we have variable specific dissimilarities is important for lovo. But, the output of mdist is the aggregation so I am not sure whether mentioning one part of the calculation explicitly here....I explicitly added point 3. But without the terminology (association-aware and response-aware). I think these can come later. I do think it is important to mention association-based here as it is something not present in other packages.

## Page 3

### 1. Highlight

> **Marked text:** Association-aware distances

**Michel's comment:**

But I think this is the first time we mention this. We had (or have) something about using other variables for individual variable distance. But I am not sure if this is the place, In fact, we may want to dedicate a bit more space, or at least a clearly separated part/section to this.

### 2. Text

**Michel's comment:**

I think that perhaps this complete paragraph should be moved elsewhere. I do not think it fits the introduction. In fact, perhaps elsewhere we need a dedicated section on association-based and response-aware. 
On the other hand, if we want to say something about the package being able to deal with this (and in a sense, the third (fourth previously) point about being able to re-use on new observations is a bit linked....Still: if you think about it, the association-aware and response aware part is in the general formulation (perhaps the response part a bit less clear). So: Not sure if we need to make a big deal of this now. (Instead of simply having it implicitly in the sentence on pre-processing)....

### 3. Highlight

> **Marked text:** Scale- and type-aware

**Michel's comment:**

new terminology

### 4. Highlight

> **Marked text:** pairwise benchmarking of candidate distances,

**Michel's comment:**

What is pairwise benchmarking? And: don;t we in fact have comparisons methods/visualizations that allow for much more than pairs?
How about: 
lovo diganostics as well as several evaluation and comparion tools. In particular, integration with tidymodels (ref) is facilitated allowing appraisal of distances based on subsequent distance-based methods like cluster analysis and k-nn classification.

### 5. Highlight

> **Marked text:** Several R packages provide distances for homogeneous numerical, categorical, or binary data.

**Michel's comment:**

Should we create a dedicated section on "other software"? Seems like this could work like that. Perhaps the only problem is that this is probably not complete. (And only R). But, we can probably mention this: Not exhaustive and only R....I mean, shouldn't really be a problem I think.

### 6. Highlight

> **Marked text:** clustering, ordination, and nearest-neighbour met hods.

**Michel's comment:**

distance-based data analysis methods

### 7. StrikeOut

> **Marked text:** manydist builds on rather than replaces these tools, computing comparable per-block dissimilarities and adding a commensurability layer they do not provide.

**Action:** Delete the marked text.

## Page 4

### 1. Highlight

> **Marked text:** matrices,

**Michel's comment:**

data?

### 2. Highlight

> **Marked text:** extensible

**Michel's comment:**

?

### 3. Highlight

> **Marked text:** provides frequency-weighted similarities.

**Michel's comment:**

What is this? (Also some of the previous descriptions (sparse transaction data?) are a bit obscure to me. Wasn't nomclust the package by Sulc? Doesn't that have all the Boriah distances in there?

### 4. Highlight

> **Marked text:** cluster::daisy(),

**Michel's comment:**

but daisy has also all regular numerical distances. Shouldn't we make that clear? (And: I assume that this is unlike, for example, the gower package that probably only has gower. Or not?

### 5. Highlight

> **Marked text:** of the same underlying construction

**Michel's comment:**

What do we mean?

### 6. Highlight

> **Marked text:** and its tendency to give relatively greater influence to categorical variables when many categorical attributes are present or when standardized numerical differences are comparatively small.

**Michel's comment:**

This has nothing to do with the packages. This is simply a property of gower, as we explained.

### 7. Highlight

> **Marked text:** for mixed data

**Michel's comment:**

I guess we can resolve my point after this (about daisy) by stating that packages implementing/allowing/providing options for mixed-variable distances are summarized in T2 and fall into three families.

### 8. Highlight

> **Marked text:** controls the relative influence of numerical and categorical information through explicit weighting.

**Michel's comment:**

concerns mixed distance alternatives in which the numerical and categorical distance parts are aggregated by assigning weights for these two parts.

### 9. Highlight

> **Marked text:** Gower, Podani, Wishart, Huang, Harikumar, and Ahmad

**Michel's comment:**

we probably need references. (For Podani, Wishart, Harikumar and Ahmad we do not have any yet I think)

### 10. Highlight

> **Marked text:** learns an adaptive

**Michel's comment:**

concerns mixed-data distances that estimate distances based on the data.

### 11. StrikeOut

> **Marked text:** automatically,

**Action:** Delete the marked text.

### 12. Highlight

> **Marked text:** Related to

**Michel's comment:**

In?

### 13. Highlight

> **Marked text:** a weighting

**Michel's comment:**

weights

### 14. Highlight

> **Marked text:** Both address the balance between variable types as part of a clustering procedure,

**Michel's comment:**

So the clustering part is also in kamila. Like cluspcamix, you cannot separate the two? So: we should probably add: 
as part of a clustering procedure and do not explicitly yield a separate distance matrix. 
(Or actually, use that to replace the whereas the .... sentence.)

### 15. Highlight

> **Marked text:** it differs from the commensurable distances developed here, which explicitly equalize the expected contributions of numerical and categorical variables.

**Michel's comment:**

it does not necessarily lead to commensurable distances.

### 16. StrikeOut

> **Marked text:** likewise

**Action:** Delete the marked text.

### 17. Highlight

> **Marked text:** the estimated quantities must be refit within a resampling scheme, exactly as recipes refits step_pca() or step_normalize () on each analysis set.

**Michel's comment:**

For new data, either new estimation is required or the estimated quantities must be stored and used for this.

### 18. Highlight

> **Marked text:** To our knowledge no other mixed-data distance implementation exposes a fitted representation that can be reapplied to new observations without re-estimation; Section 5 takes this up.

**Michel's comment:**

I think this is too cryptic. We do not need it, or we need to explain more clearly.

### 19. Highlight

> **Marked text:** Gower dissimilarity for numerical, nominal, ordinal and binary variables.

**Michel's comment:**

But daisy has more than gower, not?

### 20. Highlight

> **Marked text:** variants

**Michel's comment:**

in what sense? weights?

### 21. Highlight

> **Marked text:** share a consequence

**Michel's comment:**

But do the others not have this also? For example Ahmad and Dey? But also clustrd and prototypes? (If the weight gamma is based on data, certainly. Not?)

### 22. Text

**Michel's comment:**

It seems we do not visualize the three types mentioned in the text in the table. Why not?

## Page 5

### 1. Highlight

> **Marked text:** The remainder

**Michel's comment:**

I would put this before the "Several R packages" part (that I think would make a good other software section)

### 2. StrikeOut

> **Marked text:** These tools are useful, but their variable contributions can depend on scale, prevalence, imbalance, and measurement type.

**Action:** Delete the marked text.

### 3. Highlight

> **Marked text:** them

**Michel's comment:**

the existing packages by offering

### 4. Highlight

> **Marked text:** association-aware and response-aware

**Michel's comment:**

Then this needs to be explained. (I suggested to move from introduction). But, removing it is not a big problem I think. 
In fact, perhaps we can point to this after "...interfaces." 
e.g.:
Furthermore, in line with the general framework from vandevelden etal, mdist also allows for construction of distances exploiting association between variables, as well as association with a response variable.

### 5. Highlight

> **Marked text:** develops

**Michel's comment:**

briefly summarizes/revisits
Actually. I prefer: 
In Section 2, we 
(and the same for other Sections...)
Also: Section 2 provides an overview of existing packages....

### 6. Highlight

> **Marked text:** fitted preprocessing.

**Michel's comment:**

do we need fitted? preprocessing estimates

### 7. Highlight

> **Marked text:** pairwise benchmarking

**Action:** Marked without a written comment.

### 8. Highlight

> **Marked text:** embeds distances

**Michel's comment:**

In Section 5, we show how distances can be embedded in

### 9. Highlight

> **Marked text:** I

**Michel's comment:**

I doubt that I used that for number of observations. Actually, a lot of the notation (the capital Qs etc. is not the same as in JCGS). To be honest, unless a specific reason to divert, I would prefer the n's, ps and qs. (But I do not care very much; still, if we change better not to say we follow the notation)

## Page 6

### 1. Highlight

> **Marked text:** the contribution

**Michel's comment:**

distance contribution?

### 2. Highlight

> **Marked text:** Each component may depend on all elements of the two observations; the index identifies the associated contribution rather than restricting it to a single coordinate. This permits additive constructions that incorporate dependencies between variables. If every dj is a metric on XQ, their sum is also a metric.

**Michel's comment:**

I do think a lot seems very similar to what we have in JCGS. To some degree ok because referring is also a bit annoying. But, perhaps some of the details (this property of a distance being metric) is not needed here. On the other hand, perhaps it is good to keep here. (I see we do not proof, so we could refer to the paper for proof(s).

### 3. Highlight

> **Marked text:** Scale- and type-aware distances

**Michel's comment:**

Not sure about this title

## Page 7

### 1. Highlight

> **Marked text:** The formulation also permits constrained weights, such as common weights within numerical or categorical groups, when only within-type commensurability is required. Under such constrained weights the per-variable cancellation described below no longer applies, and the numerical scaling choice again affects the result.

**Michel's comment:**

But: since we do not do anything with this in the package (or do we?) we can leave it perhaps

### 2. StrikeOut

> **Marked text:** beyond commensurability

**Action:** Delete the marked text.

### 3. Highlight

> **Marked text:** Association-aware distances:

**Michel's comment:**

I find this a very interesting section. It made me rethink about some of the things we do. In fact, it made me seriously doubt the usefulness of the PCA-scaling option in combination with commensurability. 
The PCA rotation leads to independent dimensions. (By construction). But: The dimensions are also by definition not equally "important". The first cover max variance, etc. If after rotating, we make commensurable, we increase the weights of the lower (not important in terms of variance) dimensions by a potentially large factor...
The reason for doing PCA rotation is that this mimics Mahalanobis. Mahalabonis is insensitive to rotation and as such dimensions (variables) are not independent. Euclidean on Mahalanobis is the same as Euclidean on PCA. But: Manhattan differs. So: we used the PCA orientation as it provides "independent" directions. However, following this by the commensurability step leads to the "problem" I just mentioned: Unimportant dimensions get high weights, important (in terms of variance) dimensions get downweighted...
I guess using fewer dimensions would be a good way to remedy...(But: this is all not for now).

But: What we need to do is rewrite/reorganize this section. 
Start by observing that the general formulation allows for inclusion of association as the d's are a function of all x's (of the same type). Then, separate numerical and categorical. For numerical, explain that one way, related to Mahalanobis distances, involves PCA rotation. Perhaps explain why this could be a good idea (using the Palmer plots). But: perhaps we need to be a bit careful as I am not so sure about this being such a great idea....On the other hand: I think we had the option of doing PCA scores using k-dimensions? Do we want to go there?

Next, explain that for categorical variables, the association can be exploited by using it to define category dissimilarities. This is obviously very different from the PCA/Mahalanobis idea. It is probably a bit easier to justify. (But: Perhaps that is not really needed here. We just describe what is possible). 

Finally, we should introduce/formalize the response-awareness concept. In particular: this is relevant for the categorical variables I think. (We do not have it for numerical yet. Do we?). It follows quite naturally when describing association-based for categorical. In fact, we can/should refer to the Pattern Recogniton paper where we propose this. And already "tested" it.

### 4. StrikeOut

> **Marked text:** by itself

**Action:** Delete the marked text.

### 5. Highlight

> **Marked text:** When several variables encode the same underlying source of variation, that variation contributes repeatedly to the overall distance, even after the variables have been standardized or made commensurable. An association-based construction instead allows a contribution to depend on the complete block of

**Michel's comment:**

Maybe we should start here. And then, before An association-based construction. Simply note that definitions for variable specific distances allow for contributions of other variables. (I think we had a remark about that in the introduction. Here would be a good place to have that): Note that, in the general formulation variable specific distances are a function of the values of other variables (of the same type). This means that we can use dependencies/association between variables when defining variable specific distances. manydist allows for this both for numerical as categorical variables.

## Page 9

### 1. Highlight

> **Marked text:** whitening,

**Action:** Marked without a written comment.

### 2. Highlight

> **Marked text:** Figure 1:

**Michel's comment:**

I think we do not refer to this Figure....

### 3. Highlight

> **Marked text:** whitening

**Michel's comment:**

needs to be define. Or avoided.

### 4. Highlight

> **Marked text:** to remove linear redundancy.

**Michel's comment:**

Indeed, if you put it like that it sort of becomes obvious that doing PCA (whitening) and then making commensurable is maybe not a good idea....We don't remove redundancies and in fact may "explode" them...

## Page 10

### 1. Highlight

> **Marked text:** Common transformations include: • standard-deviation scaling: (x - x)/sx, which uses all observations rather than only the extremes; • range scaling: (x - as in Gower's distance, which preserves the Xmin)/(xmax - Xmin), distributional shape but, for the non-commensurable presets, is sensitive to outliers through the extremes; • robust scaling: (x - which replaces the mean and range with their XMed)/(xQ 3 - XQ 1 ), robust counterparts and so resists outliers in those same presets.

**Michel's comment:**

This should all come before the association part. Not?

## Page 11

### 1. Highlight

> **Marked text:** Indicator-based dissimilarities treat a binary indicator representation as numerical. Without scaling, Manhattan distance gives values zero or two. This fixed spacing ignores both the number of categories and their prevalence. Hennig- Liao scaling (Hennig and Liao 2013), standard-deviation scaling, and direct category-dissimilarity scaling recalibrate these indicator distances so that they become comparable, before the commensurability weights equalize the expected contributions.

**Michel's comment:**

I think we do not need this

### 2. Highlight

> **Marked text:** Association-based dissimilarities incorporate relationships among categorical variables. Let = Zk)-1ZI R klk' (ZI Zk, be the conditional distribution of variable k' given the categories of variable k. Dissimilarities can compare rows and using total variation (Ahmad and Dey r a r b 2007) or a symmetric Kullback- Leibler construction (Le and Ho 2005). Under independence every category of k induces the same condit ional distribution, the rows of R klk' coincide, and all category dissimilarities vanish; the dissimilarities grow as the categories become more distinct in what they imply about the other variables. This mirrors, for categorical variables, the way principal-component preprocessing removes linear redundancy among the numerical variables. Because the construction conditions on the remaining categorical variables, it requires at least two of them; with a single categorical variable the fitted specification falls back to simple matching.

**Michel's comment:**

I think this also follows to closely our JCGS paper. We should probably rethink a bit what is needed. And what not.

### 3. Highlight

> **Marked text:** The same construction accommodates a response. When an outcome is supplied, the conditional distributions are taken with respect to it, so R records the kly distribution of the outcome within each category of variable and two categories k, are dissimilar to the extent that they imply different outcome distributions. This is what the argument of O supplies, and it is the only route by response mdist which the outcome enters a distance in manydist. It also means a response-aware categorical dissimilarity is available with a single categorical predictor, where the predictor-only construction would fall back to matching. The numerical block has no counterpart to this. Principal-component preprocessing uses the covariance of the predictors alone, so the outcome does not enter the numerical contributions under any preset.

**Michel's comment:**

Needs to be rephrase/organized

## Page 12

### 1. StrikeOut

**Michel's comment:**

,

### 2. Highlight

> **Marked text:** A preset identifies a complete distance specification.

**Michel's comment:**

Before the preset, we should explain basic components/workgins fo mdist: what goes in (without the preset). Than later the presets are easily explained as combinations of the inputs. More advanced options (response) can also be added. (Or we immediately mention, but revert).

So: mdist(data, method_num, method_cat, commensurable). Where data is the data, method_num the preprocessing/scaling for numerical, method_cat the distance/dissimilarity definition for categorical variables and commensurable whether distance should be made commensurable or not. 

We also probably need to mention defaults. And all options (or: put that in a table).

Then, list the presets and explain. 
Then: describe the resulting MDist object.

### 3. Highlight

> **Marked text:** u_dep

**Action:** Marked without a written comment.

### 4. Highlight

> **Marked text:** carries

**Michel's comment:**

contains

### 5. Highlight

> **Marked text:** consumed by downstream

**Michel's comment:**

used or necessary for ?

### 6. Highlight

> **Marked text:** hclust () and cmdscale ().

**Michel's comment:**

bit random subset of methods....
Actually: why not stop at returns an R dist object. (Also not call it plain).

## Page 13

### 1. Highlight

> **Marked text:** in the gower, euclidean, and hl presets.

**Michel's comment:**

We need a list of presets (and a bit on why)

## Page 14

### 1. Highlight

> **Marked text:** The custom preset exposes the three choices directly through the method_num, method_cat, and commensurable arguments.

**Michel's comment:**

Yes, move forward. (And: custom is in fact not really needed, is it as it defaults to this. Right?)

### 2. Highlight

> **Marked text:** Table 3 summarizes the presets

**Michel's comment:**

We should state somewhere why some options are in here. (In particular, how some presets concern versions that do not fit the framework and hence cannot be replicated using custom)

## Page 15

### 1. Highlight

> **Marked text:** A distance is an intermediate data representation, so its sensitivity should be examined before it is supplied to a learning algorithm.

**Michel's comment:**

Bit strong statement....Would keep it more casual/objective and simply state that it could be interesting/insightful to inspect properties of certain distances. (Especially as there are so many options available)
Actually: We in fact offer diagnostics based on what happens in the next step (cluster analysis and MDS). 

I think this is an important feature of our package. So we should clearly present (and sell this). The distinction mentioned here (before analysis, after analysis: I think you already referred to this as "downstream" tasks) is good. Let's keep this. Introduce by observing this and explain how for both we propose LOVO measures. Then, we first introduce the direct (distance-based) ones (MAD, and relative MAD (actually: I am getting confused, do we have both? For the JCGS paper we did. Here I do not see this clearly). Then we do "downstream". And for that is is pretty important that we make clear what options we have and why. (The cluster ones are, I think, easier to describe and justify. The MDS one is smaller, more peculiar/specific. We need to therefore be clear about what we do. What it does. What it "assumes" (an underlying metric MDS solution?). 
After that, we can go to the comparison functions. Because: after explaining the lovo functionality, it should be easy to mention that one can do this to compare different distances. So: we explain how one can do this in practice (by creating the list of methods defined by the method_cat, method_num arguments and (possibly), commensurability and preset ones. The response aware should probably not be treated here yet. If we have a response, we enter a "new world" where performance can be assessed using the y variable. This is the next section I think. So: here we first explain how lovo_mdist wokrs. (One dist only). Than, compare_lovo_mdist. Then, we explain how we have a function to give a huge list with candidates (?)

### 2. Highlight

> **Marked text:** one predictor

**Michel's comment:**

We do not have predictors....Just a variable

### 3. Highlight

> **Marked text:** The diagnostics

**Michel's comment:**

I guess this follows on the "quantifying resulting change" part. I think we need to connect that here. (What you now call diagnostics concern quantifications of such changes)

### 4. Highlight

> **Marked text:** relative_distance

**Michel's comment:**

I think we probably need a list/itimized list with options and explanations. And probably some context as well (as to why use one or the other etc.)

### 5. Highlight

> **Marked text:** the alienation between the full and reduced configurations.

**Michel's comment:**

What configuration? Did we ever explain how this relates to/uses MDS?

### 6. Text

> **Marked text:** in the [O,

**Michel's comment:**

@

### 7. Highlight

> **Marked text:** measuring the configuration variance left unexplained

**Michel's comment:**

Is that correct? Is that what alienation coefficient really is?

### 8. Highlight

> **Marked text:** Reporting both is deliberate: a variable can shift many pairwise distances while barely moving the low-dimensional configuration those distances induce, or move few distances in a way that nonetheless reshapes the configuration. The two metrics can therefore disagree, and each answers a different question about the variable's role. Larger values indicate a larger contribution for both quantities.

**Michel's comment:**

Yes, but I wonder if we should do the two immediately together. The distance one is very simple and intuitive. The AC one involves the MDS step (which we need to explain). It is less obvious that this is always appropriate. After all, it says something about the influence on MDS solution. If MDS is not appropriate, what does it say? Also: doesn't it fit better with the cluster analysis measures? (There too: influence on a different method applied after. So: If cluster analysis is not appropriate, what do the values say?)

### 9. Highlight

> **Marked text:** benchmarking

**Action:** Marked without a written comment.

## Page 16

### 1. Highlight

> **Marked text:** Figure 2: Relative

**Michel's comment:**

Do we refer to this? I don't think so. I think these plots are really nice. We should show them. And explain why they are nice and useful.

### 2. Highlight

> **Marked text:** lovo_gower$autoplot( = metric "relative_distance", + reorder= + TRUE +) + = ggplot2:: guides ( colour "none") +

**Michel's comment:**

wouldn't:
autoplot(lovo_gower) work? Wouldn't that be nicer? (Or: plot(lovo_gower) ?)

### 3. Highlight

> **Marked text:** The two metrics disagree

**Action:** Marked without a written comment.

## Page 17

### 1. Highlight

> **Marked text:** Comparing distance specifications

**Michel's comment:**

This includes the alienation results. Also: confusing with respect for what is to come. The way I now (after reading several times) understand this is that this concerns lovo comparison. The next section (benchmarking) concerns overall comparisons. But: in unsupervised setting. (Supervised is Secition 5?) We have to be much more explicit structured concerning this.

## Page 18

### 1. Highlight

> **Marked text:** Figure 3: Stability

**Michel's comment:**

Also not referenced

### 2. Highlight

> **Marked text:** and can disagree,

**Michel's comment:**

I am not sure I agree with your use of the term disagree....

### 3. Highlight

> **Marked text:** The pairwise diagnostics

**Michel's comment:**

OK. So you compare two distance matrices. What does the pairwise add/avoid/help in understanding this? And: Shouldn't we have this earlier? (Before the previous section with relative distances?)

### 4. Highlight

> **Marked text:** changes in magnitude from changes in geometry.

**Michel's comment:**

Tricky....In any case, let's avoid. I think the distinction is between magnitude of distances (relative and mad!) and the effect on subsequent analyses. Let's make that clear. And: Let's start with the magnitudes.

### 5. Highlight

> **Marked text:** every pair of successful distances.

**Michel's comment:**

succesful?

### 6. Highlight

> **Marked text:** The MDS congruence coefficient and its corresponding alienation coefficient compare the configurations induced by the two distances.

**Michel's comment:**

Yes, and specifically for MDS (metric? non-metric? Which did we use?)

### 7. Highlight

> **Marked text:** If is specified, 0 cluster_ k benchmark_mdist also applies each requested clustering method to every distance and computes pairwise ARls. With the default cluster_k = NULL, no clustering is performed.

**Michel's comment:**

Too many things happen at the same time

### 8. Highlight

> **Marked text:** Benchmarking

**Michel's comment:**

Not sure whether this is the correct word. We do offer the function benchmark_mdist to do comparisons

### 9. Highlight

> **Marked text:** evaluates an explicit table of distance spec ifications.

**Michel's comment:**

What is that? (I mean: I know were we are going, but if you don't, this is not so straightforward)

## Page 19

### 1. Highlight

> **Marked text:** benchmark_pairs <- benchmark_comparisons(distance_benchmark)

**Action:** Marked without a written comment.

### 2. Highlight

> **Marked text:** candidate_specs <- all_dist_method_specs( mode = "presets_only", + preset = c("gower", "u_indep", "u_dep") + I> +) dplyr::mutate( + label= unname(candidate_labels[preset]) + + )

**Michel's comment:**

not a very simple example/illustration of what our package offers and can do. Too many other stuff (dplyr, mutate) is going on.

### 3. Highlight

> **Marked text:** under pairs of distances;

**Action:** Marked without a written comment.

### 4. Highlight

> **Marked text:** low-dimensional geometry,

**Action:** Marked without a written comment.

### 5. Highlight

> **Marked text:** determines

**Michel's comment:**

measures?

### 6. Highlight

> **Marked text:** whether those changes alter the partition.

**Michel's comment:**

well, particularly how much they changed

### 7. Highlight

> **Marked text:** are therefore complementary rather than alternative rankings.

**Michel's comment:**

you mean complementary to each other. They really measure different things. Perhaps we should be much stronger about that. (And have users make a choice before just running what ever)

### 8. Highlight

> **Marked text:** renders any pairwise diagnostic autoplot () as an annotated triangular heatmap.

**Michel's comment:**

where? how?

### 9. Highlight

> **Marked text:** An error in one specification is captured

**Michel's comment:**

but what is an error? How can there be errors in one specification and not in another?

### 10. Highlight

> **Marked text:** a benchmark remains a sensitivity analysis rather than a model-selection rule.

**Michel's comment:**

??

### 11. Highlight

> **Marked text:** Gower

**Michel's comment:**

distance?

### 12. Highlight

> **Marked text:** has more than twice the alienation

**Michel's comment:**

what does that mean? twice as dissimilar? twice as bad MDS fit?

### 13. Highlight

> **Marked text:** a much lower PAM ARI.

**Michel's comment:**

I thought the ARI is between the methods. So: this just means they are dissimilar. (You get different cluster solutions).

### 14. Highlight

> **Marked text:** This

**Michel's comment:**

what exactly?

## Page 20

### 1. Highlight

> **Marked text:** distance, and so the same gap from the Gower

**Michel's comment:**

?

### 2. Highlight

> **Marked text:** Taken together, the diagnostics indicate how much of an analysis rests on the choice of distance. When the measures converge - a variable is influential, or two distances agree - the downstream conclusion is robust to that choice. When they disagree, the distance is not a neutral preprocessing step but a modeling decision in its own right, one that should be carried forward and evaluated alongside the learning algorithm.

**Michel's comment:**

Not sure if we need to go into this. We are just presenting the package. Some guidance is good. But, we can be cautious. Especially as the downstream comparison is obviously not complete and multi-faceted (as we do mention ourselves).

### 3. Highlight

> **Marked text:** Distance-based learning pipelines

**Michel's comment:**

This is an important contribution of our package. We need to carefully introduce the problem. Explain when/why this is important. (We made some remarks about this already, best to reiterate here more clearly)

### 4. Highlight

> **Marked text:** Refitting the distance within resamples

**Michel's comment:**

To fast, to specific/technical

### 5. Highlight

> **Marked text:** recipe.

**Michel's comment:**

for me a recipe is really food related. Using it here requires an introduction

### 6. Highlight

> **Marked text:** The fitted step stores the training data and preprocessing parameters, so the distance is re-estimated inside each resample and reused when genuinely new observations are processed. T his re-estimation is the point: the commensurability weights, any principal-component rotation, and the category dissimilarities are quantities estimated from data, not fixed constants. Computing one distance matrix on the full data and then resampling would let these estimates see the assessment observations, so the step instead refits them on each analysis set alone. This is the pipeline-level counterpart of the square-versus-rectangular distinction of Section 3.

**Michel's comment:**

I think we can be more concise. It is important to explain the point. But I have the impression that we are repeating ourselves. By doing so, using slightly different words, we may in fact not make things more clear but confuse the reader.

I do not understnd the last sentence (This is the pipeline-level?)

### 7. Highlight

> **Marked text:** Its

**Michel's comment:**

This is about step_mdist()...very far away.

### 8. Highlight

> **Marked text:** required by the downstream task: for clustering and "pai rwise" for prediction. "distance_to_training"

**Michel's comment:**

Isn't this the first time we do prediction? Shouldn't we introduce this? Separate it from the other downstream task? (Clustering. But, can't we also have MDS downstream? Didn't we more or less consider that before?

### 9. Highlight

> **Marked text:** whereas nearest-neighbour prediction

**Michel's comment:**

only mentioned in introduction

## Page 21

### 1. Highlight

> **Marked text:** a fixed snapshot

**Action:** Marked without a written comment.
