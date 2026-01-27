## {.center}

[SC for mixed data?]{style="color:indianred;font-size:2em;"}

. . .

the distance measure definition for mixed data is key








##

[**desirable properties**]{style="color:indianred"} ^[`r Citet(bib=wp_bologna_26_bib,"vdv_jcgs")`]

::: {.callout-tip appearence="simple" icon="false"}
## Multivariate Additivity

Let $\mathbf{x}_i=\left(x_{i1}, \dots, x_{iQ}\right)$ denote a $Q-$dimensional vector. A distance function $d\left(\mathbf{x}_i,\mathbf{x}_\ell\right)$ between observations $i$ and $\ell$ is multivariate [additive]{style="color:indianred"} if

$$
  d\left(\mathbf{x}_i,\mathbf{x}_\ell\right)=\sum_{j=1}^{Q} d_j\left(\mathbf{x}_i,\mathbf{x}_\ell\right),
$$

  where $d_j\left(\mathbf{x}_i,\mathbf{x}_\ell\right)$ denotes the $j-$th variable specific distance.

:::

  - Manhattan distance satisfies the additivity property; the Euclidean distance does not





. . .

If additivity holds, by-variable distances are added together: they should be on equivalent scales


::: {.callout-note appearence="simple" icon="false"}
## Commensurability

Let ${\boldsymbol X}_i =\left(X_{i1}, \dots, X_{iQ}\right)$ denote a $Q-$dimensional random variable corresponding to an observation $i$. Furthermore, let $d_{j}$ denote the distance function corresponding to the $j-$th variable. We have [commensurability]{style="color:indianred"} if,
for all $j$, and $i \neq \ell$,

$$
  E[d_{j}({ X}_{ij}, {X}_{\ell j})] = c,
$$

  where $c$ is some constant.
:::



  ## {.center}


  If the multivariate distance function $d(\cdot,\cdot)$ satisfies [additivity]{style="color:forestgreen"} and
[commensurability]{style="color:dodgerblue"},
*ad hoc* distance functions can be used for each variable and then aggregated:

  :::{.columns}

::::{.column width=10%}
&nbsp;

then
::::

  ::::{.column width=90%}
>
  >   one can pick the appropriate $d_{j}(\cdot,\cdot)$, given the nature of $X_{j}$
  >
  >   - well suited in the mixed data case
>
  ::::

  :::

  ##

  [**mixed-data setup**]{style="color:indianred"}

::: {.callout-tip appearence="simple" icon="false"}
## a [mixed]{style="color: dodgerblue"} data set

-   $I$ observations described by $Q$ variables,  $Q_{n}$ numerical  and $Q_{c}$ categorical

-   the $I\times Q$ data matrix ${\bf X}=\left[{\bf X}_{n},{\bf X}_{c}\right]$ is column-wise partitioned
:::


  A formulation for mixed distance between observations $i$ and $\ell$:

  \begin{eqnarray}\label{genmixeddist_formula}
d\left(\mathbf{x}_i,\mathbf{x}_\ell\right)&=& \sum_{j_n=1}^{Q_n} d_{j_n}\left(\mathbf{x}^n_i,\mathbf{x}^n_\ell\right)+ \sum_{j_c=1}^{Q_c} d_{j_c}\left(\mathbf{x}^c_i,\mathbf{x}^c_\ell\right)=\\
&=& \sum_{j_n=1}^{Q_n} w_{j_n} \delta^n_{j_n}\left(\mathbf{x}^n_i,\mathbf{x}^n_\ell\right)+ \sum_{j_c=1}^{Q_c} w_{j_c}\delta^c_{j_c}\left(\mathbf{x}^c_i,\mathbf{x}^c_\ell\right)
\end{eqnarray}

:::{.columns}
::::{.column width=45%}

:::{.callout-note appearence="simple" icon="false"}
## numeric case
- $\delta^n_{j_n}$ is a function quantifying the dissimilarity between observations on the $j_n-$th numerical variable

- $w_{j_n}$ is a  weight for the $j_n-$th variable.
:::
  ::::

  ::::{.column width=45%}
:::{.callout-important appearence="simple" icon="false"}
## categorical case
dissimilarity between the categories chosen by subjects $i$ and $\ell$ for categorical variable $j_c$

  - $w_{j_c}$ is a  weight for the $j_c-$th variable
:::
  ::::
  :::

  ## {.center}

  [association-based (AB) distances?]{style="color:indianred;font-size:2em;"}

. . .

not all differences weighs the same

##

::: {.callout-note appearence="simple" icon="false"}
##  distances for categorical data

[the general (delta) framework]{style="color:indianred"}^[`r Citet(bib=wp_bologna_26_bib,"vdv_pr")`]

Let ${\bf Z}=\left[{\bf Z}_{1},{\bf Z}_{2},\ldots,{\bf Z}_{Q_c}\right]$ be the one-hot encoding ${\bf X}_{c}$

  The pair-wise distances between categorical observations are given by

$${\bf D}_{c}={\bf Z}{\bf \Delta}{\bf Z}^{\sf T}=
  \left[\begin{array}{ccc} {\bf Z}_{1} & \dots & {\bf Z}_{Q_{c}} \end{array} \right]\left[\begin{array}{ccc}
                                                                                          {\bf\Delta}_1  & & \\
                                                                                          & \ddots &\\
                                                                                          & & {\bf\Delta}_{Q_{c}} \end{array} \right] \left[ \begin{array}{c}
                                                                                                                                             {\bf Z}_{1}^{\sf T}\\
                                                                                                                                             \vdots \\
                                                                                                                                             {\bf Z}_{Q_{c}}^{\sf T}
                                                                                                                                             \end{array} \right]$$

  -   the definition of ${\bf \Delta}$ determines the distance in use

- if the off-diagonal terms of the $\Delta_{j}$'s only depend on the $j^{th}$ variable, then ${\bf D}_{c}$ is [independence-based]{style="color: indianred"}

  - if the off-diagonal terms of the $\Delta_{j}$'s depend on the $j^{th}$ variable AND on the relation of variable $j$ with the other $Q_{c}-1$ variables, then ${\bf D}_{c}$ is [association-based]{style="color: indianred"}


:::





  ##

  [**$\Delta_{j}$'s for association-based distances**]{style="color:indianred"}

the matrix co-occurrence proportions is


$$
{\bf P} =\frac{1}{I} \begin{bmatrix}
{\bf Z}_{1}^{\sf T}{\bf Z}_{1} & {\bf Z}_{1}^{\sf T}{\bf Z}_{2}&\ldots &{\bf Z}_{1}^{\sf T}{\bf Z}_{Q_{c}}\\
\vdots                         & \ddots                        &\vdots & \vdots                       \\
% {\bf Z}_{2}^{\sf T}{\bf Z}_{1} & {\bf Z}_{2}^{\sf T}{\bf Z}_{2}&\ldots &{\bf Z}_{2}^{\sf T}{\bf Z}_{Q}\\
\vdots                         & \vdots                        &\ddots & \vdots                       \\
{\bf Z}_{Q_{c}}^{\sf T}{\bf Z}_{1} & {\bf Z}_{Q_{c}}^{\sf T}{\bf Z}_{2}&\ldots &{\bf Z}_{Q_{c}}^{\sf T}{\bf Z}_{Q_{c}}\\
\end{bmatrix}
$$



- ${\bf R} = {\bf P}_{d}^{-1}\left({\bf P}-{\bf P}_{d}\right)$, with ${\bf P}_{d}=diag({\bf P})$, is a block matrix such that

- the general off-diagonal block is ${\bf R}_{ij}$ ( $q_{i}\times q_{j}$ )

- the $a^{th}$ row of ${\bf R}_{ij}$, ${\bf r}^{ij}_{a}$, is the  conditional distribution of the $j^{th}$ variable, given the $a^{th}$ category of the $i^{th}$ variable

##

[**Total variation distance (TVD): pair-wise dissimilarity between categories**]{style="color:indianred"}^[`r Cite(bib=wp_bologna_26_bib,"le2005association")`]

consider any pair of categories $a$ and $b$ from the $i^{th}$ categorical variable,
their overall dissimilarity is $\delta^{i}(a,b)$, that is
$$\delta^{i}(a,b)=\sum_{i\neq j}^{q}w_{ij}\Phi^{ij}({\bf r}^{ij}_{a},{\bf r}^{ij}_{b})$$
 upon defining  $$\Phi^{ij}({\bf r}^{ij}_{a},{\bf r}^{ij}_{b})=\frac{1}{2}\sum_{\ell=1}^{q_{j}}|{\bf r}^{ij}_{\ell a}-{\bf r}^{ij}_{\ell b}|$$
then $\Phi^{ij}()$ corresponds the [total variation distance]{style="color:dodgerblue"} between two discrete probability distributions





## {.center}

[AB spectral clustering for mixed data]{style="color:indianred;font-size:2em;"}


<!-- ## {.center} -->

<!-- **pairwise object distances (TVD-based)** -->

<!-- - For the $i^{th}$ categorical variable, compute $\delta^{i}(a,b)$ for all the possible categories pairs, and store them in the $q_{i}\times q_{i}$ matrix ${\bf \Delta}_{i}$ -->

<!-- - define the block-diagonal matrix -->

<!-- $${\bf \Delta}=\begin{bmatrix} -->
<!-- {\bf \Delta}_{1}&             & & & \\ -->
<!--                 &{\bf \Delta}_{2} & & & \\ -->
<!--                 &                 & \ddots & &  \\ -->
<!--           &                 & & & {\bf \Delta}_{q}\\ -->
<!-- \end{bmatrix}$$ -->
<!-- - the pairwise object distances matrix is -->
<!-- $${\bf D}_{cat} = {\bf Z}{\bf \Delta}{\bf Z}^{\sf T}$$ -->

##
[*naive* spectral clustering for mixed data]{style="color:indianred"}

the [SC for mixed data]{style="color:dodgerblue"}  proposed by `r Citet(bib=wp_bologna_26_bib,"mbuga21")`
is the NSW procedure for SC  applied on the  convex combination of the numerical and categorical distances:

$${\bf D}^{*} =\alpha{\bf {D}}_{n}^{*}+(1-\alpha){\bf {D}}^{*}_{c}$$
where

-  for continuous variables: ${\bf D}^{*}_{n}$ is the Euclidean distance

-  for categorical variables: ${\bf D}^{*}_{c}$ is the Hamming distance

- $\alpha$  can be tuned, its default is [$\alpha=.5$]{style="color:dodgerblue"}


##
**from**  [SC for mixed data]{style="color:dodgerblue"} **to** [AB SC for mixed data]{style="color:indianred"}

building upon the proposal by `r Citet(bib=wp_bologna_26_bib,"mbuga21")`, we replace the independence-based ${\bf D}^{*}_{n}$ and ${\bf D}^{*}_{c}$ with an association-based version. In particular:

- for continuous variables:  ${\bf D}_{n}$ is the [*commensurable* whithened Manhattan]{style="color:forestgreen"} distance (L1 equivalent of Mahalanobis distance)

- for categorical variables: ${\bf D}_{c}$ is [*commensurable* total variation]{style="color:forestgreen"} distance.

Since  additivity and commensurability hold, the considered general distance becomes:

$${\bf D} ={\bf {D}}_{n}+{\bf {D}}_{c}$$

<!-- ## -->

<!-- #### The aim: extend SC to mixed data via association-based distances -->

<!-- ::: {.callout-note  appearence="simple" icon=false} -->
<!-- ## association-based for categorical data: total variation distance (TVD) `r Citet(bib=wp_bologna_26_bib,"le2005association", before="see, e.g. , ")` -->

<!-- the weight assigned to the difference between a pair of categories [$a$]{style="color: dodgerblue"} and [$b$]{style="color: dodgerblue"} depends on the distributions of the other attributes conditional to $a$ and $b$: the higher the similarity, the lower the weight -->

<!-- $$\Phi^{ij}({\bf r}^{ij}_{a},{\bf r}^{ij}_{b})=\frac{1}{2}\sum_{\ell=1}^{q_{j}}|{\bf r}^{ij}_{\ell a}-{\bf r}^{ij}_{\ell b}|$$  -->

<!-- - ${\bf r}^{ij}_{a}$ is the $a^{th}$ row of ${\bf R}_{ij}$, $q_{i}\times q_{j}$ general off-diagonal block of ${\bf R}$ -->
<!-- - ${\bf r}^{ij}_{a}$ is the conditional distribution of the $j^{th}$ variable, given the $a^{th}$ category of the $i^{th}$ variable -->

<!-- - ${\bf R} = {\bf P}_{d}^{-1}\left({\bf P}-{\bf P}_{d}\right)$  -->
<!--   - ${\bf P}_{d}=diag({\bf P})$  -->


<!-- ::: -->



<!-- # Methodology and main results -->
