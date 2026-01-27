## {.center}

[**cluster analysis** ]{style="color:indianred"}

the general aim is to find homogeneous groups of observations


- observations from the same group are similar to each other

- observations from different groups are not similar to each other

- multiple clustering approaches exist

```{r,fig.align='center', out.width="40%", include=FALSE}
set.seed(1)

cl1 = rmvnorm(100,c(3,4),sigma= matrix(c(.05,0,0,.05),nrow=2)) %>% as_tibble() %>% mutate(true_clust="a")
cl2 = rmvnorm(75,c(2,1),sigma= matrix(c(.1,0,0,.1),nrow=2))%>% as_tibble() %>% mutate(true_clust="b")
cl3 = rmvnorm(100,c(-2,2),sigma= matrix(c(.25,.2,.2,.25),nrow=2))%>% as_tibble() %>% mutate(true_clust="c")
cl4 = rmvnorm(75,c(-4,-2),sigma= matrix(c(1,0,0,.01),nrow=2))%>% as_tibble() %>% mutate(true_clust="d")
cl5 = rmvnorm(75,c(-5,4),sigma= matrix(c(.1,-.05,-.05,.1),nrow=2))%>% as_tibble() %>% mutate(true_clust="e")

easy_clusters =  rbind(cl1,cl2,cl3,cl4,cl5) %>%
  rename(x=V1,y=V2) %>%
  mutate(kmeans_clust = kmeans(tibble(x,y),centers = 5,
                               nstart = 100)$cluster,
         kmeans_clust = fct(letters[kmeans_clust]),
         hiera_clust = cutree(hclust(dist(tibble(x,y)),
                                     method = "single"),k=5),
         hiera_clust = fct(letters[hiera_clust]),
         gmix_clust = Mclust(tibble(x,y),G = 5)$classification,
         gmix_clust = fct(letters[gmix_clust]),
         dbs_clust = fpc::dbscan(tibble(x,y),eps=1)$cluster,
         dbs_clust = fct(letters[dbs_clust+1]),
         spect_clust = core_spectral(Dist = dist(tibble(x,y)) %>% as.matrix(), K = 5)$labels,
         spect_clust = fct(letters[spect_clust])
  )

easy_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```

## {.center}

[**clustering approaches**]{style="color:indianred"}

- If a cluster structure does underlie the observations at hand, any approach will work

```{r,fig.align='center'}
easy_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```

## {.center}

[**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}
<!-- partitioning^[MacQueen (1967).] -->
  partitioning^[`r Citet(bib=cladag_25_bib,"MacQueen")`]
:::

  :::{.column width="80%"}
```{r,fig.align='center'}
easy_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = kmeans_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```
:::

  ::::


  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}

agglomerative^[`r Citet(bib=cladag_25_bib,"kaufman2009finding")`]

:::


  :::{.column width="80%"}

```{r,fig.align='center'}
easy_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = hiera_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```

:::

  ::::


  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}
model-based^[`r Citet(bib=cladag_25_bib,"mclachlan1988mixture")`]
:::



  :::{.column width="80%"}

```{r,fig.align='center'}
easy_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = gmix_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```

:::

  ::::


  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}
density-based^[`r Citet(bib=cladag_25_bib,"ester1996density", before="DBSCAN, ")`]
:::

  :::{.column width="80%"}
```{r,fig.align='center'}
easy_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = dbs_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```

:::

  ::::


  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

:::: {.columns}

::: {.column width="20%"}
spectral^[`r Citet(bib=cladag_25_bib,"ng2001spectral", before="NJW, ")`]
:::


  ::: {.column width="80%"}

```{r,fig.align='center'}
easy_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = spect_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```

:::

  ::::


  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

If a cluster structure does underlie the observations at hand, any approach will work

```{r}
multi_s_clusters = multishapes %>% as_tibble() %>% filter(shape!=5) %>%
  mutate(true_clust = fct(letters[shape]),
         true_clust = fct_recode(.f = true_clust,e ="f"),
         kmeans_clust = kmeans(tibble(x,y),
                               centers = 5,nstart=100)$cluster,
         kmeans_clust = fct(letters[kmeans_clust]),
         hiera_clust = cutree(hclust(dist(tibble(x,y))),k=5),
         hiera_clust = fct(letters[hiera_clust]),
         gmix_clust = Mclust(tibble(x,y),G = 5)$classification,
         gmix_clust = fct(letters[gmix_clust]),
         dbs_clust = fpc::dbscan(tibble(x,y),eps=.2)$cluster,
         dbs_clust = fct(letters[dbs_clust]),
         spect_clust = core_spectral(Dist = dist(tibble(x,y)) %>% as.matrix(), K = 5)$labels,
         spect_clust = fct(letters[spect_clust])
  )

multi_s_clusters = multi_s_clusters %>%
  mutate(spect_scores = core_spectral(Dist = dist(tibble(x,y)) %>%
                                        as.matrix(), K = 5)$spectral_scores)

multi_s_clusters = cbind(multi_s_clusters %>% dplyr::select(-spect_scores),
                         as.matrix(multi_s_clusters %>% dplyr::select(spect_scores))
) %>% as_tibble() %>%  clean_names()


```

&nbsp;

::::{.columns}

:::{.column width="20%"}

&nbsp;

&nbsp;

&nbsp;

::: {.fragment}
[OR NOT?]{style="color:indianred"}
:::

  :::


  :::{.column width="80%"}

```{r,fig.align='center'}
multi_s_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```

:::

  ::::


  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}

partitioning
(Kmeans)

:::

  :::{.column width="80%"}
```{r,fig.align='center'}
multi_s_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = kmeans_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```
:::

  ::::


  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}

agglomerative

:::

  :::{.column width="80%"}


```{r,fig.align='center'}
multi_s_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = hiera_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```
:::

  ::::

  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}

model-based
(GMM)

:::

  :::{.column width="80%"}

```{r,fig.align='center'}
multi_s_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = gmix_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```
:::

  ::::

  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}

density-based
(dbscan)

:::

  :::{.column width="80%"}

```{r,fig.align='center'}
multi_s_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = dbs_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)
```
:::

  ::::



  ## {.center}

  [**clustering approaches**]{style="color:indianred"}

::::{.columns}

:::{.column width="20%"}
spectral
:::


  :::{.column width="80%"}
```{r,fig.align='center'}
multi_s_clusters %>%
  ggplot(aes(x=x,y=y,shape = true_clust,color = spect_clust)) +
  theme_void() + geom_point(size=3.5,alpha=.75)

```
:::

  ::::
