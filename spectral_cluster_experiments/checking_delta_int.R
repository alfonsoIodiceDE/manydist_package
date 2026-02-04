rm(list=ls())

#| echo: FALSE
#| message: FALSE
#| warning: FALSE



packages <- c(
  "aricode", "cluster", "fpc", "manydist", "patchwork", "tidyverse",
  "tidymodels", "varhandle", "vegan", "kdml", "ggplot2", "viridis",
  "dplyr", "entropy", "kernlab","data.table","fastDummies", "mclust","movMF", "Matrix","mlbench", "philentropy", "conflicted",
  "kableExtra"
)

# Install any missing packages from CRAN
missing <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(missing) > 0) {
  install.packages(missing)
}

# Load all packages quietly
invisible(
  lapply(packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  })
)

devtools::load_all()

# Optional: Resolve common function conflicts
tidymodels_prefer()
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("distance", "philentropy")
conflicts_prefer(purrr::transpose)
conflicts_prefer(base::as.matrix)

n_a = 1200
n_b = 800

library(movMF)


non_convex_data_a_cl_1_moo = rmovMF(n_a,  c(0, 2))
non_convex_data_a_cl_2_moo = rmovMF(n_b,  c(0, -2))

non_convex_data_a_moo = rbind(
  tibble(V1=non_convex_data_a_cl_1_moo[,1],V2=non_convex_data_a_cl_1_moo[,2]) %>%
    mutate(V1=(V1+.5)+rnorm(n(),sd=.1),
           V2=(V2)+rnorm(n(),sd=.1),
           category="a"),
  tibble(V1=non_convex_data_a_cl_2_moo[,1],V2=non_convex_data_a_cl_2_moo[,2]) %>%
    mutate(V1=(V1-.5)+rnorm(n(),sd=.1),
           V2=(V2)+rnorm(n(),sd=.1),
           category="b")
) %>% mutate(overlap="no_overlap")%>% dplyr::select(V1,V2,category,overlap)





toy_structure = tibble(replicate=1:25) |>
  mutate(seed=replicate,
         toy_data=map(.x=seed,~{
           set.seed(.x)
           non_convex_data_a_moo %>%
             mutate(
               cat_extra = sample(c("C","D"), n(), replace = TRUE),
               cat_signal = case_when(
                 category == "a" ~ "A",   # keep category "a" → well-separated group A
                 category == "b" ~ "B"    # keep category "b" → well-separated group B
               ),
               switching = rbinom(n(), 1, 0.3),  # 30% chance to switch to C or D
               cat_signal = ifelse(switching == 1, cat_extra, cat_signal),
               cat_signal = factor(cat_signal, levels = c("A","B","C","D")),
               cat_noise3 = sample(c("noise_1","noise_2","noise_3"), n(), replace = TRUE),
               cat_noise5 = sample(paste0("noise_",1:5), n(), replace = TRUE),
               cat_noise7 = sample(paste0("noise_",1:7), n(), replace = TRUE),
               cat_noise9 = sample(paste0("noise_",1:9), n(), replace = TRUE),
               cat_noise11 = sample(paste0("noise_",1:11), n(), replace = TRUE),
               truth = case_when(
                 cat_signal == "A" ~ "K1",
                 cat_signal == "B" ~ "K2",
                 cat_signal %in% c("C","D") ~ "K3"
               ) |> factor(levels = c("K1","K2","K3")),
               flip_K = rbinom(n(), 1, 0.05),  # 5% chance to flip
               truth = if_else(
                 flip_K == 1,
                 sample(c("K1","K2","K3"), n(), replace = TRUE),  # random reassignment
                 truth
               )
             ) |> select(V1, V2, cat_signal, starts_with("cat_noise"),truth)
         }
         )
  ) |> mutate(
        dfc = map(toy_data,~.x |> select(-truth) |> mutate(across(where(is.character), as.factor))),
        truth = map(toy_data,~.x$truth),
        myk = map_int(truth, ~length(unique(.x)))) |>
  select(-seed)




# |> select(-seed) |> crossing(methods) |>
#   mutate(
#     dfc = map(toy_data,~.x |> select(-truth) |> mutate(across(where(is.character), as.factor))),
#     truth = map(toy_data,~.x$truth),
#     myk = map_int(truth, ~length(unique(.x)))) |> select(-seed) |> crossing(methods) |>
#   mutate(
#     dfc = map(toy_data,~.x |> select(-truth) |> mutate(across(where(is.character), as.factor))),
#     truth = map(toy_data,~.x$truth),
#     myk = map_int(truth, ~length(unique(.x))),
#     distance_matrix = map2(.x=dfc,.y=methods,
#                            ~distance_by_method(.x, method=ifelse(.y=="u_dep_int","u_dep",.y),
#                                                interaction = ifelse(.y=="u_dep_int",TRUE,FALSE)) |> as.matrix()),
#     clustering = map2(.x=distance_matrix, .y=myk, ~spectral_from_dist(as.matrix(.x), k=.y)),
#     ARI = map2_dbl(clustering, truth, adjustedRandIndex)
#     # ,
#     # clustering_pam = map2(.x=distance_matrix, .y=myk, ~cluster::pam(as.dist(.x), k=.y)$clustering),
#     # ARI_pam = map2_dbl(clustering_pam, truth, adjustedRandIndex)
#   )

#
# toy_structure |> select(replicate,methods,ARI) |>
#   ggplot(aes(x=methods,y=ARI,fill=methods)) + geom_boxplot() +
#   theme_bw() + theme(legend.position = "none")



df1 <- toy_structure$dfc[[1]]
truth<- toy_structure$truth[[1]]
myk <- toy_structure$myk[[1]]

myD=stats::dist(df1 |> select(where(is.numeric)),method="manhattan") |> as.matrix()

delta_int_knn(myD,labels=truth,pi_nn = 0.1,score="ba",decision="prior_corrected")
delta_knn_ba(myD,labels=truth,pi_nn = 0.1,decision="prior_corrected")

udep_int = mdist(x=df1,preset="custom", distance_cont = "manhattan",
                 distance_cat = "tot_var_dist",
                 commensurable = FALSE,
                 scaling_cont = 'pc_scores',interaction = TRUE,prop_nn = .1,score="logloss")

udep_int_ba = mdist(x=df1,preset="custom", distance_cont = "manhattan",
                 distance_cat = "tot_var_dist",
                 commensurable = FALSE,
                 scaling_cont = 'pc_scores',interaction = TRUE,prop_nn = .1,score="ba",decision="prior_corrected")

udep_no_int = mdist(x=df1,preset="custom", distance_cont = "manhattan",
                 distance_cat = "tot_var_dist",
                 commensurable = FALSE,
                 scaling_cont = 'pc_scores',interaction = FALSE)

D_int = udep_int$distance |> as.matrix()
D_int_ba = udep_int_ba$distance |> as.matrix()
D_no_int = udep_no_int$distance |> as.matrix()


naive = mdist(x=df1,preset="euclidean_onehot")
D_naive = naive$distance |> as.matrix()
gow = mdist(x=df1,preset="gower")
D_gow = gow$distance |> as.matrix()
mod_gow = mdist(x=df1,preset="mod_gower")
D_mod_gow = mod_gow$distance |> as.matrix()

cl_int = spectral_from_dist(D_int, k=myk,affinity_method = "gaussian")
cl_int_ba = spectral_from_dist(D_int_ba, k=myk,affinity_method = "gaussian")
cl_no_int = spectral_from_dist(D_no_int, k=myk,affinity_method = "gaussian")
cl_naive = spectral_from_dist(D_naive, k=myk,affinity_method = "gaussian")
cl_gow = spectral_from_dist(D_gow, k=myk,affinity_method = "gaussian")
cl_mod_gow = spectral_from_dist(D_mod_gow, k=myk,affinity_method = "gaussian")

adjustedRandIndex(cl_int, truth)
adjustedRandIndex(cl_int_ba, truth)
adjustedRandIndex(cl_no_int, truth)
adjustedRandIndex(cl_naive, truth)
adjustedRandIndex(cl_gow, truth)
adjustedRandIndex(cl_mod_gow, truth)
