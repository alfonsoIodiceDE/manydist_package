set.seed(123)
library(mvtnorm)
library(dplyr)
library(aricode)
library(cluster)    
library(mvtnorm)
library(dplyr)
library(manydist)
library(tidyverse)
library(FPDclustering)
source("R/delta_knn_ba.R")
pengs = palmerpenguins::penguins |> na.omit()

data(Star)

# load(file="Xkp.RData")
# load(file="XkpEasy.rdata")
# load(file="XkpMedium.rdata")
# 
dfc = XkpMedium |> as_tibble()
dfc = pengs |> as_tibble() |> select(-species,-year)

star = Star |> select(-Type) |> mutate(across(where(is.character), as.factor))



source("R/distance_by_method_cri.R")
# source("R/distance_by_method.R")
source("R/spectral_from_dist.R")

D_int = distance_by_method(star,method="u_dep",interaction=TRUE)

D_no_int = distance_by_method(dfc,method="u_dep",interaction =FALSE)
# D_gudmm = distance_by_method(dfc,method="gudmm",interaction=FALSE) |> as.matrix()
D_naive= distance_by_method(dfc,method="naive",interaction=FALSE)
D_gow = distance_by_method(dfc,method="gower",interaction=FALSE)
# D_dkss = distance_by_method(dfc,method="dkss",interaction=FALSE)
D_mod_gow = distance_by_method(dfc,method="mod_gower",interaction=FALSE) |> as.matrix()


source("R/spectral_from_dist.R")

myk=3

cl_int = spectral_from_dist(D_int, k=myk,affinity_method = "gaussian")
cl_no_int = spectral_from_dist(D_no_int, k=myk,affinity_method = "gaussian")
# cl_gudmm = spectral_from_dist(D_gudmm, k=myk,affinity_method = "gaussian")
cl_naive = spectral_from_dist(D_naive, k=myk,affinity_method = "gaussian")
cl_gow = spectral_from_dist(D_gow, k=myk,affinity_method = "gaussian")
cl_mod_gow = spectral_from_dist(D_mod_gow, k=myk,affinity_method = "gaussian")
# cl_dkss = spectral_from_dist(D_dkss, k=myk)


truth = c(rep(1,200),rep(2,300),rep(3,400))
adjustedRandIndex(cl_int, truth)
adjustedRandIndex(cl_no_int, truth)
# adjustedRandIndex(cl_gudmm, truth)
adjustedRandIndex(cl_gow, truth)
adjustedRandIndex(cl_mod_gow, truth)
adjustedRandIndex(cl_naive, truth)
# adjustedRandIndex(cl_dkss, truth)

# describe how data is being generated and under what conditions the interaction term is worth
dist_list=list(ab_int=D_int,ab_no_int=D_no_int,naive=D_naive,gower=D_gow)
pairs(dfc[,1:4],col=truth)
save(dist_list,file="dist_list_toy.RData")
