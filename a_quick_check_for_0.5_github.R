# options(repos = c(CRAN = "https://cloud.r-project.org"))
#
# install.packages(c("parallelDist","klaR","arules","FD","StatMatch","clustMixType"))
#
# remotes::install_github("alfonsoIodiceDE/manydist_package", subdir = "manydist")


if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("palmerpenguins", quietly = TRUE)) install.packages("palmerpenguins")
if (!requireNamespace("cluster", quietly = TRUE)) install.packages("cluster")
if (!requireNamespace("recipes", quietly = TRUE)) install.packages("recipes")
if (!requireNamespace("Rfast", quietly = TRUE)) install.packages("Rfast")

library(devtools)
library(manydist)
# Adjust the path to your package root as needed
# devtools::load_all("../manydist")
library(tidyverse)
library(tidymodels)
library(palmerpenguins)
set.seed(123)


peng <- palmerpenguins::penguins |>
  dplyr::select(species, island, bill_length_mm, bill_depth_mm,
                flipper_length_mm, body_mass_g, sex) |>
  tidyr::drop_na()

cont_peng = peng |> select(where(is.numeric))
cat_peng = peng |> select(where(is.factor))

peng_split=initial_split(peng,prop = .7,strata = species)
tr_peng=training(peng_split)
ts_peng=testing(peng_split)

udep_dist = mdist(x=tr_peng,new_data=ts_peng,preset="unbiased_dependent")
udep_dist
udep_dist$summary()



gow_dist = mdist(x=tr_peng,new_data=ts_peng,preset="gower")
gow_dist

gow_dist$distance |> as.matrix()
gow_dist$summary()
