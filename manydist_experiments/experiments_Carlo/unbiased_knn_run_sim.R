## packages
packages <- c("mclust", "class", "tidyverse", "tidymodels", "data.table", "dplyr", "MASS", "rsample", "cluster", "clusterGeneration", "parallelDist","klaR","arules","FD","StatMatch","clustMixType", "kdml", "remotes")

to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

install.packages("manydist_0.5.0.tar.gz", repos = NULL, type = "source")
library(manydist)

## source functions
source("functions_unbiased_knn.R")

source("unbiased_knn_4_4_4_4.R")
# source("unbiased_knn_6_2_6_2.R")
# source("unbiased_knn_2_6_2_6.R")
# source("unbiased_knn_6_2_4_4.R")
# source("unbiased_knn_4_4_6_2.R")
source("unbiased_knn_8_0_8_0.R")
source("unbiased_knn_8_0_4_4.R")
source("unbiased_knn_4_4_8_0.R")
source("unbiased_knn_0_8_4_4.R")
source("unbiased_knn_4_4_0_8.R")

load("unbiased_knn_4_4_4_4.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

load("unbiased_knn_6_2_6_2.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

load("unbiased_knn_2_6_2_6.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

load("unbiased_knn_6_2_4_4.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

load("unbiased_knn_4_4_6_2.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

load("unbiased_knn_8_0_8_0.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

load("unbiased_knn_8_0_4_4.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

load("unbiased_knn_4_4_8_0.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)


load("unbiased_knn_0_8_4_4.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)

load("unbiased_knn_4_4_0_8.RData")
ariplot <- ari_summary
variant <- paste(numsignal, "num. signal +", numnoise, "num. noise +",
                 catsignal, "cat. signal +", catnoise, "cat. noise", sep = " ")

ggplot(ariplot, aes(x = method, y = ari, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(x = "Distance Method", y = "Adjusted Rand Index (PAM)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none") +
  ggtitle(variant)


