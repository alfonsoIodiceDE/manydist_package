packages <- c("mclust", "class", "manydist", "tidyverse", "tidymodels", "data.table", "dplyr","MASS", "rsample")
to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

source("functions_unbiased_knn.R")

source("unbiased_knn_balanced.R")
source("unbiased_knn_continuous.R")
source("unbiased_knn_categorical.R")
source("unbiased_knn_balanced_noise.R")
source("unbiased_knn_continuous_noise.R")
source("unbiased_knn_categorical_noise.R")

load("../knn_experiments/knn_unbiased/unbiased_knn_balanced.RData")
ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.4, 1)) +
  theme(legend.position = "none")

load("../knn_experiments/knn_unbiased/unbiased_knn_continuous.RData")
ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.4, 1)) +
  theme(legend.position = "none")

load("../knn_experiments/knn_unbiased/unbiased_knn_categorical.RData")
ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.4, 1)) +
  theme(legend.position = "none")

load("../knn_experiments/knn_unbiased/unbiased_knn_balanced_noise.RData")
ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.4, 1)) +
  theme(legend.position = "none")

load("../knn_experiments/knn_unbiased/unbiased_knn_continuous_noise.RData")
ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.4, 1)) +
  theme(legend.position = "none")

load("../knn_experiments/knn_unbiased/unbiased_knn_categorical_noise.RData")
ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.4, 1)) +
  theme(legend.position = "none")

