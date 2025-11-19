packages <- c("mclust", "class", "manydist", "tidyverse", "tidymodels", "data.table", "dplyr","MASS", "rsample")
to_install <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

source("functions.R")

source("sim_balanced.R")
source("sim_balanced_noise.R")
source("sim_numerical.R")
source("sim_numerical_noise.R")
source("sim_categorical.R")
source("sim_categorical_noise.R")

load("sim_balanced.RData")
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

load("sim_balanced_noise.RData")
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

load("sim_numerical.RData")
ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.6, 1)) +
  theme(legend.position = "none")

load("sim_numerical_noise.RData")
ggplot(acc_summary, aes(x = method, y = accuracy, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    x = "Distance Method",
    y = "Accuracy"
  ) +
  scale_y_continuous(limits = c(0.6, 1)) +
  theme(legend.position = "none")

load("sim_categorical.RData")
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

load("sim_categorical_noise.RData")
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
