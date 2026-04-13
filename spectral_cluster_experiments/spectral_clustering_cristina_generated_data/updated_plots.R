# ============================================================
# Manuscript figures
#   - Figure 2: illustrative mixed-data interaction plot
#   - fig_sim1: Simulation 1 ARI boxplots (no interaction)
#   - fig_sim2: Simulation 2 ARI boxplot (interaction only)
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(colorspace)
library(ggplot2)
library(mvtnorm)
library(GGally)

# ----------------------------
# Paths
# ----------------------------
data_dir <- "../spectral_cluster_experiments/data"
fig_dir  <- "../spectral_cluster_experiments/figures"

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ============================================================
# PART A — FIGURE 2
# ============================================================

set.seed(123)

pnum <- 4
pcat <- 2
p <- pnum + pcat
n <- 500

# -------- Cluster 1
n1 <- round(0.2 * n)
sig <- matrix(0.5, pnum, pnum)
diag(sig) <- 1
Xn <- rmvnorm(n1, rep(0, pnum), sigma = sig)

pbin <- pcat
betaBin <- matrix(c(1, rep(2, p)), pbin, (p + 1), byrow = TRUE)

XX <- cbind(1, Xn)
Xcat <- matrix(0, n1, pcat)
XX <- cbind(XX, Xcat)

for (i in seq_len(n1)) {
  for (j in seq_len(pbin)) {
    Xcat[i, j] <- rbinom(
      1, 1,
      p = 1 / (1 + exp(betaBin[j, 1:(j + pnum - 1)] %*% XX[i, 1:(j + pnum - 1)]))
    )
  }
}
X1 <- cbind(Xn, Xcat)

# -------- Cluster 2
n2 <- round(0.3 * n)
sig <- matrix(0.8, pnum, pnum)
diag(sig) <- 1
Xn <- rmvnorm(n2, rep(0, pnum), sigma = sig)

betaBin <- matrix(c(1, rep(-2, p)), pbin, (p + 1), byrow = TRUE)

XX <- cbind(1, Xn)
Xcat <- matrix(0, n2, pcat)
XX <- cbind(XX, Xcat)

for (i in seq_len(n2)) {
  for (j in seq_len(pbin)) {
    Xcat[i, j] <- rbinom(
      1, 1,
      p = 1 / (1 + exp(betaBin[j, 1:(j + pnum - 1)] %*% XX[i, 1:(j + pnum - 1)]))
    )
  }
}
X2 <- cbind(Xn + 4, Xcat)

# -------- Cluster 3
n3 <- n - n1 - n2

# For pnum = 4, -0.5 is not PSD; use -0.3 in the illustrative generator
sig <- matrix(-0.3, pnum, pnum)
diag(sig) <- 1
Xn <- rmvnorm(n3, rep(0, pnum), sigma = sig)

betaBin <- matrix(c(1, rep(2, p)), pbin, (p + 1), byrow = TRUE)

XX <- cbind(1, Xn)
Xcat <- matrix(0, n3, pcat)
XX <- cbind(XX, Xcat)

for (i in seq_len(n3)) {
  for (j in seq_len(pbin)) {
    Xcat[i, j] <- rbinom(
      1, 1,
      p = 1 / (1 + exp(betaBin[j, 1:(j + pnum - 1)] %*% XX[i, 1:(j + pnum - 1)]))
    )
  }
}
Xn <- sweep(Xn, 2, rep(c(-1, 5), length.out = pnum), FUN = "+")
X3 <- cbind(Xn, Xcat)

# -------- Combine
X <- rbind(X1, X2, X3)
truth <- c(rep(1, n1), rep(2, n2), rep(3, n3))

Xkp <- as.data.frame(X)
colnames(Xkp) <- paste0("V", 1:6)
for (i in 1:pcat) {
  Xkp[[pnum + i]] <- factor(Xkp[[pnum + i]])
}
Xkp$cluster <- factor(truth, labels = c("C1", "C2", "C3"))


# ============================================================
# Figure 1 — Scatterplot matrix of the illustrative dataset
# ============================================================
fig1plot <- ggpairs(
  data = Xkp,
  columns = 1:4,
  mapping = aes(color = cluster, fill = cluster, shape = cluster),
  upper = list(continuous = wrap("points", alpha = 0.7, size = 1.2)),
  lower = list(continuous = wrap("points", alpha = 0.7, size = 1.2)),
  diag  = list(
    continuous = wrap(
      "densityDiag",
      mapping = aes(fill = cluster),
      alpha = 0.5
    )
  )
) +
  scale_color_manual(values = c("indianred", "dodgerblue3", "darkseagreen4")) +
  scale_fill_manual(values = c("indianred", "dodgerblue3", "darkseagreen4")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(fig_dir, "pairsToy.pdf"),
  plot = fig1plot,
  width = 7,
  height = 7,
  units = "in"
)

# -------- Build Figure 2 data
profiles <- bind_rows(
  Xkp %>%
    group_by(cluster, V5) %>%
    summarise(across(V1:V4, mean), .groups = "drop") %>%
    pivot_longer(V1:V4, names_to = "variable", values_to = "value") %>%
    mutate(category = as.character(V5), cat_var = "V5") %>%
    select(cluster, cat_var, category, variable, value),

  Xkp %>%
    group_by(cluster, V6) %>%
    summarise(across(V1:V4, mean), .groups = "drop") %>%
    pivot_longer(V1:V4, names_to = "variable", values_to = "value") %>%
    mutate(category = as.character(V6), cat_var = "V6") %>%
    select(cluster, cat_var, category, variable, value)
)

diff_df <- profiles %>%
  pivot_wider(names_from = category, values_from = value) %>%
  mutate(diff = `1` - `0`) %>%
  filter(!is.na(diff))

fig2plot <- ggplot(diff_df, aes(x = variable, y = diff, fill = diff > 0)) +
  geom_col(width = 0.5, show.legend = FALSE, alpha = 0.75) +
  scale_fill_manual(values = c("indianred", "dodgerblue")) +
  facet_grid(cluster ~ cat_var) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Continuous variables",
    y = "Conditional mean difference (1 - 0)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(fig_dir, "figure2.pdf"),
  plot = fig2plot,
  width = 7,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = file.path(fig_dir, "figure2.png"),
  plot = fig2plot,
  width = 7,
  height = 4.8,
  units = "in",
  dpi = 300
)

# ============================================================
# PART B — SIMULATION ARI FIGURES
# ============================================================

files <- c(
  "ARIM_N1E.rdata",
  "ARIM_N2E.rdata",
  "ARIM_N1_1020.rdata",
  "ARIM_N2_1020.rdata",
  "ARIM_N1_2010.rdata",
  "ARIM_N2_2010.rdata",
  "ARIM_noise_N1.rdata",
  "ARIM_noise_N2.rdata",
  "ARIM_x_N1.rdata",
  "ARIM_x_N2.rdata"
)

scenario_order_raw <- c(
  "balanced", "balanced",
  "more_cat", "more_cat",
  "more_cont", "more_cont",
  "balanced_noise", "balanced_noise",
  "interaction", "interaction"
)

scenario_labels <- c(
  balanced       = "a) 15 con, 15 cat",
  balanced_noise = "b) 10 con, 10 cat + noise",
  more_cat       = "c) 10 con, 20 cat",
  more_cont      = "d) 20 con, 10 cat",
  interaction    = "e) interaction"
)

method_names_raw  <- c("Udep Int.", "Udep", "Gower", "mod G.", "naive")
method_names_plot <- c("ab_dis_int", "ab_dis", "gower", "mod_gower", "naive")

load_rdata_object <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e[[ls(e)[1]]]
}

build_ari_table <- function(ari_list) {
  sc_scenario_params <- tibble(
    scenario = c("more_cat", "more_cont", "balanced", "balanced_noise", "interaction"),
    cont_n   = c(10, 20, 15, 10, 3),
    cat_n    = c(20, 10, 15, 10, 3),
    noise    = c(0, 0, 0, 10, 0)
  ) |>
    crossing(obs_n = c(500, 1000))

  tibble(
    scenario = names(ari_list),
    ari = ari_list
  ) |>
    group_by(scenario) |>
    mutate(obs_n = c(500, 1000)[row_number()]) |>
    ungroup() |>
    mutate(
      ari_long = map(ari, \(x) {
        as_tibble(x, .name_repair = "minimal") |>
          setNames(method_names_raw) |>
          mutate(rep = row_number()) |>
          pivot_longer(
            cols = -rep,
            names_to = "measure",
            values_to = "ari"
          )
      })
    ) |>
    select(-ari) |>
    unnest(ari_long) |>
    left_join(sc_scenario_params, by = c("scenario", "obs_n")) |>
    mutate(
      measure = recode(
        measure,
        "Udep"      = "ab_dis",
        "Udep Int." = "ab_dis_int",
        "naive"     = "naive",
        "Gower"     = "gower",
        "mod G."    = "mod_gower"
      )
    )
}

make_plot_data <- function(ari_tbl) {
  ari_tbl |>
    mutate(
      scenario = factor(
        scenario,
        levels = names(scenario_labels),
        labels = unname(scenario_labels)
      ),
      obs_n = factor(obs_n, levels = c(500, 1000)),
      measure = factor(measure, levels = method_names_plot),
      method_n = paste0(measure, "_", obs_n),
      slide_group = case_when(
        scenario %in% c("a) 15 con, 15 cat", "b) 10 con, 10 cat + noise") ~ "Balanced / Noise",
        scenario %in% c("c) 10 con, 20 cat", "d) 20 con, 10 cat")         ~ "Unbalanced",
        TRUE                                                               ~ "Interaction"
      )
    )
}

make_method_palette <- function(methods, obs_levels = c(500, 1000)) {
  base_cols <- qualitative_hcl(length(methods))
  names(base_cols) <- methods

  cols <- c(
    setNames(lighten(base_cols, 0.3), paste0(methods, "_", obs_levels[1])),
    setNames(darken(base_cols,  0.3), paste0(methods, "_", obs_levels[2]))
  )

  breaks <- as.vector(rbind(
    paste0(methods, "_", obs_levels[1]),
    paste0(methods, "_", obs_levels[2])
  ))

  list(values = cols, breaks = breaks)
}

plot_ari_box <- function(dat, method_cols, legend_breaks, facet = TRUE) {
  p <- ggplot(dat, aes(x = measure, y = ari)) +
    geom_boxplot(
      aes(fill = method_n, group = interaction(measure, obs_n)),
      position = position_dodge2(width = 0.8, preserve = "single"),
      width = 0.7,
      outlier.alpha = 0.4
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(
      values = method_cols,
      breaks = legend_breaks,
      labels = gsub("_(500|1000)$", " (n=\\1)", legend_breaks)
    ) +
    labs(x = NULL, y = "ARI", fill = "Method") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "top"
    )

  if (facet) {
    p + facet_wrap(~ scenario, ncol = 2, scales = "free_x")
  } else {
    p
  }
}

ari_list <- setNames(
  lapply(file.path(data_dir, files), load_rdata_object),
  tools::file_path_sans_ext(basename(files))
)

names(ari_list) <- scenario_order_raw

ari_tbl  <- build_ari_table(ari_list)
ari_plot <- make_plot_data(ari_tbl)

pal <- make_method_palette(method_names_plot)
method_cols   <- pal$values
legend_breaks <- pal$breaks

fig_sim1 <- plot_ari_box(
  dat = filter(ari_plot, slide_group != "Interaction"),
  method_cols = method_cols,
  legend_breaks = legend_breaks,
  facet = TRUE
)

fig_sim2 <- plot_ari_box(
  dat = filter(ari_plot, slide_group == "Interaction"),
  method_cols = method_cols,
  legend_breaks = legend_breaks,
  facet = FALSE
)

ggsave(
  filename = file.path(fig_dir, "fig_sim1.pdf"),
  plot = fig_sim1,
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  filename = file.path(fig_dir, "fig_sim2.pdf"),
  plot = fig_sim2,
  width = 10,
  height = 8,
  units = "in"
)


# ============================================================
# Illustrative interaction dataset + scatterplot matrix
# ============================================================

library(tidyverse)
library(sampling)
library(SimDesign)
library(GGally)
library(ggplot2)

conflicted::conflicts_prefer(SimDesign::rmvnorm)

generate_interaction_example <- function(n = 1000, pnum = 6, seed = 123) {
  set.seed(seed)

  c1 <- rmvnorm(n, rep(0, pnum))

  cat1 <- rep(1, n)
  sel <- which(srswor(n * 0.6, n) == 1)
  cat1[sel[1:(n * 0.2)]] <- 2
  cat1[sel[(n * 0.2 + 1):(n * 0.6)]] <- 3
  cat1 <- factor(cat1)

  cat2 <- rep(1, n)
  cat2[which(srswor(n * 0.2, n) == 1)] <- 2
  cat2 <- factor(cat2)

  cat4 <- rep(1, n)
  sel <- which(srswor(n * 0.6, n) == 1)
  cat4[sel[1:(n * 0.3)]] <- 2
  cat4[sel[(n * 0.3 + 1):(n * 0.6)]] <- 3
  cat4 <- factor(cat4)

  cat3 <- rep(1, n)
  cat3[which(srswor(n * 0.4, n) == 1)] <- 2
  cat3 <- factor(cat3)

  pcat <- 4
  beta <- matrix(c(0.8, 0.5, -0.4, 0.8, -0.4, 0.5), 3, 2, byrow = TRUE)

  truth <- rep(3, nrow(c1))

  for (i in 1:3) {
    for (j in 1:2) {
      idx <- which(cat1 == i & cat2 == j)

      for (h in (pnum / 2):1) {
        c1[idx, h] <- apply(c1[idx, (h + 1):pnum, drop = FALSE] * beta[i, j], 1, sum)
        if (beta[i, j] == 0.5) {
          c1[idx, h] <- c1[idx, h] + 8
        }
      }

      if (beta[i, j] == 0.8) {
        truth[idx] <- 1
      } else if (beta[i, j] == -0.4) {
        truth[idx] <- 2
      }
    }
  }

  Xkp <- data.frame(cbind(c1, cat1, cat2, cat3, cat4))
  names(Xkp) <- c(paste0("V", 1:pnum), paste0("C", 1:pcat))

  for (i in 1:pnum) {
    Xkp[[i]] <- as.vector(Xkp[[i]])
  }

  for (i in (pnum + 1):(pnum + pcat)) {
    Xkp[[i]] <- factor(unlist(Xkp[[i]]))
  }

  Xkp$cluster <- factor(truth, labels = c("C1", "C2", "C3"))

  Xkp
}

# Generate one illustrative dataset
df_int <- generate_interaction_example(n = 1000, pnum = 6, seed = 123)

# Full scatterplot matrix of continuous variables
fig_interaction_pairs <- ggpairs(
  data = df_int,
  columns = paste0("V", 1:6),
  mapping = aes(color = cluster, shape = cluster),
  upper = list(continuous = wrap("points", alpha = 0.5, size = 0.9)),
  lower = list(continuous = wrap("points", alpha = 0.5, size = 0.9)),
  diag  = list(
    continuous = wrap(
      "densityDiag",
      mapping = aes(fill = cluster),
      alpha = 0.4
    )
  )
  ) +
  scale_color_manual(values = c("indianred", "dodgerblue3", "darkseagreen4")) +
  scale_fill_manual(values = c("indianred", "dodgerblue3", "darkseagreen4")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 9),
    legend.title = element_blank()
  )
# fig_interaction_pairs

ggsave(
  filename = "../spectral_cluster_experiments/figures/interaction_pairs.pdf",
  plot = fig_interaction_pairs,
  width = 8.5,
  height = 8.5,
  units = "in"
)
