suppressPackageStartupMessages({
  library(clue)
  library(ggplot2)
  library(manydist)
  library(mclust)
  library(patchwork)
})

# -------------------------------------------------------------------------
# Interaction-only half-moon design
#
# Continuous block: two moons.
# Signal factor: three unequal arc bands (A/B/C) repeated on both moons.
# Nuisance factor: exactly independent of the signal factor.
# Target: the six combinations moon x band.
#
# Consequently:
#   * the continuous block has only the two-moon structure;
#   * the categorical TVD for `band` is zero because `band` and `nuisance`
#     are independent;
#   * the six target groups are visible only after the continuous-categorical
#     interaction is included.
# -------------------------------------------------------------------------

find_script_path <- function() {
  # Rscript supplies --file=...; source() records the current file in `ofile`.
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", script_arg[[1]])))
  }

  source_files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) "" else as.character(frame$ofile[[1]])
    },
    character(1)
  )
  source_files <- source_files[nzchar(source_files)]
  if (length(source_files) > 0L) {
    return(normalizePath(tail(source_files, 1L)))
  }

  # When code is sent to the console with RStudio's Run button, neither of the
  # mechanisms above is populated. In that case, ask RStudio for the active
  # editor path without making rstudioapi a required dependency.
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    editor_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(...) ""
    )
    if (nzchar(editor_path) && file.exists(editor_path)) {
      return(normalizePath(editor_path))
    }
  }

  # Final fallbacks cover execution from the talk root, its R/ directory, or
  # the parent manydist_package directory.
  script_name <- "interaction_moons_demo.R"
  candidates <- c(
    file.path(getwd(), script_name),
    file.path(getwd(), "R", script_name),
    file.path(getwd(), "RSS26_talk", "R", script_name)
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) > 0L) {
    return(normalizePath(candidates[[1]]))
  }

  stop(
    "Cannot locate interaction_moons_demo.R. ",
    "Open the saved file in RStudio, or set the working directory to RSS26_talk."
  )
}

script_path <- find_script_path()
talk_dir <- dirname(dirname(script_path))
out_dir <- file.path(talk_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260728)
n_per_moon <- 360L
noise_sd <- 0.045
cut_1 <- 0.10
cut_2 <- 0.55

# A fixed grid makes the band sizes exact:
# A = 10%, B = 45%, C = 45% within each moon.
theta_upper <- ((seq_len(n_per_moon) - 0.5) / n_per_moon) * pi
theta_lower <- ((seq_len(n_per_moon) - 0.5) / n_per_moon) * pi

upper <- data.frame(
  V1 = cos(theta_upper) + rnorm(n_per_moon, sd = noise_sd),
  V2 = sin(theta_upper) + rnorm(n_per_moon, sd = noise_sd),
  arc = theta_upper / pi,
  moon = "upper"
)
lower <- data.frame(
  V1 = 1 - cos(theta_lower) + rnorm(n_per_moon, sd = noise_sd),
  V2 = 0.50 - sin(theta_lower) + rnorm(n_per_moon, sd = noise_sd),
  arc = theta_lower / pi,
  moon = "lower"
)
dat <- rbind(upper, lower)

dat$band <- cut(
  dat$arc,
  breaks = c(-Inf, cut_1, cut_2, Inf),
  labels = c("A", "B", "C")
)

# Balance N1/N2/N3 exactly inside every band, then randomize their locations.
dat$nuisance <- NA_character_
for (lev in levels(dat$band)) {
  idx <- which(dat$band == lev)
  dat$nuisance[idx] <- sample(rep(c("N1", "N2", "N3"), length.out = length(idx)))
}
dat$nuisance <- factor(dat$nuisance, levels = c("N1", "N2", "N3"))
dat$moon <- factor(dat$moon, levels = c("upper", "lower"))
dat$truth <- interaction(dat$moon, dat$band, sep = " x ", drop = TRUE)

x <- dat[c("V1", "V2", "band", "nuisance")]

# Keep commensurable = FALSE here. With an exactly null TVD block there is
# nothing to rescale; normalizing a near-zero empirical TVD to unit mean would
# manufacture a categorical main effect and obscure the interaction-only case.
d_no_int <- mdist(
  x,
  preset = "custom",
  method_cat = "tvd",
  method_num = "pc_scores",
  commensurable = FALSE,
  interaction = FALSE
)
d_with_int <- mdist(
  x,
  preset = "custom",
  method_cat = "tvd",
  method_num = "pc_scores",
  commensurable = FALSE,
  interaction = TRUE,
  prop_nn = 0.05,
  score = "ba",
  decision = "prior_corrected"
)

k <- nlevels(dat$truth)
set.seed(101)
cl_no_int <- spectral_from_dist(
  as.matrix(d_no_int$distance),
  k = k,
  affinity_method = "selftune"
)
set.seed(101)
cl_with_int <- spectral_from_dist(
  as.matrix(d_with_int$distance),
  k = k,
  affinity_method = "selftune"
)

ari_no_int <- adjustedRandIndex(cl_no_int, dat$truth)
ari_with_int <- adjustedRandIndex(cl_with_int, dat$truth)

# Extract the two level-dissimilarity matrices used in the argument.
d_num <- manydist:::ndist(
  dat[c("V1", "V2")],
  method = "manhattan",
  commensurable = FALSE,
  scaling = "pc_scores"
) |>
  as.matrix()

delta_int <- manydist:::delta_int_knn(
  D = d_num,
  labels = dat$band,
  pi_nn = 0.05,
  score = "ba",
  decision = "prior_corrected"
)

cat_fit <- manydist:::cdist(
  dat[c("band", "nuisance")],
  method = "tvd",
  commensurable = FALSE
)
delta_tvd <- as.matrix(cat_fit$delta)[seq_len(nlevels(dat$band)),
                                      seq_len(nlevels(dat$band))]
dimnames(delta_tvd) <- list(levels(dat$band), levels(dat$band))

metrics <- data.frame(
  distance = c("without delta_int", "with delta_int"),
  ARI = c(ari_no_int, ari_with_int)
)
write.csv(
  metrics,
  file.path(out_dir, "interaction_moons_metrics.csv"),
  row.names = FALSE
)

# Match arbitrary spectral labels to truth labels for consistent plot colours.
align_to_truth <- function(pred, truth) {
  tab <- table(factor(pred), truth)
  assignment <- as.integer(solve_LSAP(tab, maximum = TRUE))
  factor(levels(truth)[assignment[pred]], levels = levels(truth))
}
dat$cl_no_int <- align_to_truth(cl_no_int, dat$truth)
dat$cl_with_int <- align_to_truth(cl_with_int, dat$truth)

band_cols <- c(A = "#E76F51", B = "#2A9D8F", C = "#3A86FF")
truth_cols <- c(
  "upper x A" = "#9B2226",
  "lower x A" = "#EE9B00",
  "upper x B" = "#0A9396",
  "lower x B" = "#94D2BD",
  "upper x C" = "#3A0CA3",
  "lower x C" = "#4CC9F0"
)

base_theme <- theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11)
  )

p_geometry <- ggplot(dat, aes(V1, V2)) +
  geom_point(colour = "#5B7F86", size = 2.1, alpha = 0.86) +
  coord_equal() +
  labs(
    title = "The numerical geometry",
    subtitle = "Two connected half-moons"
  ) +
  base_theme +
  theme(legend.position = "none")

p_band <- ggplot(dat, aes(V1, V2, colour = band, shape = nuisance)) +
  geom_point(size = 3, alpha = 0.88) +
  scale_colour_manual(values = band_cols) +
  scale_shape_manual(values = c(N1 = 16, N2 = 17, N3 = 15)) +
  coord_equal() +
  labs(
    title = "Both categorical variables mapped",
    subtitle = "Band forms unequal arcs; nuisance is spatially random",
    colour = "Signal band",
    shape = "Nuisance"
  ) +
  guides(
    colour = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  base_theme

p_truth <- ggplot(dat, aes(V1, V2, colour = truth)) +
  geom_point(size = 3.0, alpha = 0.88) +
  scale_colour_manual(values = truth_cols) +
  coord_equal() +
  labs(
    title = "Target: moon x band",
    subtitle = "Six interaction-defined groups"
  ) +
  base_theme +
  theme(legend.position = "none")

p_no <- ggplot(dat, aes(V1, V2, colour = cl_no_int)) +
  geom_point(size = 3.0, alpha = 0.88) +
  scale_colour_manual(values = truth_cols) +
  coord_equal() +
  labs(
    title = "Spectral without interaction",
    subtitle = sprintf("ARI = %.3f", ari_no_int)
  ) +
  base_theme +
  theme(legend.position = "none")

p_yes <- ggplot(dat, aes(V1, V2, colour = cl_with_int)) +
  geom_point(size = 3.0, alpha = 0.88) +
  scale_colour_manual(values = truth_cols) +
  coord_equal() +
  labs(
    title = "Spectral with interaction",
    subtitle = sprintf("ARI = %.3f", ari_with_int)
  ) +
  base_theme +
  theme(legend.position = "none")

main_plot <- (p_band | p_truth) / (p_no | p_yes) +
  plot_annotation(
    title = expression(Delta[int] * " reveals the cross-type partition"),
    theme = theme(plot.title = element_text(face = "bold", size = 21))
  )

ggsave(
  file.path(out_dir, "interaction_moons_demo.png"),
  main_plot,
  width = 13,
  height = 10,
  dpi = 220,
  bg = "white"
)

clean_for_slide <- function(p, keep_legend = FALSE) {
  p +
    labs(title = NULL, subtitle = NULL) +
    theme(
      plot.margin = margin(6, 6, 6, 6),
      legend.position = if (keep_legend) "bottom" else "none"
    )
}

individual_plots <- list(
  interaction_moons_geometry = clean_for_slide(p_geometry),
  interaction_moons_categories = clean_for_slide(p_band, keep_legend = TRUE),
  interaction_moons_truth = clean_for_slide(p_truth),
  interaction_moons_no_interaction = clean_for_slide(p_no),
  interaction_moons_with_interaction = clean_for_slide(p_yes)
)

for (nm in names(individual_plots)) {
  ggsave(
    file.path(out_dir, paste0(nm, ".png")),
    individual_plots[[nm]],
    width = 8,
    height = 5.2,
    dpi = 220,
    bg = "white"
  )
}

matrix_long <- function(mat, component) {
  out <- as.data.frame(as.table(mat))
  names(out) <- c("from", "to", "value")
  out$component <- component
  out
}
delta_df <- rbind(
  matrix_long(delta_tvd, "Categorical main effect: delta_tvd"),
  matrix_long(delta_int, "Cross-type interaction: delta_int")
)
delta_df$to <- factor(delta_df$to, levels = levels(dat$band))
delta_df$from <- factor(delta_df$from, levels = rev(levels(dat$band)))

delta_plot <- ggplot(delta_df, aes(to, from, fill = value)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 5) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#D1495B",
    limits = c(0, 1)
  ) +
  coord_equal() +
  facet_wrap(~component, nrow = 1) +
  labs(
    title = "Why the interaction changes the graph",
    subtitle = "The nuisance factor makes the TVD block null; the numerical geometry separates the band levels.",
    x = NULL,
    y = NULL,
    fill = "dissimilarity"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "interaction_moons_deltas.png"),
  delta_plot,
  width = 11,
  height = 5,
  dpi = 220,
  bg = "white"
)

ggsave(
  file.path(out_dir, "interaction_moons_deltas_clean.png"),
  delta_plot +
    labs(title = NULL, subtitle = NULL) +
    theme(plot.margin = margin(6, 6, 6, 6)),
  width = 10,
  height = 4.4,
  dpi = 220,
  bg = "white"
)

cat(sprintf("ARI without delta_int: %.3f\n", ari_no_int))
cat(sprintf("ARI with delta_int:    %.3f\n", ari_with_int))
cat(sprintf("Plots and metrics written to: %s\n", out_dir))
cat("delta_tvd for band:\n")
print(round(delta_tvd, 3))
cat("delta_int for band:\n")
print(round(delta_int, 3))
