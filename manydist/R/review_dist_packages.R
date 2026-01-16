# ============================================================
# Check, install (if needed), and load R packages for
# distance computation competitors to manydist
# ============================================================

competitor_pkgs <- c(
  # --- Continuous ---
  "proxy",          # many distance metrics + custom ones
  "parallelDist",   # multithreaded distance computation
  "Rfast",          # fast C++ distances for large matrices

  # --- Categorical / binary ---
  "klaR",           # kmodes -> simple matching distance
  "arules",         # dissimilarities on transactions (Jaccard, cosine, etc.)

  # --- Mixed (numeric + categorical) ---
  "cluster",        # daisy(metric = "gower")
  "FD",             # gowdis()
  "StatMatch",      # gower.dist()
  "gower",          # fast C implementation of Gower
  "clustMixType"    # kproto (Gower variant and lambda-weighted)
)

# Identify missing packages
missing_pkgs <- setdiff(competitor_pkgs, rownames(installed.packages()))

# Install missing ones if needed
if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:\n")
  cat(paste0("  - ", missing_pkgs, collapse = "\n"), "\n")
  install.packages(missing_pkgs, dependencies = TRUE)
} else {
  cat("OK: All competitor packages are already installed.\n")
}

# Load all packages quietly
suppressPackageStartupMessages(
  lapply(competitor_pkgs, require, character.only = TRUE)
)

cat("\nCompetitor packages loaded successfully:\n")
print(competitor_pkgs)
