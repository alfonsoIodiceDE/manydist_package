# ============================================================================
# Comparison Test Suite: mdist vs mdist_evolved
# ============================================================================

library(palmerpenguins)
library(dplyr)
library(tidyr)

# Load data
data("penguins")
penguins_complete <- penguins %>% drop_na()

set.seed(123)
train_idx <- sample(1:nrow(penguins_complete), 0.7 * nrow(penguins_complete))
train <- penguins_complete[train_idx, ]
test <- penguins_complete[-train_idx, ]

# Prepare data subsets
mixed_data <- train
cat_only <- train %>% select(species, island, sex)
cont_only <- train %>% select(bill_length_mm, bill_depth_mm, 
                              flipper_length_mm, body_mass_g)

cat("=== COMPARISON TEST SUITE: mdist vs mdist_evolved ===\n\n")

# ============================================================================
# Helper Functions
# ============================================================================

compare_distances <- function(dist1, dist2, test_name, tolerance = 1e-10) {
  # Extract distance matrices
  mat1 <- as.matrix(dist1)
  mat2 <- as.matrix(dist2)
  
  # Check dimensions
  dims_match <- all(dim(mat1) == dim(mat2))
  
  # Check values
  if (dims_match) {
    max_diff <- max(abs(mat1 - mat2))
    mean_diff <- mean(abs(mat1 - mat2))
    all_close <- all(abs(mat1 - mat2) < tolerance)
  } else {
    max_diff <- NA
    mean_diff <- NA
    all_close <- FALSE
  }
  
  # Print results
  cat(sprintf("%-50s: ", test_name))
  if (all_close && dims_match) {
    cat("✓ PASS\n")
  } else {
    cat("✗ FAIL\n")
    cat(sprintf("  Dimensions match: %s (%s vs %s)\n", 
                dims_match, 
                paste(dim(mat1), collapse = "x"),
                paste(dim(mat2), collapse = "x")))
    if (dims_match) {
      cat(sprintf("  Max difference: %.2e\n", max_diff))
      cat(sprintf("  Mean difference: %.2e\n", mean_diff))
    }
  }
  
  return(list(
    test_name = test_name,
    pass = all_close && dims_match,
    dims_match = dims_match,
    max_diff = max_diff,
    mean_diff = mean_diff
  ))
}

# ============================================================================
# SECTION 1: PRESET COMPARISONS - MIXED DATA
# ============================================================================

cat("\n=== SECTION 1: PRESET COMPARISONS - MIXED DATA ===\n\n")

results <- list()

# Test 1.1: Gower preset - non-commensurable
cat("Test 1.1: ")
old_1_1 <- mdist(
  x = mixed_data,
  preset = "gower",
  commensurable = FALSE
)
new_1_1 <- mdist_evolved(
  x = mixed_data,
  preset = "gower",
  commensurable = FALSE
)
results$test_1_1 <- compare_distances(old_1_1, new_1_1$distance, 
                                      "Gower (non-commensurable) - mixed")

# Test 1.2: Gower preset - commensurable
cat("Test 1.2: ")
old_1_2 <- mdist(
  x = mixed_data,
  preset = "gower",
  commensurable = TRUE
)
new_1_2 <- mdist_evolved(
  x = mixed_data,
  preset = "gower",
  commensurable = TRUE
)
results$test_1_2 <- compare_distances(old_1_2, new_1_2$distance, 
                                      "Gower (commensurable) - mixed")

# Test 1.3: Unbiased dependent preset
cat("Test 1.3: ")
old_1_3 <- mdist(
  x = mixed_data,
  preset = "unbiased_dependent"
)
new_1_3 <- mdist_evolved(
  x = mixed_data,
  preset = "unbiased_dependent"
)
results$test_1_3 <- compare_distances(old_1_3, new_1_3$distance, 
                                      "Unbiased dependent - mixed")

# Test 1.4: Euclidean one-hot preset
cat("Test 1.4: ")
old_1_4 <- mdist(
  x = mixed_data,
  preset = "euclidean_onehot"
)
new_1_4 <- mdist_evolved(
  x = mixed_data,
  preset = "euclidean_onehot"
)
results$test_1_4 <- compare_distances(old_1_4, new_1_4$distance, 
                                      "Euclidean one-hot - mixed")

# ============================================================================
# SECTION 2: CUSTOM PRESET COMPARISONS - MIXED DATA
# ============================================================================

cat("\n=== SECTION 2: CUSTOM PRESET - MIXED DATA ===\n\n")

# Test 2.1: Manhattan + Matching
cat("Test 2.1: ")
old_2_1 <- mdist(
  x = mixed_data,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = FALSE,
  scaling_cont = "none"
)
new_2_1 <- mdist_evolved(
  x = mixed_data,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = FALSE,
  scaling_cont = "none"
)
results$test_2_1 <- compare_distances(old_2_1, new_2_1$distance, 
                                      "Manhattan + Matching (non-comm)")

# Test 2.2: Euclidean + Total Variation (commensurable)
cat("Test 2.2: ")
old_2_2 <- mdist(
  x = mixed_data,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "tot_var_dist",
  commensurable = TRUE,
  scaling_cont = "std"
)
new_2_2 <- mdist_evolved(
  x = mixed_data,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "tot_var_dist",
  commensurable = TRUE,
  scaling_cont = "std"
)
results$test_2_2 <- compare_distances(old_2_2, new_2_2$distance, 
                                      "Euclidean + TotVar (comm, std)")

# Test 2.3: Manhattan with PC scores
cat("Test 2.3: ")
old_2_3 <- mdist(
  x = mixed_data,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  scaling_cont = "pc_scores",
  ncomp = 2, 
  commensurable = FALSE
)
new_2_3 <- mdist_evolved(
  x = mixed_data,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  scaling_cont = "pc_scores",
  ncomp = 2
)
results$test_2_3 <- compare_distances(old_2_3, new_2_3$distance, 
                                      "Manhattan + Matching (PC scores, ncomp=2)")

# Test 2.4: Euclidean + HLeucl special combination
cat("Test 2.4: ")
old_2_4 <- mdist(
  x = mixed_data,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "HLeucl",
  scaling_cont = "std",
  commensurable = FALSE
)
new_2_4 <- mdist_evolved(
  x = mixed_data,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "HLeucl",
  scaling_cont = "std"
)
results$test_2_4 <- compare_distances(old_2_4, new_2_4$distance, 
                                      "Euclidean + HLeucl (special sqrt)")

# Test 2.5: With threshold parameter
cat("Test 2.5: ")
old_2_5 <- mdist(
  x = mixed_data,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "matching",
  scaling_cont = "pc_scores",
  threshold = 0.95,
  commensurable = FALSE
)
new_2_5 <- mdist_evolved(
  x = mixed_data,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "matching",
  scaling_cont = "pc_scores",
  threshold = 0.95
)
results$test_2_5 <- compare_distances(old_2_5, new_2_5$distance, 
                                      "Euclidean + Matching (threshold=0.95)")

# ============================================================================
# SECTION 3: CATEGORICAL DATA ONLY
# ============================================================================

cat("\n=== SECTION 3: CATEGORICAL DATA ONLY ===\n\n")

# Test 3.1: Gower - categorical only
cat("Test 3.1: ")
old_3_1 <- mdist(
  x = cat_only,
  preset = "gower"
)
new_3_1 <- mdist_evolved(
  x = cat_only,
  preset = "gower"
)
results$test_3_1 <- compare_distances(old_3_1, new_3_1$distance, 
                                      "Gower - categorical only")

# Test 3.2: Unbiased dependent - categorical only
cat("Test 3.2: ")
old_3_2 <- mdist(
  x = cat_only,
  preset = "unbiased_dependent"
)
new_3_2 <- mdist_evolved(
  x = cat_only,
  preset = "unbiased_dependent"
)
results$test_3_2 <- compare_distances(old_3_2, new_3_2$distance, 
                                      "Unbiased dependent - categorical only")

# Test 3.3: Custom matching - categorical only
cat("Test 3.3: ")
old_3_3 <- mdist(
  x = cat_only,
  preset = "custom",
  distance_cat = "matching",
  commensurable = FALSE
)
new_3_3 <- mdist_evolved(
  x = cat_only,
  preset = "custom",
  distance_cat = "matching",
  commensurable = FALSE
)
results$test_3_3 <- compare_distances(old_3_3, new_3_3$distance, 
                                      "Custom matching - categorical only")

# Test 3.4: Custom total variation - categorical only
cat("Test 3.4: ")
old_3_4 <- mdist(
  x = cat_only,
  preset = "custom",
  distance_cat = "tot_var_dist",
  commensurable = TRUE
)
new_3_4 <- mdist_evolved(
  x = cat_only,
  preset = "custom",
  distance_cat = "tot_var_dist",
  commensurable = TRUE
)
results$test_3_4 <- compare_distances(old_3_4, new_3_4$distance, 
                                      "Custom tot_var_dist - categorical only")

# ============================================================================
# SECTION 4: CONTINUOUS DATA ONLY
# ============================================================================

cat("\n=== SECTION 4: CONTINUOUS DATA ONLY ===\n\n")

# Test 4.1: Gower - continuous only
cat("Test 4.1: ")
old_4_1 <- mdist(
  x = cont_only,
  preset = "gower"
)
new_4_1 <- mdist_evolved(
  x = cont_only,
  preset = "gower"
)
results$test_4_1 <- compare_distances(old_4_1, new_4_1$distance, 
                                      "Gower - continuous only")

# Test 4.2: Unbiased dependent - continuous only
cat("Test 4.2: ")
old_4_2 <- mdist(
  x = cont_only,
  preset = "unbiased_dependent"
)
new_4_2 <- mdist_evolved(
  x = cont_only,
  preset = "unbiased_dependent"
)
results$test_4_2 <- compare_distances(old_4_2, new_4_2$distance, 
                                      "Unbiased dependent - continuous only")

# Test 4.3: Euclidean one-hot - continuous only
cat("Test 4.3: ")
old_4_3 <- mdist(
  x = cont_only,
  preset = "euclidean_onehot"
)
new_4_3 <- mdist_evolved(
  x = cont_only,
  preset = "euclidean_onehot"
)
results$test_4_3 <- compare_distances(old_4_3, new_4_3$distance, 
                                      "Euclidean one-hot - continuous only")

# Test 4.4: Custom Manhattan - continuous only
cat("Test 4.4: ")
old_4_4 <- mdist(
  x = cont_only,
  preset = "custom",
  distance_cont = "manhattan",
  scaling_cont = "none",
  commensurable = FALSE
)
new_4_4 <- mdist_evolved(
  x = cont_only,
  preset = "custom",
  distance_cont = "manhattan",
  scaling_cont = "none"
)
results$test_4_4 <- compare_distances(old_4_4, new_4_4$distance, 
                                      "Custom Manhattan - continuous only")

# Test 4.5: Custom Euclidean with standardization
cat("Test 4.5: ")
old_4_5 <- mdist(
  x = cont_only,
  preset = "custom",
  distance_cont = "euclidean",
  scaling_cont = "std",
  commensurable = FALSE
)
new_4_5 <- mdist_evolved(
  x = cont_only,
  preset = "custom",
  distance_cont = "euclidean",
  scaling_cont = "std"
)
results$test_4_5 <- compare_distances(old_4_5, new_4_5$distance, 
                                      "Custom Euclidean (std) - continuous only")

# Test 4.6: Custom with PC scores
cat("Test 4.6: ")
old_4_6 <- mdist(
  x = cont_only,
  preset = "custom",
  distance_cont = "manhattan",
  scaling_cont = "pc_scores",
  ncomp = 2, 
  commensurable = FALSE
)
new_4_6 <- mdist_evolved(
  x = cont_only,
  preset = "custom",
  distance_cont = "manhattan",
  scaling_cont = "pc_scores",
  ncomp = 2
)
results$test_4_6 <- compare_distances(old_4_6, new_4_6$distance, 
                                      "Custom Manhattan (PC scores, ncomp=2)")

# ============================================================================
# SECTION 5: VALIDATION SET COMPARISONS
# ============================================================================

cat("\n=== SECTION 5: VALIDATION SET (TRAIN-TO-TEST) ===\n\n")

# Test 5.1: Gower with validation set
cat("Test 5.1: ")
old_5_1 <- mdist(
  x = train,
  validate_x = test,
  preset = "gower",
  commensurable = FALSE
)
new_5_1 <- mdist_evolved(
  x = train,
  validate_x = test,
  preset = "gower",
  commensurable = FALSE
)
results$test_5_1 <- compare_distances(old_5_1, new_5_1$distance, 
                                      "Gower with validation (non-comm)")

# Test 5.2: Gower with validation set - commensurable
cat("Test 5.2: ")
old_5_2 <- mdist(
  x = train,
  validate_x = test,
  preset = "gower",
  commensurable = TRUE
)
new_5_2 <- mdist_evolved(
  x = train,
  validate_x = test,
  preset = "gower",
  commensurable = TRUE
)
results$test_5_2 <- compare_distances(old_5_2, new_5_2$distance, 
                                      "Gower with validation (comm)")

# Test 5.3: Unbiased dependent with validation
cat("Test 5.3: ")
old_5_3 <- mdist(
  x = train,
  validate_x = test,
  preset = "unbiased_dependent"
)
new_5_3 <- mdist_evolved(
  x = train,
  validate_x = test,
  preset = "unbiased_dependent"
)
results$test_5_3 <- compare_distances(old_5_3, new_5_3$distance, 
                                      "Unbiased dependent with validation")

# Test 5.4: Euclidean one-hot with validation
cat("Test 5.4: ")
old_5_4 <- mdist(
  x = train,
  validate_x = test,
  preset = "euclidean_onehot"
)
new_5_4 <- mdist_evolved(
  x = train,
  validate_x = test,
  preset = "euclidean_onehot"
)
results$test_5_4 <- compare_distances(old_5_4, new_5_4$distance, 
                                      "Euclidean one-hot with validation")

# Test 5.5: Custom with validation and response
cat("Test 5.5: ")
old_5_5 <- mdist(
  x = train,
  validate_x = test,
  response = "species",
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "tot_var_dist",
  scaling_cont = "std",
  commensurable = TRUE
)
new_5_5 <- mdist_evolved(
  x = train,
  validate_x = test,
  response = "species",
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "tot_var_dist",
  scaling_cont = "std",
  commensurable = TRUE
)
results$test_5_5 <- compare_distances(old_5_5, new_5_5$distance, 
                                      "Custom with validation + response")

# ============================================================================
# SECTION 6: COMMENSURABLE TESTS
# ============================================================================

cat("\n=== SECTION 6: COMMENSURABLE PARAMETER TESTS ===\n\n")

# Test 6.1: Non-commensurable
cat("Test 6.1: ")
old_6_1 <- mdist(
  x = mixed_data,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = FALSE,
  scaling_cont = "std"
)
new_6_1 <- mdist_evolved(
  x = mixed_data,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = FALSE,
  scaling_cont = "std"
)
results$test_6_1 <- compare_distances(old_6_1, new_6_1$distance, 
                                      "Commensurable = FALSE")

# Test 6.2: Commensurable
cat("Test 6.2: ")
old_6_2 <- mdist(
  x = mixed_data,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = TRUE,
  scaling_cont = "std"
)
new_6_2 <- mdist_evolved(
  x = mixed_data,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = TRUE,
  scaling_cont = "std"
)
results$test_6_2 <- compare_distances(old_6_2, new_6_2$distance, 
                                      "Commensurable = TRUE")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== SUMMARY ===\n\n")

# Convert results to data frame
summary_df <- data.frame(
  test = sapply(results, function(x) x$test_name),
  pass = sapply(results, function(x) x$pass),
  dims_match = sapply(results, function(x) x$dims_match),
  max_diff = sapply(results, function(x) ifelse(is.na(x$max_diff), NA, x$max_diff)),
  mean_diff = sapply(results, function(x) ifelse(is.na(x$mean_diff), NA, x$mean_diff)),
  row.names = NULL,
  stringsAsFactors = FALSE
)

# Count passes and fails
n_pass <- sum(summary_df$pass)
n_fail <- sum(!summary_df$pass)
n_total <- nrow(summary_df)

cat(sprintf("Total tests: %d\n", n_total))
cat(sprintf("Passed: %d (%.1f%%)\n", n_pass, 100 * n_pass / n_total))
cat(sprintf("Failed: %d (%.1f%%)\n", n_fail, 100 * n_fail / n_total))

if (n_fail > 0) {
  cat("\n=== FAILED TESTS ===\n\n")
  failed_tests <- summary_df[!summary_df$pass, ]
  for (i in 1:nrow(failed_tests)) {
    cat(sprintf("%s:\n", failed_tests$test[i]))
    cat(sprintf("  Dimensions match: %s\n", failed_tests$dims_match[i]))
    if (!is.na(failed_tests$max_diff[i])) {
      cat(sprintf("  Max difference: %.2e\n", failed_tests$max_diff[i]))
      cat(sprintf("  Mean difference: %.2e\n", failed_tests$mean_diff[i]))
    }
    cat("\n")
  }
}

# Print full summary table
cat("\n=== DETAILED RESULTS ===\n\n")
print(summary_df, row.names = FALSE)

# Save results
saveRDS(list(
  summary = summary_df,
  detailed_results = results,
  timestamp = Sys.time()
), file = "comparison_results.rds")

cat("\n=== COMPARISON COMPLETE ===\n")
cat("Results saved to: comparison_results.rds\n")
