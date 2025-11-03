# ============================================================================
# Comprehensive Testing Suite for mdist_evolved Function
# Using Palmer Penguins Dataset
# ============================================================================

library(palmerpenguins)
library(dplyr)
library(tidyr)

# Load the complete dataset
data("penguins")

# Prepare different versions of the data for testing
penguins_complete <- penguins %>% drop_na()

# Create train/test split for validation testing
set.seed(123)
train_idx <- sample(1:nrow(penguins_complete), size = 0.7 * nrow(penguins_complete))
penguins_train <- penguins_complete[train_idx, ]
penguins_test <- penguins_complete[-train_idx, ]

# ============================================================================
# SECTION 1: MIXED DATA TESTS (Categorical + Continuous)
# ============================================================================

cat("\n=== SECTION 1: MIXED DATA (Categorical + Continuous) ===\n")

# Test 1.1: Gower preset - non-commensurable
cat("\nTest 1.1: Gower preset (non-commensurable)\n")
test_1_1 <- mdist_evolved(
  x = penguins_train,
  preset = "gower",
  commensurable = FALSE
)
print(summary(as.vector(test_1_1$distance)))

# Test 1.2: Gower preset - commensurable
cat("\nTest 1.2: Gower preset (commensurable)\n")
test_1_2 <- mdist_evolved(
  x = penguins_train,
  preset = "gower",
  commensurable = TRUE
)
print(summary(as.vector(test_1_2$distance)))

# Test 1.3: Gower preset with validation set
cat("\nTest 1.3: Gower preset with validation set\n")
test_1_3 <- mdist_evolved(
  x = penguins_train,
  validate_x = penguins_test,
  preset = "gower",
  commensurable = FALSE
)
print(dim(as.matrix(test_1_3$distance)))

# Test 1.4: Unbiased dependent preset
cat("\nTest 1.4: Unbiased dependent preset\n")
test_1_4 <- mdist_evolved(
  x = penguins_train,
  preset = "unbiased_dependent"
)
print(summary(as.vector(test_1_4$distance)))

# Test 1.5: Unbiased dependent with validation
cat("\nTest 1.5: Unbiased dependent with validation\n")
test_1_5 <- mdist_evolved(
  x = penguins_train,
  validate_x = penguins_test,
  preset = "unbiased_dependent"
)
print(dim(as.matrix(test_1_5$distance)))

# Test 1.6: Euclidean one-hot preset
cat("\nTest 1.6: Euclidean one-hot preset\n")
test_1_6 <- mdist_evolved(
  x = penguins_train,
  preset = "euclidean_onehot"
)
print(summary(as.vector(test_1_6$distance)))

# Test 1.7: Euclidean one-hot with validation
cat("\nTest 1.7: Euclidean one-hot with validation\n")
test_1_7 <- mdist_evolved(
  x = penguins_train,
  validate_x = penguins_test,
  preset = "euclidean_onehot"
)
print(dim(as.matrix(test_1_7$distance)))

# ============================================================================
# SECTION 2: CUSTOM PRESET TESTS - MIXED DATA
# ============================================================================

cat("\n\n=== SECTION 2: CUSTOM PRESET - MIXED DATA ===\n")

# Test 2.1: Custom - Manhattan + Matching
cat("\nTest 2.1: Custom - Manhattan + Matching\n")
test_2_1 <- mdist_evolved(
  x = penguins_train,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = FALSE,
  scaling_cont = "none"
)
print(summary(as.vector(test_2_1$distance)))

# Test 2.2: Custom - Euclidean + Total Variation
cat("\nTest 2.2: Custom - Euclidean + Total Variation\n")
test_2_2 <- mdist_evolved(
  x = penguins_train,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "tot_var_dist",
  commensurable = TRUE,
  scaling_cont = "std"
)
print(summary(as.vector(test_2_2$distance)))

# Test 2.3: Custom - Manhattan with standardization
cat("\nTest 2.3: Custom - Manhattan with standardization\n")
test_2_3 <- mdist_evolved(
  x = penguins_train,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "tot_var_dist",
  commensurable = FALSE,
  scaling_cont = "std"
)
print(summary(as.vector(test_2_3$distance)))

# Test 2.4: Custom - with PC scores scaling
cat("\nTest 2.4: Custom - PC scores scaling\n")
test_2_4 <- mdist_evolved(
  x = penguins_train,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "matching",
  scaling_cont = "pc_scores",
  ncomp = 2
)
print(summary(as.vector(test_2_4$distance)))

# Test 2.5: Custom - Euclidean + HLeucl combination
cat("\nTest 2.5: Custom - Euclidean + HLeucl\n")
test_2_5 <- mdist_evolved(
  x = penguins_train,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "HLeucl",
  commensurable = FALSE,
  scaling_cont = "std"
)
print(summary(as.vector(test_2_5$distance)))

# Test 2.6: Custom with validation and response
cat("\nTest 2.6: Custom with validation and response variable\n")
test_2_6 <- mdist_evolved(
  x = penguins_train,
  validate_x = penguins_test,
  response = "species",
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "tot_var_dist",
  commensurable = TRUE,
  scaling_cont = "pc_scores"
)
print(dim(as.matrix(test_2_6$distance)))

# ============================================================================
# SECTION 3: CATEGORICAL DATA ONLY
# ============================================================================

cat("\n\n=== SECTION 3: CATEGORICAL DATA ONLY ===\n")

# Prepare categorical-only data
penguins_cat <- penguins_complete %>% 
  select(species, island, sex)

penguins_cat_train <- penguins_cat[train_idx, ]
penguins_cat_test <- penguins_cat[-train_idx, ]

# Test 3.1: Gower with categorical only
cat("\nTest 3.1: Gower preset - categorical only\n")
test_3_1 <- mdist_evolved(
  x = penguins_cat_train,
  preset = "gower",
  commensurable = FALSE
)
print(summary(as.vector(test_3_1$distance)))

# Test 3.2: Unbiased dependent - categorical only
cat("\nTest 3.2: Unbiased dependent - categorical only\n")
test_3_2 <- mdist_evolved(
  x = penguins_cat_train,
  preset = "unbiased_dependent"
)
print(summary(as.vector(test_3_2$distance)))

# Test 3.3: Euclidean one-hot - categorical only
cat("\nTest 3.3: Euclidean one-hot - categorical only\n")
test_3_3 <- mdist_evolved(
  x = penguins_cat_train,
  preset = "euclidean_onehot"
)
print(summary(as.vector(test_3_3$distance)))

# Test 3.4: Custom matching distance - categorical only
cat("\nTest 3.4: Custom matching - categorical only\n")
test_3_4 <- mdist_evolved(
  x = penguins_cat_train,
  preset = "custom",
  distance_cat = "matching",
  commensurable = FALSE
)
print(summary(as.vector(test_3_4$distance)))

# Test 3.5: Custom total variation - categorical only with validation
cat("\nTest 3.5: Custom total variation with validation\n")
test_3_5 <- mdist_evolved(
  x = penguins_cat_train,
  validate_x = penguins_cat_test,
  preset = "custom",
  distance_cat = "tot_var_dist",
  commensurable = TRUE
)
print(dim(as.matrix(test_3_5$distance)))

# ============================================================================
# SECTION 4: CONTINUOUS DATA ONLY
# ============================================================================

cat("\n\n=== SECTION 4: CONTINUOUS DATA ONLY ===\n")

# Prepare continuous-only data
penguins_cont <- penguins_complete %>% 
  select(bill_length_mm, bill_depth_mm, flipper_length_mm, body_mass_g)

penguins_cont_train <- penguins_cont[train_idx, ]
penguins_cont_test <- penguins_cont[-train_idx, ]

# Test 4.1: Gower - continuous only
cat("\nTest 4.1: Gower preset - continuous only\n")
test_4_1 <- mdist_evolved(
  x = penguins_cont_train,
  preset = "gower",
  commensurable = FALSE
)
print(summary(as.vector(test_4_1$distance)))

# Test 4.2: Unbiased dependent - continuous only
cat("\nTest 4.2: Unbiased dependent - continuous only\n")
test_4_2 <- mdist_evolved(
  x = penguins_cont_train,
  preset = "unbiased_dependent"
)
print(summary(as.vector(test_4_2$distance)))

# Test 4.3: Euclidean one-hot - continuous only
cat("\nTest 4.3: Euclidean one-hot - continuous only\n")
test_4_3 <- mdist_evolved(
  x = penguins_cont_train,
  preset = "euclidean_onehot"
)
print(summary(as.vector(test_4_3$distance)))

# Test 4.4: Custom Manhattan - continuous only
cat("\nTest 4.4: Custom Manhattan - continuous only\n")
test_4_4 <- mdist_evolved(
  x = penguins_cont_train,
  preset = "custom",
  distance_cont = "manhattan",
  scaling_cont = "none"
)
print(summary(as.vector(test_4_4$distance)))

# Test 4.5: Custom Euclidean with standardization
cat("\nTest 4.5: Custom Euclidean with standardization\n")
test_4_5 <- mdist_evolved(
  x = penguins_cont_train,
  preset = "custom",
  distance_cont = "euclidean",
  scaling_cont = "std",
  commensurable = FALSE
)
print(summary(as.vector(test_4_5$distance)))

# Test 4.6: Custom with PC scores and ncomp
cat("\nTest 4.6: Custom PC scores with ncomp=2\n")
test_4_6 <- mdist_evolved(
  x = penguins_cont_train,
  preset = "custom",
  distance_cont = "manhattan",
  scaling_cont = "pc_scores",
  ncomp = 2
)
print(summary(as.vector(test_4_6$distance)))

# Test 4.7: Custom with threshold parameter
cat("\nTest 4.7: Custom with threshold parameter\n")
test_4_7 <- mdist_evolved(
  x = penguins_cont_train,
  preset = "custom",
  distance_cont = "euclidean",
  scaling_cont = "pc_scores",
  threshold = 0.95
)
print(summary(as.vector(test_4_7$distance)))

# Test 4.8: Custom with validation set
cat("\nTest 4.8: Custom with validation - continuous only\n")
test_4_8 <- mdist_evolved(
  x = penguins_cont_train,
  validate_x = penguins_cont_test,
  preset = "custom",
  distance_cont = "euclidean",
  scaling_cont = "std",
  commensurable = FALSE
)
print(dim(as.matrix(test_4_8$distance)))

# ============================================================================
# SECTION 5: EXTERNAL METHODS (if available)
# ============================================================================

cat("\n\n=== SECTION 5: EXTERNAL METHODS ===\n")

# Test 5.1: GUDMM method (if functions available)
cat("\nTest 5.1: GUDMM preset\n")
tryCatch({
  test_5_1 <- mdist_evolved(
    x = penguins_train,
    preset = "gudmm"
  )
  print(summary(as.vector(test_5_1$distance)))
}, error = function(e) {
  cat("GUDMM not available or error:", e$message, "\n")
})

# Test 5.2: DKSS method (if kdml package available)
cat("\nTest 5.2: DKSS preset\n")
tryCatch({
  test_5_2 <- mdist_evolved(
    x = penguins_train,
    preset = "dkss"
  )
  print(summary(as.vector(test_5_2$distance)))
}, error = function(e) {
  cat("DKSS not available or error:", e$message, "\n")
})

# Test 5.3: Modified Gower method
cat("\nTest 5.3: Modified Gower preset\n")
tryCatch({
  test_5_3 <- mdist_evolved(
    x = penguins_train,
    preset = "mod_gower"
  )
  print(summary(as.vector(test_5_3$distance)))
}, error = function(e) {
  cat("Modified Gower not available or error:", e$message, "\n")
})

# ============================================================================
# SECTION 6: EDGE CASES AND WARNINGS
# ============================================================================

cat("\n\n=== SECTION 6: EDGE CASES AND WARNINGS ===\n")

# Test 6.1: Single categorical variable with tot_var_dist (should warn)
cat("\nTest 6.1: Single categorical variable with tot_var_dist\n")
penguins_single_cat <- penguins_complete %>% select(species, bill_length_mm)
test_6_1 <- mdist_evolved(
  x = penguins_single_cat,
  preset = "custom",
  distance_cat = "tot_var_dist",
  distance_cont = "euclidean"
)

# Test 6.2: Single continuous variable with pc_scores (should warn)
cat("\nTest 6.2: Single continuous variable with pc_scores\n")
penguins_single_cont <- penguins_complete %>% select(species, island, bill_length_mm)
test_6_2 <- mdist_evolved(
  x = penguins_single_cont,
  preset = "custom",
  scaling_cont = "pc_scores",
  distance_cat = "matching"
)

# Test 6.3: Commensurable distances
cat("\nTest 6.3: Testing commensurable option\n")
test_6_3a <- mdist_evolved(
  x = penguins_train,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = TRUE,
  scaling_cont = "std"
)
test_6_3b <- mdist_evolved(
  x = penguins_train,
  preset = "custom",
  distance_cont = "manhattan",
  distance_cat = "matching",
  commensurable = FALSE,
  scaling_cont = "std"
)
cat("Commensurable TRUE:", mean(as.vector(test_6_3a$distance)), "\n")
cat("Commensurable FALSE:", mean(as.vector(test_6_3b$distance)), "\n")

# ============================================================================
# SECTION 7: RETURN OBJECT INSPECTION
# ============================================================================

cat("\n\n=== SECTION 7: RETURN OBJECT INSPECTION ===\n")

test_obj <- mdist_evolved(
  x = penguins_train,
  validate_x = penguins_test,
  preset = "custom",
  distance_cont = "euclidean",
  distance_cat = "tot_var_dist",
  scaling_cont = "std",
  commensurable = TRUE
)

cat("\nObject class:", class(test_obj), "\n")
cat("Distance matrix dimensions:", dim(as.matrix(test_obj$distance)), "\n")
cat("Preset used:", test_obj$preset, "\n")
cat("\nParameters:\n")
print(test_obj$params)

# ============================================================================
# SECTION 8: PERFORMANCE COMPARISON
# ============================================================================

cat("\n\n=== SECTION 8: PERFORMANCE COMPARISON ===\n")

cat("\nComparing computation times for different presets:\n")

time_gower <- system.time({
  mdist_evolved(x = penguins_train, preset = "gower")
})

time_unbiased <- system.time({
  mdist_evolved(x = penguins_train, preset = "unbiased_dependent")
})

time_euclidean <- system.time({
  mdist_evolved(x = penguins_train, preset = "euclidean_onehot")
})

time_custom <- system.time({
  mdist_evolved(x = penguins_train, preset = "custom", 
                distance_cont = "manhattan", distance_cat = "matching")
})

cat("Gower:", time_gower["elapsed"], "seconds\n")
cat("Unbiased dependent:", time_unbiased["elapsed"], "seconds\n")
cat("Euclidean one-hot:", time_euclidean["elapsed"], "seconds\n")
cat("Custom:", time_custom["elapsed"], "seconds\n")

# ============================================================================
# SECTION 9: DISTANCE MATRIX PROPERTIES
# ============================================================================

cat("\n\n=== SECTION 9: DISTANCE MATRIX PROPERTIES ===\n")

test_props <- mdist_evolved(
  x = penguins_train,
  preset = "gower"
)

dist_mat <- as.matrix(test_props$distance)

cat("\nChecking distance matrix properties:\n")
cat("Symmetric:", isSymmetric(dist_mat), "\n")
cat("Non-negative:", all(dist_mat >= 0), "\n")
cat("Zero diagonal:", all(diag(dist_mat) == 0), "\n")
cat("Min distance:", min(dist_mat[dist_mat > 0]), "\n")
cat("Max distance:", max(dist_mat), "\n")
cat("Mean distance:", mean(dist_mat[upper.tri(dist_mat)]), "\n")

cat("\n=== ALL TESTS COMPLETED ===\n")
