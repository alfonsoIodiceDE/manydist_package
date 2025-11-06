pkgname <- "manydist"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "manydist-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('manydist')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cdist")
### * cdist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cdist
### Title: Calculation of Pairwise Distances for Categorical Data
### Aliases: cdist

### ** Examples

library(palmerpenguins)
library(rsample)

# Prepare data with complete cases for both categorical variables and response
complete_vars <- c("species", "island", "sex", "body_mass_g")
penguins_complete <- penguins[complete.cases(penguins[, complete_vars]), ]
penguins_cat <- penguins_complete[, c("species", "island", "sex")]
response <- penguins_complete$body_mass_g

# Create training-test split
set.seed(123)
penguins_split <- initial_split(penguins_cat, prop = 0.8)
tr_penguins <- training(penguins_split)
ts_penguins <- testing(penguins_split)
response_tr <- response[penguins_split$in_id]
response_ts <- response[-penguins_split$in_id]

# Basic usage
result <- cdist(tr_penguins)

# With validation data
val_result <- cdist(x = tr_penguins, 
                   validate_x = ts_penguins,
                   method = "tot_var_dist")
                   
# ...and commensurability
val_result_COMM <- cdist(x = tr_penguins, 
                   validate_x = ts_penguins,
                   method = "tot_var_dist",
                   commensurable = TRUE)

# Supervised distance with response variable
sup_result <- cdist(x = tr_penguins, 
                   response = response_tr,
                   method = "supervised")

# Supervised with validation data
sup_val_result <- cdist(x = tr_penguins,
                       validate_x = ts_penguins,
                       response = response_tr,
                       method = "supervised")

# Commensurable distances with custom weights
comm_result <- cdist(tr_penguins,
                    commensurable = TRUE,
                    weights = c(2, 1, 1))

# Different methods per variable
multi_method <- cdist(tr_penguins,
                     method = c("matching", "goodall_3", "tot_var_dist"))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cdist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fifa_nl")
### * fifa_nl

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fifa_nl
### Title: FIFA 21 Player Data - Dutch League
### Aliases: fifa_nl
### Keywords: datasets

### ** Examples

data(fifa_nl)
summary(fifa_nl)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fifa_nl", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hello")
### * hello

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hello
### Title: Hello, World!
### Aliases: hello

### ** Examples

hello()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hello", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mdist")
### * mdist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mdist
### Title: Calculation of Pairwise Distances for Mixed-Type Data
### Aliases: mdist

### ** Examples

library(palmerpenguins)
library(rsample)

# Prepare complete data
pengmix <- palmerpenguins::penguins[complete.cases(palmerpenguins::penguins), ]

# Create training-test split
set.seed(123)
pengmix_split <- initial_split(pengmix, prop = 0.8)
tr_pengmix <- training(pengmix_split)
ts_pengmix <- testing(pengmix_split)

# Example 1: Basic usage with validation data
dist_matrix <- mdist(x = tr_pengmix, 
                    validate_x = ts_pengmix)

# Example 2: Gower preset with validation
dist_gower <- mdist(x = tr_pengmix, 
                   validate_x = ts_pengmix,
                   preset = "gower", 
                   commensurable = TRUE)

# Example 3: Euclidean one-hot preset with validation
dist_onehot <- mdist(x = tr_pengmix, 
                    validate_x = ts_pengmix,
                    preset = "euclidean_onehot")

# Example 4: Custom preset with standardization
dist_custom <- mdist(x = tr_pengmix,
                    validate_x = ts_pengmix,
                    preset = "custom",
                    distance_cont = "manhattan",
                    distance_cat = "matching",
                    commensurable = TRUE,
                    scaling_cont = "std")

# Example 5: PCA-based scaling with threshold
dist_pca <- mdist(x = tr_pengmix,
                 validate_x = ts_pengmix,
                 distance_cont = "euclidean",
                 scaling_cont = "pc_scores",
                 threshold = 0.85)

# Example 6: Categorical variables only
cat_vars <- c("species", "island", "sex")
dist_cat <- mdist(tr_pengmix[, cat_vars],
                 validate_x = ts_pengmix[, cat_vars],
                 distance_cat = "tot_var_dist")

# Example 7: Continuous variables only
num_vars <- c("bill_length_mm", "bill_depth_mm", 
              "flipper_length_mm", "body_mass_g")
dist_cont <- mdist(tr_pengmix[, num_vars],
                  validate_x = ts_pengmix[, num_vars],
                  distance_cont = "manhattan",
                  scaling_cont = "std")

# Example 8: Supervised distance with response
response_tr <- tr_pengmix$body_mass_g
dist_sup <- mdist(tr_pengmix,
                 validate_x = ts_pengmix,
                 response = response_tr,
                 distance_cat = "supervised")





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mdist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ndist")
### * ndist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ndist
### Title: Calculation of Pairwise Distances for Continuous Data
### Aliases: ndist

### ** Examples

library(palmerpenguins)
library(rsample)

penguins_cont <- palmerpenguins::penguins[, c("bill_length_mm",
"bill_depth_mm", "flipper_length_mm", "body_mass_g")]
penguins_cont <- penguins_cont[complete.cases(penguins_cont), ]

# Basic usage
dist_matrix <- ndist(penguins_cont)

# Commensurable distances with standardization
dist_matrix <- ndist(penguins_cont, 
                    commensurable = TRUE, 
                    scaling = "std")

# PCA-based dimensionality reduction
dist_matrix <- ndist(penguins_cont, 
                    scaling = "pc_scores", 
                    threshold = 0.95)

# Mahalanobis distance
dist_matrix <- ndist(penguins_cont, 
                    method = "mahalanobis")

# Weighted Euclidean distance
dist_matrix <- ndist(penguins_cont, 
                    method = "euclidean",
                    weights = c(1, 0.5, 2, 1))
                    
# Training-test split example with validation data
set.seed(123)
# Create training-test split using rsample
penguins_split <- initial_split(penguins_cont, prop = 0.8)
tr_penguins <- training(penguins_split)
ts_penguins <- testing(penguins_split)

# Basic usage with training data only
dist_matrix <- ndist(tr_penguins)

# Computing distances between test and training sets
val_dist_matrix <- ndist(x = tr_penguins, 
                        validate_x = ts_penguins,
                        method = "euclidean")

# Using validation data with standardization
val_dist_matrix_std <- ndist(x = tr_penguins,
                            validate_x = ts_penguins,
                            scaling = "std",
                            method = "manhattan")

# Validation with PCA and commensurability
val_dist_matrix_pca <- ndist(x = tr_penguins,
                            validate_x = ts_penguins,
                            scaling = "pc_scores",
                            ncomp = 2,
                            commensurable = TRUE)

# Validation with robust scaling and custom weights
val_dist_matrix_robust <- ndist(x = tr_penguins,
                               validate_x = ts_penguins,
                               scaling = "robust",
                               weights = c(1, 0.5, 2, 1))

# Mahalanobis distance with validation data
val_dist_matrix_mahal <- ndist(x = tr_penguins,
                              validate_x = ts_penguins,
                              method = "mahalanobis")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ndist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
