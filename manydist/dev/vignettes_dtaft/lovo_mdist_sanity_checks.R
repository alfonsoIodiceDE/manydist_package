library(dplyr)
library(purrr)
library(tibble)

load(file = "../manydist/vignettes/Dartpoints.RData")

#----------------------------
# 1. mdist(): quoted vs unquoted response
#----------------------------
md_sup_unquoted <- mdist(df, response = Name, distance_cat = "tvd")
md_sup_quoted   <- mdist(df, response = "Name", distance_cat = "tvd")

md_unsup_drop   <- mdist(df |> select(-Name), distance_cat = "tvd")

cat("mdist quoted vs unquoted response:\n")
print(all.equal(
  as.matrix(md_sup_unquoted$to_dist()),
  as.matrix(md_sup_quoted$to_dist())
))

#----------------------------
# 2. lovo_mdist(): quoted vs unquoted response
#----------------------------
lovo_sup_unquoted <- lovo_mdist(
  df,
  response = Name,
  distance_cat = "tvd",
  response_used = TRUE,
  commensurable = FALSE,
  keep_dist = TRUE
)

lovo_sup_quoted <- lovo_mdist(
  df,
  response = "Name",
  distance_cat = "tvd",
  response_used = TRUE,
  commensurable = FALSE,
  keep_dist = TRUE
)

cat("\nlovo_mdist quoted vs unquoted response:\n")
print(all.equal(lovo_sup_unquoted$results, lovo_sup_quoted$results))
print(all.equal(lovo_sup_unquoted$full_dist, lovo_sup_quoted$full_dist))

#----------------------------
# 3. response_used = FALSE should match explicit unsupervised call
#----------------------------
lovo_unsup_flag <- lovo_mdist(
  df,
  response = Name,
  distance_cat = "tvd",
  response_used = FALSE,
  commensurable = FALSE,
  keep_dist = TRUE
)

lovo_unsup_drop <- lovo_mdist(
  df |> select(-Name),
  distance_cat = "tvd",
  response_used = FALSE,
  commensurable = FALSE,
  keep_dist = TRUE
)

cat("\nresponse_used = FALSE equivalence check:\n")
print(all.equal(lovo_unsup_flag$results, lovo_unsup_drop$results))
print(all.equal(lovo_unsup_flag$full_dist, lovo_unsup_drop$full_dist))

#----------------------------
# 4. supervised vs unsupervised should generally differ
#----------------------------
cat("\nsupervised vs unsupervised difference check:\n")
print(all.equal(lovo_sup_unquoted$results, lovo_unsup_flag$results))
print(sum(abs(lovo_sup_unquoted$full_dist - lovo_unsup_flag$full_dist)))

#----------------------------
# 5. check clustering columns appear only when requested
#----------------------------
cat("\ncolumns without cluster_k:\n")
print(names(lovo_sup_unquoted$results))

lovo_sup_cluster <- lovo_mdist(
  df,
  response = Name,
  distance_cat = "tvd",
  response_used = TRUE,
  commensurable = FALSE,
  cluster_k = 3,
  keep_dist = FALSE
)

cat("\ncolumns with cluster_k:\n")
print(names(lovo_sup_cluster$results))

#----------------------------
# 6. summary table of checks
#----------------------------
check_tbl <- tibble(
  check = c(
    "mdist quoted = unquoted",
    "lovo quoted = unquoted (results)",
    "lovo quoted = unquoted (full_dist)",
    "response_used = FALSE matches dropped-response LOVO (results)",
    "response_used = FALSE matches dropped-response LOVO (full_dist)",
    "supervised differs from unsupervised"
  ),
  passed = c(
    isTRUE(all.equal(
      as.matrix(md_sup_unquoted$to_dist()),
      as.matrix(md_sup_quoted$to_dist())
    )),
    isTRUE(all.equal(lovo_sup_unquoted$results, lovo_sup_quoted$results)),
    isTRUE(all.equal(lovo_sup_unquoted$full_dist, lovo_sup_quoted$full_dist)),
    isTRUE(all.equal(lovo_unsup_flag$results, lovo_unsup_drop$results)),
    isTRUE(all.equal(lovo_unsup_flag$full_dist, lovo_unsup_drop$full_dist)),
    !isTRUE(all.equal(lovo_sup_unquoted$results, lovo_unsup_flag$results))
  )
)

check_tbl
