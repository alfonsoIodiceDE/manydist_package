library(targets)

source(file.path("R", "paths.R"))
source(file.path("R", "data.R"))
source(file.path("R", "specifications.R"))
source(file.path("R", "analysis.R"))
source(file.path("R", "tables.R"))
source(file.path("R", "replication.R"))

tar_option_set(
  packages = c(
    "manydist",
    "palmerpenguins",
    "tibble"
  )
)

list(
  tar_target(
    replication_manifest,
    run_replication(paper_root()),
    format = "file"
  )
)

