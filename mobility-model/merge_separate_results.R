#!/bin/Rscript

"Usage:
  merge_separate_results <location> <destination>
" -> opt_desc

library(magrittr)
library(tidyverse)

script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, "")
} else {
  root_path <- ".."

  setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd("mobility-model")

  docopt::docopt(opt_desc)
}

merged_results <- dir(script_options$location, pattern = "\\w{2}_\\d+_\\d+_mob_results\\.rds$", full.names = TRUE) %>%
  map_dfr(read_rds)

write_rds(merged_results, script_options$destination)

