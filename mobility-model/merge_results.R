#!/bin/Rscript

"Usage:
  merge_results <location> <destination> [--lite --results-file-pattern=<pattern>]

Options:
  --lite  Drop all sampling iterations
  --results-file-pattern=<pattern>  Pattern used to find results .rds files [default: _results\\.rds$]
" -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  # docopt::docopt(opt_desc, "../data/mobility/results/clean_epi3deaths/ ../data/mobility/results/lite_merged_epi3deaths.rds --lite")
  docopt::docopt(opt_desc, "/n/holyscratch01/kremer_lab/karimn/mob_results/run_65564444/ ../temp --lite")
} else {
  root_path <- ".."

  old_dir <- setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd(old_dir)

  docopt::docopt(opt_desc)
}

library(magrittr)
library(tidyverse)

source(file.path(root_path, "mobility-model", "mob_util.R"))

merged_results <- dir(script_options$location, pattern = script_options$`results-file-pattern`, full.names = TRUE) %>%
  map_dfr(~ read_rds(.x) %>% mutate(results_file = .x))

if (script_options$lite) {
  merged_results %<>%
    trim_iter_data()
}

write_rds(merged_results, script_options$destination)

