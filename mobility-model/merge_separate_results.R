#!/bin/Rscript

"Usage:
  merge_separate_results <location> <destination> <diagnostics> [--lite --run-suffix=<suffix>]

Options:
  --run-suffix=<suffix>  Run suffix used in results files [default: mob]
" -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, "../data/mobility/results/clean_epi3deaths/ ../data/mobility/results/lite_merged_epi3deaths.rds ../temp/log/diag_61500659.tsv --lite")
} else {
  root_path <- ".."

  setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd("mobility-model")

  docopt::docopt(opt_desc)
}

library(magrittr)
library(tidyverse)

diagnostics_data <- read_tsv(script_options$diagnostics, col_names = c("job_id", "country_index", "divergent_trans", "max_rhat", "min_ess_bulk", "min_ess_tail"))

merged_results <- dir(script_options$location, pattern = str_glue("\\w{2}_\\d+_\\d+_{script_options$`run-suffix`}_results\\.rds$"), full.names = TRUE) %>%
  map_dfr(~ tibble(
    results_file = .x,
    results = list(read_rds(.x))
  )) %>%
  tidyr::extract(results_file, "job_id", "\\w{2}_\\d+_(\\d+)", remove = FALSE) %>%
  left_join(diagnostics_data, by = c("job_id", "country_index")) %>%
  unnest(results)

if (script_options$lite) {
  merged_results %<>%
    trim_iter_data()
}

write_rds(merged_results, script_options$destination)

