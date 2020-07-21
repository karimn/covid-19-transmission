#!/bin/Rscript

"Usage:
  merge_separate_results <location> <destination> <diagnostics> [--lite --run-suffix=<suffix>]

Options:
  --run-suffix=<suffix>  Run suffix used in results files [default: mob]
" -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, "../data/mobility/results/clean_epi3deaths/ ../data/mobility/results/lite_merged_epi3deaths.rds ../temp/log/diag_61500659.tsv --lite --run-suffix=epi3deaths_mob")
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

diagnostics_data <- read_tsv(script_options$diagnostics,
                             col_names = c("job_id", "country_index", "country_code", "divergent_trans", "low_bfmi", "max_rhat", "min_ess_bulk", "min_ess_tail", "iter1", "iter2", "iter3", "iter4", "status", "time", "country_name"),
                             col_types = cols(job_id = col_integer(), country_index = col_integer())) %>%
  select(-status, -country_name, -max_rhat, -starts_with("min_ess"))

merged_results <- dir(script_options$location, pattern = str_glue("\\w{{2}}_\\d+_\\d+_{script_options$`run-suffix`}_results\\.rds$"), full.names = TRUE) %>%
  map_dfr(~ tibble(
    results_file = .x,
    results = list(read_rds(.x))
  )) %>%
  tidyr::extract(results_file, "job_id", "\\w{2}_(\\d+)_\\d+", remove = FALSE, convert = TRUE) %>%
  unnest(results) %>%
  left_join(diagnostics_data, by = c("job_id", "country_index", "country_code"))

if (script_options$lite) {
  merged_results %<>%
    trim_iter_data()
}

write_rds(merged_results, script_options$destination)

