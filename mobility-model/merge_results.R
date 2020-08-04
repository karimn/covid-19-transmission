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

# diagnostics_data <- read_tsv(script_options$diagnostics,
#                              col_names = c("job_id", "country_index", "country_code", "divergent_trans", "low_bfmi", "max_rhat", "min_ess_bulk", "min_ess_tail", "iter1", "iter2", "iter3", "iter4", "status", "time", "country_name"),
#                              col_types = cols(job_id = col_integer(), country_index = col_integer())) %>%
#   select(-status, -country_name, -max_rhat, -starts_with("min_ess"))

merged_results <- dir(script_options$location, pattern = script_options$`results-file-pattern`, full.names = TRUE) %>%
  # map_dfr(~ tibble(
  #   results_file = .x,
  #   results = list(read_rds(.x))
  # )) %>%
  map_dfr(~ read_rds(.x) %>% mutate(results_file = .x))
  # tidyr::extract(results_file, "job_id", "\\w{2}_(\\d+)_\\d+", remove = FALSE, convert = TRUE) %>%
  # unnest(results)
  # left_join(diagnostics_data, by = c("job_id", "country_index", "country_code"))

if (script_options$lite) {
  merged_results %<>%
    trim_iter_data()
}

write_rds(merged_results, script_options$destination)

