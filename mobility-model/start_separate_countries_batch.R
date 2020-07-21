#!/bin/Rscript

root_path <- if (interactive()) "." else ".."

source(file.path(root_path, "mobility-model", "constants.R"))

stringr::str_glue("Usage:
  start_separate_countries_batch.R (fit | prior) [<country-code> ... | --exclude-us --include-singletons] [options]

Options:
  --iter=<iterations>, -i <iterations>  Total number of iterations [default: 2000].
  --outputname=<name, -o <name>  Name to use as suffix for results [default: mob].
  --no-sbatch  Dry run.
  --exponential  Use exponential model for mobility.
  --epidemic-cutoff=<num-deaths>  Number of cumulative deaths that defines the start of an epidemic [default: {min_deaths_day_before_epidemic}]
  --raw-data-file=<raw file>  Path to raw data [default: {file.path(root_path, 'data', 'mergecleaned.csv')}]
") -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, "fit it --epidemic-cutoff=3")
} else {
  script_path <- setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd(script_path)

  docopt::docopt(opt_desc)
}

library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

script_options %<>%
  modify_at(c("epidemic-cutoff"), as.integer) %>%
  modify_at(c("country-code"), str_to_upper)

source(file.path(root_path, "mobility-model", "mob_util.R"))

subnat_data <- if (is_empty(script_options$`epidemic-cutoff`)) {
  read_rds(file.path(root_path, "data", "mobility", "cleaned_subnat_data.rds"))
} else {
  prepare_subnat_data(script_options$`raw-data-file`, script_options$`epidemic-cutoff`)
}

countries <- subnat_data %>%
    filter(has_epidemic) %>% {
      if (!is_empty(script_options$`country-code`)) {
        filter(., fct_match(country_code, script_options$`country-code`))
      } else {
        multi_region_countries <- if (!script_options$`include-singletons`) {
          semi_join(., count(., country_code) %>% filter(n > 1), by = "country_code")
        } else .

        if (script_options$`exclude-us`) {
          filter(multi_region_countries, !fct_match(country_code, "US"))
        } else multi_region_countries
      }
    } %>%
    pull(country_index) %>%
    unique()

job_country_dict <- subnat_data %>%
  filter(country_index %in% countries) %>%
  distinct(country_index, country_code, country_name)

job_id <- NULL

run_type <- if (script_options$fit) "fit" else "prior"
iter <- as.integer(script_options$iter)
mob_model_type <- if (script_options$exponential) "exponential" else "inv_logit"

batchcmd <- str_glue("sbatch --parsable --array={str_c(countries, collapse = ',')} separate_countries_slurm.sh {run_type} {script_options$outputname} {iter} {mob_model_type} {script_options$`epidemic-cutoff`}")

if (!script_options$`no-sbatch`) {

  cat("Running:", batchcmd, "\n")

  job_id <- system(batchcmd, intern = TRUE)

  cat("Submitted job", job_id, "\n")

  job_country_dict %<>%
    mutate(job_id)
} else {
  cat("Would have run: ", batchcmd, "\n")
}

write_tsv(job_country_dict, str_c(str_c("country_dict", job_id, sep = "_"), ".tsv"))

