#!/bin/Rscript

root_path <- if (interactive()) "." else ".."

source(file.path(root_path, "mobility-model", "constants.R"))

stringr::str_glue("Usage:
  start_joint_countries_batch.R (fit | prior) [<country-code> ...] [options] [--top-deaths=<num-countries> | --batch=<batch-size>]

Options:
  --iter=<iterations>, -i <iterations>  Total number of iterations [default: 2000].
  --outputname=<name, -o <name>  Name to use as suffix for results [default: mob].
  --no-sbatch  Dry run.
  --exponential  Use exponential model for mobility.
  --singletons-only  Use countries with a single sub-region only.
  --exclude-countries  Exclude countries listed.
  --top-deaths=<num-countries>  Use only the top <num-countries> in terms for number of deaths.
  --batch=<batch-size>  Group countries in batched runs.
  --epidemic-cutoff=<num-deaths>  Number of cumulative deaths that defines the start of an epidemic [default: {min_deaths_day_before_epidemic}]
  --raw-data-file=<raw file>  Path to raw data [default: {file.path(root_path, 'data', 'mergecleaned.csv')}]
") -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, "fit --singletons-only --no-sbatch --batch=10")
} else {
  script_path <- setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd(script_path)

  docopt::docopt(opt_desc)
}

library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

script_options %<>%
  modify_at(c("epidemic-cutoff", "top-deaths", "batch"), as.integer) %>%
  modify_at(c("country-code"), str_to_upper)

source(file.path(root_path, "mobility-model", "mob_util.R"))

subnat_data <- if (is_empty(script_options$`epidemic-cutoff`)) {
  read_rds(file.path(root_path, "data", "mobility", "cleaned_subnat_data.rds"))
} else {
  prepare_subnat_data(script_options$`raw-data-file`, script_options$`epidemic-cutoff`)
}

countries <- subnat_data %>%
    filter(has_epidemic) %>% {
      filtered_countries <- if (!is_empty(script_options$`country-code`) && !script_options$`exclude-countries`) {
        filter(., fct_match(country_code, script_options$`country-code`))
      } else if (script_options$`singletons-only`) {
        semi_join(., count(., country_code) %>% filter(n == 1), by = "country_code")
      } else .

      filtered_countries <- if (!is_empty(script_options$`country-code`) && script_options$`exclude-countries`) {
        filter(filtered_countries, !fct_match(country_code, script_options$`country-code`))
      } else filtered_countries

      if (!is_empty(script_options$`top-deaths`)) {
        filtered_countries %>%
          filter(row_number(-total_deaths) <= script_options$`top-deaths`)
      } else filtered_countries
    } %>%
  distinct(country_index, country_code, country_name)

if (!is_empty(script_options$batch)) {
  countries %<>%
    mutate(batch_index = ((seq(n()) - 1) %% script_options$batch) + 1)
} else {
  countries %<>%
    mutate(batch_index = 1)
}

job_id <- NULL

run_type <- if (script_options$fit) "fit" else "prior"
iter <- as.integer(script_options$iter)
mob_model_type <- if (script_options$exponential) "exponential" else "inv_logit"

batchcmd_builder <- function(country_codes) {
  str_glue("sbatch --parsable joint_countries_slurm.sh {run_type} {script_options$outputname} {iter} {mob_model_type} {script_options$`epidemic-cutoff`} {str_c(country_codes, collapse = ' ')}")
}

batchcmd <- countries %>%
  group_by(batch_index) %>%
  group_map(~ .x$country_code) %>%
  map_chr(batchcmd_builder)

if (!script_options$`no-sbatch`) {
  cat("Running:\n", str_c(batchcmd, collapse = "\n "))

  job_id <- map_chr(batchcmd, system, intern = TRUE)

  cat("\nSubmitted job(s):\n", job_id, "\n")

  countries %<>%
    group_nest(batch_index) %>%
    mutate(job_id) %>%
    unnest(data)
} else {
  cat("Would have run:\n", str_c(batchcmd, collapse = "\n "))
}

write_tsv(countries, "country_dict.tsv")

