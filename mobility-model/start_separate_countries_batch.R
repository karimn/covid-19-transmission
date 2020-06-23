#!/bin/Rscript

stringr::str_glue(
"Usage:
  start_separate_countries_batch.R [<countries> | --exclude-us] [--no-sbatch]
") -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, "")
} else {
  root_path <- ".."

  script_path <- setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd(script_path)

  docopt::docopt(opt_desc)
}

library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

subnat_data <- read_rds(file.path(root_path, "data", "mobility", "cleaned_subnat_data.rds"))

countries <- if (!is_empty(script_options$countries)) {
  script_options$countries
} else {
  subnat_data %>%
    filter(has_epidemic) %>%
    semi_join(count(., country_code) %>% filter(n > 1), by = "country_code") %>% {
      if (script_options$`exclude-us`) {
        filter(., !fct_match(country_code, "US"))
      } else .
    } %>%
    pull(country_index) %>%
    unique() %>%
    str_c(collapse = ",")
}

job_country_dict <- subnat_data %>%
  filter(country_index %in% as.integer(c(str_split(countries, ",", simplify = TRUE)))) %>%
  distinct(country_index, country_code, country_name)

job_id <- NULL

if (!script_options$`no-sbatch`) {
  batchcmd <- str_glue("sbatch --array={countries} separate_countries_slurm.sh")

  cat("Running:", batchcmd, "\n")

  job_id <- system(batchcmd, intern = TRUE)

  cat("Submitted job", job_id, "\n")

  job_country_dict %<>%
    mutate(job_id)
}

write_rds(job_country_dict, str_c(str_c("country_dict", job_id, sep = "_"), ".rds"))

