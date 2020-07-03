#!/bin/Rscript

"Usage:
  start_separate_countries_batch.R (fit | prior) [<country-code> ... | --exclude-us] [--no-sbatch --outputname=<name> --iter=<iterations> --exponential]

Options:
  --iter=<iterations>, -i <iterations>  Total number of iterations [default: 2000]
  --outputname=<name, -o <name>  Name to use as suffix for results [default: mob].
" -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, "fit it")
} else {
  root_path <- ".."

  script_path <- setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd(script_path)

  docopt::docopt(opt_desc)
}

library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)

script_options %<>%
  modify_at(c("country-code"), str_to_upper)

subnat_data <- read_rds(file.path(root_path, "data", "mobility", "cleaned_subnat_data.rds"))

countries <- subnat_data %>%
    filter(has_epidemic) %>%
    semi_join(count(., country_code) %>% filter(n > 1), by = "country_code") %>% {
      if (script_options$`exclude-us`) {
        filter(., !fct_match(country_code, "US"))
      } else if (!is_empty(script_options$`country-code`)) {
        filter(., fct_match(country_code, script_options$`country-code`))
      } else .
    } %>%
    pull(country_index) %>%
    unique() %>%
    str_c(collapse = ",")

job_country_dict <- subnat_data %>%
  filter(country_index %in% as.integer(c(str_split(countries, ",", simplify = TRUE)))) %>%
  distinct(country_index, country_code, country_name)

job_id <- NULL

if (!script_options$`no-sbatch`) {
  run_type <- if (script_options$fit) "fit" else "prior"
  iter <- as.integer(script_options$iter)
  mob_model_type <- if (script_options$exponential) "exponential" else "inv_logit"

  batchcmd <- str_glue("sbatch --parsable --array={countries} separate_countries_slurm.sh {run_type} {script_options$outputname} {iter} {mob_model_type}")

  cat("Running:", batchcmd, "\n")

  job_id <- system(batchcmd, intern = TRUE)

  cat("Submitted job", job_id, "\n")

  job_country_dict %<>%
    mutate(job_id)
}

write_rds(job_country_dict, str_c(str_c("country_dict", job_id, sep = "_"), ".rds"))

