#!/bin/Rscript

stringr::str_glue(
"Usage:
  start_separate_countries_batch.R [<countries> | --exclude-us]
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

countries <- if (!is_empty(script_options$countries)) {
  script_options$countries
} else {
  subnat_data <- read_rds(file.path(root_path, "data", "mobility", "cleaned_subnat_data.rds"))

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

batchcmd <- str_glue("sbatch --array={countries} separate_countries_slurm.sh")

cat("Running:", batchcmd, "\n")

system(batchcmd)

