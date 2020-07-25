#!/bin/Rscript

stringr::str_glue(
"Usage:
  render_multi_country_report.R <job-id>
") -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, '65038549')
} else {
  root_path <- ".."

  old_dir <- setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd(old_dir)

  docopt::docopt(opt_desc)
}

library(magrittr)
library(tidyverse)
library(rmarkdown)

source(file.path(root_path, "mobility-model", "constants.R"))
source(file.path(root_path, "mobility-model", "mob_util.R"))

results_path <- file.path(root_path, "temp")
job_id <- script_options$`job-id`
output_file <- str_glue("mobility_report_{job_id}.pdf")

rmarkdown::render(file.path(root_path, "mobility-model", "mobility_report.Rmd"),
                  output_file = output_file,
                  output_dir = results_path,
                  params = lst(
                    results_file = str_glue("multi_{job_id}_mob_results.rds"),
                    fit_file = str_glue("multi_{job_id}_mob.RData"),
                    subnat_hier = FALSE,
                    multi_singletons = TRUE,
                  ),
                  knit_root_dir = "..")

system(str_glue("evince {file.path(results_path, output_file)}"), wait = FALSE)

