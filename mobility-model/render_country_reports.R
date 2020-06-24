#!/bin/Rscript

stringr::str_glue(
"Usage:
  render_country_reports <merged-data-file> <report-dir> [--report-id=<report-id>]
") -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, 'lite_merged.rds mobility-model/country-reports/')
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

merged_data <- read_rds(file.path(root_path, "data", "mobility", "results", script_options$`merged-data-file`))

merged_data %>%
  render_country_reports(file.path(root_path, "mobility-model", "mobility_report.Rmd"),
                         results_file = script_options$`merged-data-file`,
                         script_options$`report-dir`,
                         script_options$`report-id`)
