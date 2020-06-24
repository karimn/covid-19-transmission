#!/bin/Rscript

stringr::str_glue(
"Usage:
  render_country_reports <merged-data-file> <report-dir> [<report-id>]
") -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, '')
} else {
  root_path <- ".."

  old_dir <- setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd(old_dir)

  docopt::docopt(opt_desc)
}

library(magrittr)
library(tidyverse)
library(wpp2019)

source(file.path(root_path, "mobility-model", "constants.R"))
source(file.path(root_path, "mobility-model", "mob_util.R"))

merged_data <- read_rds(script_options$`merged_data-file`)

merged_data %>%
  render_country_reports(file.path(root_path, "mobility-model", "mobility_report.Rmd"),
                         script_options$`report-dir`,
                         script_options$`report-id`)
