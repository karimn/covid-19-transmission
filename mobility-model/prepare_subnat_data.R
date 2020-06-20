#!/usr/bin/Rscript

"Usage:
  prepare_subnat_data <raw-data-file> <output-data-file>
" -> opt_desc


script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, "data/mergecleaned.csv data/mobility/cleaned_subnat_data.rds") # Add the files here if running interactively
} else {
  root_path <- ".."

  setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd("mobility-model")

  docopt::docopt(opt_desc)
}

library(magrittr)
library(tidyverse)

source(file.path(root_path, "mobility-model", "constants.R"))
source(file.path(root_path, "util.R"))

# Cleanup -----------------------------------------------------------------

subnat_data <- prepare_subnat_data(script_options$`raw-data-file`, min_deaths_day_before_epidemic)

write_rds(subnat_data, path = script_options$`output-data-file`)
