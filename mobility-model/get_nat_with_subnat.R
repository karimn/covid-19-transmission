#!/bin/Rscript

library(magrittr)
library(tidyverse)

if (interactive()) {
  root_path <- "."
} else {
  root_path <- ".."
}

subnat_data <- read_rds(file.path(root_path, "data", "mobility", "cleaned_subnat_data.rds"))

subnat_data %>%
  filter(has_epidemic) %>%
  semi_join(count(., country_code) %>% filter(n > 1), by = "country_code") %>%
  pull(country_index) %>%
  unique() %>%
  str_c(collapse = ",") %>%
  cat()
