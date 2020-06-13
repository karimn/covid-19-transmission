#!/usr/bin/Rscript

"Usage:
  prepare_subnat_data <raw-data-file> <output-data-file>
" -> opt_desc


script_options <- if (interactive()) {
  root_path <- "."

  docopt::docopt(opt_desc, "data/cleaned.csv data/mobility/cleaned_subnat_data.csv") # Add the files here if running interactively
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

# Cleanup -----------------------------------------------------------------

subnat_data_raw <- read_csv(
  script_options$`raw-data-file`,
  col_types = cols("subregion1_code" = col_character(),
                   "subregion2_code" = col_character(),
                   "subregion2_name" = col_character())
)

subnat_data <- subnat_data_raw %>%
  rename(subnat_day_index = X1) %>%
  left_join(transmute(countrycode::codelist, country_code = iso2c, countrycode_iso3n = iso3n), by = "country_code") %>% # Getting the countrycode
  rename_at(vars(matches("^(mobility|total)_")), str_replace_all, c("^mobility_" = "g_", "^total_" = "cum_")) %>%
  rename_at(vars(matches("(confirmed|deceased)$")), str_replace_all, c("confirmed$" = "cases", "deceased$" = "deaths")) %>%
  mutate_at(vars(starts_with("g_")), ~ . / 100) %>%
  mutate(
    sub_region = case_when(aggregation_level == 1 ~ subregion1_name,
                           aggregation_level == 2 ~ subregion2_name),

    all_mob_observed = select(., starts_with("g_")) %>%
      map(is.na) %>%
      transpose() %>%
      map_lgl(~ all(unlist(.x))) %>%
      not(),
  )

clean_missing_spread <- function(daily_data, first_observed_death) {
  if (is.finite(first_observed_death)) {
    daily_data %>%
      mutate(
        cum_deaths = if_else(is.na(cum_deaths) & date < first_observed_death, 0, cum_deaths),
        new_deaths = if_else(is.na(new_deaths), cum_deaths - lag(cum_deaths, default = 0), new_deaths),
        cum_deaths = if_else(is.na(cum_deaths), cumsum(new_deaths), cum_deaths),
      )
  } else return(daily_data)
}

clean_date_ranges <- function(daily_data, first_infection_seeding_day, last_effective_observed_day, is_valid) {
  if (is_valid) {
    daily_data %>%
      filter(between(date, first_infection_seeding_day, last_effective_observed_day)) %>%
      complete(date = seq.Date(from = first_infection_seeding_day, to = max(date), by = "day"))
  } else {
    return(daily_data)
  }
}

subnat_data %<>%
  group_by_at(vars(starts_with("country"), starts_with("subregion"), sub_region, key, aggregation_level, population)) %>%
  group_nest(.key = "daily_data") %>%
  mutate(
    first_observed_death = map_dbl(daily_data,
                                   ~ filter(.x, !is.na(cum_deaths) | !is.na(new_deaths)) %>%
                                     pull(date) %>%
                                     min()) %>%
      lubridate::as_date(),

    last_observed_death = map_dbl(daily_data,
                                 ~ filter(.x, !is.na(cum_deaths) | !is.na(new_deaths)) %>%
                                   pull(date) %>%
                                   max()) %>%
      lubridate::as_date(),

    daily_data = map2(daily_data, first_observed_death, clean_missing_spread),

    epidemic_start_date = map_dbl(daily_data,
                                  ~ filter(.x, cum_deaths >= min_deaths_day_before_epidemic) %>% # First day with >= 10 deaths
                                    pull(date) %>%
                                    min()) %>%
      lubridate::as_date() %>%
      add(1), # The day after

    first_infection_seeding_day = epidemic_start_date - 1 - seeding_days_before_epidemic, # 30 days before the first day with >=10 deaths

    first_mob_day = map_dbl(daily_data,
                       ~ filter(.x, all_mob_observed) %>%
                         pull(date) %>%
                         min()) %>%
      lubridate::as_date(),

    last_mob_day = map_dbl(daily_data,
                       ~ filter(.x, all_mob_observed) %>%
                         pull(date) %>%
                         max()) %>%
      lubridate::as_date(),

    last_effective_observed_day = pmin(last_observed_death, last_mob_day), # We need both deaths and mobility day for the tail end of observed data (used in likelihood)
    num_days_observed = last_effective_observed_day - first_infection_seeding_day + 1,
  ) %>%
  mutate(
    is_valid = is.finite(first_observed_death) & is.finite(first_mob_day),
    has_epidemic = is_valid & is.finite(epidemic_start_date),

    daily_data = pmap(lst(daily_data, first_infection_seeding_day, last_effective_observed_day, is_valid = has_epidemic), clean_date_ranges) %>%
      map2(first_observed_death, clean_missing_spread) %>%
      map(mutate_at, vars(starts_with("g_")), zoo::na.locf, na.rm = FALSE) %>% # Any non-leading NAs are replaced with last prior non-NA
      map(mutate_at, vars(starts_with("g_")), coalesce, 0) %>%  # Any leading NAs are replaced with zeroes (as in Vollmer et al.)
      map(mutate, average_mob = (g_grocery_and_pharmacy + g_parks + g_retail_and_recreation + g_workplaces) / 4) %>%
      map(arrange, date) %>%
      map(mutate, day_index = seq(n()))
  )

write_rds(subnat_data, path = script_options$`output-data-file`)
