#!/usr/bin/Rscript

library(magrittr)
library(tidyverse)

"Usage:
  prepare_subnat_data <raw-data-file> <output-data-file>
" -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, "") # Add the files here if running interactively
} else {
  docopt::docopt(opt_desc)
}

source("constants.R")

# Cleanup -----------------------------------------------------------------

subnat_data_raw <- read_csv(script_options$`raw-data-file`)

subnat_data <- subnat_data_raw %>%
  rename(subnat_day_index = X1) %>%
  left_join(transmute(countrycode::codelist, countrycode_string = iso3c, countrycode = iso3n), by = "countrycode_string") %>% # Getting the countrycode
  mutate_at(vars(starts_with("g_")), ~ . / 100) %>%
  mutate(
    all_mob_observed = select(., starts_with("g_")) %>%
      map(is.na) %>%
      transpose() %>%
      map_lgl(~ all(unlist(.x))) %>%
      not(),
  )

temp_pop_data <- subnat_data %>%
  distinct_at(vars(countrycode, sub_region, starts_with("pop"))) %>%
  filter(!is.na(pop_ana))

clean_missing_spread <- function(daily_data, first_observed_death) {
  if (is.finite(first_observed_death)) {
    daily_data %>%
      mutate(
        cum_deaths = if_else(date < first_observed_death, 0, cum_deaths),
        new_deaths = pmax(cum_deaths - lag(cum_deaths, default = 0), 0) # TODO remove the max() part when data cleaning is fixed
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
  select(-starts_with("pop")) %>%
  group_by_at(vars(countrycode, countrycode_string, countryname, sub_region, starts_with("pop"))) %>%
  group_nest(.key = "daily_data") %>%
  left_join(temp_pop_data, by = c("countrycode", "sub_region")) %>%
  mutate(
    first_observed_death = map_dbl(daily_data,
                                   ~ filter(.x, !is.na(cum_deaths)) %>%
                                     pull(date) %>%
                                     min()) %>%
      lubridate::as_date(),

    last_observed_death = map_dbl(daily_data,
                                 ~ filter(.x, !is.na(cum_deaths)) %>%
                                   pull(date) %>%
                                   max()) %>%
      lubridate::as_date(),

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
    is_valid = is.finite(first_observed_death) & is.finite(first_mob_day) & is.finite(epidemic_start_date),

    daily_data = pmap(lst(daily_data, first_infection_seeding_day, last_effective_observed_day, is_valid), clean_date_ranges) %>%
      map2(first_observed_death, clean_missing_spread) %>%
      map(mutate_at, vars(starts_with("g_")), zoo::na.locf, na.rm = FALSE) %>% # Any non-leading NAs are replaced with last prior non-NA
      map(mutate_at, vars(starts_with("g_")), coalesce, 0) %>%  # Any leading NAs are replaced with zeroes (as in Vollmer et al.)
      map(mutate, average_mob = (g_grocery_pharma + g_parks + g_retail_recreation + g_workplaces) / 4),
  )

write_rds(subnat_data, path = script_options$`output-data-file`)
