#!/usr/bin/Rscript

library(magrittr)
library(tidyverse)
library(rstan)
library(wpp2019)

# Population Data --------------------------------------------------------

data(pop)

age_conv <- list(
  "[0,10)" = c("0-4", "5-9"),
  "[10,20)" = c("10-14", "15-19"),
  "[20,30)" = c("20-24", "25-29"),
  "[30,40)" = c("30-34", "35-39"),
  "[40,50)" = c("40-44", "45-49"),
  "[50,60)" = c("50-54", "55-59"),
  "[60,70)" = c("60-64", "65-69"),
  "[70,80)" = c("70-74", "75-79"),
  "80+" =   c("80-84", "85-89", "90-94", "95-99", "100+")
)

# This is age stratified population used in calculating IFR. Subnational population figures are not available.

pop_data <- bind_rows(
  female = popF,
  male = popM,
  .id = "sex"
) %>%
  rename(age_group = age) %>%
  pivot_longer(cols = `1950`:`2020`, names_to = "year", values_to = "pop") %>%
  mutate(
    age_group = fct_collapse(age_group, !!!age_conv) %>%
      factor(levels = names(age_conv))
  ) %>%
  group_by(country_code, name, age_group, year) %>%
  summarise(pop = sum(pop * 1000)) %>%
  ungroup()

# Trying to get some subnational population sizes for some of the subnational entities in the main infection/deaths data.
# So far really only using Brazil, India, and China.

wb_pop_data <- read_csv(file.path("data", "population", "wb_subnational_pop.csv")) %>%
  filter(fct_match(`Series Code`, "SP.POP.TOTL")) %>%
  tidyr::extract(`Country Code`, c("countrycode_string", "sub_region"), "^(\\w{3})_[^\\.]+\\.(\\w{2})", remove = FALSE) %>%
  mutate(
    # For some subnational regions we don't have their code in the below infection/deaths data
    sub_region = if_else(fct_match(countrycode_string, c("IND", "CHN")), str_extract(`Country Name`, "(?<=,\\s).+"), sub_region),
    sub_region = str_remove(sub_region, "\\s+(Sheng|Shi)$")
  ) %>%
  select(-`Series Name`, -`Series Code`, -Level_attr, -`Country Name`, -`Country Code`) %>%
  filter(!is.na(countrycode_string), !is.na(sub_region)) %>%
  pivot_longer(-c(countrycode_string, sub_region), names_to = "year", names_pattern = "(^\\d{4})", values_to = "pop") %>%
  group_by(countrycode_string, sub_region) %>%
  filter(min_rank(year) == n()) %>% # Just the last year available
  ungroup() %>%
  select(-year)


# Subnational Mobility Data -----------------------------------------------

subnat_mob_data <- haven::read_dta(file.path("~", "Dropbox", "COVID-19", "Data Collection", "Travel", "clean", "google_trends.dta")) %>%
  select(countrycode_string, sub_region, date, starts_with("g_")) %>%
  mutate(average_mob = (g_grocery_pharma + g_parks + g_retail_recreation + g_workplaces) / 4) %>%
  group_nest(countrycode_string, sub_region, .key = "daily_data") %>%
  mutate(
    first_day_mob = map(daily_data, pull, date) %>% map_dbl(min),
    last_day_mob = map(daily_data, pull, date) %>% map_dbl(max)
  ) %>%
  mutate_at(vars(first_day_mob, last_day_mob), lubridate::as_date)

# Subnational Deaths Data -------------------------------------------------

subnat_data <- haven::read_dta(file.path("~", "Dropbox", "COVID-19", "Analysis", "Data", "subspread.dta")) %>%
  left_join(transmute(countrycode::codelist, countrycode_string = iso3c, countrycode = iso3n), by = "countrycode_string") %>%
  group_nest(countrycode, countrycode_string, sub_region, .key = "daily_data") %>%
  inner_join(subnat_mob_data, by = c("countrycode_string", "sub_region"), suffix = c("_spread", "_mob")) %>%
  mutate(
    first_day_deaths = map(daily_data_spread, pull, date) %>% map_dbl(min, na.rm = TRUE) %>% lubridate::as_date(), # TODO There shouldn't be any NA dates
    daily_data_spread = pmap(lst(daily_data_spread, first_day_mob, last_day_mob), ~ filter(..1, between(date, ..2, ..3))),
    daily_data_mob = pmap(lst(daily_data_mob, first_day_deaths), ~ filter(..1, date >= ..2)),
    daily_data = map2(daily_data_spread, daily_data_mob, right_join, by = "date")
  ) %>%
  select(-daily_data_mob, -daily_data_spread) %>%
  filter(map_int(daily_data, nrow) > 0) %>% # Only entities with any observations
  mutate(
    first_day = map(daily_data, pull, date) %>% map_dbl(min, na.rm = TRUE) %>% lubridate::as_date(),
    last_day = map(daily_data, pull, date) %>% map_dbl(max, na.rm = TRUE) %>% lubridate::as_date(),
    daily_data = pmap(lst(daily_data, first_day, last_day), ~ complete(..1, date = seq.Date(from = ..2, to = ..3, by = "day"))),
    days_observed = map_int(daily_data, nrow),
    num_missing_new_deaths = map(daily_data, filter, is.na(new_deaths)) %>% map_int(nrow),
    num_missing_mob = map(daily_data, select, starts_with("g_")) %>% map(map, is.na) %>% map(reduce, ~ .x | .y) %>% map_int(sum)
    ) %>%
  arrange(countrycode, sub_region) %>%
  left_join(wb_pop_data, by = c("countrycode_string", "sub_region")) # Only available for regions (not all) in Brazil, India, and China

# IFR ---------------------------------------------------------------------

age_metadata <- age_conv %>%
  enframe(name = "age_group", value = "old_age_group") %>%
  mutate(
    age_group = factor(age_group),
    ifr = c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3) / 100
  )

ifr_adj <- pop_data %>%
  filter(fct_match(year, "2020")) %>%
  left_join(age_metadata, by = "age_group") %>%
  select(-old_age_group) %>%
  group_by(country_code, name) %>%
  summarize(ifr = sum(ifr * pop) / sum(pop)) %>%
  ungroup()

# Time to Death -----------------------------------------------------------

quant_time_to_death <- ecdf(
  EnvStats::rgammaAlt(1e6, 5.1, 0.86) + # infection-to-onset distribution
    EnvStats::rgammaAlt(1e6, 18.8, 0.45) # onset-to-death distribution
)

time_to_death <- c(0, seq(1.5, by = 1, length.out = max(subnat_deaths_data$days_observed))) %>%
  map_dbl(quant_time_to_death) %>%
  subtract(lag(.)) %>%
  discard(is.na)

# Stan Data ---------------------------------------------------------------

use_subnat_data <- subnat_data %>%
  filter(!is.na(pop), days_observed > 0) %>%
  arrange(countrycode)

stan_data <- lst(
  fit_model = 1,

  N_national = n_distinct(use_subnat_data$countrycode),
  N_subnational = use_subnat_data %>% count(countrycode) %>% pull(n),

  days_observed = use_subnat_data$days_observed,
  days_to_impute_cases = 6,
  total_days = max(days_observed),
)
