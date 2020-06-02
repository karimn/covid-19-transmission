#!/usr/bin/Rscript

library(magrittr)
library(tidyverse)
library(rstan)
library(wpp2019)

options(mc.cores = max(1, parallel::detectCores()))
rstan_options(auto_write = TRUE)

# Constants ---------------------------------------------------------------

min_deaths_day_before_epidemic <- 10
seeding_days_before_epidemic <- 30
start_epidemic_offset <- seeding_days_before_epidemic + 1
days_seeding <- 6
days_to_forecast <- 0
mob_formula <- ~ 0 + g_residential + g_transit_stations + average_mob

# Population Data For IFR -------------------------------------------------

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

# wb_pop_data <- read_csv(file.path("data", "population", "wb_subnational_pop.csv")) %>%
#   filter(fct_match(`Series Code`, "SP.POP.TOTL")) %>%
#   tidyr::extract(`Country Code`, c("countrycode_string", "sub_region"), "^(\\w{3})_[^\\.]+\\.(\\w{2})", remove = FALSE) %>%
#   mutate(
#     # For some subnational regions we don't have their code in the below infection/deaths data
#     sub_region = if_else(fct_match(countrycode_string, c("IND", "CHN")), str_extract(`Country Name`, "(?<=,\\s).+"), sub_region),
#     sub_region = str_remove(sub_region, "\\s+(Sheng|Shi)$")
#   ) %>%
#   select(-`Series Name`, -`Series Code`, -Level_attr, -`Country Name`, -`Country Code`) %>%
#   filter(!is.na(countrycode_string), !is.na(sub_region)) %>%
#   pivot_longer(-c(countrycode_string, sub_region), names_to = "year", names_pattern = "(^\\d{4})", values_to = "pop") %>%
#   group_by(countrycode_string, sub_region) %>%
#   filter(min_rank(year) == n()) %>% # Just the last year available
#   ungroup() %>%
#   select(-year)

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
  select(-old_age_group, -name) %>%
  rename(countrycode = country_code) %>%
  group_by(countrycode) %>%
  summarize(ifr = sum(ifr * pop) / sum(pop)) %>%
  ungroup()

# Subnational Mobility Data -----------------------------------------------

# subnat_mob_data <- haven::read_dta(file.path("~", "Dropbox", "COVID-19", "Data Collection", "Travel", "clean", "google_trends.dta")) %>%
#   select(countrycode_string, sub_region, date, starts_with("g_")) %>%
#   mutate_at(vars(starts_with("g_")), ~ . / 100) %>%
#   group_nest(countrycode_string, sub_region, .key = "daily_data") %>%
#   mutate(
#     first_day_mob = map(daily_data, pull, date) %>% map_dbl(min),
#     last_day_mob = map(daily_data, pull, date) %>% map_dbl(max),
#   ) %>%
#   mutate_at(vars(first_day_mob, last_day_mob), lubridate::as_date) %>%
#   mutate(
#     daily_data = pmap(lst(daily_data, first_day_mob, last_day_mob), ~ complete(..1, date = seq.Date(from = ..2, to = ..3, by = "day"))),
#     # Try to impute mobility. Doesn't work if variable starts with NAs
#     daily_data = map(daily_data, ~ tryCatch(mutate_at(.x, vars(starts_with("g_")), zoo::na.locf), error = function(err) .x)) %>%
#       map(mutate, average_mob = (g_grocery_pharma + g_parks + g_retail_recreation + g_workplaces) / 4)
#   )

# Subnational Deaths Data -------------------------------------------------

subnat_data_raw <- read_csv(file.path("data", "subnational.csv"))

subnat_data <- subnat_data_raw %>%
  rename(subnat_day_index = X1) %>%
  left_join(transmute(countrycode::codelist, countrycode_string = iso3c, countrycode = iso3n), by = "countrycode_string") %>% # Getting the countrycode
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
        new_deaths = max(cum_deaths - lag(cum_deaths, default = 0), 0) # TODO remove the max() part when data cleaning is fixed
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

    # daily_data = map2(daily_data, first_observed_death, clean_missing_spread),

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
  ) %>%
  left_join(ifr_adj, by = "countrycode")

use_subnat_data <- subnat_data %>%
  filter(is_valid,
         fct_match(countryname, "Italy"))

# Time to Death -----------------------------------------------------------

quant_time_to_death <- ecdf(
  EnvStats::rgammaAlt(1e6, 5.1, 0.86) + # infection-to-onset distribution
    EnvStats::rgammaAlt(1e6, 18.8, 0.45) # onset-to-death distribution
)

time_to_death <- c(0, seq(1.5, by = 1, length.out = max(use_subnat_data$num_days_observed) + days_to_forecast)) %>%
  map_dbl(quant_time_to_death) %>%
  subtract(lag(.)) %>%
  discard(is.na)

# Plot timelines ----------------------------------------------------------

use_subnat_data %>%
  mutate(sub_region = fct_reorder(sub_region, epidemic_start_date)) %>%
  ggplot(aes(y = sub_region)) +
  geom_errorbar(aes(xmin = first_observed_death, xmax = last_observed_death), width = 0.25) +
  geom_linerange(aes(xmin = first_mob_day, xmax = last_mob_day), size = 7, alpha = 0.15) +
  geom_point(aes(first_infection_seeding_day)) +
  geom_point(aes(epidemic_start_date)) +
  geom_point(aes(last_effective_observed_day)) +
  ggrepel::geom_text_repel(aes(x = epidemic_start_date, label = "Epidemic"), nudge_y = 0.25, size = 3.5) +
  ggrepel::geom_text_repel(aes(x = first_infection_seeding_day, label = "First Seeding"), nudge_y = 0.25, size = 3.5) +
  ggrepel::geom_text_repel(aes(x = last_effective_observed_day, label = "Last Effective Observation"), nudge_y = 0.25, size = 3.5) +
  geom_linerange(aes(xmin = first_infection_seeding_day, xmax = epidemic_start_date - 1), linetype = "dashed") +
  scale_color_discrete("") +
  labs(x = "", y = "",
       caption = "Grey bars: range of mobility data. Black lines: range of infection/deaths data.") +
  facet_wrap(vars(countrycode_string), ncol = 1) +
  theme_minimal()

# Stan Data ---------------------------------------------------------------

# use_subnat_data %<>% slice(1:2) # Testing with one unit

stan_data <- lst(
  fit_model = 1,

  N_national = n_distinct(use_subnat_data$countrycode),
  N_subnational = use_subnat_data %>% count(countrycode) %>% pull(n) %>% as.array(),

  start_epidemic_offset,
  days_observed = as.integer(use_subnat_data$num_days_observed) %>% as.array(),
  days_to_impute_cases = days_seeding,
  days_to_forecast,
  total_days = max(days_observed),

  time_to_death = time_to_death[1:total_days],

  population = use_subnat_data$pop_ana %>%
    as.array(),
  mean_ifr = use_subnat_data$ifr %>%
    as.array(),

  design_matrix = use_subnat_data %>%
    unnest(daily_data) %>%
    modelr::model_matrix(mob_formula),

  num_coef = ncol(design_matrix),

  deaths = use_subnat_data %>%
    unnest(daily_data) %>%
    pull(new_deaths),
)

cat(str_glue("Running model with {stan_data$N_national} countries and {stan_data$N_subnational} total subnational entities."))

# Sampling ----------------------------------------------------------------

mob_model <- stan_model(file.path("mobility-model", "mobility.stan"))

mob_fit <- mob_model %>%
  sampling(
    data = stan_data,
    chains = 8,
    iter = 8000,
    control = lst(
      adapt_delta = 0.9,
    )
  )

# Save --------------------------------------------------------------------

save(mob_fit, stan_data, use_subnat_data, file = file.path("data", "mobility", "results", "mob_results.RData"))
