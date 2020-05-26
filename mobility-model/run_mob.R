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


# Subnational Deaths Data -------------------------------------------------

subnat_deaths_data <- haven::read_dta(file.path("~", "Dropbox", "COVID-19", "Analysis", "Data", "subspread.dta")) %>%
  left_join(transmute(countrycode::codelist, countrycode_string = iso3c, countrycode = iso3n), by = "countrycode_string")

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

c(0, seq(1.5, by = 1, length.out = 20)) %>% # TODO figure out  correct
  map_dbl(quant_time_to_death) %>%
  subtract(lag(.)) %>%
  discard(is.na)

