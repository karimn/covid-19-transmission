library(tidyverse)

# Consider India dataset
india_data <- readRDS("data/mobility/cleaned_subnat_data.rds") %>% 
  filter(countryname == "India") %>%
  unnest(cols = c("daily_data"))

head(india_data)

india_data %>% ggplot(aes(x = date, y= cum_deaths, color = sub_region)) + geom_line()
india_data %>% filter(cum_deaths > 0) %>% ggplot(aes(x = date, y= cum_deaths, color = sub_region)) + geom_line()


india_data %>% group_by(sub_region) %>% summarise(m = max(cum_deaths, na.rm=TRUE)) %>% pull(m) %>% table()
india_data %>% group_by(sub_region) %>% summarise(m = max(cum_deaths, na.rm=TRUE)) %>% View
# 23 states with no deaths, 5 are NAs, 11 have over 10 deaths





# Population Data For IFR -------------------------------------------------
library(wpp2019)
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





# Let's select 7 states in the south, one epicenter, 3 states with some deaths, 3 states below threshold
st <- c("Karnataka", "Telengana", "Tamil Nadu", "Kerala", "Goa", "Maharashtra", "Arunachal Pradesh")
st <- c("Karnataka", "Tamil Nadu", "Kerala", "Goa", "Maharashtra", "Arunachal Pradesh")


india_data %>% filter(sub_region %in% st) %>%
  group_by(sub_region) %>% summarise(m = max(cum_deaths, na.rm=TRUE)) 

use_subnat_data <- india_data %>% filter(sub_region %in% st) %>% filter(!is.na(cum_deaths)) %>%
  filter(date > "2020-03-01") %>%
  left_join(ifr_adj, by = "countrycode") %>%
  select(sub_region, date, average_mob, cum_deaths, new_deaths, pop_ana, ifr) %>%
  group_by(sub_region) %>% mutate(deaths = max(cum_deaths),
                                  num_days_observed = max(date) - min(date) + 1) %>% ungroup() %>%
  nest(data = c("date", "average_mob", "cum_deaths", "new_deaths")) %>%
  arrange(desc(deaths))




# Simplified data inputs
unnest(use_subnat_data) %>% 
  ggplot(aes(x=date, y=average_mob, color=sub_region)) + geom_line()
unnest(use_subnat_data) %>% 
  ggplot(aes(x=date, y=cum_deaths, color=sub_region)) + geom_line()
select(unnest(use_subnat_data), sub_region, date, average_mob) %>% 
  group_by(sub_region) %>% summarise(mean(average_mob))



# Extra inputs from Karim's script -----

# subnat_data <- readRDS("data/mobility/cleaned_subnat_data.rds")

days_to_forecast <- 0

time_to_death_draws <- EnvStats::rgammaAlt(1e6, 5.1, 0.86) + # infection-to-onset distribution
  EnvStats::rgammaAlt(1e6, 18.8, 0.45) # onset-to-death distribution

quant_time_to_death <- ecdf(time_to_death_draws)

time_to_death <- c(0, seq(1.5, by = 1, length.out = max(use_subnat_data$num_days_observed) + days_to_forecast)) %>%
  map_dbl(quant_time_to_death) %>%
  subtract(lag(.)) %>%
  discard(is.na)


# Run the model now:

days_seeding <- 6
min_deaths_day_before_epidemic <- 10
seeding_days_before_epidemic <- 30
start_epidemic_offset <- seeding_days_before_epidemic + 1

stan_data <- lst(
  # Settings
  fit_model = 1,
  hierarchical_mobility_model = 1,
  mobility_model_type = 1, 
  
  # Hyperparameters
  hyperparam_tau_beta_toplevel = 0.5,
  hyperparam_tau_beta_national_sd = 1,
  hyperparam_tau_beta_subnational_sd = 0.5,
  
  # Data
  N_national = 1,
  N_subnational = as.array(nrow(use_subnat_data), dim = c(1)),
  
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
    unnest(data) %>%
    modelr::model_matrix(new_deaths ~ average_mob),
  num_coef = 2,
  
  deaths = use_subnat_data %>%
      unnest(data) %>%
      pull(new_deaths),
  
  # WW addition
  travel_delay = 5,
  import_pr_sd = 1e-02
)


options(mc.cores = 4)
mob_model <- stan_model(file.path("mobility-model", "mobility.stan"))

mob_model_s <- stan_model("low-count/mobility-simple.stan")
mob_fit <- sampling(mob_model_s,
                    data = stan_data,
                    iter = 500,
                    chains = 4,
                    # control = lst(
                      # adapt_delta = 0.95,
                      # max_treedepth = 12
                    # ),
                    par = "mean_deaths",
                    include = FALSE)


print(mob_fit, c("imputed_cases", "import_pr"))

save(mob_fit, day_results, file = "res.Rdata")

# Results (Karim) _-----

load("res.Rdata")
source("util.R")
library(rstan)

subnat_results <- mob_fit %>%
  extract_subnat_results(c("R0"))

unnest(unnest(subnat_results))

day_results <- mob_fit %>%
  as.array(par = c("Rt", "Rt_adj")) %>%
  plyr::adply(3, diagnose) %>%
  tidyr::extract(parameters, c("parameter", "long_day_index"), ("(\\w+)\\[(\\d+)\\]"), convert = TRUE) %>%
  mutate(
    iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
    sub_region = rep(rep(use_subnat_data$sub_region, times = stan_data$days_observed), 2),
    quants = map(iter_data, quantilize, iter_value),
    mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
  ) %>%
  unnest(quants) %>%
  nest(day_data = -c(sub_region)) %>%
  mutate(
    day_data = map(day_data,
                   ~ group_by(., parameter) %>%
                     mutate(day_index = seq_along(long_day_index)) %>%
                     ungroup() %>%
                     select(-long_day_index) %>%
                     nest(param_results = -c(day_index)))
  )

use_subnat_data %<>%
  mutate(
    subnat_index = seq(n()),
    daily_data = map(data, mutate, day_index = seq(n()))
  ) %>% select(-data) %>%
  left_join(subnat_results, by = "subnat_index") %>%
  left_join(day_results, by = c("sub_region")) %>%
  mutate(
    daily_data = map2(daily_data, day_data, left_join, by = "day_index")
  ) %>%
  select(-day_data) %>%
  mutate(countrycode_string = "IND")
