#!/usr/bin/Rscript

"Usage:
  run_mob fit [<country-code> ...] [--output=<output-name> --no-partial-pooling --mobility-model-type=<model-type> --chains=<chains> --iter=<iterations> --cmdstan --old-r0 --fixed-tau-beta]
  run_mob prior [<country-code> ...] [--output=<output-name> --no-partial-pooling --mobility-model-type=<model-type> --chains=<chains> --iter=<iterations> --old-r0 --fixed-tau-beta]

Options:
  -c <chains>, --chains=<chains>  Number of chains to use [default: 4]
  -i <iterations>, --iter=<iterations>  Number of iterations [default: 2000]
  -o <output-name>, --output=<output-name>  Output name to use in file names [default: mob]
  --no-partial-pooling  Do not use a hierarchical model
  --mobility-model-type=<model-type>  Type of mobility model (one of: inv_logit, exponential) [default: inv_logit]
  --cmdstan  Use {cmdstanr} instead of {rstan}
  --old-r0  Don't use log R0, instead follow same model as Vollmer et al.
  --fixed-tau-beta  Homogenous partial pooling for all mobility model parameters as in the Vollmer et al. model.
" -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, "fit ita -i 2000")
} else {
  docopt::docopt(opt_desc)
}

library(magrittr)
library(tidyverse)
library(wpp2019)

script_options %<>%
  modify_at(c("chains", "iter"), as.integer) %>%
  modify_at("mobility-model-type", factor, levels = c("inv_logit", "exponential")) %>%
  modify_at("country-code", str_to_upper)
  # modify_at(c(), as.numeric)

if (script_options$cmdstan) {
  library(cmdstanr)
} else {
  library(rstan)
  options(mc.cores = max(1, parallel::detectCores()))
  rstan_options(auto_write = TRUE)
}

source("util.R")
source(file.path("mobility-model", "constants.R"))

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

# Alternative IFR from Vollmer et al. -------------------------------------

italian_region_rename <- list(
  "Friuli Venezia Giulia" = "Friuli Venezia-Giulia",
  "P.A. Trento" = "Provincia Autonoma Trento",
  "Valle d'Aosta" = "Valle D'Aosta"
)

ita_ifr_data <- read_csv(file.path("data", "population", "ita_ifr.csv")) %>%
  select(-X1) %>%
  rename(ifr = IFR, sub_region = state, pop_ana = total_pop) %>%
  mutate(sub_region = fct_recode(sub_region, !!!italian_region_rename))

# Subnational Deaths Data -------------------------------------------------

subnat_data <- read_rds(file.path("data", "mobility", "cleaned_subnat_data.rds"))

if (length(script_options$`country-code`) == 1 && script_options$`country-code` == "ITA") { # Using the same IFR as Vollmer et al. to check for differences in estimation
  use_subnat_data <- subnat_data %>%
    filter(is_valid,
           fct_match(countrycode_string, script_options$`country-code`)) %>%
    select(-pop_ana) %>%
    left_join(ita_ifr_data, by = "sub_region")
} else {
  subnat_data %<>%
    left_join(ifr_adj, by = "countrycode")

  use_subnat_data <- subnat_data %>%
    filter(is_valid,
           fct_match(countrycode_string, script_options$`country-code`))
}

# Time to Death -----------------------------------------------------------

time_to_death_draws <- EnvStats::rgammaAlt(1e6, 5.1, 0.86) + # infection-to-onset distribution
                       EnvStats::rgammaAlt(1e6, 18.8, 0.45) # onset-to-death distribution

quant_time_to_death <- ecdf(time_to_death_draws)

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

# TESTING TESTING
# use_subnat_data %<>% sample_n(10)

stan_data <- lst(
  # Configuration

  fit_model = if (script_options$fit) 1 else if (script_options$prior) 0 else stop("Unsupported run type."),
  hierarchical_mobility_model = !script_options$`no-partial-pooling`,
  mobility_model_type = as.integer(script_options$`mobility-model-type`), # 1: 2 * inv_logit(), 2: exp()
  use_log_R0 = !script_options$`old-r0`,
  use_fixed_tau_beta = script_options$`fixed-tau-beta`,

  # Hyperparameters

  hyperparam_tau_beta_toplevel = if (fct_match(script_options$`mobility-model-type`, "exponential")) 0.1 else 0.5,
  hyperparam_tau_beta_national_sd = if (fct_match(script_options$`mobility-model-type`, "exponential")) 0.1 else 0.5,
  hyperparam_tau_beta_subnational_sd = if (fct_match(script_options$`mobility-model-type`, "exponential")) 0.05 else 0.25,

  hyperparam_toplevel_R0_sd = 0.25,
  hyperparam_tau_national_effect_log_R0_sd = 0.1,
  hyperparam_tau_subnational_effect_log_R0_sd = 0.08,

  # Data

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
    # mutate_all(~ (.x - mean(.x)) / sd(.x)), # Standardization
    # select(1) %>% # Testing speed when all feature are (artificially) uncorrelated
    # mutate(x2 = rnorm(n()), x3 = rnorm(n())),

  num_coef = ncol(design_matrix),

  deaths = if (script_options$fit) {
    use_subnat_data %>%
      unnest(daily_data) %>%
      pull(new_deaths)
  } else {
    array(dim = 0)
  },
)

make_initializer <- function(stan_data) {
  N <- sum(stan_data$N_subnational)
  D_total <- sum(stan_data$days_observed) + stan_data$N_subnational * stan_data$days_to_forecast
  is_multinational <- stan_data$N_national > 1

  function(chain_id) {
    lst(
      beta_toplevel = rnorm(stan_data$num_coef, 0, 0.1),
      beta_national_sd = if (is_multinational) abs(rnorm(stan_data$num_coef, 0, 0.1)) else array(dim = c(0)),

      beta_national_raw = if (is_multinational)
        matrix(rnorm(stan_data$num_coef * stan_data$N_national, 0, 1), nrow = stan_data$num_coef, ncol = stan_data$N_national)
      else array(dim = c(stan_data$num_coef, 0)),

      beta_subnational_sd = if (stan_data$hierarchical_mobility_model)
        matrix(abs(rnorm(stan_data$num_coef * stan_data$N_national, 0, 0.1)), nrow = stan_data$num_coef, ncol = stan_data$N_national)
      else array(dim = c(0, stan_data$N_national)),

      beta_subnational_raw = if (stan_data$hierarchical_mobility_model)
        matrix(rnorm(stan_data$num_coef * N, 0, 1), nrow = stan_data$num_coef, ncol = N)
      else array(dim = c(0, N)),

      beta = matrix(rnorm(stan_data$num_coef * N, 0, 0.1), ncol = stan_data$num_coef, nrow = N),

      mean_deaths = rdunif(D_total, 2, 0),
      new_cases = rdunif(D_total, 2, 0),
      log_R0 = if (stan_data$use_log_R0) rnorm(N, 0, 0.1) else array(dim = 0),
      Rt = runif(D_total, 0, 0.25),
      Rt_adj = runif(D_total, 0, 0.25),

      toplevel_log_R0 = rnorm(1, 0, 0.1),
      national_effect_log_R0_raw = if (is_multinational && stan_data$use_log_R0) rnorm(stan_data$N_national, 0, 0.1) else array(dim = 0),
      national_effect_log_R0_sd = abs(rnorm(stan_data$N_national, 0, 0.1)),
      subnational_effect_log_R0_raw = if (stan_data$use_log_R0) rnorm(N, 0, 0.1) else array(dim = 0),
      subnational_effect_log_R0_sd = if (stan_data$use_log_R0) as.array(abs(rnorm(stan_data$N_national, 0, 0.075))) else array(dim = 0),

      ifr_noise = abs(rnorm(N, 0, 0.1)),
    )
  }
}

cat(str_glue("Running model with {stan_data$N_national} countries and {sum(stan_data$N_subnational)} total subnational entities.\n\n"))
use_subnat_data %>%
  count(countryname) %$% {
    cat(str_glue("\n\t{countryname}: {n} sub regions.\n\n"))
  }
cat("\n")

# Sampling ----------------------------------------------------------------

if (script_options$cmdstan) {
  mob_model <- cmdstan_model(file.path("mobility-model", "mobility.stan"), cpp_options = list(stan_threads = TRUE))
  set_num_threads(3)

  mob_fit <- mob_model$sample(
    stan_data,
    cores = script_options$chains,
    chains = script_options$chains,
    adapt_delta = 0.9,
    # init = 0,
    iter_sampling = script_options$iter %/% 2,
    iter_warmup = script_options$iter %/% 2
  )
} else {
  mob_model <- stan_model(file.path("mobility-model", "mobility.stan"))

  mob_fit <- mob_model %>%
    sampling(
      data = stan_data,
      iter = script_options$iter,
      chains = script_options$chains,
      control = lst(
        adapt_delta = 0.95,
        max_treedepth = 12
      ),
      # init = if (script_options$`random-init`) "random" else make_initializer(stan_data)
      init = make_initializer(stan_data),
      par = "mean_deaths",
      include = FALSE
    )
}

# Extract Results ----------------------------------------------------------

cat("Extracting results...")

subnat_results <- mob_fit %>%
  extract_subnat_results(c("log_R0", "subnational_effect_log_R0"))

day_results <- mob_fit %>%
  extract_day_results(c("Rt", "Rt_adj", "mobility_effect"))

beta_results <- mob_fit %>%
  extract_beta()

use_subnat_data %<>%
  mutate(
    subnat_index = seq(n()),
    daily_data = map(daily_data, mutate, day_index = seq(n()))
  ) %>%
  left_join(subnat_results, by = "subnat_index") %>%
  left_join(day_results, by = c("countrycode_string", "sub_region")) %>%
  mutate(
    daily_data = map2(daily_data, day_data, left_join, by = "day_index")
  ) %>%
  select(-day_data)

cat("done.\n")

# Save --------------------------------------------------------------------

save_file <- file.path("data", "mobility", "results", str_c(script_options$output , ".RData"))
save_results_file <- file.path("data", "mobility", "results", str_c(script_options$output , "_results.rds"))

cat(str_glue("Saving results to {save_file} and {save_results_file} ..."))
save(mob_fit, stan_data, file = save_file)
write_rds(use_subnat_data, save_results_file)
cat("done.\n")
