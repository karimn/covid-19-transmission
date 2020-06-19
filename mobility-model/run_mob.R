#!/bin/Rscript

"Usage:
  run_mob (fit | prior) [<country-code> ...] [options]

Options:
  -c <chains>, --chains=<chains>  Number of chains to use [default: 4]
  -i <iterations>, --iter=<iterations>  Total number of iterations [default: 2000]
  -w <iterations>, --warmup=<iterations>  Number of warmup iteration. By default this would be half the total number of iterations.
  -o <output-name>, --output=<output-name>  Output name to use in file names [default: mob]
  --no-partial-pooling=<which-parts>  Do not use a hierarchical model (parts: all,mob,r0)
  --mobility-model-type=<model-type>  Type of mobility model (one of: inv_logit, exponential) [default: inv_logit]
  --mobility-model=<model-formula>  Linear mobility model. Makes sure there are no spaces. Don't forget to remove the intercept from the formula.
  --hyperparam=<hyperparam-file>  Use YAML file to specify hyperparameter values
  --merge-days=<num-days>  Number of days to merge together.
  --cmdstan  Use {cmdstanr} instead of {rstan}
  --old-r0  Don't use log R0, instead follow same model as Vollmer et al.
  --fixed-tau-beta  Homogenous partial pooling for all mobility model parameters as in the Vollmer et al. model.
  --no-post-predict  Don't do posterior prediction
  --rand-sample-subnat=<sample-size>  Instead of running all of subnational units, run with a random sample. Only allowed with one country.
  --show-script-options
" -> opt_desc

script_options <- if (interactive()) {
  root_path <- "."

  # docopt::docopt(opt_desc, 'fit ar au ca pt pl -i 1000 -o ar_au_ca_pt_pl_mob_all_pooling --no-partial-pooling=all --mobility-model=~0+average_all_mob')
  # docopt::docopt(opt_desc, 'fit it -i 1000 --merge-days=2')
  # docopt::docopt(opt_desc, "fit ar au ca pt pl -i 2000 -o ar_au_ca_pt_pl_mob_r0_pooling --no-partial-pooling=r0")
  # docopt::docopt(opt_desc, "fit ar au ca pt pl -i 1000 --hyperparam=mobility-model/test_hyperparam.yaml")
  docopt::docopt(opt_desc, "fit 1 3 -i 1000 --hyperparam=mobility-model/test_hyperparam.yaml -o test_{all_country_codes}")
} else {
  root_path <- ".."

  setwd(root_path)
  source(file.path("renv", "activate.R"))
  setwd("mobility-model")

  docopt::docopt(opt_desc)
}

if (script_options$`show-script-options`) {
  print(script_options)
}

library(magrittr)
library(tidyverse)
library(wpp2019)

script_options %<>%
  modify_at(c("chains", "iter", "warmup", "rand-sample-subnat", "merge-days"), as.integer) %>%
  modify_at("mobility-model-type", factor, levels = c("inv_logit", "exponential")) %>%
  modify_at("country-code", str_to_upper)
  # modify_at(c(), as.numeric)

if (script_options$cmdstan) {
  library(cmdstanr)
} else {
  library(rstan)
  options(mc.cores = script_options$chains) # max(1, parallel::detectCores()))
  rstan_options(auto_write = TRUE)
}

time_resolution <- if (is_empty(script_options$`merge-days`)) 1 else script_options$`merge-days`

if (!is_null(script_options$`no-partial-pooling`)) {
  tryCatch(
    script_options$`no-partial-pooling` %<>%
      rlang::arg_match(values = c("all", "mob", "r0")) %>%
      factor(levels = c("all", "mob", "r0")),

    error = function(err) stop("Unexpected value for --no-partial-pooling")
  )
}

source(file.path(root_path, "util.R"))
source(file.path(root_path, "mobility-model", "constants.R"))

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
  group_by(country_code) %>%
  summarize(ifr = sum(ifr * pop) / sum(pop)) %>%
  ungroup()

# Alternative IFR from Vollmer et al. -------------------------------------

italian_region_rename <- list(
  "Friuli Venezia Giulia" = "Friuli Venezia-Giulia",
  "P.A. Trento" = "Provincia Autonoma Trento",
  "Valle d'Aosta" = "Valle D'Aosta"
)

ita_ifr_data <- read_csv(file.path(root_path, "data", "population", "ita_ifr.csv")) %>%
  select(-X1) %>%
  rename(ifr = IFR, sub_region = state, population = total_pop) %>%
  mutate(sub_region = fct_recode(sub_region, !!!italian_region_rename))

# Subnational Deaths Data -------------------------------------------------

subnat_data <- read_rds(file.path(root_path, "data", "mobility", "cleaned_subnat_data.rds"))

# if (length(script_options$`country-code`) == 1 && script_options$`country-code` == "IT") { # Using the same IFR as Vollmer et al. to check for differences in estimation
#   use_subnat_data <- subnat_data %>%
#     filter(has_epidemic,
#            fct_match(country_code, "IT")) %>%
#     select(-population) %>%
#     left_join(ita_ifr_data, by = "sub_region")
# } else {
  subnat_data %<>%
    left_join(ifr_adj, by = c("countrycode_iso3n" = "country_code"))

  if (!is_empty(script_options$`country-code`)) {
    if (all(str_detect(script_options$`country-code`, "^\\d+$"))) {
      cat("\nFiltering countries by index...")
      use_subnat_data <- subnat_data %>%
        filter(country_index %in% as.integer(script_options$`country-code`))
      cat("done.\n")
    } else {
      cat("\nFiltering countries by code...")
      use_subnat_data <- subnat_data %>%
        filter(fct_match(country_code, script_options$`country-code`))
      cat("done.\n")
    }
  } else {
    use_subnat_data <- subnat_data
  }
# }

if (any(!use_subnat_data$is_valid)) {
  warning("{nrow(filter(use_subnat_data, !is_valid))} rows found with invalid data.")
}

if (any(!use_subnat_data$has_epidemic)) {
  cat(str_glue("\n{sum(!use_subnat_data$has_epidemic)} subregions found without an epidemic. Removing them.\n"))
}

use_subnat_data %<>%
    filter(has_epidemic)

if (!is_empty(script_options$`rand-sample-subnat`)) {
  if (n_distinct(use_subnat_data$country_code) > 1) {
    stop("Only one country allowed with the --rand-sample-subnat option.")
  }

  use_subnat_data %<>%
    sample_n(script_options$`rand-sample-subnat`)
}

if (!is_empty(script_options$`merge-days`)) {
  cat(str_glue("Merging days (resolution = {time_resolution})..."))

  use_subnat_data %<>%
    mutate(
      daily_data = map(
        daily_data,
        ~ mutate(.x, day_group_index = ((day_index - 1) %/% time_resolution) + 1) %>%
          group_by(day_group_index) %>% # Using mutate instead of summarize because I need to combine mutate with mutate_at.
          mutate(
            date = first(date),
            new_deaths = sum(new_deaths),
            cum_deaths = last(cum_deaths),
            all_mob_observed = sum(all_mob_observed) == n(),
            group_size = n()
          ) %>%
          mutate_at(vars(starts_with("g_"), average_mob), mean) %>%
          slice(1) %>%
          select(day_group_index, group_size, date, new_deaths, cum_deaths, all_mob_observed, starts_with("g_"), average_mob) %>%
          ungroup() %>%
          filter(group_size == time_resolution)),

      num_days_observed = map_int(daily_data, nrow)
    )

  cat("done.\n")
}

all_country_codes <- use_subnat_data %>%
  pull(country_code) %>%
  unique() %>%
  str_to_lower() %>%
  str_c(collapse = "_")

save_file <- file.path(root_path, "data", "mobility", "results", str_c(str_glue(script_options$output), ".RData"))
save_results_file <- file.path(root_path, "data", "mobility", "results", str_c(str_glue(script_options$output), "_results.rds"))

# Time to Death -----------------------------------------------------------

# time_to_death_draws <- EnvStats::rgammaAlt(1e8, 5.1, 0.86) + # infection-to-onset distribution
#                        EnvStats::rgammaAlt(1e8, 18.8, 0.45) # onset-to-death distribution
#
# quant_time_to_death <- ecdf(time_to_death_draws)
#
# time_to_death_ecdf <- c(0, seq(1.5, by = 1, length.out = 365)) %>%
#   map_dbl(quant_time_to_death)
#
# write_rds(time_to_death_ecdf, path = file.path(root_path, "data", "mobility", "time_to_death.rds"))

time_to_death_ecdf <- read_rds(file.path(root_path, "data", "mobility", "time_to_death.rds"))

time_to_death <- time_to_death_ecdf %>%
  magrittr::extract(((seq_along(.) - 1) %% time_resolution) == 0) %>%
  subtract(lag(.)) %>%
  discard(is.na)

# Plot timelines ----------------------------------------------------------

# use_subnat_data %>%
#   mutate(sub_region = fct_reorder(sub_region, epidemic_start_date)) %>%
#   ggplot(aes(y = sub_region)) +
#   geom_errorbar(aes(xmin = first_observed_death, xmax = last_observed_death), width = 0.25) +
#   geom_linerange(aes(xmin = first_mob_day, xmax = last_mob_day), size = 7, alpha = 0.15) +
#   geom_point(aes(first_infection_seeding_day)) +
#   geom_point(aes(epidemic_start_date)) +
#   geom_point(aes(last_effective_observed_day)) +
#   ggrepel::geom_text_repel(aes(x = epidemic_start_date, label = "Epidemic"), nudge_y = 0.25, size = 3.5) +
#   ggrepel::geom_text_repel(aes(x = first_infection_seeding_day, label = "First Seeding"), nudge_y = 0.25, size = 3.5) +
#   ggrepel::geom_text_repel(aes(x = last_effective_observed_day, label = "Last Effective Observation"), nudge_y = 0.25, size = 3.5) +
#   geom_linerange(aes(xmin = first_infection_seeding_day, xmax = epidemic_start_date - 1), linetype = "dashed") +
#   scale_color_discrete("") +
#   labs(x = "", y = "",
#        caption = "Grey bars: range of mobility data. Black lines: range of infection/deaths data.") +
#   facet_wrap(vars(country_code), ncol = 1) +
#   theme_minimal()

# Stan Data ---------------------------------------------------------------

stan_data <- lst(
  # Configuration

  fit_model = if (script_options$fit) 1 else if (script_options$prior) 0 else stop("Unsupported run type."),
  hierarchical_R0_model = is_null(script_options$`no-partial-pooling`) || !fct_match(script_options$`no-partial-pooling`, c("all", "r0")),
  hierarchical_mobility_model = is_null(script_options$`no-partial-pooling`) || !fct_match(script_options$`no-partial-pooling`, c("all", "mob")),
  mobility_model_type = as.integer(script_options$`mobility-model-type`), # 1: 2 * inv_logit(), 2: exp()
  use_log_R0 = !script_options$`old-r0`,
  use_fixed_tau_beta = script_options$`fixed-tau-beta`,
  generate_post_prediction = !script_options$`no-post-predict`,
  use_transformed_param_constraints = 0,

  # Hyperparameters

  hyperparam_tau_beta_toplevel = if (fct_match(script_options$`mobility-model-type`, "exponential")) 0.1 else 0.5,
  hyperparam_tau_beta_national_sd = if (fct_match(script_options$`mobility-model-type`, "exponential")) 0.1 else 0.5,
  hyperparam_tau_beta_subnational_sd = if (fct_match(script_options$`mobility-model-type`, "exponential")) 0.05 else 0.25,

  hyperparam_toplevel_R0_sd = 0.25,
  hyperparam_tau_national_effect_log_R0_sd = 0.1,
  hyperparam_tau_subnational_effect_log_R0_sd = 0.08,

  # Data

  N_national = n_distinct(use_subnat_data$country_code),
  N_subnational = use_subnat_data %>% count(country_code) %>% pull(n) %>% as.array(),

  start_epidemic_offset = (start_epidemic_offset - 1) %/% time_resolution + 1,
  days_observed = as.integer(use_subnat_data$num_days_observed) %>% as.array(),
  days_to_impute_cases = days_seeding %/% time_resolution,
  days_to_forecast = days_to_forecast %/% time_resolution,
  total_days = max(days_observed),

  time_to_death = time_to_death[1:total_days],

  population = use_subnat_data$population %>%
    as.array(),
  mean_ifr = use_subnat_data$ifr %>%
    as.array(),

  design_matrix = use_subnat_data %>%
    unnest(daily_data) %>%
    modelr::model_matrix(if (is_null(script_options$`mobility-model`)) mob_formula else as.formula(script_options$`mobility-model`)),
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

  time_resolution,
)

if (!is_null(script_options$hyperparam)) {
  hyperparam <- yaml::yaml.load_file(script_options$hyperparam)

  cat("\nUsing hyperparameters:\n")
  iwalk(hyperparam, ~ cat(str_c("\t", .y, " = ", .x, "\n")))
  cat("\n")

  stan_data %<>%
    list_modify(!!!hyperparam)
}

# Initializer -------------------------------------------------------------

make_initializer <- function(stan_data) {
  N <- sum(stan_data$N_subnational)
  D_total <- sum(stan_data$days_observed) + stan_data$N_subnational * stan_data$days_to_forecast
  is_multinational <- stan_data$N_national > 1
  num_singleton_countries <- sum(stan_data$N_subnational == 1)

  function(chain_id) {
    lst(
      beta_toplevel = as.array(rnorm(stan_data$num_coef, 0, 0.1)),
      beta_national_sd = if (is_multinational && stan_data$hierarchical_mobility_model) as.array(abs(rnorm(stan_data$num_coef, 0, 0.1))) else array(dim = c(0)),

      beta_national_raw = if (is_multinational && stan_data$hierarchical_mobility_model)
        matrix(rnorm(stan_data$num_coef * stan_data$N_national, 0, 1), nrow = stan_data$num_coef, ncol = stan_data$N_national)
      else array(dim = c(stan_data$num_coef, 0)),

      beta_subnational_sd = if (stan_data$hierarchical_mobility_model)
        matrix(abs(rnorm(stan_data$num_coef * (stan_data$N_national - num_singleton_countries), 0, 0.1)),
               nrow = stan_data$num_coef,
               ncol = stan_data$N_national - num_singleton_countries)
      else array(dim = c(0, stan_data$N_national - num_singleton_countries)),

      beta_subnational_raw = if (stan_data$hierarchical_mobility_model)
        matrix(rnorm(stan_data$num_coef * (N - num_singleton_countries), 0, 1), nrow = stan_data$num_coef, ncol = N - num_singleton_countries)
      else array(dim = c(0, N - num_singleton_countries)),

      beta = matrix(rnorm(stan_data$num_coef * N, 0, 0.1), ncol = stan_data$num_coef, nrow = N),

      mean_deaths = rdunif(D_total, 2, 0),
      new_cases = rdunif(D_total, 2, 0),
      log_R0 = if (stan_data$use_log_R0) rnorm(N, 0, 0.1) else array(dim = 0),
      Rt = runif(D_total, 0, 0.25),
      Rt_adj = runif(D_total, 0, 0.25),

      toplevel_log_R0 = rnorm(1, 0, 0.1),
      national_effect_log_R0_raw = if (is_multinational && stan_data$use_log_R0 && stan_data$hierarchical_R0_model) rnorm(stan_data$N_national, 0, 0.1) else array(dim = 0),
      national_effect_log_R0_sd = abs(rnorm(1, 0, 0.1)),
      subnational_effect_log_R0_raw = if (stan_data$use_log_R0 && stan_data$hierarchical_R0_model) rnorm(N - num_singleton_countries, 0, 0.1) else array(dim = 0),
      subnational_effect_log_R0_sd = if (stan_data$use_log_R0 && stan_data$hierarchical_R0_model) as.array(abs(rnorm(stan_data$N_national - num_singleton_countries, 0, 0.075))) else array(dim = 0),

      ifr_noise = as.array(abs(rnorm(N, 0, 0.1))),
    )
  }
}

cat(str_glue("Running model with {stan_data$N_national} countries and {sum(stan_data$N_subnational)} total subnational entities.\n\n"))
use_subnat_data %>%
  count(country_name) %$% {
    cat(str_glue("\n\t{country_name}: {n} sub regions.\n\n"))
  }
cat("\n")

# Sampling ----------------------------------------------------------------

if (script_options$cmdstan) {
  mob_model <- cmdstan_model(file.path(root_path, "mobility-model", "mobility.stan"), cpp_options = list(stan_threads = TRUE))
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
  mob_model <- stan_model(file.path(root_path, "mobility-model", "mobility.stan"))

  sampling_args <- rlang::list2(
    mob_model,
    data = stan_data,
    iter = script_options$iter,
    chains = script_options$chains,
    control = lst(
      adapt_delta = 0.95,
      max_treedepth = 12
    ),
    # init = if (script_options$`random-init`) "random" else make_initializer(stan_data)
    init = make_initializer(stan_data),
    save_warmup = FALSE,
    # par = "mean_deaths",
    # include = FALSE
  )

  if (!is_empty(script_options$warmup)) {
    sampling_args %<>%
      update_list(
        warmup = script_options$warmup
      )
  }

  tictoc::tic("Sampling")

  mob_fit <- exec(sampling, !!!sampling_args)

  tictoc::toc()
}

# Extract Results ----------------------------------------------------------

tryCatch({
  cat("\nExtracting results...\n\n")

  all_parameters <- extract_parameters(mob_fit)

  cat("Maximum Rhat = ")
  all_parameters %>%
    select(rhat) %>%
    summarize_all(max, na.rm = TRUE) %>%
    first() %>%
    cat()
  cat("\nMinimum ESS Bulk = ")
  all_parameters %>%
    select(ess_bulk) %>%
    summarize_all(min, na.rm = TRUE) %>%
    first() %>%
    cat()
  cat("\nMinimum ESS Tail = ")
  all_parameters %>%
    select(ess_tail) %>%
    summarize_all(min, na.rm = TRUE) %>%
    first() %>%
    cat()

  cat("\n\n")

  subnat_results <- mob_fit %>%
    extract_subnat_results(c("log_R0", "national_effect_log_R0", "subnational_effect_log_R0", "imputed_cases"))

  day_param <- c("Rt", "Rt_adj", "adj_factor", "mobility_effect", "mean_deaths")

  if (!script_options$`no-post-predict` && !script_options$prior) {
    day_param %<>% c("deaths_rep")
  }

  day_results <- mob_fit %>%
    extract_day_results(day_param)

  beta_results <- mob_fit %>%
    extract_beta() %>%
    rename(beta_results = param_results)

  use_subnat_data %<>%
    mutate(
      subnat_index = seq(n()),
      daily_data = map(daily_data, mutate, day_index = seq(n()))
    ) %>%
    left_join(subnat_results, by = "subnat_index") %>%
    left_join(beta_results, by = "subnat_index") %>%
    left_join(day_results, by = c("country_code", "sub_region")) %>%
    mutate(
      daily_data = map2(daily_data, day_data, left_join, by = "day_index")
    ) %>%
    select(-day_data)

  cat("done.\n")

  cat(str_glue("Saving results to {save_results_file} ..."))
  write_rds(use_subnat_data, save_results_file)
  cat("done.\n")
},
error = function(err) {
  cat("error encountered while extracting results.\n")
},
finally = {
  cat(str_glue("Saving fit to {save_file} ..."))
  save(mob_fit, stan_data, file = save_file)
  cat("done.\n")
})

# Save --------------------------------------------------------------------

