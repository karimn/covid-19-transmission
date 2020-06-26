library(magrittr)
library(tidyverse)
library(rstan)
library(bayesplot)
library(cowplot)

source(file.path("util.R"))
source(file.path("mobility-model", "mob_util.R"))

load("data/mobility/results/my_62090526_83_mob.RData")

check_hmc_diagnostics(mob_fit)

my_posterior <- as.array(mob_fit)

sub_iter <- sample(nrow(my_posterior), 500)

my_posterior %<>% magrittr::extract(sub_iter, ,)

my_lp <- log_posterior(mob_fit) %>%
  filter(Iteration %in% sub_iter)
my_np <- nuts_params(mob_fit) %>%
  filter(Iteration %in% sub_iter)

my_all_parameters <- extract_parameters(mob_fit) %>%
  select(-iter_data)

load("data/mobility/results/co_62090526_27_mob.RData")

check_hmc_diagnostics(mob_fit)

co_lp <- log_posterior(mob_fit)
co_np <- nuts_params(mob_fit)
co_posterior <- as.array(mob_fit)

# co_all_parameters <- extract_parameters(mob_fit) %>%
#   select(-iter_data)

color_scheme_set("darkgray")

# mcmc_parcoord(my_posterior, np = my_np, pars = c("overdisp_deaths", "tau_impute_cases", "mean_deaths[1075]"))
mcmc_parcoord(my_posterior, np = my_np, pars = c("overdisp_deaths", "mean_deaths[1075]"))

plot_grid(
  mcmc_parcoord(co_posterior, np = co_np, pars = c("overdisp_deaths", "tau_impute_cases", "mean_deaths[1075]"))
)

mcmc_pairs(my_posterior, np = my_np, pars = c("overdisp_deaths", "mean_deaths[1075]"),
           transformations = lst(overdisp_deaths = "log", "mean_deaths[1075]" = "log"),
           off_diag_args = list(size = 0.75))

mcmc_pairs(my_posterior, np = my_np, pars = c("overdisp_deaths",
                                              "mean_deaths[1075]",
                                              "toplevel_log_R0",
                                              "subnational_effect_log_R0_raw[1]", "subnational_effect_log_R0[1]",
                                              "subnational_effect_log_R0_raw[2]", "subnational_effect_log_R0[2]",
                                              #"subnational_effect_log_R0_raw[3]", "subnational_effect_log_R0[3]",
                                              "subnational_effect_log_R0_sd[1]",
                                              # "imputed_cases[1]",
                                              "imputed_cases[2]"),
                                              # "beta_toplevel[1]", "beta_subnational_raw[1,2]", "beta_subnational_sd[1,1]"),
           transformations = lst("overdisp_deaths" = "log",
                                 "subnational_effect_log_R0_sd[1]" = "log",
                                 "mean_deaths[1075]" = "log",
                                 # "beta_subnational_sd[1,1]" = "log",
                                 # "imputed_cases[1]" = "log",
                                 "imputed_cases[2]" = "log"),
           off_diag_args = list(size = 0.75))

