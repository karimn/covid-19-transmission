library(magrittr)
library(tidyverse)
library(rstan)
library(bayesplot)
library(cowplot)

source(file.path("util.R"))
source(file.path("mobility-model", "mob_util.R"))

# load("data/mobility/results/pk_63034482_97_mob.RData")
# load("data/mobility/results/ca_63034482_23_mob.RData")
# load("data/mobility/results/ca_63036887_23_mob.RData")
load("data/mobility/results/pk_63036887_97_mob.RData")

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

# load("data/mobility/results/co_62090526_27_mob.RData")
#
# check_hmc_diagnostics(mob_fit)
#
# co_lp <- log_posterior(mob_fit)
# co_np <- nuts_params(mob_fit)
# co_posterior <- as.array(mob_fit)

# co_all_parameters <- extract_parameters(mob_fit) %>%
#   select(-iter_data)

color_scheme_set("darkgray")

mcmc_parcoord(my_posterior, np = my_np, pars = c("overdisp_deaths", "mean_deaths[10]", # "imputed_cases[2]",
                                                 "toplevel_log_R0", "subnational_log_R0[1]", "subnational_log_R0_sd[1]",
                                                 "beta_toplevel[1]", "beta_toplevel[2]", "beta_toplevel[3]",
                                                 "beta_subnational[2,1]", "beta_subnational_sd[2,1]"))
              # alpha = 0.1, np_style = parcoord_style_np(div_alpha = 1, div_size = 0.5))

mcmc_pairs(
  my_posterior, np = my_np,
  pars = c("overdisp_deaths", "mean_deaths[10]", "imputed_cases[2]",
           "toplevel_log_R0", "subnational_log_R0[1]", "subnational_log_R0_sd[1]",
           "beta_toplevel[2]", "beta_subnational[2,1]", "beta_subnational_sd[2,1]",
           "ifr_noise[1]",
           "trend_lambda[1]", "toplevel_trend_kappa"),
  transformations = lst(overdisp_deaths = "log", "mean_deaths[10]" = "log", "subnational_log_R0_sd[1]" = "log",
                        "beta_subnational_sd[2,1]" = "log", "imputed_cases[2]" = "log",
                        "trend_lambda[1]" = "log", "toplevel_trend_kappa" = function(kappa) log(-kappa),
                        "ifr_noise[1]" = "log"),
  off_diag_args = list(size = 0.75)
  # np_style = pairs_style_np(div_size = 2, div_shape = 16, td_alpha = 0.1)
)

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

