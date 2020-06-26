functions {
  int[] add_array(int[] arr, int to_add) {
    int arr_size = num_elements(arr);
    int result[arr_size] = arr;

    for (arr_index in 1:arr_size) {
      result[arr_index] += to_add;
    }

    return result;
  }

  real neg_binomial_partial_sum(int[] deaths, int start, int end, vector mean_deaths, real overdisp_deaths) {
    return neg_binomial_2_lpmf(deaths | mean_deaths[start:end], overdisp_deaths);
    // return neg_binomial_2_lpmf(deaths | mean_deaths[start:end], overdisp_deaths[day_subnat_idx[start:end]]);
  }
}

data {
  int<lower = 0, upper = 1> fit_model; // Fit vs prior-predict
  int<lower = 0, upper = 1> hierarchical_R0_model;
  int<lower = 0, upper = 1> hierarchical_mobility_model;
  int<lower = 1, upper = 2> mobility_model_type;
  int<lower = 0, upper = 1> use_log_R0;
  int<lower = 0, upper = 1> use_fixed_tau_beta; // Use the same SD for all mobility effects -- homogenous partial pooling.
  int<lower = 0, upper = 1> generate_prediction;
  int<lower = 0, upper = 1> use_transformed_param_constraints;
  int<lower = 0, upper = 1> use_parametric_trend;

  int<lower = 1> N_national; // Number of countries
  int<lower = 1> N_subnational[N_national]; // Number of subnational entities for each country

  int<lower = 1> start_epidemic_offset;

  int<lower = 1> days_observed[sum(N_subnational)]; // Days observed per national/subnational entity can be different
  int<lower = 1> first_case_day_index[sum(N_subnational)];
  int<lower = 1, upper = min(days_observed)> days_to_impute_cases; // Same for all national/subnational entities. We can change that to vary if needed.
  int<lower = 0> days_to_forecast;

  vector<lower = 0>[max(days_observed) + days_to_forecast] time_to_death; // This is pi in the Vollmer model. Needs to be pre-computed empirically using ecdf(); it's the sum of two gamma distributions

  int<lower = 0> num_coef; // Number of coefficients for linear model
  matrix[sum(days_observed) + sum(N_subnational) * days_to_forecast, num_coef] design_matrix; // E.g., mobility

  int<lower = 0> population[sum(N_subnational)];
  vector<lower = 0>[sum(N_subnational)] mean_ifr; // Not sure we'd have this at the subnational level. If we don't we change this to the national level.
  int<lower = 0> deaths[fit_model ? sum(days_observed) : 0];

  int<lower = 1> time_resolution;

  // Hyperparameters

  real<lower = 0> hyperparam_tau_beta_toplevel;
  real<lower = 0> hyperparam_tau_beta_national_sd;
  real<lower = 0> hyperparam_tau_beta_subnational_sd;

  real<lower = 0> hyperparam_tau_national_effect_log_R0_sd;
  real<lower = 0> hyperparam_tau_subnational_effect_log_R0_sd;

  real<lower = 0> hyperparam_toplevel_R0_sd;
}

transformed data {
  int MOBILITY_MODEL_INV_LOGIT = 1;
  int MOBILITY_MODEL_EXPONENTIAL = 2;

  int N = sum(N_subnational);
  int D = sum(days_observed);
  int total_days[N] = add_array(days_observed, days_to_forecast);
  int D_total = sum(total_days);
  int max_days_observed = max(days_observed) + days_to_forecast;
  int is_multinational = N_national > 1;

  int subnat_national_id[N];

  int num_singleton_countries = 0; // Number of countries with only one subnational entity

  real toplevel_R0_mean = 3.28;

  real toplevel_log_R0_sd = sqrt(log((hyperparam_toplevel_R0_sd^2 / toplevel_R0_mean^2) + 1));
  real toplevel_log_R0_mean = log(toplevel_R0_mean) - (toplevel_log_R0_sd^2 / 2);

  int num_likelihood_days[N] = add_array(days_observed, 1 - start_epidemic_offset);
  int total_num_likelihood_days = sum(num_likelihood_days);
  int likelihood_day_idx[total_num_likelihood_days];
  int day_subnat_idx[total_num_likelihood_days];

  vector<lower = 0>[max_days_observed] rev_time_to_death;
  vector<lower = 0>[max_days_observed] rev_gen_factor;

  vector[D_total] trend_day_index;

  {
    real gen_factor_alpha = 1 / (0.62^2);              // alpha = 1 / tau^2;
    real gen_factor_beta = gen_factor_alpha / 6.5;     // beta = 1 / (tau^2 * mu)

    // Need to discretize
    vector[max_days_observed] gen_factor;
    vector[max_days_observed] gamma_cdfs;

    gamma_cdfs[1] = gamma_cdf(1.5 * time_resolution, gen_factor_alpha, gen_factor_beta);
    gen_factor[1] = gamma_cdfs[1];

    // Discretization
    for (day_index in 2:max_days_observed) {
      gamma_cdfs[day_index] = gamma_cdf(time_resolution * day_index + 0.5, gen_factor_alpha, gen_factor_beta);
      gen_factor[day_index] = gamma_cdfs[day_index] - gamma_cdfs[day_index - 1];
    }

    for (day_index in 1:max_days_observed) {
      rev_time_to_death[day_index] = time_to_death[max_days_observed - day_index + 1];
      rev_gen_factor[day_index] = gen_factor[max_days_observed - day_index + 1];
    }
  }

  {
    int lkh_pos = 1;

    for (subnat_index in 1:N) {
      for (lkh_day_index in 1:num_likelihood_days[subnat_index]) {
        likelihood_day_idx[lkh_pos] = start_epidemic_offset + lkh_day_index - 1;
        day_subnat_idx[lkh_pos] = subnat_index;
        lkh_pos += 1;
      }
    }
  }

  {
    int subnat_pos = 1;
    int day_pos = 1;

    for (country_index in 1:N_national) {
      int subnat_end = subnat_pos + N_subnational[country_index] - 1;

      if (N_subnational[country_index] == 1) {
        num_singleton_countries += 1;
      }

      subnat_national_id[subnat_pos:subnat_end] = rep_array(country_index, N_subnational[country_index]);

      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_subnat_pos = subnat_pos + subnat_index - 1;
        int num_days = days_observed[curr_subnat_pos] + days_to_forecast;
        int day_end = day_pos + num_days - 1;

        for (day_index in 1:num_days) {
          int curr_day_pos = day_pos + day_index - 1;

          trend_day_index[curr_day_pos] = first_case_day_index[curr_subnat_pos] - day_index;
        }

        day_pos = day_end + 1;
      }

      subnat_pos = subnat_end + 1;
    }
  }
}

parameters {
  // vector<lower = 0>[N] overdisp_deaths;
  real<lower = 0> overdisp_deaths;

  real<lower = 0> tau_impute_cases;
  vector<lower = 0>[N] imputed_cases;

  // Linear model parameters
  vector[num_coef] beta_toplevel;

  vector<lower = 0>[is_multinational && hierarchical_mobility_model ? (use_fixed_tau_beta ? 1 : num_coef) : 0] beta_national_sd; // Separate SD for each parameters. Vollmer et al. use same for all parameters
  matrix[num_coef, is_multinational && hierarchical_mobility_model ? N_national : 0] beta_national_raw; // Uncentered for now to avoid divergence

  matrix<lower = 0>[hierarchical_mobility_model ? (use_fixed_tau_beta ? 1 : num_coef) : 0, N_national - num_singleton_countries] beta_subnational_sd;
  matrix[hierarchical_mobility_model ? num_coef : 0, N - num_singleton_countries] beta_subnational_raw;

  vector<lower = 0>[N] ifr_noise;

  vector<lower = 0>[use_log_R0 ? 0 : N] original_R0;
  real<lower = 0> original_R0_sd;

  real toplevel_log_R0;
  vector[is_multinational && use_log_R0 && hierarchical_R0_model ? N_national : 0] national_effect_log_R0_raw;
  real<lower = 0> national_effect_log_R0_sd;
  vector[use_log_R0 && hierarchical_R0_model ? N - num_singleton_countries : 0] subnational_effect_log_R0_raw;
  vector<lower = 0>[use_log_R0 && hierarchical_R0_model ? N_national - num_singleton_countries : 0] subnational_effect_log_R0_sd;

  vector<lower = 0, upper = 1>[use_parametric_trend ? N : 0] trend_lambda;
  real toplevel_trend_kappa;
  vector<lower = 0>[use_parametric_trend ? N_national : 0] trend_kappa_subnational_sd;
  vector[use_parametric_trend ? N : 0] trend_kappa_subnational_raw;
}

transformed parameters {
  vector[N] log_R0 = use_log_R0 ? rep_vector(toplevel_log_R0, N) : log(original_R0);
  vector[use_log_R0 && hierarchical_R0_model ? N_national : 0] national_effect_log_R0;
  vector[use_log_R0 && hierarchical_R0_model ? N - num_singleton_countries : 0] subnational_effect_log_R0;

  // vector<lower = (use_transformed_param_constraints ? 0 : negative_infinity())>[D_total] mean_deaths; // Not a matrix; this could be a ragged data structure
  vector<lower = 0>[D_total] mean_deaths = rep_vector(0, D_total); // Not a matrix; this could be a ragged data structure
  vector<lower = (use_transformed_param_constraints ? 0 : negative_infinity())>[D_total] Rt = rep_vector(0, D_total);
  vector<lower = (use_transformed_param_constraints ? 0 : negative_infinity())>[D_total] Rt_adj = Rt;
  row_vector<lower = (use_transformed_param_constraints ? 0 : negative_infinity())>[D_total] new_cases = rep_row_vector(0, D_total);

  vector[D_total] adj_factor = rep_vector(1, D_total);

  matrix[num_coef, N] beta = rep_matrix(beta_toplevel, N);

  vector[D_total] mobility_effect = rep_vector(1, D_total);

  vector[N] ifr = mean_ifr .* ifr_noise;

  vector[use_parametric_trend ? N : 0] trend_kappa;
  vector[D_total] log_trend = rep_vector(0, D_total);

  if (use_parametric_trend) {
    trend_kappa = toplevel_trend_kappa + trend_kappa_subnational_raw .* trend_kappa_subnational_sd[subnat_national_id];
  }

  if (use_log_R0 && hierarchical_R0_model) {
    national_effect_log_R0 = rep_vector(0, N_national);
    subnational_effect_log_R0 = rep_vector(0, N - num_singleton_countries);
  }

  {
    int full_subnat_pos = 1; // The "full" pointers do not exclude subnational entities in singleton countries (have only one subnational entity)
    int subnat_pos = 1;
    int days_pos = 1;

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int full_subnat_end = full_subnat_pos + num_subnat - 1;
      int subnat_end = subnat_pos + (num_subnat > 1 ? num_subnat : 0) - 1; // Drop the single subnational entity

      if (is_multinational) {
        if (hierarchical_mobility_model) {
          if (use_fixed_tau_beta) {
            beta[, full_subnat_pos:full_subnat_end] += rep_matrix(beta_national_raw[, country_index] * beta_national_sd[1], num_subnat);
          } else {
            beta[, full_subnat_pos:full_subnat_end] += rep_matrix(beta_national_raw[, country_index] .* beta_national_sd, num_subnat);
          }
        }

        if (use_log_R0 && hierarchical_R0_model) {
          national_effect_log_R0[country_index] = national_effect_log_R0_raw[country_index] * national_effect_log_R0_sd;
          log_R0[full_subnat_pos:full_subnat_end] += national_effect_log_R0[country_index];
        }
      }

      if (use_log_R0 && hierarchical_R0_model && num_subnat > 1) {
        subnational_effect_log_R0[subnat_pos:subnat_end] = subnational_effect_log_R0_raw[subnat_pos:subnat_end] * subnational_effect_log_R0_sd[country_index];
        log_R0[full_subnat_pos:full_subnat_end] += subnational_effect_log_R0[subnat_pos:subnat_end];
      }

      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_full_subnat_pos = full_subnat_pos + subnat_index - 1;
        int curr_subnat_pos = subnat_pos + subnat_index - 1;
        int num_days = days_observed[curr_full_subnat_pos] + days_to_forecast;
        int days_end = days_pos + num_days - 1;

        vector[total_days[curr_full_subnat_pos]] cumulative_cases = rep_vector(0, total_days[curr_full_subnat_pos]);

        if (use_parametric_trend) {
          log_trend[days_pos:days_end] = log(1 - trend_lambda[curr_full_subnat_pos]) + log_inv_logit(trend_kappa[curr_full_subnat_pos] * trend_day_index[days_pos:days_end]);
        }

        if (hierarchical_mobility_model && num_subnat > 1) {
          if (use_fixed_tau_beta) {
            beta[, curr_full_subnat_pos] += beta_subnational_raw[, curr_subnat_pos] * beta_subnational_sd[1, country_index];
          } else {
            beta[, curr_full_subnat_pos] += beta_subnational_raw[, curr_subnat_pos] .* beta_subnational_sd[, country_index];
          }
        }

        new_cases[days_pos:(days_pos + days_to_impute_cases - 1)] = rep_row_vector(imputed_cases[curr_full_subnat_pos] * time_resolution, days_to_impute_cases);

        if (mobility_model_type == MOBILITY_MODEL_INV_LOGIT) {
          mobility_effect[days_pos:days_end] = 2 * inv_logit(design_matrix[days_pos:days_end] * beta[, curr_full_subnat_pos] + log_trend[days_pos:days_end]);

          Rt[days_pos:days_end] = exp(log_R0[curr_full_subnat_pos]) * mobility_effect[days_pos:days_end];
        } else {
          mobility_effect[days_pos:days_end] = exp(design_matrix[days_pos:days_end] * beta[, curr_full_subnat_pos] + log_trend[days_pos:days_end]);

          Rt[days_pos:days_end] = exp(log_R0[curr_full_subnat_pos] + design_matrix[days_pos:days_end] * beta[, curr_full_subnat_pos] + log_trend[days_pos:days_end]);
        }

        for (day_index in 1:total_days[curr_full_subnat_pos]) {
          int curr_day_pos = days_pos + day_index - 1;

          if (day_index > 1) {
            adj_factor[curr_day_pos] = 1 - (cumulative_cases[day_index - 1] / population[curr_full_subnat_pos]);

            if (day_index > days_to_impute_cases) {
              Rt_adj[curr_day_pos] = adj_factor[curr_day_pos] * Rt[curr_day_pos];

              new_cases[curr_day_pos] = Rt_adj[curr_day_pos] * new_cases[days_pos:(curr_day_pos - 1)] * tail(rev_gen_factor, day_index - 1);
            } else {
              Rt_adj[curr_day_pos] = Rt[curr_day_pos];
            }
          } else {
            Rt_adj[curr_day_pos] = Rt[curr_day_pos];
          }

          cumulative_cases[day_index] = new_cases[curr_day_pos] + (day_index > 1 ? cumulative_cases[day_index - 1] : 0);

          if (day_index > 1) {
            mean_deaths[curr_day_pos] = ifr[curr_full_subnat_pos] * new_cases[days_pos:(curr_day_pos - 1)] * tail(rev_time_to_death, day_index - 1);
          } else {
            mean_deaths[curr_day_pos] = 1e-15 * new_cases[curr_day_pos];
          }
        }

        days_pos = days_end + 1;
      }

      full_subnat_pos = full_subnat_end + 1;
      subnat_pos = subnat_end + 1;
    }
  }
}

model {
  // These are not model as hierarchical over countries yet.
  overdisp_deaths ~ normal(0, 5);
  ifr_noise ~ normal(1, 0.1);

  tau_impute_cases ~ exponential(0.03);
  imputed_cases ~ exponential(1 / tau_impute_cases);

  beta_toplevel ~ normal(0, hyperparam_tau_beta_toplevel);

  if (hierarchical_mobility_model) {
    if (is_multinational) {
      beta_national_sd ~ normal(0, hyperparam_tau_beta_national_sd);
      to_vector(beta_national_raw) ~ std_normal();
    }

    to_vector(beta_subnational_sd) ~ normal(0, hyperparam_tau_beta_subnational_sd);
    to_vector(beta_subnational_raw) ~ std_normal();
  }

  toplevel_log_R0 ~ normal(toplevel_log_R0_mean, toplevel_log_R0_sd);
  national_effect_log_R0_sd ~ normal(0, hyperparam_tau_national_effect_log_R0_sd);

  original_R0_sd ~ normal(0, 0.5);

  if (use_log_R0 && hierarchical_R0_model) {
    if (is_multinational) {
      national_effect_log_R0_raw ~ std_normal();
    }

    subnational_effect_log_R0_raw ~ std_normal();
    subnational_effect_log_R0_sd ~ normal(0, hyperparam_tau_subnational_effect_log_R0_sd);
  } else {
    original_R0 ~ normal(toplevel_R0_mean, original_R0_sd);
  }

  toplevel_trend_kappa ~ std_normal();

  if (use_parametric_trend) {
    trend_lambda ~ beta(3, 1);
    trend_kappa_subnational_sd ~ std_normal();
    trend_kappa_subnational_raw ~ std_normal();
  }

  if (fit_model) {
    // target += reduce_sum(neg_binomial_partial_sum, deaths, 1, mean_deaths, day_subnat_idx, overdisp_deaths);
    // target += reduce_sum(neg_binomial_partial_sum, deaths, 1, mean_deaths, overdisp_deaths);

    int subnat_pos = 1;
    int days_pos = 1;

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int subnat_end = subnat_pos + num_subnat - 1;

      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_subnat_pos = subnat_pos + subnat_index - 1;
        int days_end = days_pos + days_observed[curr_subnat_pos] - 1;

        deaths[(days_pos + start_epidemic_offset - 1):days_end] ~ neg_binomial_2(mean_deaths[(days_pos + start_epidemic_offset - 1):days_end],
                                                                                 overdisp_deaths);
                                                                                 // overdisp_deaths[curr_subnat_pos]);

        days_pos = days_end + days_to_forecast + 1;
      }

      subnat_pos = subnat_end + 1;
    }
  }
}

generated quantities {
  vector<lower = 0>[generate_prediction ? D : 0] deaths_rep;
  vector<lower = 0>[generate_prediction ? D : 0] cum_deaths_rep;

  vector[D_total] trend = exp(log_trend);

  if (generate_prediction) {
    int subnat_pos = 1;
    int days_pos = 1;

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int subnat_end = subnat_pos + num_subnat - 1;

      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_subnat_pos = subnat_pos + subnat_index - 1;
        int days_end = days_pos + days_observed[curr_subnat_pos] - 1;

        deaths_rep[days_pos:days_end] = to_vector(neg_binomial_2_rng(mean_deaths[days_pos:days_end], overdisp_deaths));
        cum_deaths_rep[days_pos:days_end] = cumulative_sum(deaths_rep[days_pos:days_end]);

        days_pos = days_end + days_to_forecast + 1;
      }

      subnat_pos = subnat_end + 1;
    }
  }
}

