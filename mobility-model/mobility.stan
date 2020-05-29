functions {
  int[] add_array(int[] arr, int to_add) {
    int arr_size = num_elements(arr);
    int result[arr_size] = arr;

    for (arr_index in 1:arr_size) {
      result[arr_index] += to_add;
    }

    return result;
  }
}

data {
  int<lower = 0, upper = 1> fit_model; // Fit vs prior-predict

  int<lower = 1> N_national; // Number of countries
  int<lower = 1> N_subnational[N_national]; // Number of subnational entities for each country

  int<lower = 1> start_epidemic_offset;

  int<lower = 1> days_observed[sum(N_subnational)]; // Days observed per national/subnational entity can be different
  int<lower = 1, upper = min(days_observed)> days_to_impute_cases; // Same for all national/subnational entities. We can change that to vary if needed.
  int<lower = 0> days_to_forecast;

  vector<lower = 0>[max(days_observed) + days_to_forecast] time_to_death; // This is pi in the Vollmer model. Needs to be pre-computed empirically using ecdf(); it's the sum of two gamma distributions

  int<lower = 0> num_coef; // Number of coefficients for linear model
  matrix[sum(days_observed) + sum(N_subnational) * days_to_forecast, num_coef] design_matrix; // E.g., mobility

  int<lower = 0> population[sum(N_subnational)];
  vector<lower = 0>[sum(N_subnational)] mean_ifr; // Not sure we'd have this at the subnational level. If we don't we change this to the national level.
  int<lower = 0> deaths[fit_model ? sum(days_observed) : 0];
}

transformed data {
  int N = sum(N_subnational);
  int D = sum(days_observed);
  int total_days[N] = add_array(days_observed, days_to_forecast);
  int D_total = sum(total_days);
  int max_days_observed = max(days_observed);

  real gen_factor_alpha = 1 / (0.62^2);          // alpha = 1 / tau^2;
  real gen_factor_beta = gen_factor_alpha / 6.5;     // beta = 1 / (tau^2 * mu)

  // Need to discretize
  vector<lower = 0>[max_days_observed] gen_factor;

  {
    vector[max_days_observed] gamma_cdfs;

    gamma_cdfs[1] = gamma_cdf(1.5, gen_factor_alpha, gen_factor_beta);
    gen_factor[1] = gamma_cdfs[1];

    // Discretization
    for (day_index in 2:max_days_observed) {
      gamma_cdfs[day_index] = gamma_cdf(day_index + 0.5, gen_factor_alpha, gen_factor_beta);
      gen_factor[day_index] = gamma_cdfs[day_index] - gamma_cdfs[day_index - 1];
    }
  }
}

parameters {
  vector<lower = 0>[N] overdisp_deaths;

  real<lower = 0> tau_impute_cases;
  vector<lower = 0>[N] imputed_cases;

  // Linear model parameters
  vector[num_coef] beta_toplevel;

  vector<lower = 0>[num_coef] beta_national_sd;
  matrix[num_coef, N_national] beta_national_raw; // Uncentered for now to avoid divergence

  matrix<lower = 0>[num_coef, N_national] beta_subnational_sd;
  matrix[num_coef, N] beta_subnational_raw;

  vector<lower = 0>[N] ifr_noise;

  vector<lower = 0>[N] R0;
  real<lower = 0> R0_sd;
}

transformed parameters {
  vector<lower = 0>[D_total] mean_deaths; // Not a matrix; this could be a ragged data structure
  vector<lower = 0>[D_total] Rt = rep_vector(0, D_total);
  row_vector<lower = 0>[D_total] new_cases = rep_row_vector(0, D_total);

  {
    int subnat_pos = 1;
    int days_pos = 1;

    vector[N] ifr = mean_ifr .* ifr_noise;

    matrix[num_coef, N] beta = rep_matrix(beta_toplevel, N);

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int subnat_end = subnat_pos + num_subnat - 1;

      beta[, subnat_pos:subnat_end] += rep_matrix(beta_national_raw[, country_index] .* beta_national_sd, num_subnat);

      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_subnat_pos = subnat_pos + subnat_index - 1;
        int days_end = days_pos + days_observed[curr_subnat_pos] - 1;

        vector[total_days[curr_subnat_pos]] cumulative_cases;

        // print("subnat_index = ", subnat_index);

        beta[, curr_subnat_pos] += beta_subnational_raw[, curr_subnat_pos] .* beta_subnational_sd[, country_index];

        new_cases[days_pos:(days_pos + days_to_impute_cases - 1)] = rep_row_vector(imputed_cases[subnat_index], days_to_impute_cases);

        Rt[days_pos:days_end] = 2 * R0[curr_subnat_pos] * inv_logit(design_matrix[days_pos:days_end] * beta[, curr_subnat_pos]);

        for (day_index in 1:total_days[curr_subnat_pos]) {
          int curr_day_pos = days_pos + day_index - 1;

          if (day_index > days_to_impute_cases) {
            real adjust_factor = 1 - (cumulative_cases[day_index - 1] / population[curr_subnat_pos]);

            // print("adjust_factor[", day_index , "] = ", adjust_factor);
            // print("Rt[", day_index, "] = ", Rt[curr_day_pos]);
            // print("new_cases[:] = ", new_cases[days_pos:(curr_day_pos - 1)]);
            // print("gen_factor[:] = ", gen_factor[1:(day_index - 1)]);

            new_cases[curr_day_pos] = adjust_factor * Rt[curr_day_pos] * new_cases[days_pos:(curr_day_pos - 1)] * gen_factor[1:(day_index - 1)];
          }

          // print("new_cases[", day_index, "] = ", new_cases[curr_day_pos]);

          cumulative_cases[day_index] = new_cases[curr_day_pos] + (day_index > 1 ? cumulative_cases[day_index - 1] : 0);

          if (day_index > 1) {
            mean_deaths[curr_day_pos] = ifr[curr_subnat_pos] * new_cases[days_pos:(curr_day_pos - 1)] * time_to_death[1:(day_index - 1)];
          } else {
            mean_deaths[curr_day_pos] = 1e-15 * new_cases[curr_day_pos];
          }

          if (mean_deaths[curr_day_pos] < 0) {
            // print("mean_deaths[", curr_day_pos, "] = ", mean_deaths[curr_day_pos]);
            // print("day_index = ", day_index);
            // print("new_cases = ", new_cases[days_pos:(curr_day_pos - 1)]);
          }
        }

        days_pos = days_end + 1;
      }

      subnat_pos = subnat_end + 1;
    }
  }
}

model {
  // These are not model as hierarchical over countries yet.
  overdisp_deaths ~ normal(0, 5);
  ifr_noise ~ normal(1, 0.1);

  tau_impute_cases ~ exponential(0.03);
  imputed_cases ~ exponential(1 / tau_impute_cases); // TODO imputed cases should vary by subnational unit

  beta_toplevel ~ normal(0, 0.5);
  beta_national_sd ~ normal(0, 0.5);
  to_vector(beta_national_raw) ~ std_normal();
  to_vector(beta_subnational_sd) ~ normal(0, 0.5);
  to_vector(beta_subnational_raw) ~ std_normal();

  R0_sd ~ normal(0, 0.5);
  R0 ~ normal(3.28, R0_sd); // Again, not properly heirarchical.

  if (fit_model) {
    int subnat_pos = 1;
    int days_pos = 1;

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int subnat_end = subnat_pos + num_subnat - 1;

      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_subnat_pos = subnat_pos + subnat_index - 1;
        int days_end = days_pos + days_observed[curr_subnat_pos] - 1;

        deaths[(days_pos + start_epidemic_offset - 1):days_end] ~ neg_binomial_2(mean_deaths[(days_pos + start_epidemic_offset - 1):days_end], overdisp_deaths[curr_subnat_pos]);

        days_pos = days_end + 1;
      }

      subnat_pos = subnat_end + 1;
    }
  }
}

generated quantities {
  // Forecasting
}
