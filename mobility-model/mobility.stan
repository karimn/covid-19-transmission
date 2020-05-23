functions {

}

data {
  int<lower = 0, upper = 1> fit_model; // Fit vs prior-predict

  int<lower = 1> N_national; // Number of countries
  int<lower = 1> N_subnational[N_national]; // Number of subnational entities for each country

  int<lower = 1> days_observed[sum(N_subnational)]; // Days observed per subnational entity can be different
  int<lower = 1, upper = min(days_observed)> days_to_impute_cases; // Same for all national/subnational entities. We can change that to vary if needed.
  int<lower = max(days_observed)> total_days; // days observed + days to predict to

  int<lower = 0> population[sum(N_subnational)];

  int<lower = 0> deaths[sum(days_observed)];

  vector<lower = 0>[sum(N_subnational)] mean_ifr; // Not sure we'd have this at the subnational level. If we don't we change this to the national level.

  vector<lower = 0>[max(days_observed)] time_to_death; // This is pi in the Vollmer model. Needs to be pre-computed using ecdf(); it's the sum of two gamma distributions

  int<lower = 0> num_coef; // Number of coefficients for linear model

  matrix[total_days, num_coef] design_matrix[sum(N_subnational)];
}

transformed data {
  int N = sum(N_subnational);
  int D = sum(days_observed);
  int D_total = N * total_days;
  int max_days_observed = max(days_observed);

  real gen_factor_alpha = 6.5;
  real gen_factor_beta = 0.62;

  // Need to discretize
  vector<lower = 0>[max_days_observed] gen_factor;

  gen_factor[1] = gamma_cdf(1.5, gen_factor_alpha, gen_factor_beta);

  // Discretization
  for (day_index in 2:max_days_observed) {
    gen_factor[day_index] = gamma_cdf(day_index + 0.5, gen_factor_alpha, gen_factor_beta) - gen_factor[day_index - 1];
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
  matrix<lower = 0>[N_national, num_coef] beta_subnational_sd;
  matrix[num_coef, N] beta_subnational_raw;

  vector<lower = 0>[N] ifr_noise;

  vector<lower = 0>[N] R0;
  real<lower = 0> R0_sd;
}

transformed parameters {
  matrix<lower = 0>[total_days, N] mean_deaths; // Not a matrix; this could be a ragged data structure
  matrix<lower = 0>[total_days, N] Rt = rep_matrix(0, total_days, N);
  matrix<lower = 0>[N, total_days] new_cases = rep_matrix(0, N, total_days);

  {
    int subnat_pos = 1;

    vector[N] ifr = mean_ifr .* ifr_noise;

    matrix[num_coef, N] beta = rep_matrix(beta_toplevel, N);

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int subnat_end = subnat_pos + num_subnat - 1;

      beta[, subnat_pos:subnat_end] += rep_matrix(beta_national_raw[, country_index] .* beta_national_sd, num_subnat);

      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_subnat_pos = subnat_pos + subnat_index - 1;

        vector[total_days] cumulative_cases;

        beta[, curr_subnat_pos] += beta_subnational_raw[, curr_subnat_pos] .* beta_subnational_sd[, country_index];

        new_cases[curr_subnat_pos] = rep_row_vector(imputed_cases[subnat_index], days_to_impute_cases);

        Rt[, curr_subnat_pos] = 2 * R0[curr_subnat_pos] * inv_logit(design_matrix[curr_subnat_pos] * beta[, curr_subnat_pos]);

        for (day_index in 1:total_days) {
          if (day_index > days_to_impute_cases) {
            real adjust_factor = 1 - (cumulative_cases[day_index - 1] / population[curr_subnat_pos]);

            new_cases[curr_subnat_pos, day_index] = adjust_factor * Rt[day_index, curr_subnat_pos] * (new_cases[curr_subnat_pos, 1:(day_index - 1)] * gen_factor[1:(day_index - 1)]);
          }

          cumulative_cases[day_index] = new_cases[curr_subnat_pos, day_index] + (day_index > 1 ? cumulative_cases[day_index - 1] : 0);
          mean_deaths[day_index, curr_subnat_pos] = ifr[curr_subnat_pos] * (new_cases[curr_subnat_pos, 1:(day_index - 1)] * time_to_death[1:(day_index - 1)]);
        }
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
  imputed_cases ~ exponential(1 / tau_impute_cases);

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

        deaths[days_pos:days_end] ~ neg_binomial_2(mean_deaths[1:days_observed[curr_subnat_pos], curr_subnat_pos], overdisp_deaths[curr_subnat_pos]);

        days_pos = days_end + 1;
      }

      subnat_pos = subnat_end + 1;
    }
  }
}

generated quantities {

}
