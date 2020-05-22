functions {

}

data {
  int<lower = 1> N_national; // Number of countries
  int<lower = 1> N_subnational[N_national]; // Number of subnational entities for each country

  int<lower = 1> days_observed[sum(N_subnational)]; // Days observed per subnational entity can be different

  int<lower = 0> deaths[sum(days_observed)];

  vector<lower = 0>[sum(N_subnational)] mean_ifr; // Not sure we'd have this at the subnational level. If we don't we change this to the national level.
}

transformed data {
  int N = sum(N_subnational);
  int D = sum(days_observed);
  int max_days_observed = max(days_observed);

  real gen_factor_alpha = 6.5;
  real gen_factor_beta = 0.62;

  // Need to discretize these two
  vector<lower = 0>[max_days_observed] time_to_death; // Ugh, this is a sum of two gamma distributions
  vector<lower = 0>[max_days_observed] gen_factor;

  gen_factor[1] = gamma_cdf(1.5, gen_factor_alpha, gen_factor_beta);

  // Discretization
  for (day_index in 2:max_days_observed) {
    gen_factor[day_index] = gamma_cdf(day_index + 0.5, gen_factor_alpha, gen_factor_beta) - gen_factor[day_index - 1];
  }
}

parameters {
  vector<lower = 0>[D] mean_deaths; // Not a matrix; this could be a ragged data structure

  vector<lower = 0>[N] overdisp_deaths;

  vector<lower = 0>[N] ifr_noise;

  vector<lower = 0>[N] R0;
  real<lower = 0> R0_sd;
}

transformed parameters {
  vector<lower = 0>[N] ifr = mean_ifr .* ifr_noise;
  vector[D] R_predictors; // TODO Linear mobility model
  vector<lower = 0>[D] Rt;

  // int<lower = 0> new_cases[sum(days_observed)];

  {
    int subnat_pos = 1;
    int days_pos = 1;

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int subnat_end = subnat_pos + num_subnat - 1;

      for (subnat_index in 1:N_subnational[country_index]) {
        int days_end = days_pos + days_observed[subnat_pos + subnat_index - 1] - 1;

        Rt[days_pos:days_end] = 2 * R0[subnat_pos + subnat_index - 1] * inv_logit(R_predictors[days_pos:days_end]);

        // TODO: infection process
        // TODO: mean deaths process

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

  R0_sd ~ normal(0, 0.5);
  R0 ~ normal(3.28, R0_sd); // Again, not properly heirarchical.

  // TODO: Negative binomial likelihood for observed deathns
}

generated quantities {
}
