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
  int<lower = 1> N_national; // Number of countries
  int<lower = 1> N_subnational[N_national]; // Number of subnational entities for each country

  int<lower = 1> start_epidemic_offset;

  // Days observed per national/subnational entity can be different
  int<lower = 1> days_observed[sum(N_subnational)]; 
  // Same for all national/subnational entities. We can change that to vary if needed.
  int<lower = 1, upper = min(days_observed)> days_to_impute_cases; 
  int<lower = 0> days_to_forecast;

  //For now let's try only one region where cases can originate per country
  int<lower = 1> travel_delay;
  real<lower=0> import_pr_sd;

  // This is pi in the Vollmer model. Needs to be pre-computed empirically using ecdf(); 
  //it's the sum of two gamma distributions
  vector<lower = 0>[max(days_observed) + days_to_forecast] time_to_death; 

  int<lower = 0> num_coef; // Number of coefficients for linear model
  // E.g., mobility
  matrix[sum(days_observed) + sum(N_subnational) * days_to_forecast, num_coef] design_matrix; 

  int<lower = 0> population[sum(N_subnational)];
  vector<lower = 0>[sum(N_subnational)] mean_ifr; 
  int<lower = 0> deaths[sum(days_observed)];

  // Hyperparameters
  real hyperparam_tau_beta_toplevel;
  real hyperparam_tau_beta_national_sd;
  real hyperparam_tau_beta_subnational_sd;
}

transformed data {
  int N = sum(N_subnational);
  int D = sum(days_observed);
  int total_days[N] = add_array(days_observed, days_to_forecast);
  int D_total = sum(total_days);
  int max_days_observed = max(days_observed) + days_to_forecast;
  int is_multinational = N_national > 1;

  real toplevel_R0_mean = 3.28;
  real toplevel_R0_sd = 0.5;

  real toplevel_log_R0_sd = sqrt(log((toplevel_R0_sd^2 / toplevel_R0_mean^2) + 1));
  real toplevel_log_R0_mean = log(toplevel_R0_mean) - (toplevel_log_R0_sd^2 / 2);

  int num_likelihood_days[N] = add_array(days_observed, 1 - start_epidemic_offset);
  int total_num_likelihood_days = sum(num_likelihood_days);
  int likelihood_day_idx[total_num_likelihood_days];
  int day_subnat_idx[total_num_likelihood_days];

  vector<lower = 0>[max_days_observed] rev_time_to_death;
  vector<lower = 0>[max_days_observed] rev_gen_factor;

  {
    real gen_factor_alpha = 1 / (0.62^2);              // alpha = 1 / tau^2;
    real gen_factor_beta = gen_factor_alpha / 6.5;     // beta = 1 / (tau^2 * mu)

    // Need to discretize
    vector[max_days_observed] gen_factor;
    vector[max_days_observed] gamma_cdfs;

    gamma_cdfs[1] = gamma_cdf(1.5, gen_factor_alpha, gen_factor_beta);
    gen_factor[1] = gamma_cdfs[1];

    // Discretization
    for (day_index in 2:max_days_observed) {
      gamma_cdfs[day_index] = gamma_cdf(day_index + 0.5, gen_factor_alpha, gen_factor_beta);
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
}

parameters {
  // vector<lower = 0>[N] overdisp_deaths;
  real<lower = 0> overdisp_deaths;

  real<lower = 0> tau_impute_cases;
  vector<lower = 0>[N] imputed_cases;

  // Linear model parameters
  vector[num_coef] beta_toplevel;

  // Separate SD for each parameters. Vollmer et al. use same for all parameters
  vector<lower = 0>[is_multinational ? num_coef : 0] beta_national_sd; 
  // Uncentered for now to avoid divergence
  matrix[num_coef, is_multinational ? N_national : 0] beta_national_raw; 

  matrix<lower = 0>[num_coef, N_national] beta_subnational_sd;
  matrix[num_coef, N] beta_subnational_raw;

  vector<lower = 0>[N] ifr_noise;

  vector<lower = 0>[N] R0_raw;
  real<lower = 0> R0_sd;
  
  vector<lower = 0, upper = 1>[N] import_pr;
}

transformed parameters {
  vector<lower = 0>[N] R0 = 3.28 + R0_raw * R0_sd;
  // Not a matrix; this could be a ragged data structure
  vector[D_total] mean_deaths = rep_vector(0, D_total); 
  vector[D_total] Rt = rep_vector(0, D_total);
  vector[D_total] Rt_adj = Rt;
  row_vector[D_total] new_cases = rep_row_vector(0, D_total);

  matrix[num_coef, N] beta = rep_matrix(beta_toplevel, N);

  {
    int subnat_pos = 1;
    int days_pos = 1;

    vector[N] ifr = mean_ifr .* ifr_noise;

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int subnat_end = subnat_pos + num_subnat - 1;
      real center_cases[total_days[subnat_pos]];

      if (is_multinational) {
        beta[, subnat_pos:subnat_end] += rep_matrix(beta_national_raw[, country_index] .* beta_national_sd, num_subnat);
      }

      
      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_subnat_pos = subnat_pos + subnat_index - 1;
        int days_end = days_pos + days_observed[curr_subnat_pos] + days_to_forecast - 1;

        vector[total_days[curr_subnat_pos]] cumulative_cases = rep_vector(0, total_days[curr_subnat_pos]);

        vector[total_days[curr_subnat_pos]] mobility_effect = rep_vector(0, total_days[curr_subnat_pos]);

        beta[, curr_subnat_pos] += beta_subnational_raw[, curr_subnat_pos] .* beta_subnational_sd[, country_index];
        
          new_cases[days_pos:(days_pos + days_to_impute_cases - 1)] = 
            rep_row_vector(imputed_cases[subnat_index], days_to_impute_cases);
        if(subnat_index == 1){
            for(i in 1:days_to_impute_cases)
              center_cases[i] = new_cases[days_pos + i - 1];
        }
        // Set to 0 for now
        // else
          // new_cases[days_pos:(days_pos + days_to_impute_cases - 1)] = 
            // rep_row_vector(0, days_to_impute_cases);

        mobility_effect = 2 * inv_logit(design_matrix[days_pos:days_end] * beta[, curr_subnat_pos]);

        Rt[days_pos:days_end] = R0[curr_subnat_pos] * mobility_effect;

        for (day_index in 1:total_days[curr_subnat_pos]) {
          int curr_day_pos = days_pos + day_index - 1;
          if (day_index > days_to_impute_cases) {
            real adjust_factor = 1 - (cumulative_cases[day_index - 1] / population[curr_subnat_pos]);
            Rt_adj[curr_day_pos] = adjust_factor * Rt[curr_day_pos];
            
            new_cases[curr_day_pos] = 
              Rt_adj[curr_day_pos] * new_cases[days_pos:(curr_day_pos - 1)] * tail(rev_gen_factor, day_index - 1);
            if(subnat_index == 1)
              center_cases[day_index] = new_cases[curr_day_pos];
            else
              new_cases[curr_day_pos] = new_cases[curr_day_pos] + 
                                        import_pr[subnat_index]*center_cases[day_index-travel_delay]*
                                        (Rt_adj[curr_day_pos]/R0[curr_subnat_pos]); //fix this
                                        
            // print("subnat", subnat_index, "center_cases[", day_index, "] = ", center_cases[day_index]);
          } else {
            Rt_adj[curr_day_pos] = Rt[curr_day_pos];
          }
          // print("center_cases[", day_index, "] = ", new_cases[curr_day_pos]);
          
          cumulative_cases[day_index] = 
            new_cases[curr_day_pos] + (day_index > 1 ? cumulative_cases[day_index - 1] : 0);
          if (day_index > 1) {
            mean_deaths[curr_day_pos] = 
              ifr[curr_subnat_pos] * new_cases[days_pos:(curr_day_pos - 1)] * tail(rev_time_to_death, day_index - 1);
          } else {
            mean_deaths[curr_day_pos] = 1e-15 * new_cases[curr_day_pos];
          }
          if (mean_deaths[curr_day_pos] < 0) {
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
  imputed_cases ~ exponential(1 / tau_impute_cases);

  beta_toplevel ~ normal(0, hyperparam_tau_beta_toplevel);
  
  import_pr ~ normal(0, import_pr_sd);

  if (is_multinational) {
    beta_national_sd ~ normal(0, hyperparam_tau_beta_national_sd);
    to_vector(beta_national_raw) ~ std_normal();
  }

  to_vector(beta_subnational_sd) ~ normal(0, hyperparam_tau_beta_subnational_sd);
  to_vector(beta_subnational_raw) ~ std_normal();

  R0_sd ~ normal(0, 0.5);
  R0_raw ~ std_normal(); // Again, not properly heirarchical.

  {
    int subnat_pos = 1;
    int days_pos = 1;

    for (country_index in 1:N_national) {
      int num_subnat = N_subnational[country_index];
      int subnat_end = subnat_pos + num_subnat - 1;

      for (subnat_index in 1:N_subnational[country_index]) {
        int curr_subnat_pos = subnat_pos + subnat_index - 1;
        int days_end = days_pos + days_observed[curr_subnat_pos] - 1;

        deaths[(days_pos + start_epidemic_offset - 1):days_end] ~ neg_binomial_2(
          mean_deaths[(days_pos + start_epidemic_offset - 1):days_end], overdisp_deaths);

        days_pos = days_end + 1;
      }

      subnat_pos = subnat_end + days_to_forecast + 1;
    }
  }
}

