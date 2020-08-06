prepare_subnat_data <- function(raw_data_file, min_deaths) {
  subnat_data_raw <- read_csv(
    raw_data_file,
    col_types = cols("new_recovered" = col_integer(),
                     "total_recovered" = col_integer(),
                     "subregion1_name" = col_character(),
                     "subregion1_code" = col_character(),
                     "subregion2_code" = col_character(),
                     "subregion2_name" = col_character())
  )

  subnat_data <- subnat_data_raw %>%
    left_join(transmute(countrycode::codelist, country_code = iso2c, countrycode_iso3n = iso3n), by = "country_code") %>% # Getting the countrycode
    rename_at(vars(matches("^(mobility|total)_")), str_replace_all, c("^mobility_" = "g_", "^total_" = "cum_")) %>%
    rename_at(vars(matches("(confirmed|deceased)$")), str_replace_all, c("confirmed$" = "cases", "deceased$" = "deaths")) %>%
    mutate_at(vars(starts_with("g_")), ~ . / 100) %>%
    mutate(
      sub_region = case_when(aggregation_level == 1 ~ subregion1_name,
                             aggregation_level == 2 ~ subregion2_name),

      all_mob_observed = select(., starts_with("g_")) %>%
        map(is.na) %>%
        transpose() %>%
        map_lgl(~ all(unlist(.x))) %>%
        not(),
    )

  clean_missing_spread <- function(daily_data, first_observed_death) {
    if (is.finite(first_observed_death)) {
      daily_data %>%
        mutate(
          new_deaths = coalesce(new_deaths, 0),
          cum_deaths = if_else(is.na(cum_deaths) & date < first_observed_death, 0, cum_deaths),
          new_deaths = coalesce(new_deaths, cum_deaths - lag(cum_deaths, default = 0)),
          cum_deaths = coalesce(cum_deaths, cumsum(new_deaths)),

          new_cases = coalesce(new_cases, 0),
          cum_cases = zoo::na.locf0(cum_cases) %>% coalesce(0),
          stringency_index = coalesce(stringency_index,0)
        )
    } else return(daily_data)
  }

  clean_date_ranges <- function(daily_data, first_infection_seeding_day, last_effective_observed_day, is_valid) {
    if (is_valid) {
      daily_data %>%
        filter(between(date, first_infection_seeding_day, last_effective_observed_day)) %>%
        complete(date = seq.Date(from = first_infection_seeding_day, to = max(date), by = "day"))
    } else {
      return(daily_data)
    }
  }

  subnat_data %<>%
    group_by_at(vars(starts_with("country"), starts_with("subregion"), sub_region, key, aggregation_level, population)) %>%
    group_nest(.key = "daily_data") %>%
    mutate(
      first_observed_death = map_dbl(daily_data,
                                     ~ filter(.x, !is.na(cum_deaths) | !is.na(new_deaths)) %>%
                                       pull(date) %>%
                                       min()) %>%
        lubridate::as_date(),

      last_observed_death = map_dbl(daily_data,
                                   ~ filter(.x, !is.na(cum_deaths) | !is.na(new_deaths)) %>%
                                     pull(date) %>%
                                     max()) %>%
        lubridate::as_date(),

      first_case_day = map_dbl(daily_data,
                         ~ filter(.x, cum_cases > 0) %>%
                           pull(date) %>%
                           min()) %>%
        lubridate::as_date(),

      first_mob_day = map_dbl(daily_data,
                         ~ filter(.x, all_mob_observed) %>%
                           pull(date) %>%
                           min()) %>%
        lubridate::as_date(),

      last_mob_day = map_dbl(daily_data,
                         ~ filter(.x, all_mob_observed) %>%
                           pull(date) %>%
                           max()) %>%
        lubridate::as_date(),

      is_valid = is.finite(first_observed_death) & is.finite(first_mob_day) & !is.na(population),

      daily_data = map2(daily_data, first_observed_death, clean_missing_spread),

      epidemic_start_date = map_dbl(daily_data,
                                    ~ filter(.x, cum_deaths >= min_deaths) %>% # First day with >= 10 deaths
                                      pull(date) %>%
                                      min()) %>%
        lubridate::as_date() %>%
        add(1), # The day after

      first_infection_seeding_day = epidemic_start_date - 1 - seeding_days_before_epidemic, # 30 days before the first day with >=10 deaths

      last_effective_observed_day = pmin(last_observed_death, last_mob_day), # We need both deaths and mobility day for the tail end of observed data (used in likelihood)
      num_days_observed = last_effective_observed_day - first_infection_seeding_day + 1,
    ) %>%
    mutate(
      has_epidemic = is_valid & is.finite(epidemic_start_date),

      daily_data = pmap(lst(daily_data, first_infection_seeding_day, last_effective_observed_day, is_valid = has_epidemic), clean_date_ranges) %>%
        map2(first_observed_death, clean_missing_spread) %>%
        map(mutate_at, vars(starts_with("g_")), zoo::na.locf, na.rm = FALSE) %>% # Any non-leading NAs are replaced with last prior non-NA
        map(mutate_at, vars(starts_with("g_")), coalesce, 0) %>%  # Any leading NAs are replaced with zeroes (as in Vollmer et al.)
        map(mutate, average_mob = (g_grocery_and_pharmacy + g_parks + g_retail_and_recreation + g_workplaces) / 4,
                    average_all_mob = (g_transit_stations + g_grocery_and_pharmacy + g_parks + g_retail_and_recreation + g_workplaces) / 5) %>%
        map(arrange, date) %>%
        map(mutate, day_index = seq(n())),

      first_case_day_index = map_dbl(daily_data,
                         ~ filter(.x, cum_cases > 0) %>%
                           pull(day_index) %>%
                           min()),

      total_deaths = map(daily_data, pull, new_deaths) %>% map(as.integer) %>% map_int(sum),

      country_index = group_indices(., country_code) # For SLURM runs on cluster server
    )

  return(subnat_data)
}

extract_toplevel_results <- function(fit, par, exp_logs = TRUE) {
  fit %>%
    extract_parameters(par = par) %>%
    tidyr::extract(parameters, c("parameter", "coef_index"), ("(\\w+)(?:\\[(\\d+)\\])?"), convert = TRUE) %>%
    bind_rows(
      filter(., str_detect(parameter, "log")) %>%
        mutate(iter_data = map(iter_data, exp),
               parameter = str_remove(parameter, "log_?"))
    ) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants)
}

extract_hyper_results <- function(fit, par, exp_logs = TRUE) {
  fit %>%
    extract_parameters(par = par) %>%
    tidyr::extract(parameters, c("parameter"), ("(\\w+)"), convert = TRUE) %>%
    bind_rows(
      filter(., str_detect(parameter, "log")) %>%
        mutate(iter_data = map(iter_data, exp),
               parameter = str_remove(parameter, "log_?"))
    ) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants)
}

extract_nat_results <- function(fit, par, exp_logs = TRUE) {
  fit %>%
    extract_parameters(par = par) %>%
    tidyr::extract(parameters, c("parameter", "run_country_index"), ("(\\w+)\\[(\\d+)\\]"), convert = TRUE) %>%
    bind_rows(
      filter(., str_detect(parameter, "log")) %>%
        mutate(iter_data = map(iter_data, exp),
               parameter = str_remove(parameter, "log_?"))
    ) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants)
}

extract_subnat_results <- function(fit, par, exp_logs = TRUE) {
  fit %>%
    extract_parameters(par = par) %>%
    tidyr::extract(parameters, c("parameter", "subnat_index"), ("(\\w+)\\[(\\d+)\\]"), convert = TRUE) %>%
    bind_rows(
      filter(., str_detect(parameter, "log")) %>%
        mutate(iter_data = map(iter_data, exp),
               parameter = str_remove(parameter, "log_?"))
    ) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants)
}

extract_day_results <- function(fit, par, use_subnat_data) {
  days_observed <- use_subnat_data %$% map_int(daily_data, nrow)

  fit %>%
    extract_parameters(par = par) %>%
    tidyr::extract(parameters, c("parameter", "long_day_index"), ("(\\w+)\\[(\\d+)\\]"), convert = TRUE) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      country_code = rep(rep(use_subnat_data$country_code, times = days_observed), length(par)),
      sub_region = rep(rep(use_subnat_data$sub_region, times = days_observed), length(par)),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants) %>%
    nest(day_data = -c(country_code, sub_region)) %>%
    mutate(
      day_data = map(day_data,
                     ~ group_by(., parameter) %>%
                       mutate(day_index = seq_along(long_day_index)) %>%
                       ungroup() %>%
                       select(-long_day_index) %>%
                       nest(param_results = -c(day_index)))
    )
}

extract_hyper_beta <- function(fit) {
  fit %>%
    extract_parameters(par = c("beta_toplevel", "beta_national_sd")) %>%
    tidyr::extract(parameters, c("parameter", "coef_index"), ("(\\w+)\\[(\\d+)\\]"), convert = TRUE) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants)
}

extract_nat_beta <- function(fit) {
  fit %>%
    extract_parameters(par = "beta_subnational_sd") %>%
    tidyr::extract(parameters, c("parameter", "coef_index", "run_country_index"), ("(\\w+)\\[(\\d+),(\\d+)\\]"), convert = TRUE) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants)
}

extract_subnat_beta <- function(fit) {
  fit %>%
    extract_parameters(par = "beta") %>%
    tidyr::extract(parameters, c("parameter", "coef_index", "subnat_index"), ("(\\w+)\\[(\\d+),(\\d+)\\]"), convert = TRUE) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants)
}

extract_results <- function(use_subnat_data, mob_fit, hyper_params, nat_params, subnat_params, day_params) {
  hyper_results <- mob_fit %>%
    extract_hyper_results(hyper_params) %>%
    bind_rows(
      extract_hyper_beta(mob_fit)
    )

  multi_national <- n_distinct(use_subnat_data$country_code) > 1
  multi_subregion <- use_subnat_data %>%
    count(country_code) %>%
    pull(n) %>%
    max() %>%
    is_greater_than(1)

  nat_beta_results <- if (multi_national && multi_subregion) {
    mob_fit %>%
      extract_nat_beta()
  }

  nat_results <- mob_fit %>%
    extract_nat_results(nat_params) %>%
    bind_rows(nat_beta_results) %>%
    nest(param_results = -run_country_index)

  beta_results <- mob_fit %>%
    extract_subnat_beta()

  subnat_results <- mob_fit %>%
    extract_subnat_results(subnat_params) %>%
    bind_rows(beta_results) %>%
    nest(param_results = -subnat_index)

  day_results <- mob_fit %>%
    extract_day_results(day_params, use_subnat_data)

  use_subnat_data %>%
    mutate(
      subnat_index = seq(n()),
    ) %>%
    select(-any_of("param_results")) %>%
    left_join(subnat_results, by = "subnat_index") %>%
    left_join(day_results, by = c("country_code", "sub_region")) %>%
    mutate(
      daily_data = map2(daily_data, day_data,
                        ~ select(.x, -any_of("param_results")) %>%
                          left_join(.y, by = "day_index"))
    ) %>%
    select(-day_data) %>%
    nest(country_data = -c(country_index, country_code, country_name, countrycode_iso3n)) %>%
    mutate(run_country_index = seq(n())) %>%
    select(-any_of("param_results")) %>%
    left_join(nat_results, by = "run_country_index") %>%
    select(-run_country_index) %>%
    nest(run_data = everything()) %>%
    mutate(
      param_results = list(hyper_results)
    )
}

plot_subnat_ci <- function(results, par, beta = FALSE, violin_density = FALSE) {
  dodge_width <- if (violin_density) 1 else 0.5

  plot_obj <- results %>% {
    if (beta) {
      unnest(., beta_results) %>%
        mutate(parameter = str_c("$\\beta_", coef_index, "$"))
    } else {
      unnest(., param_results) %>%
        filter(fct_match(parameter, par))
    }
  } %>%
    mutate(sub_region = fct_reorder(sub_region, per_0.5) %>% fct_explicit_na(na_level = "")) %>%
    ggplot(aes(y = sub_region))

  if (violin_density) {
    plot_obj <- plot_obj + geom_violin(aes(x = iter_value, group = parameter),
                                       position = position_dodge(width = dodge_width),
                                       data = . %>% unnest(iter_data))
  }

  plot_obj <- plot_obj +
    geom_pointrange(aes(x = per_0.5, xmin = per_0.1, xmax = per_0.9, color = parameter),
                    fatten = 1.5, position = position_dodge(width = dodge_width), show.legend = beta || length(par) > 1) +
    geom_point(aes(x = mean, color = parameter), size = 1.5, shape = 5, position = position_dodge(width = dodge_width), show.legend = beta || length(par) > 1) +
    scale_color_discrete("", labels = if (beta) TeX else waiver()) +
    labs(x = "", y = "") +
    facet_wrap(vars(country_name), ncol = 2, scales = "free_y")

  return(plot_obj)
}

plot_population <- function(results, at_level = "sub_region") {
  at_level_sym <- sym(at_level)

  results %>%
    mutate(!!at_level_sym := fct_reorder(!!at_level_sym, population)) %>%
    ggplot(aes(population, !!at_level_sym)) +
    geom_col(alpha = 0.5) +
    scale_x_continuous("", labels = scales::label_number(scale = 1/1000, suffix = "K")) +
    labs(title = "Population", y = "")
}

plot_day_data <- function(results, par, use_date = FALSE, use_bar = FALSE, use_free_y_scale = FALSE, facet_by = "sub_region") {
  time_aes <- if (use_date) aes(x = date) else aes(x = day_index)
  y_aes <- if (length(par) > 1) aes(y = value, color = name) else aes(y = value)

  plot_obj <- results %>%
    select(country_code, all_of(facet_by), daily_data, first_infection_seeding_day, epidemic_start_date) %>%
    unnest(daily_data) %>%
    select(day_index, date, all_of(c(par, facet_by)), first_infection_seeding_day, epidemic_start_date) %>%
    pivot_longer(all_of(par)) %>%
    ggplot(time_aes) +
    geom_vline(aes(xintercept = first_infection_seeding_day), linetype = "dotted", data = . %>% distinct_at(vars(all_of(facet_by), first_infection_seeding_day))) +
    geom_vline(aes(xintercept = epidemic_start_date), linetype = "dotted", data = . %>% distinct_at(vars(all_of(facet_by), epidemic_start_date)))

  plot_obj <- if (use_bar) {
    plot_obj + geom_col(y_aes, color = "black", alpha = 0.5, size = 0.25)
  } else {
    plot_obj + geom_line(y_aes)
  }

  plot_obj +
    scale_color_discrete("") +
    labs(x = "", y = "",
         caption = "Vertical dotted lines represent the first seeding day and the epidemic start date.") +
    facet_wrap(facet_by, ncol = 3, strip.position = "left", scales = if (use_free_y_scale) "free_y" else "fixed") +
    theme(
      # strip.placement = "outside",
      # strip.text = element_text(angle = 0),
      axis.text.x = if (!use_date) element_blank()
    )
}

plot_day_ci <- function(results, par, use_date = FALSE, facet_by = "sub_region") {
  time_aes <- if (use_date) aes(x = date) else aes(x = day_index)

  results %>%
    select(country_code, all_of(facet_by), daily_data, first_infection_seeding_day, epidemic_start_date) %>%
    unnest(daily_data) %>%
    select(country_code, all_of(facet_by), day_index, date, param_results, first_infection_seeding_day, epidemic_start_date) %>%
    unnest(param_results) %>%
    filter(fct_match(parameter, par)) %>%
    ggplot(time_aes) +
    geom_vline(aes(xintercept = first_infection_seeding_day), linetype = "dotted", data = . %>% distinct_at(vars(all_of(facet_by), first_infection_seeding_day))) +
    geom_vline(aes(xintercept = epidemic_start_date), linetype = "dotted", data = . %>% distinct_at(vars(all_of(facet_by), epidemic_start_date))) +
    geom_line(aes(y = per_0.5, color = parameter), show.legend = length(par) > 1) +
    geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, group = parameter, fill = parameter), alpha = 0.25) +
    scale_color_discrete("") +
    scale_fill_discrete("") +
    labs(x = "", y = "",
         caption = "Vertical dotted lines represent the first seeding day and the epidemic start date.
                    Ribbons represent the 80% credible intervals.") +
    facet_wrap(facet_by, ncol = 3) + #, strip.position = "left") +
    theme(
      # strip.placement = "outside",
      # strip.text.y.left = element_text(angle = 0),
      axis.text.x = if (!use_date) element_blank()
    )
}

plot_last_day_ci <- function(results) {
  results %>%
    select(country_code, sub_region, daily_data) %>%
    mutate(daily_data = map(daily_data, filter, min_rank(date) == n())) %>%
    unnest(daily_data) %>%
    select(country_code, sub_region, day_index, param_results) %>%
    unnest(param_results) %>%
    filter(fct_match(parameter, "Rt_adj")) %>%
    mutate(sub_region = fct_reorder(sub_region, per_0.5)) %>%
    ggplot(aes(y = sub_region)) +
    geom_pointrange(aes(x = per_0.5, xmin = per_0.1, xmax = per_0.9), fatten = 1.5) +
    geom_point(aes(x = mean), size = 1.5, shape = 5) +
    labs(x = "", y = "") +
    facet_wrap(vars(country_code), ncol = 1)
}

plot_pred_vs_obs <- function(results, pred, obs, y_axis_lab = "", facet_by = "sub_region") {
  results %>%
    select(all_of(facet_by), daily_data, first_infection_seeding_day, epidemic_start_date) %>%
    unnest(daily_data) %>%
    select(day_index, date, all_of(c(facet_by, obs)), param_results, first_infection_seeding_day, epidemic_start_date) %>%
    unnest(param_results) %>%
    filter(fct_match(parameter, pred)) %>%
    ggplot(aes(x = date)) +
    geom_vline(aes(xintercept = first_infection_seeding_day), linetype = "dotted", data = . %>% distinct_at(vars(all_of(facet_by), first_infection_seeding_day))) +
    geom_vline(aes(xintercept = epidemic_start_date), linetype = "dotted", data = . %>% distinct_at(vars(all_of(facet_by), epidemic_start_date))) +
    geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9), alpha = 0.15) +
    geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75), alpha = 0.25) +
    geom_line(aes(y = per_0.5, color = "predicted median"), size = 1) +
    geom_line(aes(y = !!sym(obs), color = "observed"), size = 1) +
    scale_color_discrete("") +
    labs(x = "", y = y_axis_lab,
         caption = "Solid black line: observed new deaths. Grey ribbon: posterior predicted new deaths.
                    Vertical dotted lines represent the first seeding day and the epidemic start date.") +
    facet_wrap(facet_by, scales = "free_y", ncol = 3) +
    NULL
}

plot_post_deaths <- function(results, facet_by = "sub_region") {
  plot_pred_vs_obs(results, "deaths_rep", "new_deaths", "New Deaths", facet_by)
}

render_country_reports <- function(results,
                                   results_file = "lite_merged.rds",
                                   report_template = file.path("mobility-model", "mobility_report.Rmd"),
                                   reports_dir = file.path("mobility-model", "country-reports"),
                                   reports_id = NULL, ...) {
  extra_param <- rlang::list2(...)

  results %>%
    pull(country_code) %>%
    unique() %>%
    walk(~ {
      rmarkdown::render(report_template,
                        output_file = str_c(str_to_lower(.x), reports_id, "report.pdf", sep = "_"),
                        output_dir = reports_dir,
                        params = lst(country_code = .x, results_file, !!!extra_param),
                        knit_root_dir = "..")
    })

}

trim_iter_data <- function(results) {
  results %>%
    mutate(
      param_results = map(param_results, select, -iter_data),
      run_data = map(
        run_data,
        mutate,
        param_results = map(param_results, select, -iter_data),
        country_data = map(
          country_data,
          mutate,
          param_results = map(param_results, select, -iter_data),
          daily_data = map(daily_data, mutate, param_results = map(param_results, select, - iter_data))
        )
      )
    )
}
