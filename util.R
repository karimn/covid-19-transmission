quantilize <- function(iter_data, var, quants = c(.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)) {
  iter_data %>%
    pull({{ var }}) %>%
    quantile(probs = quants, names = FALSE) %>%
    enframe(name = NULL, value = "est") %>%
    mutate(per = quants) %>%
    pivot_wider(names_from = per, values_from = est, names_prefix = "per_")
}

diagnose <- function(cell, no_sim_diag = FALSE) {
  tibble(iter_data = list(cell)) %>% {
    if (!no_sim_diag) {
      mutate(.,
             ess_bulk = ess_bulk(cell),
             ess_tail = ess_tail(cell),
             rhat = Rhat(cell))
    } else .
  }
}

extract_subnat_results <- function(fit, par, exp_logs = TRUE) {
  fit %>%
    as.array(par = par) %>%
    plyr::adply(3, diagnose) %>%
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
    unnest(quants) %>%
    nest(param_results = -subnat_index)
}

extract_day_results <- function(fit, par) {
  fit %>%
    as.array(par = par) %>%
    plyr::adply(3, diagnose) %>%
    tidyr::extract(parameters, c("parameter", "long_day_index"), ("(\\w+)\\[(\\d+)\\]"), convert = TRUE) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      countrycode_string = rep(rep(use_subnat_data$countrycode_string, times = stan_data$days_observed), length(par)),
      sub_region = rep(rep(use_subnat_data$sub_region, times = stan_data$days_observed), length(par)),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants) %>%
    nest(day_data = -c(countrycode_string, sub_region)) %>%
    mutate(
      day_data = map(day_data,
                     ~ group_by(., parameter) %>%
                       mutate(day_index = seq_along(long_day_index)) %>%
                       ungroup() %>%
                       select(-long_day_index) %>%
                       nest(param_results = -c(day_index)))
    )
}

extract_beta <- function(fit) {
  fit %>%
    as.array(par = "beta") %>%
    plyr::adply(3, diagnose) %>%
    tidyr::extract(parameters, c("coef_index", "subnat_index"), ("\\[(\\d+),(\\d+)\\]"), convert = TRUE) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_value = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      quants = map(iter_data, quantilize, iter_value),
      mean = map(iter_data, pull, iter_value) %>% map_dbl(mean),
    ) %>%
    unnest(quants) %>%
    nest(param_results = -subnat_index)

}

plot_subnat_ci <- function(results) {
  results %>%
    unnest(param_results) %>%
    filter(fct_match(parameter, "R0")) %>%
    mutate(sub_region = fct_reorder(sub_region, per_0.5)) %>%
    ggplot(aes(y = sub_region)) +
    geom_pointrange(aes(x = per_0.5, xmin = per_0.1, xmax = per_0.9), fatten = 1.5) +
    geom_point(aes(x = mean), size = 1.5, shape = 5) +
    labs(x = "", y = "") +
    facet_wrap(vars(countrycode_string), ncol = 1)
}

plot_day_ci <- function(results) {
  results %>%
    select(countrycode_string, sub_region, daily_data) %>%
    unnest(daily_data) %>%
    select(countrycode_string, sub_region, day_index, param_results) %>%
    unnest(param_results) %>%
    filter(fct_match(parameter, "Rt_adj")) %>%
    ggplot(aes(x = day_index)) +
    geom_line(aes(y = per_0.5)) +
    geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9), alpha = 0.25) +
    labs(x = "", y = "") +
    facet_wrap(vars(sub_region), ncol = 2, strip.position = "left") +
    theme(strip.placement = "outside", strip.text = element_text(angle = 0))
}

plot_last_day_ci <- function(results) {
  results %>%
    select(countrycode_string, sub_region, daily_data) %>%
    mutate(daily_data = map(daily_data, filter, min_rank(date) == n())) %>%
    unnest(daily_data) %>%
    select(countrycode_string, sub_region, day_index, param_results) %>%
    unnest(param_results) %>%
    filter(fct_match(parameter, "Rt_adj")) %>%
    mutate(sub_region = fct_reorder(sub_region, per_0.5)) %>%
    ggplot(aes(y = sub_region)) +
    geom_pointrange(aes(x = per_0.5, xmin = per_0.1, xmax = per_0.9), fatten = 1.5) +
    geom_point(aes(x = mean), size = 1.5, shape = 5) +
    labs(x = "", y = "") +
    facet_wrap(vars(countrycode_string), ncol = 1)
}
