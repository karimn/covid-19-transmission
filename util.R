quantilize <- function(iter_data, var, quants = c(.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)) {
  iter_data %>%
    pull({{ var }}) %>%
    quantile(probs = quants, names = FALSE) %>%
    enframe(name = NULL, value = "est") %>%
    mutate(per = quants) %>%
    pivot_wider(names_from = per, values_from = est, names_prefix = "per_")
}

plot_subnat_ci <- function(results) {
  results %>%
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
    ggplot(aes(x = day_index)) +
    geom_line(aes(y = per_0.5)) +
    geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9), alpha = 0.25) +
    labs(x = "", y = "") +
    facet_wrap(vars(sub_region), ncol = 2, strip.position = "left") +
    theme(strip.placement = "outside", strip.text = element_text(angle = 0))
}
