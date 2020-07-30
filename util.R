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

extract_parameters <- function(fit, ...) {
  fit %>%
    as.array(...) %>%
    plyr::adply(3, diagnose) %>%
    as_tibble()
}
