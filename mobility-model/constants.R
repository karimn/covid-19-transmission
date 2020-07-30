min_deaths_day_before_epidemic <- 3
seeding_days_before_epidemic <- 30
start_epidemic_offset <- seeding_days_before_epidemic + 1
days_seeding <- 6
days_to_forecast <- 0
mob_formula <- ~ 0 + g_residential + g_transit_stations + average_mob