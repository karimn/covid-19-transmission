source(file.path("ifr", "data_prep.R"))

ifr      <- c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100
#IFR needs upper and lower bounds (and double-checking)

# pbc = pop by country (WPP 2019)
head(pbc_spread)

countries

# Assuming same attack rates in all age groups
ifr_adj <- apply(pbc_spread, 1, function(x) sum(as.numeric(x)*ifr)/sum(as.numeric(x)))
names(ifr_adj) <- names(countries)

# Can also do this:
# ar_vec <- rep(1, 9)
# apply(pbc_spread, 1, function(x) sum(as.numeric(x)*ifr)/sum(as.numeric(x)))

hist(ifr_adj)

