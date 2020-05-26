library(tidyverse)
library(socialmixr)
default_cm <- contact_matrix(polymod, age.limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80))$matrix


library(wpp2019)
data(pop)

age_group_names <- rownames(default_cm)

pop_by_country <- rbind(data.frame(popF[c("country_code", "name", "age", "2020")], sex = "F"),
                        data.frame(popM[c("country_code", "name", "age", "2020")], sex = "M")) %>%
  setNames(c("code", "country", "age", "pop", "sex")) %>%
  mutate(age = forcats::fct_collapse(age, "[0,10)" = c("0-4", "5-9"),
                                     "[10,20)" = c("10-14", "15-19"),
                                     "[20,30)" = c("20-24", "25-29"),
                                     "[30,40)" = c("30-34", "35-39"),
                                     "[40,50)" = c("40-44", "45-49"),
                                     "[50,60)" = c("50-54", "55-59"),
                                     "[60,70)" = c("60-64", "65-69"),
                                     "[70,80)" = c("70-74", "75-79"),
                                     "80+" =   c("80-84", "85-89", "90-94", "95-99", "100+"))) %>%
  mutate(age = factor(age, levels = age_group_names)) %>%
  group_by(code, country, age) %>%
  summarise(total = sum(pop*1000)) %>%
  ungroup()

pbc_spread <- pop_by_country %>% spread(age, total) %>% ungroup()
countries <- as.character(pbc_spread$code)
names(countries) <- pbc_spread$country
pbc_spread <- pbc_spread %>% select(-country) %>% column_to_rownames("code")






# Johns Hopkins dataset on counts and deaths -----

path <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
filepaths <- list(
  "D"    = paste0(path, "time_series_covid19_deaths_global.csv"),
  "cases"     = paste0(path, "time_series_covid19_confirmed_global.csv"),
  "R" = paste0(path, "time_series_covid19_recovered_global.csv")
)
csv_df <- lapply(filepaths, function(cpath) {
  read_csv(cpath) %>%
    gather(date, value, -`Province/State`, -`Country/Region`, -Lat, -Long) %>%
    setNames(c("province", "country", "lat", "long", "date", "value"))
}) %>%
  bind_rows(.id = "variable")


cases_csv_clean <- csv_df %>%
  mutate(country = fct_recode(country, "Viet Nam" = "Vietnam",
                              "United States of America" = "US",
                              "Venezuela (Bolivarian Republic of)" = "Venezuela",
                              "United Republic of Tanzania" = "Tanzania")) %>%
  mutate(date = as.Date(date, format="%m/%d/%y")) %>%
  group_by(variable, country, date) %>%
  summarise(value = sum(value)) %>%
  rename(time = date) %>%
  mutate(pop = rowSums(pbc_spread[countries[as.character(country)],]))


# Exploration graphs
cases_csv_clean %>%
  mutate(value = (value+1)/pop) %>%
  ggplot(aes(x=time, y=value, group=country)) + geom_line() + scale_y_log10() +
  facet_wrap(~variable, scales = "free")


# Save all -----
# save(countries, pbc_spread, default_cm, cases_csv_clean, file = "data/demo_inputs.Rdata")
