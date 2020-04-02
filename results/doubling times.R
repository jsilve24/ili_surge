library(tidyverse)
library(driver)
library(lubridate)
library(lme4)

# Load state abbreviations
abbr <- read_csv("../../data/population/state_abbreviations.csv", col_names=FALSE)
abbr <- select(abbr, X1, X4)
tmp <- abbr$X1
names(tmp) <- abbr$X4
abbr <- tmp

# New York Times Dataset
us_confirmed_tidy <- read_csv("../../data/covid-19-data/us-states.csv") %>% 
  mutate(day = as.numeric(date - ymd("2020-01-01"))) %>% 
  mutate(day = day-min(day))


tmp <- us_confirmed_tidy %>% 
  split(.$state) %>% 
  map(~dplyr::filter(.x, cases > 5))
tmp <- tmp %>% 
  map(~glm(cases~day, family="poisson", data=.x)) %>% 
  map(summary) %>% 
  map(coefficients) %>% 
  map(~.x[2,1:2]) %>% 
  enframe() %>%
  mutate(r = map(value, ~.x[1]), 
         std.error = map(value, ~.x[2])) %>% 
  select(-value) %>% 
  unnest() %>% 
  mutate(doubling_time = log(2)/r) %>% 
  mutate(name = factor(name, levels=name[order(doubling_time)])) %>% 
  mutate(doubling_time_std.error = log(2)*abs(doubling_time)*std.error) # https://en.wikipedia.org/wiki/Propagation_of_uncertainty
median_doubling_time <- median(tmp$doubling_time)

tmp %>% 
  ggplot(aes(x=name, y=doubling_time)) +
  geom_point() +
  geom_errorbar(aes(ymin=doubling_time-doubling_time_std.error, ymax = doubling_time + doubling_time_std.error)) +
  geom_hline(yintercept=3, color="blue") +
  geom_text(x=50, y=3.2, label="3 Days", color="blue") +
  geom_hline(yintercept=6.1, color="blue") +
  geom_text(x=50, y=6.3, label="6.1 Days", color="blue") +
  geom_hline(yintercept=median_doubling_time, color="red")+
  geom_text(x=45, y=2.5, label="Median Doubling Time", color="red") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        axis.title.x=element_blank()) +
  ylab("Doubling Time")
ggsave("doubling_times_confirmed_states.pdf", height=5, width=7, units="in")


# Across the US
tmp <- us_confirmed_tidy %>% 
  group_by(day) %>% 
  mutate(cases = sum(cases)) %>% 
  filter(cases > 5)
fit <- glm(cases~day, family="poisson", data=tmp) %>% 
  summary %>% 
  coefficients 
# Overall US estimate
log(2)/fit[2,1]
# Overall US Std error
log(2)*(abs(log(2)/fit[2,1]))*fit[2,2]



# US Deaths ---------------------------------------------------------------


tmp <- us_confirmed_tidy %>% 
  filter(deaths > 5) %>% 
  split(.$state) 
tmp <- tmp %>% 
  keep(~nrow(.x) > 3) %>% 
  map(~glm(deaths~day, family="poisson", data=.x)) %>% 
  map(summary) %>% 
  map(coefficients) %>% 
  map(~.x[2,1:2]) %>% 
  enframe() %>%
  mutate(r = map(value, ~.x[1]), 
         std.error = map(value, ~.x[2])) %>% 
  select(-value) %>% 
  unnest() %>% 
  mutate(doubling_time = log(2)/r) %>% 
  mutate(name = factor(name, levels=name[order(doubling_time)])) %>% 
  mutate(doubling_time_std.error = log(2)*abs(doubling_time)*std.error) # https://en.wikipedia.org/wiki/Propagation_of_uncertainty
median_doubling_time <- median(tmp$doubling_time)

tmp %>% 
  ggplot(aes(x=name, y=doubling_time)) +
  geom_point() +
  geom_errorbar(aes(ymin=doubling_time-doubling_time_std.error, ymax = doubling_time + doubling_time_std.error)) +
  geom_hline(yintercept=3, color="blue") +
  geom_text(x=2.5, y=3.1, label="3 Days", color="blue") +
  geom_hline(yintercept=6.1, color="blue") +
  geom_text(x=2.5, y=5.9, label="6.1 Days", color="blue") +
  geom_hline(yintercept=median_doubling_time, color="red")+
  geom_text(x=2.5, y=2.55, label="Median Doubling Time", color="red") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        axis.title.x=element_blank()) +
  ylab("Doubling Time")
ggsave("doubling_times_deaths_states.pdf", height=5, width=7, units="in")


# Across the US
tmp <- us_confirmed_tidy %>% 
  group_by(day) %>% 
  mutate(deaths = sum(deaths)) %>% 
  filter(deaths > 5)
fit <- glm(cases~day, family="poisson", data=tmp) %>% 
  summary %>% 
  coefficients 
# Overall US estimate
log(2)/fit[2,1]
# Overall US Std error
log(2)*(abs(log(2)/fit[2,1]))*fit[2,2]



# Global Confirmed/Deaths -----------------------------------------------------------

global <- read_csv("~/Research/covid19/data/CSSEGIS/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

global


