library(tidyverse)
library(driver)
library(lubridate)
library(ggpubr)

set.seed(39502)

# import data -------------------------------------------------------------

read_nyc_syndromic <- function(fn){
  read_csv(fn) %>% 
    select(-`Select Metric`, -`Data note 1`, -`Dim2Name`, -`Ind1Name`) %>% 
    rename(count = X9, 
           borough = Dim1Value, 
           age_group = Dim2Value) %>% 
    mutate(Date = mdy(Date), 
           dofy = yday(Date), 
           age_group=factor(age_group, levels=c("Ages 0-4 years", 
                                                "Ages 5-17 years", 
                                                "Ages 18-64 years", 
                                                "Ages 65+ years", 
                                                "All age groups")))
}

vomiting <- read_nyc_syndromic("data/new_york_state/2020-04-07_nyc_vomiting.csv")
iliadmit <- read_tsv("data/new_york_state/ili_daily_counts.txt")

vomiting <- vomiting %>% 
  filter(age_group == "All age groups") %>% 
  filter(borough=="Citywide") %>% 
  rename(vomiting=count, 
         date = Date) %>% 
  select(date, vomiting)

ili_admit <- iliadmit %>% 
  rename(ili=total) %>% 
  select(date, ili, ili_admit) 

# analysis ----------------------------------------------------------------

br <- ymd("2020-01-01")

# compare to baseline -- vomiting
vomiting_baseline <- vomiting %>% 
  mutate(current = date >=br, 
         dofy = yday(date)) %>% 
  group_by(dofy, current) %>% 
  summarise(vomiting = mean(vomiting)) %>% 
  spread(current, vomiting) %>% 
  rename(baseline = `FALSE`, 
         current = `TRUE`) %>% 
  mutate(vomiting_delta = current/baseline) %>% 
  mutate(date = ymd("2019-12-31") + days(dofy)) %>% 
  filter(date <= ymd("2020-03-28")) %>% 
  ungroup() %>% 
  select(-dofy) 

iliadmit_baseline <- ili_admit %>% 
  mutate(current = date >=br, 
         dofy = yday(date)) %>% 
  mutate(admit_rate = ili_admit/ili) %>% 
  group_by(dofy, current) %>% 
  summarise(admit_rate = mean(admit_rate)) %>% 
  spread(current, admit_rate) %>% 
  rename(baseline = `FALSE`, 
         current = `TRUE`) %>% 
  mutate(admit_rate_delta = current/baseline) %>% 
  mutate(date = ymd("2019-12-31") + days(dofy)) %>% 
  filter(date <= ymd("2020-03-28")) %>% 
  ungroup() %>% 
  select(-dofy)


combined_baseline <- vomiting_baseline %>% 
  full_join(iliadmit_baseline, by=c("date"),  suffix=c(".vomiting", ".iliadmit")) %>% 
  mutate(admit_rate_delta_div_vomiting_delta = (1/admit_rate_delta)/vomiting_delta)

p1 <- combined_baseline %>% 
  select(date, current.vomiting, baseline.vomiting) %>% 
  gather(param, val, -date) %>% 
  ggplot(aes(x=date, y=val)) +
  geom_line(aes(color=param)) +
  theme_bw()+
  scale_x_date(date_labels = "%b %d", date_breaks="2 weeks") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank())

p2 <- combined_baseline %>% 
  select(date, current.iliadmit, baseline.iliadmit) %>% 
  gather(param, val, -date) %>% 
  ggplot(aes(x=date, y=val)) +
  geom_line(aes(color=param)) +
  theme_bw()+
  scale_x_date(date_labels = "%b %d", date_breaks="2 weeks") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5), 
        axis.title.y=element_blank(), 
        axis.title.x=element_blank())

p3 <- combined_baseline %>% 
  select(date, vomiting_delta, admit_rate_delta, admit_rate_delta_div_vomiting_delta) %>% 
  gather(param, val, -date) %>% 
  ggplot(aes(x=date, y=val))+
  geom_line(aes(color=param)) +
  theme_bw()+
  scale_x_date(date_labels = "%b %d", date_breaks="2 weeks") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5), 
        axis.title.y=element_blank(), 
        axis.title.x=element_blank()) 

p <- ggpubr::ggarrange(p1, p2, p3, nrow=3, align="hv")
ggsave("figures/patient_behavior_nyced.pdf", height=8, width=8, units="in")
