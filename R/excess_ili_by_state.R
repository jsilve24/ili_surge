library(tidyverse)
library(stray)
library(driver)
library(MMWRweek)
library(lubridate)
library(ggrepel)
library(padr)
library(ggpubr)

set.seed(599291)

# Here we lump NYC in with NY state as I could not find primary care provider number for the city in isolation. 
# I also define my own variable "flu_week" which starts on CDC Week 30 and goes to the end of the year -- this was just a 
# way to get peak Influenza season in the middle of the year rather than spread over two calendar years. 

# load ILINet -------------------------------------------------------------

ili <- read_csv("data/ILINet/week13/ILINet.csv", skip=1, na="X")
ili_flu_clinical <- read_csv("data/ILINet/week13/WHO_NREVSS_Clinical_Labs.csv", skip=1, na="X")
ili_flu_combined <- read_csv("data/ILINet/week13/WHO_NREVSS_Combined_prior_to_2015_16.csv", skip=1, na="X")


# Add flu year/week variables

ili <- ili %>% filter(REGION %in% c("New York City", "New York")) %>% 
  group_by(YEAR, WEEK) %>% 
  summarize(ILITOTAL = sum(ILITOTAL), 
            `NUM. OF PROVIDERS` = sum(`NUM. OF PROVIDERS`), 
            `TOTAL PATIENTS` = sum(`TOTAL PATIENTS`)) %>% 
  ungroup() %>% 
  mutate(REGION = "New York") %>% 
  bind_rows(filter(ili, !(REGION %in% c("New York City", "New York"))), .)

ili <- ili %>% 
  mutate(flu_year = ifelse(WEEK < 30, YEAR-1, YEAR)) %>%
  group_by(flu_year, REGION) %>%
  mutate(flu_week_indicator = ifelse(WEEK < 30, TRUE, FALSE), 
         flu_week = ifelse(flu_week_indicator, 0, WEEK-29), 
         flu_week = ifelse(flu_week_indicator, max(flu_week)+WEEK, flu_week)) %>% 
  ungroup() %>% 
  arrange(flu_year, flu_week, REGION) %>% 
  mutate(flu_year_week_combined = paste0(flu_year, "_", flu_week))  %>% 
  mutate(CDC_date = MMWRweek2Date(YEAR, WEEK))

flu_year_week <- ili %>% 
  pull(flu_year_week_combined) %>% 
  unique()
tmp <- seq_along(flu_year_week)
names(tmp) <- flu_year_week
flu_year_week <- tmp
rm(tmp)
CDC_date <- ili %>% 
  pull(CDC_date) %>% 
  unique()


ili <- ili %>% 
  mutate(flu_week = flu_year_week[flu_year_week_combined])

# Add flu year/week variables
ili_flu_clinical <- ili_flu_clinical %>% 
  filter(REGION != "New York City") %>% 
  mutate(flu_year = ifelse(WEEK < 30, YEAR-1, YEAR)) %>%
  group_by(flu_year, REGION) %>%
  mutate(flu_week_indicator = ifelse(WEEK < 30, TRUE, FALSE), 
         flu_week = ifelse(flu_week_indicator, 0, WEEK-29), 
         flu_week = ifelse(flu_week_indicator, max(flu_week)+WEEK, flu_week)) %>% 
  ungroup() %>% 
  arrange(flu_year, flu_week, REGION) %>% 
  mutate(flu_year_week_combined = paste0(flu_year, "_", flu_week))  %>% 
  mutate(flu_week = flu_year_week[flu_year_week_combined])

ili_flu_combined <- ili_flu_combined %>% 
  filter(REGION != "New York City") %>% 
  mutate(flu_year = ifelse(WEEK < 30, YEAR-1, YEAR)) %>%
  group_by(flu_year, REGION) %>%
  mutate(flu_week_indicator = ifelse(WEEK < 30, TRUE, FALSE), 
         flu_week = ifelse(flu_week_indicator, 0, WEEK-29), 
         flu_week = ifelse(flu_week_indicator, max(flu_week)+WEEK, flu_week)) %>% 
  ungroup() %>% 
  arrange(flu_year, flu_week, REGION) %>% 
  mutate(flu_year_week_combined = paste0(flu_year, "_", flu_week))  %>% 
  mutate(flu_week = flu_year_week[flu_year_week_combined])


# Load population data ----------------------------------------------------

pop <- read_csv("data/population/data.csv")
tmp <- pop$Pop
names(tmp) <- pop$State
pop <- tmp
rm(tmp)

# Prep data for analysis --------------------------------------------------

# Subset to just states in ILI
regions <- unique(ili$REGION)
regions_to_model <- regions[-which(regions%in% c("Virgin Islands", "Puerto Rico", 
                                                 "Commonwealth of the Northern Mariana Islands", 
                                                 "Florida"))] 

ili_select <- ili %>% 
  ungroup() %>% 
  filter(REGION %in% regions_to_model) 

ili_flu_select_clinical <- ili_flu_clinical %>% 
  ungroup() %>% 
  filter(REGION %in% regions_to_model)

ili_flu_select_combined <- ili_flu_combined %>% 
  ungroup() %>% 
  filter(REGION %in% regions_to_model)
ili_flu_select <- bind_rows(ili_flu_select_clinical, ili_flu_select_combined)

Y <- matrix(NA, length(regions_to_model), max(ili_select$flu_week))
Y_flu <- matrix(NA, length(regions_to_model), max(ili_flu_select$flu_week))
Y_flu_total <- matrix(NA, length(regions_to_model), max(ili_flu_select$flu_week))
Y_total <- matrix(NA, length(regions_to_model), max(ili_select$flu_week))
rownames(Y) <- rownames (Y_flu) <- rownames(Y_flu_total) <- rownames(Y_total) <- regions_to_model
for (i in 1:nrow(ili_select)){
  Y[ili_select$REGION[i], ili_select$flu_week[i]] <- 
    ili_select[["ILITOTAL"]][i]
}
for (i in 1:nrow(ili_select)){
  Y_total[ili_select$REGION[i], ili_select$flu_week[i]] <- 
    ili_select[["TOTAL PATIENTS"]][i]
}
for (i in 1:nrow(ili_flu_select)){
  x <- (ili_flu_select[["PERCENT POSITIVE"]][i]/100)*ili_flu_select[["TOTAL SPECIMENS"]][i]
  Y_flu[ili_flu_select$REGION[i], ili_flu_select$flu_week[i]] <- x
}
for (i in 1:nrow(ili_flu_select)){
  Y_flu_total[ili_flu_select$REGION[i], ili_flu_select$flu_week[i]] <- 
    ili_flu_select[["TOTAL SPECIMENS"]][i]
}


# Calculate Factor to Get to State Level ----------------------------------

# Assumptions/model  -- all variables refer to single week
#  let w = Total Excess ILI cases for a given state
#  let v = inferred excess ili proportion in study in a given state measured by d provided
#  let d = number of physicians (assumed primary care) in study in a given state
#  let D100K = number of primary care physicians in a given state per 100K people
#  let D = number of primary care physicians ina given state
#  let p = population of a given state
#  let n = number of patients seen by d physicians in a given state
#  
#  D = D100K*p/1e5
#  ## n/d = average number of patients seen by each primary care provider
#  w = v * (n/d) * D
#  
#  Below I calculate the factors scale_factor <- (n/d)*D for each state each week.
#  I assume the number of primary care proviers in each state over the past year 
#  is constant. 
#  
#  Note: from https://tinyurl.com/wrljc78
#  DC has 453 total primary care physicians 

D100K <- read_csv("data/primary_care_drs/primary_care_per_100K.csv")
D100K <- deframe(D100K)
missing <- setdiff(names(pop), names(D100K))
D100K <- c(D100K, deframe(data.frame(missing, NA)))
D100K <- D100K[names(pop)]
D <- D100K*pop/1e5
#D["District of Columbia"] <- 453

### Old Scale Factors based on number of providers reported by CDC in ILINet
#scale_factor <- matrix(NA, length(regions_to_model), max(ili_select$flu_week))
# scale_factor <- array(NA, dim=c(length(regions_to_model), max(ili_select$flu_week), 4000))
# rownames(scale_factor) <- regions_to_model
# for (i in 1:nrow(ili_select)){
#   tmp <- rpois(4000, ili_select[["TOTAL PATIENTS"]][i] / ili_select[["NUM. OF PROVIDERS"]][i]) * D[ili_select$REGION[i]]
#   scale_factor[ili_select$REGION[i], ili_select$flu_week[i],] <- tmp
# }
# 
# scale_factor %>%  apply(3, function(X) colSums(X, na.rm=TRUE)) %>% colMeans %>% summary

### New Scale Factors
physician_activity_rate = 0.55
scale_factor <- array(NA, dim=c(length(regions_to_model), max(ili_select$flu_week), 4000))
rownames(scale_factor) <- regions_to_model
for (i in 1:nrow(ili_select)){
  tmp <- rep(20.2, 4000)*D[ili_select$REGION[i]]*5*physician_activity_rate # Assume 5 day work weeks
  scale_factor[ili_select$REGION[i], ili_select$flu_week[i],] <- tmp
}
scale_factor %>%  apply(3, function(X) colSums(X, na.rm=TRUE)) %>% colMeans %>% summary


# 1.3 billion doctors visits in 49 weeks of the year. 
# 55% of licenced physicisans working -- alignment with number of physicians in the US
# number of patients seen by physicians per day


# impute missing flu data -------------------------------------------------


mean_impute_local <- function(Y){
  for (i in 1:nrow(Y)){
    for (j in 2:(ncol(Y)-1)){
      if (is.na(Y[i,j]) & !is.na(Y[i,j-1]) & !is.na(Y[i,j+1])){
        Y[i,j] <- mean(Y[i,j-1], Y[i,j+1])
      }
    }
  }
  return(Y)
}

mean_impute_neighbors <- function(Y, state, neighbors){
  missing <- which(is.na(Y[state,]))
  if (length(neighbors)==1){
    Y[state,missing] <- Y[neighbors,missing]
  } else {
    Y[state,missing] <- colMeans(Y[neighbors,missing,drop=F],na.rm=TRUE)
  }
  return(Y)
}

imputer <- function(Y) {
  Y %>% 
    mean_impute_local() %>% 
    mean_impute_neighbors("New Hampshire", c("Vermont", "Maine", "Massachusetts")) %>% 
    mean_impute_neighbors("New Jersey", c("Pennsylvania", "New York")) %>% 
    mean_impute_neighbors("District of Columbia", c("Maryland", "Virginia")) %>% 
    mean_impute_neighbors("Rhode Island", "Massachusetts") #%>%
  # mean_impute_neighbors("New York City", c("New York", "Pennsylvania", "Massachusetts")) %>% 
  # mean_impute_neighbors("Alaska", "Washington") %>% 
  # mean_impute_neighbors("Nevada", c("California", "Arizona", "Utah", "Oregon")) %>% 
  # mean_impute_neighbors("Wyoming", c("Idaho", "Montana", "Colorado", "South Dakota", "Nebraska")) %>% 
  # mean_impute_neighbors("Idaho", c("Montana", "Oregon", "Washington", "Wyoming", "Utah")) %>% 
  # mean_impute_neighbors("Vermont", c("New Hampshire", "New York", "New York City", "Massachusetts")) %>% 
  # mean_impute_neighbors('Maine', c("New Hampshire", "Massachusetts"))  %>% 
  # mean_impute_neighbors('Arkansas', c("Oklahoma", "Missouri", "Tennessee", "Mississippi", "Louisiana")) %>%
  # mean_impute_neighbors("Delaware", c("Maryland", "Pennsylvania", "New Jersey")) %>% 
  # mean_impute_neighbors("Mississippi", c("Alabama", "Tennessee", "Arkansas", "Louisiana")) %>% 
  # mean_impute_neighbors("Tennessee", c("North Carolina", "Kentucky", "Mississippi", "Alabama", "Georgia")) %>% 
  # mean_impute_neighbors("New Mexico", c("Colorado", "Texas", "Arizona")) %>% 
  # mean_impute_neighbors("North Carolina", c("South Carolina", "Tennessee", "Virginia", "Georgia")) %>% 
  # mean_impute_neighbors("Oklahoma", c("Texas", "Kansas", "Arkansas", "Colorado")) %>%
  # mean_impute_neighbors("Wisconsin", c("Minnesota", "Iowa", "Michigan", "Illinois"))
  
}

flu_factor <- (1-Y_flu/Y_flu_total)
flu_factor <- imputer(flu_factor)


# Brief Explore -----------------------------------------------------------

ili_flu_select %>% 
  mutate(flu_week = as.numeric(unlist(map(str_split(flu_year_week_combined,"_"), ~.x[2])))) %>% 
  ggplot(aes(x=flu_week, color=factor(flu_year), y=`PERCENT POSITIVE`))+
  geom_line() +
  facet_wrap(~REGION)+
  theme_minimal()


# Plot Non-Influenza ILI Arranged By Year (in proportions to compare easily across states)
(flu_factor*(Y/Y_total)) %>% 
  gather_array(val, REGION, date) %>% 
  mutate(REGION = regions_to_model[REGION]) %>% 
  mutate(flu_year_week = names(flu_year_week)[date]) %>% 
  mutate(flu_week = as.numeric(unlist(map(str_split(flu_year_week, "_"), ~.x[2]))), 
         flu_year = as.numeric(unlist(map(str_split(flu_year_week, "_"), ~.x[1])))) %>% 
  filter(!is.na(val)) %>% 
  mutate(Season = as.character(paste(flu_year, "-", flu_year+1))) %>% 
  ggplot(aes(x=flu_week, color=Season, y=val)) +
  geom_point(alpha=0.4,size=.5) +
  facet_wrap(~REGION, scales="free_y")+
  theme_minimal() +
  xlab("Week") +
  ylab("Proportion Patient Encounters with Non-Influenza ILI") +
  scale_x_continuous(breaks=(1:53)[seq(1, 53, by=5)], labels=c(30:53,1:29)[seq(1, 53, by=5)]) +
  theme(axis.text.x=element_text(size=5, angle=90, hjust=1), 
        axis.text.y=element_text(size=5),
        strip.text.x = element_text(size=5), 
        axis.title.x = element_text(size=7), 
        axis.title.y=element_text(size=7), 
        legend.position="bottom", 
        legend.text=element_text(size=5), 
        legend.title=element_text(size=7))
ggsave("figures/non_influenza_ili_eda.pdf", height=7, width=7, units="in")

Y_nonflu <- flu_factor*Y

# Main Analysis -----------------------------------------------------------

split_point <- which(names(flu_year_week) =="2019_1")
train_idxs <-  1:(split_point-1)
test_idxs <- (split_point):length(flu_year_week)
flu_year_week_list <- str_split(names(flu_year_week), "_")
flu_year <- as.numeric(unlist(map(flu_year_week_list, ~.x[1])))
flu_week <- as.numeric(unlist(map(flu_year_week_list, ~.x[2])))
X <- rbind(flu_week)
X_test <- rbind(seq(1,52,by=.25))
X_index_for_test <- (X_test %% 1 ==0) & (X_test <= length(test_idxs))

Gamma <- function(X) SE(X, sigma=1, rho=3)
Y_predict_summary <- list()
Y_ncov_summary <- list()
Y_ncov_full_scaled <- list() # these variable names are getting worse every day (sorry :)
Y_ncov_full_scaled_summary <- list() # these variable names are getting worse every day (sorry :)
Y_predict_full <- list()

for (r in regions_to_model){
  Y_train <- rbind(Y_nonflu[r,train_idxs,drop=F], 
                   Y_total[r,train_idxs,drop=F] - Y_nonflu[r,train_idxs,drop=F])
  Y_test <- (Y_nonflu[r,test_idxs]/Y_total[r,test_idxs])
  if (all(is.na(Y_train))) next
  if (all(is.na(Y_test))) next
  not_na_idx_test_tmp <- which(!is.na(Y_test))
  not_na_idxs_test <- which( (X_test[1,] >= min(not_na_idx_test_tmp)) &  (X_test[1,] <= max(not_na_idx_test_tmp)) )
  rm(not_na_idx_test_tmp)
  not_na_idxs_train <- which(!is.na(Y_train[1,]))
  
  Y_train <- Y_train[,not_na_idxs_train]
  rownames(Y_train) <- c("ILINotFlu", "NotILIorFlu")
  X_train <- X[,train_idxs,drop=F]
  X_train <- X_train[,not_na_idxs_train,drop=F]
  upsilon <- 1
  Theta <- function(X) return(matrix(Logit(.001), 1,ncol(X)))
  Xi <- matrix(1) 
  fit <- basset(Y_train, X_train, upsilon, Gamma=Gamma, Xi=Xi, Theta=Theta, n_samples=4000)
  X_test_tmp <- X_test[,not_na_idxs_test,drop=F]
  Y_predict <- predict(fit, X_test_tmp, summary=FALSE, response="Eta")
  Y_predict <- invLogit(Y_predict)
  X_index_for_test_comparison <- (X_test_tmp %% 1 ==0) & (X_test_tmp <= length(test_idxs))
  foo <- which(!is.na(Y_test))
  Y_test_tmp <- Y_test[min(foo):max(foo)]
  Y_ncov <- -sweep(Y_predict[,X_index_for_test_comparison,,drop=F], c(2,3), Y_test_tmp, FUN=`-`)
  Y_ncov_summary[[r]] <- Y_ncov %>% 
    gather_array(val, coord, date, iter) %>%
    mutate(date = (min(foo):max(foo))[date]) %>% 
    filter(!is.na(val)) %>% 
    mutate(positive = val > 0) %>% 
    group_by(coord, date) %>% 
    summarise(p2.5 = quantile(val, prob=0.025), 
              p25 = quantile(val, prob=0.25),
              p50 = quantile(val, prob=0.50),
              p75 = quantile(val, prob=0.75),
              p97.5 = quantile(val, prob=0.975),
              mean = mean(val), 
              p.positive = sum(positive)/n())
  Y_ncov_full_scaled[[r]] <- (Y_ncov[1,,]*scale_factor[r,test_idxs[min(foo):max(foo)],]) %>%  #sweep(Y_ncov, c(1), scale_factor[r,test_idxs[min(foo):max(foo)]], FUN=`*`) %>% 
    gather_array(val, date, iter) %>%
    mutate(date = (min(foo):max(foo))[date]) %>% 
    filter(!is.na(val))
  Y_ncov_full_scaled_summary[[r]] <- Y_ncov_full_scaled[[r]] %>% 
    mutate(positive = val>0) %>% 
    group_by(date) %>% 
    summarise(p2.5 = quantile(val, prob=0.025), 
              p25 = quantile(val, prob=0.25),
              p50 = quantile(val, prob=0.50),
              p75 = quantile(val, prob=0.75),
              p97.5 = quantile(val, prob=0.975),
              mean = mean(val), 
              p.positive = sum(positive)/n())
  
  Y_predict_full[[r]] <- Y_predict %>% 
    gather_array(val, coord, date, iter) %>% 
    mutate(date = X_test_tmp[1,date]) %>% 
    filter(!is.na(val))
  
  Y_predict_summary[[r]] <- Y_predict_full[[r]] %>% 
    group_by(coord, date) %>% 
    summarise_posterior(val) %>% 
    filter(coord==1)
}

Y_test_tidy <- (Y_nonflu[,test_idxs]/Y_total[,test_idxs]) %>% 
  gather_array(median, REGION, week) %>% 
  mutate(REGION = rownames(Y_nonflu)[REGION], 
         week = flu_week[test_idxs[week]], 
         mean=median) 

Y_train_tidy <- ili_select %>%
  mutate(week = as.numeric(unlist(map(str_split(flu_year_week_combined,"_"), ~.x[2]))),
         mean = `%UNWEIGHTED ILI`/100) 

Y_nonflu_train_tidy <- (Y_nonflu/Y_total)[,train_idxs] %>% 
  gather_array(mean, REGION, date) %>% 
  mutate(week = as.numeric(unlist(      map(str_split(names(flu_year_week[train_idxs]), "_"), ~.x[2])     ))[date]) %>%
  mutate(flu_year = as.numeric(unlist(      map(str_split(names(flu_year_week[train_idxs]), "_"), ~.x[1])     ))[date]) %>%  
  mutate(REGION = regions_to_model[REGION])


br <- which(CDC_date[test_idxs] == ymd("2020-02-02"))
Y_ncov_summary_tidy <- Y_ncov_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  filter(week >= br) %>% 
  ungroup()

# Get statistics for New York
Y_ncov_summary_tidy %>% 
  filter(REGION=="New York") %>% 
  mutate(date = CDC_date[test_idxs][week]) %>% 
  filter(date==ymd("2020-03-22"))

Y_ncov_summary_tidy  %>% 
  select(REGION, week, p2.5:p97.5, mean, p.positive) %>% 
  mutate(week = CDC_date[test_idxs][week]) %>% 
  write_csv("results/ncov_signal_extraction.csv")


Y_ncov_full_scaled_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  filter(week >= br) %>% 
  ungroup() %>% 
  mutate(week = CDC_date[test_idxs][week]) %>% 
  write_csv("results/number_excess_ili_cases.csv")
rm(foo, tmp)  


last_week <- MMWRweek(CDC_date[test_idxs][length(test_idxs)])
ntest <- length(test_idxs)
breaks <- 1:53
labels <- c(CDC_date[test_idxs],
            MMWRweek2Date(rep(2020, 53-ntest),
                          (last_week$MMWRweek+1):(last_week$MMWRweek+53-ntest))) 
breaks <- breaks[seq(1,53,by=5)]
labels <- labels[seq(1,53,by=5)]
labels <- format(labels, format="%b %d")
Y_predict_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  mutate(week = date) %>%
  filter(!(REGION %in% c("Maine", "Vermont") & week < 10)) %>% 
  ggplot(aes(x=week, y=mean)) +
  #geom_ribbon(aes(ymin=p2.5, ymax=p97.5), color="grey", alpha=0.3) +
  #geom_ribbon(aes(ymin=p25, ymax=p75), color="grey", alpha=0.5) +
  geom_path(data = Y_nonflu_train_tidy, aes(group=flu_year), alpha=0.2, lwd=.5) +
  #geom_line(lwd=.5)+
  geom_line(data=Y_test_tidy, color="red", lwd=.5) +
  geom_ribbon(data=Y_ncov_summary_tidy, aes(ymin=p25, ymax=p75), fill="blue", alpha=0.4) +
  geom_ribbon(data=Y_ncov_summary_tidy, aes(ymin=p2.5, ymax=p97.5), fill="blue", alpha=0.4) +
  geom_line(data=Y_ncov_summary_tidy, color="blue", alpha=.7) +
  #geom_line(data=Y_flu_test_tidy, color="green") +
  facet_wrap(~REGION, scales="free_y", ncol=5) +
  theme_bw() +
  ylab("Non-Influenza ILI Proportion") + 
  scale_x_continuous(breaks=breaks, labels=labels) +
  # theme(axis.text.x=element_text(angle=90, hjust=1), 
  #       axis.title.x = element_blank()) +
  theme(axis.text.x=element_text(size=5, angle=90, hjust=1), 
        axis.text.y=element_text(size=5),
        strip.text.x = element_text(size=5), 
        axis.title.x = element_blank(), 
        axis.title.y=element_text(size=7), 
        legend.position="bottom", 
        legend.text=element_text(size=5), 
        legend.title=element_text(size=7), 
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(0, NA))
  #ylim(c(0, NA))
ggsave("figures/ncov_signal_extraction.pdf", height=9, width=7, units="in")


# Make Washington Specific Graphic
Y_predict_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  filter(REGION == "Washington") %>% 
  mutate(week = date) %>%
  ggplot(aes(x=week, y=mean)) +
  geom_line(data =  filter(Y_nonflu_train_tidy, REGION == "Washington"), aes(group=flu_year), alpha=0.3) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), color="grey", alpha=0.3) +
  geom_ribbon(aes(ymin=p25, ymax=p75), color="grey", alpha=0.5) +
  geom_line()+
  geom_line(data= filter(Y_test_tidy, REGION == "Washington"), color="red") +
  geom_linerange(data= filter(Y_ncov_summary_tidy, REGION == "Washington"), aes(ymin=p25, ymax=p75), color="blue", alpha=0.5) +
  geom_point(data= filter(Y_ncov_summary_tidy, REGION == "Washington"), color="blue", size=2) +
  #geom_line(data=Y_flu_test_tidy, color="green") +
  theme_bw() +
  ylab("Non-Flu ILI Proportion") + 
  scale_x_continuous(breaks=breaks, labels=labels) +
  theme(axis.text.x=element_text(angle=90, hjust=1), 
        axis.title.x = element_blank())
ggsave("figures/ncov_signal_extraction_washington.pdf", height=5, width=7, units="in")


Y_predict_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  filter(REGION == "New York") %>% 
  mutate(week = date) %>%
  ggplot(aes(x=week, y=mean)) +
  geom_line(data =  filter(Y_nonflu_train_tidy, REGION == "New York"), aes(group=flu_year), alpha=0.3) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), color="grey", alpha=0.3) +
  geom_ribbon(aes(ymin=p25, ymax=p75), color="grey", alpha=0.5) +
  geom_line()+
  geom_line(data= filter(Y_test_tidy, REGION == "New York"), color="red") +
  geom_linerange(data= filter(Y_ncov_summary_tidy, REGION == "New York"), aes(ymin=p25, ymax=p75), color="blue", alpha=0.5) +
  geom_point(data= filter(Y_ncov_summary_tidy, REGION == "New York"), color="blue", size=2) +
  #geom_line(data=Y_flu_test_tidy, color="green") +
  theme_bw() +
  ylab("Non-Flu ILI Proportion") + 
  scale_x_continuous(breaks=breaks, labels=labels) +
  theme(axis.text.x=element_text(angle=90, hjust=1), 
        axis.title.x = element_blank()) +
  ggtitle("New York State")
ggsave("figures/ncov_signal_extraction_new_york.png", height=5, width=7, units="in")


tmp <- Y_ncov_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  ungroup()

focus_fn <- function(x, filter_date=TRUE, inv_filter_date=FALSE) {
  br <- which(CDC_date[test_idxs] == ymd("2020-03-08"))
  focus <- c("New York", "Louisiana", "Washington", "Iowa")
  x <- x %>% 
    filter(REGION %in% focus) %>% 
    mutate(REGION = factor(REGION, levels = focus)) 
  if (filter_date) return(filter(x, week >= br))
  if (inv_filter_date) return(filter(x, week <= br))
  return(x)
}
Y_predict_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  mutate(week = date) %>%
  focus_fn(filter_date=FALSE) %>% 
  ggplot(aes(x=week, y=mean)) +
  geom_line(data =  focus_fn(Y_nonflu_train_tidy, FALSE), aes(group=flu_year), alpha=0.3) +
  #geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
  #geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.7) +
  geom_line(alpha=0.8, color="black")+
  geom_line(data= focus_fn(Y_test_tidy, FALSE), color="red") +
  geom_ribbon(data=focus_fn(tmp, FALSE, TRUE), aes(ymin=p25, ymax=p75), fill="blue", alpha=0.2) +
  geom_ribbon(data=focus_fn(tmp, FALSE, TRUE), aes(ymin=p2.5, ymax=p97.5), fill="blue", alpha=0.2) +
  geom_line(data=focus_fn(tmp, FALSE, TRUE), color="blue", alpha=.3) +
  geom_ribbon(data=focus_fn(tmp), aes(ymin=p25, ymax=p75), fill="blue", alpha=0.6) +
  geom_ribbon(data=focus_fn(tmp), aes(ymin=p2.5, ymax=p97.5), fill="blue", alpha=0.6) +
  geom_line(data=focus_fn(tmp), color="blue", alpha=.8) +
  #geom_line(data=Y_flu_test_tidy, color="green") +
  facet_wrap(~REGION) +
  theme_bw() +
  ylab("Non-Influenza ILI Proportion") + 
  scale_x_continuous(breaks=breaks, labels=labels) +
  theme(axis.text.x=element_text(size=7, angle=90, hjust=1), 
        axis.text.y=element_text(size=7),
        strip.text.x = element_text(size=9), 
        axis.title.x = element_blank(), 
        axis.title.y=element_text(size=9), 
        panel.grid.minor = element_blank())  +
  coord_cartesian(ylim=c(0, NA))
ggsave("figures/ncov_signal_extraction_4_states.pdf", height=4, width=7.3, units="in")


# Compare to results without extracting flu -------------------------------

Y_ncov_full_scaled_wflu_summary <- list() # these variable names are getting worse every day (sorry :)

for (r in regions_to_model){
  Y_train <- rbind(Y[r,train_idxs,drop=F], 
                   Y_total[r,train_idxs,drop=F] - Y[r,train_idxs,drop=F])
  Y_test <- (Y[r,test_idxs]/Y_total[r,test_idxs])
  if (all(is.na(Y_train))) next
  if (all(is.na(Y_test))) next
  not_na_idx_test_tmp <- which(!is.na(Y_test))
  not_na_idxs_test <- which( (X_test[1,] >= min(not_na_idx_test_tmp)) &  (X_test[1,] <= max(not_na_idx_test_tmp)) )
  rm(not_na_idx_test_tmp)
  not_na_idxs_train <- which(!is.na(Y_train[1,]))
  
  Y_train <- Y_train[,not_na_idxs_train]
  rownames(Y_train) <- c("ILI", "NotILI")
  X_train <- X[,train_idxs,drop=F]
  X_train <- X_train[,not_na_idxs_train,drop=F]
  upsilon <- 1
  Theta <- function(X) return(matrix(Logit(.001), 1,ncol(X)))
  Xi <- matrix(1) 
  fit <- basset(Y_train, X_train, upsilon, Gamma=Gamma, Xi=Xi, Theta=Theta, n_samples=4000)
  X_test_tmp <- X_test[,not_na_idxs_test,drop=F]
  Y_predict <- predict(fit, X_test_tmp, summary=FALSE, response="Eta")
  Y_predict <- invLogit(Y_predict)
  X_index_for_test_comparison <- (X_test_tmp %% 1 ==0) & (X_test_tmp <= length(test_idxs))
  foo <- which(!is.na(Y_test))
  Y_test_tmp <- Y_test[min(foo):max(foo)]
  Y_ncov <- -sweep(Y_predict[,X_index_for_test_comparison,,drop=F], c(2,3), Y_test_tmp, FUN=`-`)
  tmp <- (Y_ncov[1,,]*scale_factor[r,test_idxs[min(foo):max(foo)],]) %>% 
    gather_array(val, date, iter) %>%
    mutate(date = (min(foo):max(foo))[date]) %>% 
    filter(!is.na(val))
  Y_ncov_full_scaled_wflu_summary[[r]] <- tmp %>% 
    mutate(positive = val>0) %>% 
    group_by(date) %>% 
    summarise(p2.5 = quantile(val, prob=0.025), 
              p25 = quantile(val, prob=0.25),
              p50 = quantile(val, prob=0.50),
              p75 = quantile(val, prob=0.75),
              p97.5 = quantile(val, prob=0.975),
              mean = mean(val), 
              p.positive = sum(positive)/n())
}

br <- which(CDC_date[test_idxs] == ymd("2020-03-08"))

noflu <- Y_ncov_full_scaled_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  filter(week >= br) %>% 
  ungroup() %>% 
  mutate(week = CDC_date[test_idxs][week]) %>% 
  select(REGION, p2.5:p97.5, week) %>% 
  mutate(p2.5 = pmax(0, p2.5))

flu <- Y_ncov_full_scaled_wflu_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  filter(week >= br) %>% 
  ungroup() %>% 
  mutate(week = CDC_date[test_idxs][week])  %>% 
  mutate(p2.5 = pmax(0, p2.5))

full_join(noflu, flu, by=c("REGION", "week"), suffix=c(".noflu", ".flu")) %>% 
  ggplot(aes(x=p50.noflu, y=p50.flu)) +
  geom_segment(x=0, y=0, xend=3e+5, yend=3e+5, color="red") +
  geom_point() +
  geom_linerange(aes(xmin=p2.5.noflu, xmax=p97.5.noflu), alpha=0.7) +
  geom_linerange(aes(ymin=p2.5.flu, ymax=p97.5.flu), alpha=0.7) +
  #coord_cartesian(xlim=c(0, NA), ylim=c(0, NA)) +
  theme_bw() +
  ylab("Non-Seasonal ILI") +
  xlab("Non-Seasonal Non-Influenza ILI")
ggsave("figures/sensitivity_withflu.pdf", height=5, width=5, units="in")


# Compare to results with "simple model" -------------------------------

Y_ncov_full_scaled_simple_summary <- list() # these variable names are getting worse every day (sorry :)

for (r in regions_to_model){
  Y_train <- rbind(Y[r,train_idxs,drop=F], 
                   Y_total[r,train_idxs,drop=F] - Y[r,train_idxs,drop=F])
  Y_test <- (Y[r,test_idxs]/Y_total[r,test_idxs])
  if (all(is.na(Y_train))) next
  if (all(is.na(Y_test))) next
  not_na_idx_test_tmp <- which(!is.na(Y_test))
  not_na_idxs_test <- which( (X_test[1,] >= min(not_na_idx_test_tmp)) &  (X_test[1,] <= max(not_na_idx_test_tmp)) )
  rm(not_na_idx_test_tmp)
  not_na_idxs_train <- which(!is.na(Y_train[1,]))
  
  Y_train <- Y_train[,not_na_idxs_train]
  rownames(Y_train) <- c("ILI", "NotILI")
  X_train <- X[,train_idxs,drop=F]
  X_train <- X_train[,not_na_idxs_train,drop=F]
  upsilon <- 1
  Theta <- function(X) return(matrix(Logit(.001), 1,ncol(X)))
  Xi <- matrix(1) 
  fit <- basset(Y_train, X_train, upsilon, Gamma=Gamma, Xi=Xi, Theta=Theta, n_samples=4000)
  X_test_tmp <- X_test[,not_na_idxs_test,drop=F]
  Y_predict <- predict(fit, X_test_tmp, summary=FALSE, response="Eta")
  Y_predict <- invLogit(Y_predict)
  X_index_for_test_comparison <- (X_test_tmp %% 1 ==0) & (X_test_tmp <= length(test_idxs))
  foo <- which(!is.na(Y_test))
  Y_test_tmp <- Y_test[min(foo):max(foo)]
  Y_ncov <- -sweep(Y_predict[,X_index_for_test_comparison,,drop=F], c(2,3), Y_test_tmp, FUN=`-`)
  tmp <- (Y_ncov[1,,]*scale_factor[r,test_idxs[min(foo):max(foo)],]) %>% 
    gather_array(val, date, iter) %>%
    mutate(date = (min(foo):max(foo))[date]) %>% 
    filter(!is.na(val))
  Y_ncov_full_scaled_simple_summary[[r]] <- tmp %>% 
    mutate(positive = val>0) %>% 
    group_by(date) %>% 
    summarise(p2.5 = quantile(val, prob=0.025), 
              p25 = quantile(val, prob=0.25),
              p50 = quantile(val, prob=0.50),
              p75 = quantile(val, prob=0.75),
              p97.5 = quantile(val, prob=0.975),
              mean = mean(val), 
              p.positive = sum(positive)/n())
}

br <- which(CDC_date[test_idxs] == ymd("2020-03-08"))

noflu <- Y_ncov_full_scaled_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  filter(week >= br) %>% 
  ungroup() %>% 
  mutate(week = CDC_date[test_idxs][week]) %>% 
  select(REGION, p2.5:p97.5, week) %>% 
  mutate(p2.5 = pmax(0, p2.5))

flu <- Y_ncov_full_scaled_wflu_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  filter(week >= br) %>% 
  ungroup() %>% 
  mutate(week = CDC_date[test_idxs][week])  %>% 
  mutate(p2.5 = pmax(0, p2.5))

full_join(noflu, flu, by=c("REGION", "week"), suffix=c(".noflu", ".flu")) %>% 
  ggplot(aes(x=p50.noflu, y=p50.flu)) +
  geom_segment(x=0, y=0, xend=3e+5, yend=3e+5, color="red") +
  geom_point() +
  geom_linerange(aes(xmin=p2.5.noflu, xmax=p97.5.noflu), alpha=0.7) +
  geom_linerange(aes(ymin=p2.5.flu, ymax=p97.5.flu), alpha=0.7) +
  #coord_cartesian(xlim=c(0, NA), ylim=c(0, NA)) +
  theme_bw() +
  ylab("Non-Seasonal ILI") +
  xlab("Non-Seasonal Non-Influenza ILI")
ggsave("figures/sensitivity_withflu.pdf", height=5, width=5, units="in")



# Plot how this excess ILI compares to confirmed cases  -------------------

# Load state abbreviations
abbr <- read_csv("data/population/state_abbreviations.csv", col_names=FALSE)
abbr <- select(abbr, X1, X4)
tmp <- abbr$X1
names(tmp) <- abbr$X4
abbr <- tmp

# New York Times Dataset
us_confirmed_tidy <- read_csv("data/covid-19-data/us-states.csv") %>% 
  mutate(State=state) %>% 
  filter(state %in% regions_to_model)

us_confirmed_tidy %>% 
  filter(state=="New York") %>% 
  filter(date == ymd("2020-03-28")) %>% 
  pull(cases) -> foo 
  (foo)/19450000*100



us_new_confirmed_tidy <- us_confirmed_tidy %>%
  pad(group="State") %>%
  group_by(State) %>% 
  arrange(date) %>% 
  mutate(NewCasesDay = cases - lag(cases)) %>% 
  ungroup() %>% 
  mutate(CDC_week = MMWRweek(date)$MMWRweek) %>% 
  group_by(State, CDC_week) %>% 
  summarize(NewCasesWeek = sum(NewCasesDay, na.rm=TRUE)) %>% 
  mutate(date = MMWRweek2Date(rep(2020, length(CDC_week)), CDC_week)) %>% 
  mutate(REGION=State)

tmp <- Y_ncov_full_scaled_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  #filter(p.positive > .95) %>% 
  #mutate(week = X_test[1,X_index_for_test][date]) %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  # filter(last_consecutive_bool(week, gap_allowance = 1)) %>% 
  # filter(contains_val_bool(week, c(33:50))) %>% 
  # filter(upswing_bool(mean)) %>% 
  filter(week >= br) %>% 
  ungroup() %>% 
  mutate(date = CDC_date[test_idxs][week]) %>% 
  left_join(us_new_confirmed_tidy, by=c("REGION", "date"))



# Plot correlation with covid counts
tmp %>% 
  ggplot(aes(y=p50, x=NewCasesWeek)) +
  geom_linerange(aes(ymin=p2.5, ymax=p97.5),  size=2, alpha=0.5, color="grey") +
  geom_linerange(aes(ymin=p25, ymax=p75),  size=2, alpha=0.8, color="grey") +
  geom_point() +
  geom_smooth(method="lm") +
  scale_x_log10() +
  scale_y_log10()


p1 <- Y_ncov_full_scaled_summary %>% 
  bind_rows(.,.id="REGION") %>% 
  ungroup() %>% 
  #filter(p.positive > .95) %>% 
  #mutate(week = X_test[1,X_index_for_test][date]) %>% 
  mutate(week = date) %>% 
  group_by(REGION) %>% 
  arrange(week) %>% 
  #filter(last_consecutive_bool(week, gap_allowance = 1)) %>% 
  #filter(contains_val_bool(week, c(33:50))) %>% 
  #filter(upswing_bool(mean)) %>% 
  ungroup() %>% 
  mutate(date = CDC_date[test_idxs][week]) %>%
  filter(date >= ymd("2020-03-1")) %>% 
  left_join(us_new_confirmed_tidy, by=c("REGION", "date")) %>% 
  mutate(Population = pop[REGION]) %>%
  ggplot(aes(y=p50/Population, x=NewCasesWeek/Population)) +
  # geom_linerange(aes(ymin=p2.5/population, ymax=p97.5/population),  size=2, alpha=0.5, color="grey") +
  # geom_linerange(aes(ymin=p25/population, ymax=p75/population),  size=2, alpha=0.8, color="grey") +
  geom_point(color="darkgrey", alpha=0.6) +
  stat_cor(show.legend = FALSE, size=3) +
  geom_smooth(method="lm", show.legend = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~date) +
  theme_bw() +
  guides(label=FALSE) +
  xlab("New COVID Cases (Per Capita)") +
  ylab("Non-Influenza ILI in Excess\nof Seasonal Norms (Per Capita)")
ggsave("figures/correlation_covid_excess_ili_by_week.pdf", plot=p1, height=5, width=10)
#ggsave("correlation_covid_excess_ili_by_week.png", height=5, width=10)

# Bring in ILI Admission Data from NY EDs
ili_admit <- read_tsv("data/new_york_state/ili_daily_counts.txt")
baz <- ili_admit %>% 
  rename(`Total ILI` = total, 
         `ILI Admissions` = ili_admit)
foo <- select(baz, `date`, admit_rate) %>% rename(`ILI Admission Percent` = admit_rate)
p2 <- baz %>% 
  select(-admit_rate) %>% 
  gather(param, count, -date) %>% 
  mutate(param = factor(param, levels = c("Total ILI", "ILI Admissions"))) %>% 
  ggplot(aes(x=date, y=count)) +
  geom_bar(aes(fill=param), stat="identity") +
  geom_line(data = foo, aes(y=`ILI Admission Percent`*2000/100),color="black", size=1.5) +
  theme_bw(  ) +
  ylab("Visits") +
  scale_y_continuous(limits = c(0, 2000), sec.axis = sec_axis(~.*100/2000, name = "ILI Admissions / Total ILI")) +
  scale_x_date(date_breaks="2 days", limits = c(ymd('2020-02-15'), NA), date_labels = "%b %d") +
  theme(axis.title.x=element_blank(), 
        axis.text.x = element_text(angle=90, hjust=1), 
        legend.title = element_blank(), 
        legend.position = c(.2, .8), 
        legend.text = element_text(size=10), 
        legend.key.size = unit(.7, "cm"), 
        legend.background=element_blank()) +
  scale_fill_brewer(palette="Set1")
ggsave("figures/ili_admission_rate.pdf", plot=p2, height=4, width=7)
rm(baz)

p <- ggarrange(p1, p2, nrow=2, labels = c("A", "B"), heights=c(1, .7))
ggsave("figures/ili_admissions_and_correlation.pdf", plot=p,  height=6, width=6, units="in")

foo <- tmp %>% 
  filter(!is.na(NewCasesWeek)) 
cor.test(foo$p50, foo$NewCasesWeek)
var(log(foo$p50)/log(foo$NewCasesWeek+0.1))
rm(foo)

tmp %>% 
  mutate(p2.5 = NewCasesWeek/p2.5, 
         p25 = NewCasesWeek/p25, 
         p50 = NewCasesWeek/p50, 
         p75 = NewCasesWeek/p75, 
         p97.5 = NewCasesWeek/p97.5) %>% 
  select(-mean) %>% 
  arrange(REGION, date) %>% 
  filter(p50 !=0) %>% 
  write_csv("results/new_confirmed_divided_by_excess_ili.csv")

tmp %>% 
  mutate(p2.5 = NewCasesWeek/p2.5, 
         p25 = NewCasesWeek/p25, 
         p50 = NewCasesWeek/p50, 
         p75 = NewCasesWeek/p75, 
         p97.5 = NewCasesWeek/p97.5) %>% 
  select(-mean) %>% 
  arrange(REGION, date) %>% 
  filter(p50 !=0) %>% 
  pull(p50) %>% 
  mean()

tmp %>% 
  filter(date >= ymd("2020-03-01")) %>% 
  mutate(State = REGION) %>% 
  ggplot(aes(y=NewCasesWeek/p50, x=date, color=State, fill=State)) +
  #geom_ribbon(aes(ymin=NewCasesWeek/p2.5, ymax=NewCasesWeek/p97.5), linetype=0, alpha=0.1)+
  geom_ribbon(aes(ymin=NewCasesWeek/p25, ymax=NewCasesWeek/p75), linetype=0, alpha=0.3)+
  geom_path() +
  geom_point() +
  #ylab("Probability of SARS-CoV2 detection given patient has ILI due to SARS-CoV2") +
  ylab("Syndromic Case Detection Rate") +
  theme_bw() +
  scale_x_date(date_breaks="week", limits = c(ymd("2020-03-7"), max(CDC_date[test_idxs]))) +
  theme(axis.text.x = element_text(angle=90, hjust=1), 
        axis.title.x = element_blank(), 
        legend.position = "bottom",
        axis.title.y=element_text(size=8)) +
  scale_y_log10(labels = function(x) sprintf("%g", x))
ggsave("figures/excess_ili_vs_case_counts.pdf", height=7, width=7, units="in")



tmp %>% 
  mutate(State = REGION) %>% 
  ggplot(aes(y=p50, x=date))+
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), color="grey", alpha=0.5) +
  geom_ribbon(aes(ymin=p25, ymax=p75), color="grey", alpha=0.8) +
  geom_line()+
  geom_point() +
  facet_wrap(~State, scales="free_y") +
  theme_bw() +
  scale_x_date(date_breaks="3 weeks", limits = c(min(CDC_date[test_idxs]), max(CDC_date[test_idxs]))) +
  theme(axis.text.x = element_text(angle=90, hjust=1), 
        axis.title.x = element_blank()) +
  ylab("Non-Flu ILI Cases In Excess of Prior Seasonal Variation")
ggsave("figures/excess_ili.pdf", height=15, width=20, units="in")


tmp <- tmp %>% 
  mutate(State = REGION) %>% 
  filter(date == max(date)) 
so <- tmp %>% arrange(NewCasesWeek/mean) %>% pull(State)

tmp %>% 
  mutate(State = factor(State, levels=so)) %>% 
  ggplot(aes(y=NewCasesWeek/mean, x=State)) +
  geom_linerange(aes(ymin=NewCasesWeek/p2.5, ymax=NewCasesWeek/p97.5), size=2, alpha=0.5, color="grey")+
  geom_linerange(aes(ymin=NewCasesWeek/p25, ymax=NewCasesWeek/p75), size=2, alpha=1, color="darkgrey") +
  geom_point() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1), 
        axis.title.x=element_blank()) +
  ylab("Syndromic Case Detection Rate") +
  scale_y_log10(labels = function(x) sprintf("%g", x))
ggsave("figures/excess_ili_vs_case_counts_latest.pdf", height=7, width=7, units="in")


# Try accounting for changing patient behavior ----------------------------

scaling_admissions <-  ili_admit %>% 
  select(date, admit_rate) %>% 
  filter(date >= ymd("2020-03-01")) %>% 
  mutate(week = floor(as.integer(date - ymd("2020-03-08"))/7) + 1) %>% 
  group_by(week) %>% 
  mutate(mean_admit_rate = exp(mean(log(admit_rate)))) %>% 
  filter(row_number()==1) %>% 
  filter(date <= max(CDC_date)) %>%
  ungroup() %>% 
  select(date, mean_admit_rate) %>% 
  deframe(  )
scaling_admissions <- scaling_admissions[-1]/scaling_admissions[1]

write_csv(enframe(scaling_admissions, "date", "scaling_factor"), path="results/scaling_admission_rates.csv")

# Compare with Alex's Simlulations ----------------------------------------

# Quantile match a beta distribution to mizumoto estimates
quantile_match_beta <- function(pars){ 
  pars <- exp(pars)
  qs <- qbeta(c(.025, .975), pars[1], pars[2]) 
  abs(qs[1] - .155) + abs(qs[2] - .202) + abs(pars[1]/(pars[1]+pars[2])-.179)
}
res <- optim(c(1,1), quantile_match_beta)


# calculate subclinical factor with uncertainty bounds
delta_b <- rbeta(100000, exp(res$par[1]), exp(res$par[2]))
delta_c <- runif(100000, 0.39, 0.41)
scrate <- (1-0.6)*(1-delta_b)
quantile(scrate, probs=c(0.025, .5, 0.975))


tmp <- Y_ncov_full_scaled %>%
  map(~mutate(.x, date = CDC_date[test_idxs][date])) %>%
  bind_rows(.id="State") %>%
  filter(date >= ymd("2020-03-08")) %>%
  group_by(date, iter) %>%
  summarize(us_val = sum(val),
            delta_b = rbeta(1, exp(res$par[1]), exp(res$par[2])), 
            delta_c = .60) %>% # 60% don't show up to doctor
  ungroup()


# Number excess per state per week 
Y_ncov_full_scaled %>%
  map(~mutate(.x, date = CDC_date[test_idxs][date])) %>%
  bind_rows(.id="State") %>%
  filter(date >= ymd("2020-03-08")) %>% 
  group_by(State, date) %>% 
  summarize(mean=mean(val)) %>% 
  ungroup() %>% 
  mutate(pop = pop[State]) %>% 
  mutate(greater = pop < mean)

# Prevalance estimates by state
foo <- Y_ncov_full_scaled %>%
  map(~mutate(.x, date = CDC_date[test_idxs][date])) %>%
  bind_rows(.id="State") %>%
  filter(date >= ymd("2020-03-08")) %>%
  group_by(State, iter) %>%
  summarize(val = sum(val),
            delta_b = rbeta(1, exp(res$par[1]), exp(res$par[2])), 
            delta_c = .60) %>% # 60% don't show up to doctor
  ungroup() %>% 
  mutate(val = val/((1-delta_b)*(1-delta_c) )) %>% 
  mutate(pop = pop[State]) %>% 
  mutate(val = val/pop*100) %>% 
  group_by(State) %>% 
  summarise_posterior(val) 
foo <- us_new_confirmed_tidy %>% 
  filter(date >=ymd("2020-03-08")) %>% 
  group_by(State) %>% 
  summarise(cases = sum(NewCasesWeek)) %>% 
  mutate(pop = pop[State]) %>% 
  mutate(cases = cases/pop*100) %>% 
  right_join(foo, by="State")
foo$State <- factor(foo$State, levels=foo$State[order(foo$p50, decreasing = FALSE)])
p <- foo %>% 
  mutate(p2.5 = pmax(p2.5, cases)) %>% 
  select(State, p2.5, p50, p97.5, cases) %>% 
  ggplot(aes(x=p50, y=State))+
  geom_point() +
  geom_linerange(aes(xmin=p2.5, xmax=p97.5)) +
  geom_point(aes(x=cases), color="darkgrey", alpha=0.6) +
  #coord_cartesian(xlim=c(0, NA)) +
  theme_bw() +
  scale_x_sqrt() +
  xlab("Percentage of State Population Infected\nBetween March 8 to March 28, 2020") +
  theme(axis.title.y=element_blank()) 
ggsave("figures/PrevalanceEstimatesStates.pdf", plot=p, height=8, width=5, units="in")
  
p + ggtitle("ILI based SARS-CoV-2 Prevalence Estimates", "estimates from confirmed case counts in grey")
ggsave("figures/PrevalanceEstimatesStates_twitter.png", height=8, width=5, units="in")

foo %>% 
  filter(State=="New York")


# Number of excess ILI cases
tmp %>% group_by(date) %>% summarise_posterior(us_val) %>% 
  write_csv("results/US_total_weekly_excess_ili_no_subclinical.csv")

# Number of new infections, with subclinical
tmp %>% 
  mutate(us_val = us_val/((1-delta_b)*(1-delta_c) )) %>% 
  group_by(date) %>% 
  summarise_posterior(us_val) %>% 
  write_csv("results/US_total_weekly_excess_ili_with_subclinical.csv")

# Number of new infections, with subclinical
tmp %>% 
  mutate(us_val = us_val/((1-delta_b)*(1-delta_c) )) %>% 
  left_join(mutate(enframe(scaling_admissions, "date", "sf"), date=ymd(date)), by="date") %>% 
  mutate(us_val = us_val * sf) %>% 
  group_by(date) %>% 
  summarise_posterior(us_val) %>% 
  write_csv("results/US_total_weekly_excess_ili_with_subclinical_with_behavior.csv")

# Excess ILI over those three weeks
tmp %>% 
  group_by(iter) %>% 
  summarise(us_val = sum(us_val)) %>% 
  ungroup() %>% 
  summarise_posterior(us_val)

# Number of newly infected
tmp %>% mutate(us_val = us_val/((1-delta_b)*(1-delta_c))) %>% group_by(date) %>% summarise_posterior(us_val) %>% 
  pull(p50) %>% sum()

# Total number infected 
tmp %>% mutate(us_val = us_val/((1-delta_b)*(1-delta_c))) %>% 
  group_by(iter) %>% 
  summarise(us_val = sum(us_val)) %>% 
  summarise_posterior(us_val)

total_new_cases <- sum(filter(us_new_confirmed_tidy, CDC_week==14)[["NewCasesWeek"]])
tmp %>% 
  mutate(detection_rate = total_new_cases/us_val) %>% 
  summarise_posterior(detection_rate)






# working with alex's simulations -----------------------------------------

sims <- read_csv("results/weekly_I_across_replicates_4d_lag.csv") %>% 
  select(-X1) %>% 
  filter(date >= ymd("2020-03-08"))



tmp <- Y_ncov_full_scaled %>%
  map(~mutate(.x, date = CDC_date[test_idxs][date])) %>%
  bind_rows(.id="State") %>%
  filter(date >= ymd("2020-03-08")) %>%
  group_by(date, iter) %>%
  summarize(us_val = sum(val),
            delta_b = rbeta(1, exp(res$par[1]), exp(res$par[2]))) %>%
  ungroup()

tmp <- tmp %>% 
  mutate(us_val = us_val/(1-delta_b))
ili_sims <- array(0, c(length(unique(tmp$date)), 4000))
rownames(ili_sims) <- as.character(unique(tmp$date))
for (i in 1:nrow(tmp)){
  ili_sims[as.character(tmp$date[i]), tmp$iter[i]] <- tmp$us_val[i]
}
alex_sims <- array(0, c(length(unique(tmp$date)), max(sims$replicate)))
sims <- filter(sims, date %in% unique(tmp$date))
rownames(alex_sims) <- as.character(unique(tmp$date))
for (i in 1:max(sims$replicate)){
  alex_sims[as.character(sims$date[i]), sims$replicate[i]] <- sims$weekly_I[i]
}

d1 <- approxfun(density(rlnorm(4000, log(ili_sims[1,]), .2)))
#d2 <- approxfun(density(rlnorm(4000, log(ili_sims[2,]), .2)))
#d3 <- approxfun(density(rlnorm(4000, log(ili_sims[3,]), .2)))

sim_probs <- rep(0, max(sims$replicate))
for (i in 1:ncol(alex_sims)){
  sim_probs[i] <- d1(alex_sims[1,i])#*d2(alex_sims[2,i])*d3(alex_sims[3,i])
}
any(!is.na(sim_probs))
sim_probs[is.na(sim_probs)] <- 0

non_zero_sims <- which(sim_probs>0)
non_zero_prob <- miniclo(sim_probs)
data.frame("replicate" = 1:length(sim_probs), "probability" =non_zero_prob[1,]) %>% 
  write_csv("results/distribution_over_replicates.csv")


d <- rownames(ili_sims)
ili_sims %>% 
  gather_array(weekly_I, date, iter) %>% 
  mutate(date = d[date]) %>% 
  write_csv("results/posterior_samples_us_weekly_I.csv")

