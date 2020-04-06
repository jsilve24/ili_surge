rm(list=ls())
gc()
library(data.table)
library(magrittr)
library(ggplot2)
library(viridis)
library(ggpubr)


# Deaths ------------------------------------------------------------------

X <- read.csv('data/covid-19-data/covid-19-data/us-states.csv',stringsAsFactors = F) %>% as.data.table

X[,deaths:=c(deaths[1],diff(deaths)),by=state]


### for each state, we'll analyze growth rate in deaths starting on the first date with non-zero deaths.
X[,date:=as.Date(date)]
setkey(X,state,date)
X[,cumulative_deaths:=cumsum(deaths),by=state]
Y <- X[,list(deaths=sum(deaths)),by=date]
X[,total_deaths:=sum(deaths),by=state]
X <- X[total_deaths>20]


## for the nation, we'll analyze growth rates with deaths after March 4 - the period of regular exp growth

ggplot(Y,aes(date,deaths))+
  geom_point()+
  geom_line()+
  scale_y_continuous(trans='log')

setkey(Y,date)
Y <- Y[date>as.Date('2020-03-04') & date<as.Date('2020-04-01')]
fit_national <- glm(deaths~date,family=poisson,data=Y)
gr_national <- coef(fit_national)[2]
se_national <- sqrt(vcov(fit_national)[4])



# SEIR models vs. ILI -------------------------------------------------------------

US <- readRDS('results/US_seir_forecasts_unif_weekly.Rd')
US[,replicate:=factor(replicate)]
setkey(US,replicate,date)
ILI <- read.csv('results/US_total_weekly_excess_ili_no_mizumoto.csv',stringsAsFactors = F) %>% 
  as.data.table
ILI[,date:=as.Date(date)]

ILI[,surge:=date>=as.Date('2020-03-07')]

ili <- ILI[surge==TRUE]
US[,likelihood:=dnorm(GrowthRate,mean=gr_national,sd=se_national)]

reps <- US[,list(GrowthRate=unique(GrowthRate),
                 likelihood=unique(likelihood),
                 peak_date=date[which.max(I)]),by=replicate]



n_plotting=1e3
reps_for_plotting <- sample(reps$replicate,size=n_plotting,prob=reps$likelihood)

ILI[,likelihood:=1]
ILI[,replicate:=15000]
Y[,likelihood:=1]
Y[,replicate:=15001]
cols=magma(5)
ggplot(US[replicate %in% reps_for_plotting & date<as.Date('2020-04-07')],
       aes(date,weekly_I,by=factor(replicate),alpha=likelihood))+
  geom_line()+
  geom_point(data=Y,aes(date,deaths),alpha=1,pch=8)+
  geom_point(data=ILI[surge==TRUE],aes(date,mean),alpha=1,pch=18,cex=3,col=cols[3])+
  scale_y_continuous(trans='log',name='Number of People',breaks=10^(0:9))+
  ggtitle('SEIR vs. ILI')+
  scale_alpha_continuous(range=c(0.01,0.4))
