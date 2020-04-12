rm(list=ls())
gc()
library(data.table)
library(magrittr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(zoo)
source('R/SEIR_functions.R')

####### NOTE: In this script, we compare ILI to SEIR curves matchign the growth rate of deaths.

####### For march 1, we use the SEIR curve estimate of ILI NOT the empirical ILI
####### The reason: March 1's ILI shows low correlation w/ COVID
####### thus ILI including & prior to first week of march are thrown out

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
Y <- Y[date>as.Date('2020-03-04')]
fit_national <- glm(deaths~date,family=poisson,data=Y)
glm(fit_national)
# plot(fit_national) ## points 29-31 have high levereage - they will be removed
Y <- Y[1:28,]

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

US[,likelihood:=dnorm(GrowthRate,mean=gr_national,sd=se_national)]

reps <- US[,list(GrowthRate=unique(GrowthRate),
                 likelihood=unique(likelihood),
                 peak_date=date[which.max(I)]),by=replicate]

ili <- ILI[surge==TRUE,c('date','mean')]
ili <- rbind(US[date==as.Date('2020-03-01'),
                list(date=unique(date),
                     mean=weighted.mean(weekly_I,w=likelihood))],ili)
ili <- rbind(ili,ili)
ili[1:(.N/2),scaled:=FALSE]
ili[is.na(scaled),scaled:=TRUE]

scaling_factors <- read.csv('results/scaling_admission_rates.csv') %>% as.data.table
scaling_factors[,date:=as.Date(date)]
setkey(scaling_factors,date)
setkey(ili,date,scaled)
ili <- scaling_factors[ili]
ili[is.na(scaling_factor),scaling_factor:=1]

ili[scaled==TRUE,mean:=mean*scaling_factor]

n_plotting=1e3
reps_for_plotting <- sample(reps$replicate,size=n_plotting,prob=reps$likelihood)

ili[,likelihood:=1]
ili[,replicate:=15000]
Y[,likelihood:=1]
Y[,replicate:=15001]
cols=magma(5)
# ggplot(US[replicate %in% reps_for_plotting & date<as.Date('2020-04-07')],
#        aes(date,weekly_I/7,by=factor(replicate),alpha=likelihood))+
#   geom_line()+
#   geom_point(data=Y,aes(date,deaths),alpha=1,pch=8,col='red')+
#   geom_point(data=ili[scaled==FALSE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=18,cex=6,col=cols[2])+
#   geom_line(data=ili[scaled==FALSE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[2])+
#   geom_point(data=ili[scaled==TRUE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=18,cex=6,col=cols[4])+
#   geom_line(data=ili[scaled==TRUE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[4])+
#   scale_y_continuous(trans='log',name='Number of People',breaks=10^(0:9))+
#   ggtitle('SEIR vs. ILI')+
#   scale_alpha_continuous(range=c(0.01,0.2))+
#   theme_bw(base_size = 25)+
#   theme(legend.position = 'none')
# ggsave('figures/death_rate_excess_ili_SEIR_matching.png',height=8,width=8,units='in')

N_week <- 7*3.27e8

# ggplot(US[replicate %in% reps_for_plotting & date<as.Date('2020-04-07')],
#        aes(date,weekly_I/7,by=factor(replicate),alpha=likelihood))+
#   geom_line()+
#   geom_point(data=Y,aes(date,deaths),alpha=1,pch=8,col='red',cex=3)+
#   geom_point(data=ili[scaled==TRUE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=18,cex=6,col=cols[3])+
#   geom_line(data=ili[scaled==TRUE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[3])+
#   scale_y_continuous(trans='log',name='Number of People',breaks=10^(0:9))+
#   ggtitle('SEIR vs. ILI')+
#   scale_alpha_continuous(range=c(0.01,0.2))+
#   theme_bw(base_size = 25)+
#   theme(legend.position = 'none')
# ggsave('figures/death_rate_excess_ili_SEIR_matching_scaled_only.png',height=8,width=8,units='in')



# Simplified figure -------------------------------------------------------
## running one iteration to find cfr
r <- as.numeric(gr_national)
seir <- US_SEIRD(r,cfr=0.005)

setkey(seir,date)
setkey(Y,date)
yy <- seir[,c('date','new_deaths')][Y]
b=glm(deaths~new_deaths+0,data=yy)$coefficients


cfr <- as.numeric(0.005/b)
seir <- US_SEIRD(r,cfr=cfr)
seir[,weekly_I:=rollapply(I,FUN=sum,w=7,align='left',fill=NA)]
ili <- ili[date>as.Date('2020-03-01')]
# ili[date==as.Date('2020-03-01'),mean:=seir[date==as.Date('2020-03-01')]$weekly_I]  ### 
setkey(seir,date)
setkey(ili,date)

ili$weekly_I <- NULL
ili <- ili[seir[,c('date','weekly_I')]]
ili_sc <- ili[scaled==TRUE]
names(ili_sc)[3] <- 'mean_scaled'
i2 <- cbind(ili[scaled==FALSE,c('date','mean')],ili_sc[,c('mean_scaled','weekly_I')])


ggplot(seir[date<=as.Date('2020-04-15')],aes(date,weekly_I/7))+
  geom_line(lwd=2)+
  geom_line(data=seir[date<=max(Y$date)],aes(date,new_deaths))+
  geom_ribbon(data=ili[scaled==TRUE],aes(ymin=mean/7,ymax=weekly_I/7),alpha=0.4,fill=cols[4])+
  geom_ribbon(data=i2,aes(ymin=mean/7,ymax=mean_scaled/7),fill=rgb(0,0.1,0.5),alpha=0.8)+
  geom_point(data=Y,aes(date,deaths),alpha=1,pch=16,col='red',cex=3)+
  geom_point(data=ili[scaled==FALSE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=18,cex=6,col=cols[2])+
  geom_line(data=ili[scaled==FALSE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[2])+
  geom_point(data=ili[scaled==TRUE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=18,cex=6,col=cols[3])+
  geom_line(data=ili[scaled==TRUE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[3])+
  scale_y_continuous(trans='log',name='Number of People',breaks=10^(0:9),limits=c(1,3e8))+
  ggtitle('SEIR vs. ILI')+
  scale_alpha_continuous(range=c(0.01,0.2))+
  theme_bw(base_size = 25)+
  theme(legend.position = 'none')

seir[,replicate:=15002]
seir[,likelihood:=1]
i2[,replicate:=15003]
i2[,likelihood:=1]

ggplot(US[replicate %in% reps_for_plotting & date<as.Date('2020-04-07')],
       aes(date,weekly_I/7,by=factor(replicate),alpha=likelihood))+
  geom_line()+
  geom_line(data=seir[date<as.Date('2020-04-07')],aes(date,weekly_I/7),lwd=2,alpha=1,col='red',lty=2)+
  geom_line(data=seir[date<=max(Y$date)],aes(date,new_deaths),alpha=1)+
  geom_ribbon(data=ili[scaled==TRUE],aes(ymin=mean/7,ymax=weekly_I/7),alpha=0.8,fill=cols[4])+
  geom_ribbon(data=i2,aes(ymin=mean/7,ymax=mean_scaled/7),fill=rgb(0,0.1,0.5),alpha=0.8)+
  geom_point(data=Y,aes(date,deaths),alpha=1,pch=4,cex=3)+
  geom_point(data=ili[scaled==FALSE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=19,cex=6,col=cols[2])+
  geom_line(data=ili[scaled==FALSE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[2])+
  geom_point(data=ili[scaled==TRUE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=18,cex=6,col=cols[3])+
  geom_line(data=ili[scaled==TRUE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[3])+
  scale_y_continuous(trans='log',name='Number of People',breaks=10^(0:9),limits=c(1,3e8))+
  ggtitle('SEIR vs. ILI')+
  scale_alpha_continuous(range=c(0.01,0.2))+
  theme_bw(base_size = 25)+
  theme(legend.position = 'none')

ggsave('figures/death_rate_excess_ili_SEIR_matching_all.png',height=8,width=8,units='in')


excess_ili <- ili[date==as.Date('2020-03-08') & scaled==FALSE,mean]

clinical_rate <- excess_ili/US[replicate %in% reps_for_plotting & date==as.Date('2020-03-08')]$weekly_I
clinical_rate <- clinical_rate[clinical_rate<1]
mean_cr <- mean(clinical_rate)
mean_cr

data.frame('cr'=clinical_rate[clinical_rate<1]) %>%
  ggplot(aes(cr))+
  geom_histogram(alpha=0.3,col='black')+
  scale_y_continuous(trans='log',breaks=10^(0:10),name='Number of Simulations')+
  geom_vline(xintercept = mean_cr,lwd=2)+
  ggtitle('Implied Clinical rate')+
  scale_x_continuous(name='Clinical Rate')+
  theme_bw(base_size = 25)

ggsave('figures/clinical_rate_estimation_US_death_gr.png',height=8,width=8)


# Scenario 2: Italy deaths ---------------------------------------------------

Y <- read.csv('data/Italy/dpc-covid19-ita-andamento-nazionale.csv',stringsAsFactors = F) %>%
  as.data.table
colnames(Y)[c(1,11)] <- c('date','deaths')
Y[,date:=as.Date(date)]
Y[,deaths:=c(deaths[1],diff(deaths))]
ggplot(Y,aes(date,deaths))+
  geom_point()+
  scale_y_continuous(trans='log')
## will use up to March 12 data
Y <- Y[date<=as.Date('2020-03-12')]

fit_italia <- glm(deaths~date,family=poisson,data=Y[date<=as.Date('2020-03-12')])
gr_italia <- fit_italia$coefficients[2]
se_italia <- sqrt(vcov(fit_italia)[4])
US[,likelihood:=dnorm(GrowthRate,mean=gr_italia,sd=se_italia)]
r <- as.numeric(fit_italia$coefficients[2])


seir <- US_SEIRD(r,cfr=0.005)

setkey(seir,date)
setkey(Y,date)
yy <- seir[,c('date','new_deaths')][Y]
b=glm(deaths~new_deaths+0,data=yy)$coefficients


cfr <- as.numeric(0.005*b)
seir <- US_SEIRD(r,cfr=cfr)
seir[,weekly_I:=rollapply(I,FUN=sum,w=7,align='left',fill=NA)]


ili[date==as.Date('2020-03-01'),mean:=seir[date==as.Date('2020-03-01')]$weekly_I]  ### 
setkey(seir,date)
setkey(ili,date)
ili$weekly_I <- NULL
ili <- ili[seir[,c('date','weekly_I')]]
ili_sc <- ili[scaled==TRUE]
names(ili_sc)[3] <- 'mean_scaled'
i2 <- cbind(ili[scaled==FALSE,c('date','mean')],ili_sc[,c('mean_scaled','weekly_I')])


seir[,replicate:=15002]
seir[,likelihood:=1]
i2[,replicate:=15003]
i2[,likelihood:=1]
Y[,replicate:=15000]
Y[,likelihood:=1]


reps <- US[,list(GrowthRate=unique(GrowthRate),
                 likelihood=unique(likelihood),
                 peak_date=date[which.max(I)]),by=replicate]
reps_for_plotting <- sample(reps$replicate,size=n_plotting,prob=reps$likelihood)

ggplot(US[replicate %in% reps_for_plotting & date<as.Date('2020-04-07')],
       aes(date,weekly_I/7,by=factor(replicate),alpha=likelihood))+
  geom_line()+
  geom_line(data=seir[date<as.Date('2020-04-07')],aes(date,weekly_I/7),lty=2,lwd=2,alpha=1,col='red')+
  geom_line(data=seir[date<=max(Y$date)],aes(date,new_deaths),alpha=1)+
  geom_point(data=Y,aes(date,deaths),alpha=1,pch=4,cex=3)+
  geom_ribbon(data=ili[date>as.Date('2020-03-01') & scaled==TRUE],aes(ymin=mean/7,ymax=weekly_I/7),alpha=0.8,fill=cols[4])+
  geom_ribbon(data=i2[date>as.Date('2020-03-01')],aes(ymin=mean/7,ymax=mean_scaled/7),fill=rgb(0,0.1,0.5),alpha=0.8)+
  geom_point(data=ili[date>as.Date('2020-03-01') & scaled==FALSE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=18,cex=6,col=cols[2])+
  geom_line(data=ili[date>as.Date('2020-03-01') & scaled==FALSE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[2])+
  geom_point(data=ili[date>as.Date('2020-03-01') & scaled==TRUE & date>as.Date('2020-03-01')],aes(date,mean/7),alpha=1,pch=18,cex=6,col=cols[3])+
  geom_line(data=ili[date>as.Date('2020-03-01') & scaled==TRUE],aes(date,mean/7),alpha=1,lwd=1.5,col=cols[3])+
  scale_y_continuous(trans='log',name='Number of People',breaks=10^(0:9),limits=c(1,3e8))+
  ggtitle('SEIR vs. ILI, Italian growth rate')+
  scale_alpha_continuous(range=c(0.01,0.2))+
  theme_bw(base_size = 25)+
  theme(legend.position = 'none')

ggsave('figures/death_rate_excess_ili_SEIR_matching_all_italy.png',height=8,width=8,units='in')


excess_ili <- ili[date==as.Date('2020-03-08') & scaled==FALSE,mean]

clinical_rate <- excess_ili/US[replicate %in% reps_for_plotting & date==as.Date('2020-03-08')]$weekly_I
clinical_rate <- clinical_rate[clinical_rate<1]
mean_cr <- mean(clinical_rate)
mean_cr

data.frame('cr'=clinical_rate[clinical_rate<1]) %>%
  ggplot(aes(cr))+
  geom_histogram(alpha=0.3,col='black')+
  scale_y_continuous(trans='log',breaks=10^(0:10),name='Number of Simulations')+
  geom_vline(xintercept = mean_cr,lwd=2)+
  ggtitle('Implied Clinical rate')+
  scale_x_continuous(name='Clinical Rate')+
  theme_bw(base_size = 25)
ggsave('figures/clinical_rate_estimation_Italy_death_gr.png',height=8,width=8)

