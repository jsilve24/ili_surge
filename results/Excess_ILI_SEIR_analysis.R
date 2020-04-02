rm(list=ls())
gc()
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(zoo)

US_seir_forecasts <- readRDS('RDS/US_seir_forecasts_unif_2_7.Rd')
setkey(US_seir_forecasts,replicate,date)
## We need to rollapply a sum for every week to capture weekly ILI
US_seir_forecasts[,weekly_I:=rollapply(I,sum,fill=NA,align='left',w=7),by=replicate]
saveRDS(US_seir_forecasts,'RDS/US_seir_forecasts_unif_weekly.Rd')

ILI <- read.csv('data/number_excess_ili_cases.csv',stringsAsFactors = F) %>% as.data.table
ILI[,week:=as.Date(week)]

US_seir_forecasts[date %in% setdiff(unique(ILI$week),NA)][!is.na(weekly_I),c('replicate','date','weekly_I')] %>%
  write.csv('data/weekly_I_across_replicates.csv')

US_seir_forecasts[,date:=shift(date,4,type='lead'),by=replicate]
US_seir_forecasts[date %in% setdiff(unique(ILI$week),NA)][!is.na(weekly_I),c('replicate','date','weekly_I')] %>%
  write.csv('data/weekly_I_across_replicates_4d_lag.csv')

rm(list=ls())
gc()

# distribution over replicates --------------------------------------------



# Calculating Clinical Rates -----------------------------------------------------
clinical_rate_calculator <- function(week='latest',method='gam',
                                     start_date='2020-01-15',time_onset_to_doc=0){
  ILI <- read.csv('data/posterior_samples_us_weekly_I.csv',stringsAsFactors = F) %>% as.data.table
  ILI[,week:=as.Date(date)]
  
  if (week=='latest'){
    wk<- max(ILI$week,na.rm=T)
  } else {
    wk <- week
  }
  
  if (class(start_date)!='Date'){
    start_date <- as.Date(start_date)
  }
  # 
  # if (statistic=='mean'){
  #   Excess_ILI <- ILI[week==wk,(Excess_ILI=sum(avg,na.rm=T))]
  # } else { #median
  #   Excess_ILI <- ILI[week==wk,(Excess_ILI=sum(p50,na.rm=T))]
  # }
  Excess_ILI <- mean(ILI[week==wk]$weekly_I)
  
  US_seir_forecasts <- readRDS('RDS/US_seir_forecasts_unif_weekly.Rd')
  setkey(US_seir_forecasts,replicate,date)
  date_shift = as.numeric(min(US_seir_forecasts$date)-start_date)
  US_seir_forecasts[,date:=shift(date,type='lead',n=date_shift),by=replicate]
  
  X = US_seir_forecasts[date==wk-time_onset_to_doc,
                        list(clinical_rate=Excess_ILI/weekly_I),by=GrowthRate]
  
  if (method=='loess'){
    fit <- loess(log(clinical_rate)~log(GrowthRate),data=X)
    X[,predicted_clinical_rate:=exp(fit$fitted)]
  } else {
    fit <- mgcv::gam(log(clinical_rate)~s(log(GrowthRate)),data=X)
    X[,predicted_clinical_rate:=exp(fit$fitted.values)]
  }
  X[,possible:=clinical_rate<1]
  X[,DoublingTime:=log(2)/GrowthRate]
  setkey(X,GrowthRate)
  return(list('Data'=X,'fit'=fit))
}

CR <- clinical_rate_calculator(time_onset_to_doc=4)
ggplot(CR$Data,aes(DoublingTime,clinical_rate,alpha=factor(possible)))+
  geom_point()+
  scale_alpha_manual(values=c(0.03,0.15))+
  scale_y_continuous(trans='log',breaks=10^(-5:5))+
  scale_x_continuous(trans='log',breaks=1:7)+
  geom_line(aes(DoublingTime,predicted_clinical_rate),lwd=2,alpha=1,col=rgb(0,.2,0.8))+
  theme(legend.position = 'none')+
  ggtitle('Estimated clinical rate, 4 day lag')
ggsave('figures/Clinical_rate_v_Doubling_time_4d_lag.png',height=8,width=8,units='in')

lag_times <- c(0,4,8)
for (lag in lag_times){
  dum <- clinical_rate_calculator(time_onset_to_doc=lag)
  if (lag==0){
    CR <- dum$Data
    CR[,delay_to_doc:=lag]
  } else {
    dum$Data[,delay_to_doc:=lag]
    CR <- rbind(CR,dum$Data)
  }
}

label_figs <- function(strings)  paste(strings,'days')




g1=ggplot(CR[DoublingTime<4],aes(DoublingTime,clinical_rate,alpha=factor(possible),by=delay_to_doc))+
  geom_point()+
  scale_alpha_manual(values=c(0.01,0.15))+
  scale_y_continuous(trans='log',breaks=10^(-5:3),limits=c(6e-3,1e3))+
  scale_x_continuous(trans='log',breaks=1:7)+
  geom_line(aes(DoublingTime,predicted_clinical_rate),lwd=2,alpha=1,col=rgb(0,.2,0.8))+
  theme_bw(base_size=25)+
  theme(legend.position = 'none')+
  geom_vline(xintercept = 3.5,lty=1,lwd=2,col='black',alpha=1)+
  geom_vline(xintercept = 1.91,lty=2,lwd=2,col='black',alpha=1)+
  facet_wrap(.~delay_to_doc,labeller=labeller(delay_to_doc=label_figs),ncol=length(lag_times))
g1
# ggsave('figures/Clinical_rate_v_Doubling_time_by_delay.png',height=8,width=14,units='in')
ggsave('figures/Clinical_rate_v_Doubling_time_by_delay_transparent.png',height=8,width=14,units='in',bg='transparent')



# Comparing ILI to US forecasts -------------------------------------------

##### load ILI data
ILI <- read.csv('data/number_excess_ili_cases.csv',stringsAsFactors = F) %>% as.data.table
ILI[,week:=as.Date(week)] 

##### Merge probabilities with US SEIR forecasts
weights <- read.csv('data/distribution_over_replicates.csv') %>% as.data.table
US <- readRDS('RDS/US_seir_forecasts_unif_weekly.Rd')
setkey(weights,replicate)
setkey(US,replicate)
US <- US[weights]
#### compute doubling time
US[,DoublingTime:=log(2)/GrowthRate]

##### posterior samples - will use to visualize mean ILI vs. Forecasts
ili <- read.csv('data/posterior_samples_us_weekly_I.csv') %>% as.data.table
ili[,date:=as.Date(date)-4]
ili <- ili[,list(weekly_I=mean(weekly_I),
          replicate=unique(as.numeric(date))+14000,
          DoublingTime=2.9,
          probability=max(US$probability)),by=date]



##### function to convert probabilities to weights for transparency setting alpha=weight
convertP<-function(x){
  x[x==0]<-min(x[x>0])/10
  y <- log(x)
  y <- (y-min(y))/(max(y)-min(y))
  return(y)
}
US[,weight:=convertP(probability)]
ili[,weight:=convertP(probability)]


########## subsampling replicates for plotting
### use this if many probabilities are = 0
n_reps_for_plotting <- 3e3
n_nonzero <- length(unique(US[probability>0,replicate]))
replicates <- c(unique(US[probability>0,replicate]),
                sample(unique(US[probability==0]$replicate),
                       n_reps_for_plotting-n_nonzero))

r_set <- unique(US$GrowthRate)
r_slow <- log(2)/3.5
r_fast <- log(2)/1.91

rep_slow <- unique(US$replicate)[order(abs(r_slow-r_set))]
rep_fast <- unique(US$replicate)[order(abs(r_fast-r_set))]
## finding a fast replicate: which curve, with unknown scaling factor, 
## comes closest to interpolating observed ILI?


### use this if most probabilities > 0
# probs <- US[,list(prob=unique(probability)),by=replicate]
# replicates <- sample(probs$replicate,1e3,prob = probs$prob,replace=F)

########### plotting forecasts vs. ILI
g2=ggplot(US[replicate %in% replicates & date<=as.Date('2020-03-15')],
          aes(date,weekly_I,color=DoublingTime,by=factor(replicate),alpha=weight))+
  geom_line()+
  ggtitle('COVID Forecasts vs. ILI')+
  scale_y_continuous(trans='log',breaks=10^(1:9),name='Past Week Infectious-Days')+
  # geom_hline(yintercept = ili$weekly_I,lty=3,col='red',lwd=2)+
  theme_bw(base_size = 25)+
  theme(legend.position = 'none')+
  geom_point(data=ili,col='red',pch=16,cex=8)+
  geom_line(data=US[replicate==rep_slow[1] & date<=as.Date('2020-03-15')],lty=1,lwd=2,col='black',alpha=1)+
  geom_line(data=US[replicate==rep_fast[1] & date<=as.Date('2020-03-15')],lty=2,lwd=2,col='black',alpha=1)+
  scale_alpha_continuous(range=c(0.01,1))
g2
ggsave('figures/Epidemic_curve_estimation.png',height=8,width=8,unit='in',bg='transparent')



# Forecast plotting -------------------------------------------------------

g3=ggplot(US[replicate %in% replicates & date>as.Date('2020-02-28')],
          aes(date,I,by=factor(replicate),alpha=weight))+
  geom_line()+
  ggtitle('US SEIR forecasts')+
  theme(legend.position = c(0.9,0.85))+
  scale_alpha_continuous(range=c(0.01,0.5))+
  scale_y_continuous('Number Infected')+
  theme_bw(base_size=25)+
  theme(legend.position='none')
g3
ggsave('figures/SEIR_forecasts.png',height=10,width=14,units='in')



# Plot forecasts + clinical rate estimation -------------------------------
ggarrange(g3,g1,nrow=2)
ggsave('figures/Clinical_rate_estimation.png',height=12,width=10,units='in')



# Clinical Rate of avg growth rate ----------------------------------------

r=US[,weighted.mean(GrowthRate,probability)] ## probability-weighted avg growth rate
r_state_median <- log(2)/2.455
r_deaths <- log(2)/2.52
CR <- clinical_rate_calculator(time_onset_to_doc = 4)  ## assume 4-day lag from onset of infectiousness to doc visit
CR$fit %>% predict(newdata=data.frame('GrowthRate'=r_state_median)) %>% exp   ## clinical rate for state growths
CR$fit %>% predict(newdata=data.frame('GrowthRate'=r_deaths)) %>% exp   ## clinical rate for rate of growth of COVID deaths
