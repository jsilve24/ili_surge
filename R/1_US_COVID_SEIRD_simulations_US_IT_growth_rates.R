rm(list=ls())
gc()
library(deSolve)
library(parallel)
library(data.table)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(truncnorm)
library(adaptivetau)
library(mgcv)
library(zoo)




# growth rate estimation --------------------------------------------------

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
Y <- Y[date>as.Date('2020-03-05') & date<=as.Date('2020-04-05')]
fit_national <- glm(deaths~date,family=poisson,data=Y)
glm(fit_national)
plot(fit_national) ## points 29-31 have high levereage - they will be removed
Y <- Y[1:28,]

fit_national <- glm(deaths~date,family=poisson,data=Y)
gr_national <- coef(fit_national)[2]
gr_national
# date 
# 0.2299108 

sd_national <- 0.01

# US_SEIR sim -------------------------------------------------------------


US_SEIRD_tl <- function(r=log(2)/3,cfr=0.004,S0=3.27e8,start_date=as.Date('2020-01-15'),
                     lag_onset_to_death=20,tf=200){
  ####### Functions
  seirDRates <- function(x,params,t){
    return(c(params$lambda,
             params$beta*x['S']*x['I'],
             params$p*x['E'],
             params$n*x['I'],
             params$m*x['S'],
             params$m*x['E'],
             params$m*x['I'],
             params$m_inf*x['I'],
             params$m*x['R']))
  }
  transitions <- list('Birth'=c('S'=1),
                      'Exposure'=c('S'=-1,'E'=1),
                      'Infection'=c('E'=-1,'I'=1),
                      'Recovery'=c('I'=-1,'R'=1),
                      'Sdeath'=c('S'=-1),
                      'Edeath'=c('E'=-1),
                      'Ideath'=c('I'=-1),
                      'IDdeath'=c('I'=-1,'D'=1),
                      'Rdeath'=c('R'=-1))
  ########## Simulation #################
  
  
  init.values <- c('S'=S0,'E'=0,'I'=1,'R'=0,'D'=0)
  times <- seq(0, 200, by = 0.01)
  m=8.685/100000/365
  n=1/10
  p=1/3
  
  m_inf <- cfr*n/(1-cfr)
  alpha=m_inf/m
  m_inf=m_inf-m
  
  
  # alpha=1.1
  c=(r+n+alpha*m)/p
  beta=c*(r+n+p)/S0
  lambda=11.8/1000/365
  
  params <- as.list(c('lambda'=lambda,
                      'm'=m,
                      'm_inf'=m_inf,
                      'beta'=beta,
                      'n'=n,
                      'p'=p))
  
  out <- ssa.adaptivetau(init.values,transitions,seirDRates,params,tf) %>% as.data.table
  while(max(out$I)<1e3){ ## ensures non-extinction
    out <- ssa.adaptivetau(init.values,transitions,seirDRates,params,tf) %>% as.data.table
  }
  out[,day:=ceiling(time)]
  
  out <- out[,list(S=S[.N],
                   E=E[.N],
                   I=I[.N],
                   R=R[.N],
                   D=D[.N]),by=day]
  out[,date:=start_date+day]
  out[,new_infections:=beta*S*I]
  out[,D:=shift(D,lag_onset_to_death)]
  out[is.na(D),D:=0]
  out[,new_deaths:=c(0,diff(D))]
  return(out)
}

# US death rate sims -----------------------------------------

set.seed(1)
n_sims <- 2e3
n_cores <- 7
r_set <- rnorm(n_sims,gr_national,sd_national)
cl <- makeCluster(n_cores)
clusterExport(cl,varlist=c('US_SEIRD_tl'))
clusterEvalQ(cl,{library(data.table)
  library(deSolve)
  library(magrittr)
  library(adaptivetau)})

US_seir_forecasts <- parLapply(cl,r_set,US_SEIRD_tl)
stopCluster(cl)
rm('cl')

for (i in 1:length(US_seir_forecasts)){
  US_seir_forecasts[[i]][,replicate:=i]
  US_seir_forecasts[[i]][,GrowthRate:=r_set[i]]
}

US_seir_forecasts <- rbindlist(US_seir_forecasts)


date_map <- data.table('day'=0:200,
                       'date'=seq(as.Date('2020-01-15'),as.Date('2020-08-02'),by='day'))
setkey(date_map,day)
setkey(US_seir_forecasts,day)
US_seir_forecasts <- US_seir_forecasts[date_map]


US_seir_forecasts[,weekly_I:=rollapply(I,sum,fill=NA,align='left',w=7),by=replicate]
saveRDS(US_seir_forecasts,file='results/US_seir_forecasts_USgr.Rd')
write.csv(US_seir_forecasts,file='results/US_seir_forecasts_USgr.csv')


# Italian growth rates ------------------------------------------------------------------------


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
gr_italy <- as.numeric(fit_italia$coefficients[2])


# Simulations -------------------------------------------------------------


set.seed(1)
n_sims <- 2e3
n_cores <- 7
r_set <- rnorm(n_sims,gr_italy,sd_national)
cl <- makeCluster(n_cores)
clusterExport(cl,'US_SEIRD_tl')
clusterEvalQ(cl,{library(data.table)
  library(deSolve)
  library(magrittr)
  library(adaptivetau)})

US_seir_forecasts <- parLapply(cl,r_set,US_SEIRD_tl)
stopCluster(cl)
rm('cl')

for (i in 1:length(US_seir_forecasts)){
  US_seir_forecasts[[i]][,replicate:=i]
  US_seir_forecasts[[i]][,GrowthRate:=r_set[i]]
}

US_seir_forecasts <- rbindlist(US_seir_forecasts)


date_map <- data.table('day'=0:200,
                       'date'=seq(as.Date('2020-01-15'),as.Date('2020-08-02'),by='day'))
setkey(date_map,day)
setkey(US_seir_forecasts,day)
US_seir_forecasts <- US_seir_forecasts[date_map]


US_seir_forecasts[,weekly_I:=rollapply(I,sum,fill=NA,align='left',w=7),by=replicate]
saveRDS(US_seir_forecasts,file='results/US_seir_forecasts_ITgr.Rd')
write.csv(US_seir_forecasts,file='results/US_seir_forecasts_ITgr.csv')



# Uniform distribution sims -----------------------------------------------



set.seed(1)
n_sims <- 2e3
n_cores <- 7

max_doubling_time <- 4
min_doubling_time <- 1.9
r_max <- log(2)/min_doubling_time
r_min <- log(2)/max_doubling_time
r_set <- r_min+runif(n_sims)*(r_max-r_min)

cl <- makeCluster(n_cores)
clusterExport(cl,'US_SEIRD_tl')
clusterEvalQ(cl,{library(data.table)
  library(deSolve)
  library(magrittr)
  library(adaptivetau)})

US_seir_forecasts <- parLapply(cl,r_set,US_SEIRD_tl)
stopCluster(cl)
rm('cl')

for (i in 1:length(US_seir_forecasts)){
  US_seir_forecasts[[i]][,replicate:=i]
  US_seir_forecasts[[i]][,GrowthRate:=r_set[i]]
}

US_seir_forecasts <- rbindlist(US_seir_forecasts)


date_map <- data.table('day'=0:200,
                       'date'=seq(as.Date('2020-01-15'),as.Date('2020-08-02'),by='day'))
setkey(date_map,day)
setkey(US_seir_forecasts,day)
US_seir_forecasts <- US_seir_forecasts[date_map]

US_seir_forecasts[,weekly_I:=rollapply(I,sum,fill=NA,align='left',w=7),by=replicate]
saveRDS(US_seir_forecasts,file='results/US_seir_forecasts_unifgr.Rd')
write.csv(US_seir_forecasts,file='results/US_seir_forecasts_unifgr.csv')