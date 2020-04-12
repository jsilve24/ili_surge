rm(list=ls())
gc()
library(deSolve)
library(parallel)
library(data.table)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(truncnorm)
library(mgcv)

# US_SEIR sim -------------------------------------------------------------


SEIR <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dS=zeta-beta*S*I-omega_b*S
    dE=beta*S*I-p*E-omega_b*E
    dI=p*E-w*I-n*I
    dR=n*I-omega_b*R
    
    list(c(dS, dE, dI, dR))
  })
}

GillespiEpi <- function(r,incubation_period=3,recovery_time=10,max_val=1e2,maxT=50,
                        S0=3.27e8,omega_b=8.685/100000/365,omega_i=1.1*omega_b){
  
  c=(r+1/10+omega_i) ## ratio of I/E that ensures exponential growth at rate r
  n=1/recovery_time
  #expected growth in d(I+E)/dt=(beta-gamma)*I
  # beta=c*(r+1/incubation_period+omega_i)/S0
  p=1/incubation_period
  c=(r+n+omega_i)/p
  beta=c*(r+n+p)/S0
  zeta=omega_b*S0
  
  
  z <- data.table('E'=0,'I'=1,'R'=0,'time'=0)
  
  
  while((z[.N]$E+z[.N]$I)<max_val & z[.N]$time<maxT){
    
    propensities <- c('Exposure'=beta*z[.N,I]*S0,
                 'Infection'=p*z[.N]$E,
                 'Recovery'=n*z[.N]$I)
    
    if (z[.N]$I==1 & z[.N]$E==0){
      dt <- rexp(1,sum(propensities[1:2]))
      event <- sample(names(propensities)[1:2],size=1,prob=propensities[1:2])
    } else {
      dt <- rexp(1,sum(propensities))
      event <- sample(names(propensities),size=1,prob=propensities)
    }
    if (event=='Exposure'){
      z <- rbind(z,z[.N]+c(1,0,0,dt))
    } else if (event=='Infection'){
      z <- rbind(z,z[.N]+c(-1,1,0,dt))
    } else {
      z <- rbind(z,z[.N]+c(0,-1,1,dt))
    }
  }
  return(z)
}

US_SEIR <- function(r){
  zz <- GillespiEpi(r=r)
  zz[,S:=3.27e8-I]
  zz[,day:=ceiling(time)]
  t0=zz[.N,time]
  I0=zz[.N,I]
  E0=zz[.N,E]
  S0=zz[.N,S]
  state <- c('S'=S0,'E'=E0,'I'=I0,'R'=0)
  times <- seq(t0, 200, by = 0.01)
  alpha=1.1
  n=1/10
  # c=(r+n+omega_i)/p ## ratio of I/E that ensures exponential growth at rate r
  
  #expected growth in d(I+E)/dt=(beta-gamma)*I
  # beta=c*(r+1/incubation_period+omega_i)/S0
  p=1/3
  # c=(r+n)/p
  omega_b <- 8.685/100000/365 ## baseline mortality rate
  omega_i <- 1.1*omega_b      ## infected/symptomatic mortality rate
  c=(r+n+omega_i)/p
  beta=c*(r+n+p)/S0
  zeta=11.8/1000/365
  parameters <- c('omega_b'=omega_b,
                  'zeta'=zeta,
                  'n'=n,
                  'p'=p,
                  'beta'=beta,
                  'w'=omega_i)
  
  out <- ode(y = state, times = times, func = SEIR, parms = parameters) %>% as.data.table
  out[,day:=ceiling(time)]
  
  zz <- rbind(zz[,c('time','day','S','E','I','R')],
              out[,c('time','day','S','E','I','R')])
  zz <- zz[,list(S=S[.N],
                 E=E[.N],
                 I=I[.N],
                 R=R[.N]),by=day]
  return(zz)
}


# Uniform growth rate simulations -----------------------------------------

set.seed(1)
n_sims <- 1.4e4
n_cores <- 7
max_doubling_time <- 7
min_doubling_time <- 1.5
r_max <- log(2)/min_doubling_time
r_min <- log(2)/max_doubling_time
r_set <- r_min+runif(n_sims)*(r_max-r_min)
cl <- makeCluster(n_cores)
clusterExport(cl,varlist=c('SEIR','GillespiEpi','US_SEIR'))
clusterEvalQ(cl,{library(data.table)
  library(deSolve)
  library(magrittr)})

US_seir_forecasts_unif <- parLapply(cl,r_set,US_SEIR)
stopCluster(cl)
rm('cl')

for (i in 1:length(US_seir_forecasts_unif)){
  US_seir_forecasts_unif[[i]][,replicate:=i]
  US_seir_forecasts_unif[[i]][,GrowthRate:=r_set[i]]
}

US_seir_forecasts_unif <- rbindlist(US_seir_forecasts_unif)


date_map <- data.table('day'=0:200,
                       'date'=seq(as.Date('2020-01-15'),as.Date('2020-08-02'),by='day'))
setkey(date_map,day)
setkey(US_seir_forecasts_unif,day)
US_seir_forecasts_unif <- US_seir_forecasts_unif[date_map]

saveRDS(US_seir_forecasts_unif,file='RDS/US_seir_forecasts_unif_2_7.Rd')
write.csv(US_seir_forecasts_unif,file='RDS/US_seir_forecasts_unif_2_7.csv')
