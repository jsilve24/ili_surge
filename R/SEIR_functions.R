library(data.table)
library(deSolve)
library(magrittr)
SEIRD <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dS=lambda-beta*S*I-m*S
    dE=beta*S*I-p*E-m*E
    dI=p*E-w*I-n*I
    dR=n*I-m*R
    dD=w*I

    list(c(dS, dE, dI, dR, dD))
  })
}

US_SEIRD <- function(r,cfr=0.004,S0=3.27e8,start_date=as.Date('2020-01-15'),
                     lag_onset_to_death=20){

  state <- c('S'=S0,'E'=0,'I'=1,'R'=0,'D'=0)
  times <- seq(0, 200, by = 0.01)
  m=8.685/100000/365
  n=1/9
  p=1/3

  m_inf <- cfr*n/(1-cfr)
  alpha <- m_inf/m

  # alpha=1.1
  c=(r+n+alpha*m)/p
  beta=c*(r+n+p)/S0
  lambda=11.8/1000/365
  parameters <- c('m'=m,
                  'lambda'=lambda,
                  'n'=n,
                  'p'=p,
                  'beta'=beta,
                  'w'=alpha*m)

  out <- ode(y = state, times = times, func = SEIRD, parms = parameters) %>% as.data.table
  out[,day:=ceiling(time)]

  out <- out[,list(S=S[.N],
                   E=E[.N],
                   I=I[.N],
                   R=R[.N],
                   D=D[.N]),by=day]
  out[,date:=seq(start_date,start_date+200,by='day')]
  out[,new_infections:=beta*S*I]
  out[,D:=shift(D,lag_onset_to_death)]
  out[is.na(D),D:=0]
  out[,new_deaths:=c(0,diff(D))]
  return(out)
}


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
  n=1/9
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
