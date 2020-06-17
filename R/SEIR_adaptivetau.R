# install.packages('adaptivetau')
library(adaptivetau)
library(magrittr)
library(data.table)
library(ggplot2)

S0=3.24e8
r=log(2)/3
cfr=0.005
tf=200
m=8.685/100000/365
n=1/10
p=1/3

m_inf <- cfr*n/(1-cfr)
alpha <- m_inf/m
c=(r+n+alpha*m)/p
beta=c*(r+n+p)/S0
lambda=11.8/1000/365

init.values = c('S'=S0,'E'=0,'I'=1,'R'=0)

params <- as.list(c('lambda'=lambda,
            'm'=m,
            'm_inf'=m_inf,
            'beta'=beta,
            'n'=n,
            'p'=p))

transitions <- list('Birth'=c('S'=1),
                    'Exposure'=c('S'=-1,'E'=1),
                    'Infection'=c('E'=-1,'I'=1),
                    'Recovery'=c('I'=-1,'R'=1),
                    'Sdeath'=c('S'=-1),
                    'Edeath'=c('E'=-1),
                    'Ideath'=c('S'=-1),
                    'Rdeath'=c('R'=-1))

seirRates <- function(x,params,t){
  return(c(params$lambda,
           params$beta*x['S']*x['I'],
           params$p*x['E'],
           params$n*x['I'],
           params$m*x['S'],
           params$m*x['E'],
           params$m_inf*x['I'],
           params$m*x['R']))
}

r=ssa.adaptivetau(init.values, transitions, seirRates, params, tf=200) %>% as.data.table



ggplot(r,aes(time,I))+
  geom_line()
