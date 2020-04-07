library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)


# load & format data ------------------------------------------------------

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


ggplot(Y,aes(date,deaths))+
  geom_point()+
  geom_line()+
  geom_smooth(method='glm',method.args=list('family'=poisson))

### test: is doubling time national < 3.5?
gr_bound <- log(2)/3.5
pnorm(gr_bound,mean = gr_national,sd=se_national)
# 5.274348e-21

gr_national
# 0.2344975
log(2)/gr_national
# 2.955883 