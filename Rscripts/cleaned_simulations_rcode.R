
### bull trout integrated population model simulations ###

# 5.29.25
# madaline cochrane


library(egg)
library(ggplot2)
library(tidyr)
library(dplyr)
library(runjags)
library(mcmcplots)
library(egg)
library(tidyhydat)

# set working directory

# read in datafile for final bull trout ipm from 'cleaned_bt_ipm_rcode.R' script
full.mod2 <-readRDS(file="R/full.mod.new.ipm") # new one

# start here ##
#### simulating ipm bull trout model ####
txtstring2 <- '
# likelihood
data{ 
  
  for (t in 1:(nYears-4)){
    # Ricker model for recruitment
    alpha[t] <- exp(alpha_sr + beta.low.flow*low.flow.recruitment[t])
                      
    # Ricker model for recruitment
    log_predR[t+4] <- log(b[t]) + log(alpha[t]) - beta_sr*b[t] 
    
    b[t]<-ifelse(biomass[t]<0,0.0001,biomass[t])
    biomass[t]<-weight[t]*Nadults[t]*mean.fish.per.redd # only females spawn
    
    #Process Likelihood
    logR_Obs[t+4] ~ dnorm(log_predR[t+4], tauR_process) 
    yr1subs[t+4] <- exp(logR_Obs[t+4]) 
    
    Nsubadults[t+4] <- yr1subs[t+4] + (1-transition.prob)*S.sa*Nsubadults[t+3] 
  }
  
  for(t in 1:(nYears-1)){
    S.a.step1[t] <- a.int + beta.dens*(NumbLake[t]/10000) + beta.kok*kok[t] + beta.reg.change*ca.reg.changes[t] # scaling raw abundance
    mean.S.a[t] <- 1/(1+exp(-S.a.step1[t]))
    S.a[t] ~ dnorm(mean.S.a[t], tau.S.a)T(0,1) # yr-specific survival includes process error 
    
    NA[t]<-ifelse(bt.harv[t]>Nadults[t],0,Nadults[t]) # slightly modified from original (another trick to keep N>0)
    Nadults[t+1] <- (NA[t]-bt.harv[t])*S.a[t] + transition.prob*S.sa*Nsubadults[t]
    
    lambdaA[t] <-Nadults[t+1]/Nadults[t] # change in adult population size over time
    lambda[t] <-Nsubadults[t+1]/Nsubadults[t]
    NumbLake[t]<-Nsubadults[t]+Nadults[t] # using total bt population to understand density dependence
  }
  
  Nadults[1] <- 2898 # 2020 value
  Nsubadults[1] <- 1366 # 2020 value
  Nsubadults[2] <- 1457 # 2021 value
  Nsubadults[3] <- 2433 # 2022 value
  Nsubadults[4] <- 2700 # 2023 value
  
} 
model{
  fake <- 0
}
'





####################################################
### run simulations for figure in ms with 95% cis


### need to bring in some data first
# read in bull trout gillnet catch data by year
bull.raw.net.dat<-read.csv("Data/BullTroutCatchRawData.csv",header=TRUE) 

bull.raw.net.dat$Date <- # change date from character to posixct
  as.POSIXct(bull.raw.net.dat$Date, tz = "", format="%m/%d/%Y") 

bull.raw.net.dat <- bull.raw.net.dat %>%
  filter(NetType == 2)%>% # only use type 2 data for bull trout
  filter(Season == "Spring") # only using spring nets

# want to summarize adults by net and year
adult.bull.net.sum<-bull.raw.net.dat%>%
  mutate(Stage = if_else(Length_mm >= 500,"A","SA"))%>%  # adults are >= 500 mm
  filter(Stage == "A")%>%
  group_by(Section, YYYY,NetNumb) %>%
  summarise(n = n())
#head(adult.bull.net.sum)

# weight estimates by year from adult bull trout cpue data
weight.yr <- bull.raw.net.dat %>% 
  group_by(YYYY)%>%
  mutate(Stage = if_else(Length_mm >= 500,"A","SA"))%>% 
  filter(Stage == "A")%>%
  mutate(wt_lbs = Weight_g*0.00220462)%>% # converging to lbs
  summarise(wt=mean(wt_lbs,na.rm=T))%>%
  rename(Season = YYYY)

# interpolating missing years 
weight.yr <- weight.yr %>%
  add_row(Season = 1977)%>% # no data on weight these years
  add_row(Season = 1979)%>%
  add_row(Season = 1983)%>%
  add_row(Season = 1984)%>%
  add_row(Season = 1985)%>%
  arrange(Season)

#Interpolate missing years if only a 1 year gap
weight.yr$wt[which(is.na(weight.yr$wt))] <- lapply(which(is.na(weight.yr$wt)),FUN=function(x){
  mean(c(weight.yr$wt[x-1],weight.yr$wt[x+1]))
}) %>% unlist()

# 3 yr gap from 1983 - 1985 so interpolating manually (3.790477+3.857534)/2
weight.yr$wt[8]<- 3.824006
weight.yr$wt[9]<- 3.824006
weight.yr$wt[10]<- 3.824006


### get discharge data from bull river
bull.river.q<-hy_daily_flows(station_number = "08NG002", 
                             start_date = "1976-01-01", 
                             end_date = "2022-12-31") # only archived to end of 2022 (as of this analysis 7.3.24)
# value = daily mean discharge in cms

bull.river.q$Date <- # change date from character to posixct
  as.POSIXct(bull.river.q$Date, tz = "", format="%Y-%m-%d") 

# exported 2023 (un-archived) bull river data (real-time) to complete dataset
bull.river.q.23<-read.csv("Data/discharge/bullriver_2023.csv",header=TRUE) 

# Value is discharge in cms
bull.river.q.23$Date <- # change date from character to posixct
  as.POSIXct(bull.river.q.23$Date, tz = "", format="%m/%d/%Y %H:%M") 

bull.river.q.23 <- bull.river.q.23 %>%
  mutate(Date2 = lubridate::date(Date))
bull.river.q.23 <- bull.river.q.23 %>% 
  group_by(Date2)%>%
  summarise(Value2 = mean(Value))%>%
  rename(Date = Date2)%>%
  rename(Value = Value2)

# combine historical discharge data with 2023 data 
bull.river.q2<-full_join(bull.river.q,bull.river.q.23,by=c("Date","Value"))

bull.river.q2 <- bull.river.q2 %>%
  arrange(Date)%>% tidyr::drop_na(Date)


#### calculating the minimum 14-day summer low flow for bull river
# will need to lag to explain recruitment 
rollingQ2 <- 
  zoo::rollapplyr(bull.river.q2$Value,  14, mean,  align='right') # 14-day mean flow from 14 days prior
# need to add 13 NAs to rollingQ2 column (at beginning of dataset)
rollingQ2 <- c( rep(NA,13),rollingQ2)
bull.river.q2$rollingQ2 <- rollingQ2 

bull.river.low <- bull.river.q2 %>%
  mutate(month = lubridate::month(Date))%>%
  mutate(Year = lubridate::year(Date))%>%
  filter(month == 7 | month == 8 | month == 9)%>% # only want summer months
  group_by(Year)%>%
  summarise(min(rollingQ2))%>%
  dplyr::rename(summerLowflow = `min(rollingQ2)`)%>%
  mutate(scale_summerLowflow = scale(summerLowflow)) # scaling discharge to make comparable to other covariates


### read in kokanee data
predict <-readRDS(file="R/kokanee_biomass")
# for 3-yr running average of kokanee we need to give the first 2 yrs (1980, 81) the same values as 1983
predict$running_accessible_biomass[1]<-predict$running_accessible_biomass[3] # giving first and second values the same as third
predict$running_accessible_biomass[2]<-predict$running_accessible_biomass[3]

predict <- predict %>%
  mutate(running_biomass_scale = scale(running_accessible_biomass)) # scaling kokanee biomass to make comptabile with discharge, pool elevation data




########################################################################
########################################################################
## comparing harvest changes across different simulations 

# first simulation 
# co-variates = 500 harvest, mean kokanee availability, canadian harvest not allowed (ie catch and release only)
results <- as.data.frame(c())
iterations<-100 

### takes some time to run all of these simulations 
for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 # DATA
                 bt.harv = rep(500,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(0,54), # mean kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54), # mean flow
                 ca.reg.changes = rep(0,54), # catch-and-release in canada
                 mean.fish.per.redd = 0.5,
                 
                 # PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant 
  
  # parameters monitored 
  params <- c("Nadults","yr1subs","lambdaA") # add lambda to get median lambda over last 20 yrs to see if population is declining 
  # was Nsubadults instead of yr1subs
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
} 

# Results
results.sa<-results[,55:104] # new sub-adults - starts at year 5 (1984)
results.lambdaA<-results[,105:157] # lambda
results<-results[,1:54] # adults

results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
head(results.tibs)
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.500.harvest <- results.tibs %>% 
  rowwise() %>% # sums across all iterations/rows
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>%
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))%>%
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))

# view to make sure everything is working
ggplot(data=sums.500.harvest,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()+
  geom_hline(yintercept=2907.8176) # current population

# sub-adults 
results.tibs<-t(results.sa) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2073,1))
head(results.tibs)
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.500.harvest.sa <- results.tibs %>% 
  rowwise() %>% # sums across all iterations/rows
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(mean = mean(c_across(1:100)))

# view to make sure everything is working
ggplot(data=sums.500.harvest.sa,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Sub-Adult Population Size")+
  theme_minimal()


# view SR curve
sums.500.harvest.sa <- sums.500.harvest.sa%>%dplyr::select(yrs,mean,low,high)
head(sums.500.harvest.sa)
wt.last<-weight.yr$wt[48]
sums.500.harvest2 <- sums.500.harvest%>%dplyr::select(yrs,mean,low,high)%>%
  mutate(yrs = yrs-4)%>%
  mutate(biomass=wt.last*mean*0.5)%>%
  rename(mean.adult = mean)%>%
  rename(low.adult = low)%>%
  rename(high.adults = high)

sr<-left_join(sums.500.harvest.sa,sums.500.harvest2,by=c("yrs"))
head(sr)

seq.biomass<-seq(0,13000,1000)
a<-exp(full.mod2$BUGSoutput$mean$alpha_sr)
b<-full.mod2$BUGSoutput$mean$beta_sr
model_recruitment1 = exp(log(seq.biomass) + log(a) - b*seq.biomass)
df<-data.frame(seq.biomass,model_recruitment1)
head(df)
  
ggplot(data=df, aes(x=seq.biomass,y=model_recruitment1))+geom_line()+
  geom_point(data=sr, aes(x=biomass,y=mean))
  
### lambda
results.tibs<-t(results.lambdaA) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2025,2077,1))
head(results.tibs)
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

lambda.sims <- results.tibs %>% 
  rowwise() %>% # sums across all iterations/rows
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(mean = mean(c_across(1:100)))

ggplot(data=lambda.sims, aes(x=yrs,y=mean))+geom_line()




#############################################################################################
# 300 harvest, mean kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(300,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(0,54), # mean
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.300.harvest <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))

ggplot(data=sums.300.harvest,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()



#############################################################################################
# 0 harvest, mean kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(0,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(0,54), # mean
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.0.harvest <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))

ggplot(data=sums.0.harvest,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()+
  geom_hline(yintercept = 2982)

yrs<-seq(2024,2077,1)
# combining first 3 simulations to plot together
sim.res1<-tibble(yrs,low=sums.500.harvest$low,mean=sums.500.harvest$mean,high=sums.500.harvest$high,
                 twentyfive=sums.500.harvest$twentyfive,median=sums.500.harvest$median,seventyfive=sums.500.harvest$seventyfive,
                 Harvest="500",Group="Normal Flow",Kokanee="Mean",CaRegs="Catch Release") 

sim.res2<-tibble(yrs,low=sums.300.harvest$low,mean=sums.300.harvest$mean,high=sums.300.harvest$high,
                 twentyfive=sums.300.harvest$twentyfive,median=sums.300.harvest$median,seventyfive=sums.300.harvest$seventyfive,
                 Harvest="300",Group="Normal Flow",Kokanee="Mean",CaRegs="Catch Release") 
sim.res3<-tibble(yrs,low=sums.0.harvest$low,mean=sums.0.harvest$mean,high=sums.0.harvest$high,
                 twentyfive=sums.0.harvest$twentyfive,median=sums.0.harvest$median,seventyfive=sums.0.harvest$seventyfive,
                 Harvest="0",Group="Normal Flow",Kokanee="Mean",CaRegs="Catch Release") 
#sim.res<-full_join(sim.res,sim.res1)
#sim.res<-full_join(sim.res,sim.res2)
sim.res<-full_join(sim.res1,sim.res2)
sim.res<-full_join(sim.res,sim.res3)

cbPalette <- c("#009E73", "#56B4E9", "#D55E00") # setting color palette for ggplot

ggplot(data=sim.res, aes(x=yrs,y=mean,group=Harvest,fill=Harvest,color=Harvest))+ #facet_wrap(~Kokanee)+
  geom_ribbon(aes(ymin=low,ymax=high,group=Harvest,fill=Harvest,color=Harvest),alpha=0.2,color=NA)+
  geom_line(size=1.5)+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()




#################################################################################
# 500 harvest, max kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(500,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(max(predict$running_biomass_scale),54), # max kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.500.harvest.maxkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))


sums.500.harvest.maxkokanee

ggplot(data=sums.500.harvest.maxkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()



#############################################################################################
# 300 harvest, max kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(300,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(max(predict$running_biomass_scale),54), # max kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) 
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.300.harvest.maxkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))

ggplot(data=sums.300.harvest.maxkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()



#############################################################################################
# 0 harvest, max kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(0,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(max(predict$running_biomass_scale),54), # max kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) 
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.0.harvest.maxkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))

ggplot(data=sums.0.harvest.maxkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()

# combining next 3 simulations into a dataframe
sim.res4<-tibble(yrs,low=sums.500.harvest.maxkokanee$low,mean=sums.500.harvest.maxkokanee$mean,high=sums.500.harvest.maxkokanee$high,
                 twentyfive=sums.500.harvest.maxkokanee$twentyfive,median=sums.500.harvest.maxkokanee$median,seventyfive=sums.500.harvest.maxkokanee$seventyfive,
                 Harvest="500",Group="Normal Flow",Kokanee="Max",CaRegs="Catch Release") 

sim.res5<-tibble(yrs,low=sums.300.harvest.maxkokanee$low,mean=sums.300.harvest.maxkokanee$mean,high=sums.300.harvest.maxkokanee$high,
                 twentyfive=sums.300.harvest.maxkokanee$twentyfive,median=sums.300.harvest.maxkokanee$median,seventyfive=sums.300.harvest.maxkokanee$seventyfive,
                 Harvest="300",Group="Normal Flow",Kokanee="Max",CaRegs="Catch Release") 
sim.res6<-tibble(yrs,low=sums.0.harvest.maxkokanee$low,mean=sums.0.harvest.maxkokanee$mean,high=sums.0.harvest.maxkokanee$high,
                 twentyfive=sums.0.harvest.maxkokanee$twentyfive,median=sums.0.harvest.maxkokanee$median,seventyfive=sums.0.harvest.maxkokanee$seventyfive,
                 Harvest="0",Group="Normal Flow",Kokanee="Max",CaRegs="Catch Release") 
sim.res<-full_join(sim.res,sim.res4)
sim.res<-full_join(sim.res,sim.res5)
sim.res<-full_join(sim.res,sim.res6)

cbPalette <- c("#009E73", "#56B4E9", "#D55E00","purple","black","blue","grey") 

ggplot(data=sim.res, aes(x=yrs,y=mean,group=Harvest,fill=Harvest,color=Harvest))+facet_wrap(~Kokanee)+
  geom_ribbon(aes(ymin=low,ymax=high,group=Harvest,fill=Harvest,color=Harvest),alpha=0.2,color=NA)+
  facet_wrap(~Kokanee)+
  geom_line(size=1.5)+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()





#################################################################################
# 500 harvest, min kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(500,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(min(predict$running_biomass_scale),54), # min kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) 
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.500.harvest.minkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))

ggplot(data=sums.500.harvest.minkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()



#############################################################################################
# 300 harvest, min kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(300,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(min(predict$running_biomass_scale),54), # min kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.300.harvest.minkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))

ggplot(data=sums.300.harvest.minkokanee,aes(x=yrs,y=median))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()



#############################################################################################
# 0 harvest, min kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(0,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(min(predict$running_biomass_scale),54), # min kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.0.harvest.minkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))


ggplot(data=sums.0.harvest.minkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()

# combining simulation results
sim.res7<-tibble(yrs,low=sums.500.harvest.minkokanee$low,mean=sums.500.harvest.minkokanee$mean,high=sums.500.harvest.minkokanee$high,
                 twentyfive=sums.500.harvest.minkokanee$twentyfive,median=sums.500.harvest.minkokanee$median,seventyfive=sums.500.harvest.minkokanee$seventyfive,
                 Harvest="500",Group="Normal Flow",Kokanee="Min",CaRegs="Catch Release") 
sim.res8<-tibble(yrs,low=sums.300.harvest.minkokanee$low,mean=sums.300.harvest.minkokanee$mean,high=sums.300.harvest.minkokanee$high,
                 twentyfive=sums.300.harvest.minkokanee$twentyfive,median=sums.300.harvest.minkokanee$median,seventyfive=sums.300.harvest.minkokanee$seventyfive,
                 Harvest="300",Group="Normal Flow",Kokanee="Min",CaRegs="Catch Release") 
sim.res9<-tibble(yrs,low=sums.0.harvest.minkokanee$low,mean=sums.0.harvest.minkokanee$mean,high=sums.0.harvest.minkokanee$high,
                 twentyfive=sums.0.harvest.minkokanee$twentyfive,median=sums.0.harvest.minkokanee$median,seventyfive=sums.0.harvest.minkokanee$seventyfive,
                 Harvest="0",Group="Normal Flow",Kokanee="Min",CaRegs="Catch Release") 
sim.res<-full_join(sim.res,sim.res7)
sim.res<-full_join(sim.res,sim.res8)
sim.res<-full_join(sim.res,sim.res9)

#cbPalette <- c( "#56B4E9","#999999" ,"#D55E00") #"#009E73",

# change order of factors
sim.res$Kokanee <- factor(sim.res$Kokanee, levels=c("Min", "Mean", "Max"))
#sim.res$Harvest <- as.numeric(sim.res$Harvest)

ggplot(data=sim.res, aes(x=yrs,y=mean,group=Harvest,fill=Harvest,color=Harvest))+facet_wrap(~Kokanee)+
  geom_ribbon(aes(ymin=low,ymax=high,group=Harvest,fill=Harvest,color=Harvest),alpha=0.4,color=NA)+
  facet_wrap(~Kokanee, labeller = labeller(Kokanee = 
                                             c("Min" = "Min Kokanee",
                                               "Mean" = "Mean Kokanee",
                                               "Max" = "Max Kokanee")))+
  geom_line(size=1)+
  #scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()+
  theme(
    #plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    #strip.background = element_blank(),
    #strip.text.x = element_text(size=0),
    #strip.text.y = element_blank(),
    legend.position = "bottom")


ggplot(data=sim.res, aes(x=yrs,y=mean,group=Harvest,fill=Harvest,color=Harvest))+
  geom_ribbon(aes(ymin=low,ymax=high,group=Harvest,fill=Harvest,color=Harvest),alpha=0.7)+ #,color=NA
  facet_wrap(~Kokanee, labeller = labeller(Kokanee = 
                                             c("Min" = "Min Kokanee",
                                               "Mean" = "Mean Kokanee",
                                               "Max" = "Max Kokanee")))+
  #geom_line(size=1.5)+
  #scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()+
  theme(
    #plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    #strip.background = element_blank(),
    #strip.text.x = element_text(size=0),
    #strip.text.y = element_blank(),
    legend.position = "bottom")
#ggsave("Results/summary/simulation.harvest.kokanee.png",width = 6, height =4)







# another simulation 
# co-variates = 400 harvest, mean kokanee availability, canadian harvest not allowed (ie catch and release only)
results <- as.data.frame(c())
iterations<-100 

### takes some time to run all of these simulations 
for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 # DATA
                 bt.harv = rep(400,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(0,54),# mean
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54), # mean flow
                 ca.reg.changes = rep(0,54), # catch-and-release in canada
                 mean.fish.per.redd = 0.5,
                 
                 # PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant 
  
  # parameters monitored 
  params <- c("Nadults","yr1subs","lambdaA") # add lambda to get median lambda over last 20 yrs to see if population is declining 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
} 

# Results
results.sa<-results[,55:104] # new sub-adults - starts at year 5 (1984)
results.lambdaA<-results[,105:157] # lambda
results<-results[,1:54] # adults

results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
head(results.tibs)
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.400.harvest <- results.tibs %>% 
  rowwise() %>% # sums across all iterations/rows
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>%
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))%>%
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))

# view to make sure everything is working
ggplot(data=sums.400.harvest,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()

# sub-adults 
results.tibs<-t(results.sa) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2073,1))
head(results.tibs)
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.400.harvest.sa <- results.tibs %>% 
  rowwise() %>% # sums across all iterations/rows
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(mean = mean(c_across(1:100)))

# view to make sure everything is working
ggplot(data=sums.400.harvest.sa,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Sub-Adult Population Size")+
  theme_minimal()


# view SR curve
sums.400.harvest.sa <- sums.400.harvest.sa%>%dplyr::select(yrs,mean,low,high)
head(sums.400.harvest.sa)
wt.last<-weight.yr$wt[48]
sums.400.harvest2 <- sums.400.harvest%>%dplyr::select(yrs,mean,low,high)%>%
  mutate(yrs = yrs-4)%>%
  mutate(biomass=wt.last*mean*0.5)%>%
  rename(mean.adult = mean)%>%
  rename(low.adult = low)%>%
  rename(high.adults = high)

sr<-left_join(sums.400.harvest.sa,sums.400.harvest2,by=c("yrs"))
head(sr)

seq.biomass<-seq(0,13000,1000)
a<-exp(full.mod2$BUGSoutput$mean$alpha_sr)
b<-full.mod2$BUGSoutput$mean$beta_sr
model_recruitment1 = exp(log(seq.biomass) + log(a) - b*seq.biomass)
df<-data.frame(seq.biomass,model_recruitment1)
head(df)

ggplot(data=df, aes(x=seq.biomass,y=model_recruitment1))+geom_line()+
  geom_point(data=sr, aes(x=biomass,y=mean))

### lambda
results.tibs<-t(results.lambdaA) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2025,2077,1))
head(results.tibs)
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

lambda.sims <- results.tibs %>% 
  rowwise() %>% # sums across all iterations/rows
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(mean = mean(c_across(1:100)))

ggplot(data=lambda.sims, aes(x=yrs,y=mean))+geom_line()





#############################################################################################
# 200 harvest, mean kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(200,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(0,54),# mean
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.200.harvest <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))%>% # 
  mutate(high = quantile(c_across(1:100),probs=0.975))

ggplot(data=sums.200.harvest,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()



#############################################################################################
# 100 harvest, mean kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(100,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(0,54),# mean
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.100.harvest <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))


ggplot(data=sums.100.harvest,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()




#############################################################################################
# 600 harvest, mean kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(600,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(0,54),# mean
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.600.harvest <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))


ggplot(data=sums.600.harvest,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()








# combining simulations to plot together
yrs=seq(2024,2077,1)
sim.res28<-tibble(yrs,low=sums.400.harvest$low,mean=sums.400.harvest$mean,high=sums.400.harvest$high,
                  twentyfive=sums.400.harvest$twentyfive,median=sums.400.harvest$median,seventyfive=sums.400.harvest$seventyfive,
                  Harvest="400",Group="Normal Flow",Kokanee="Mean",CaRegs="Catch Release") 

sim.res29<-tibble(yrs,low=sums.200.harvest$low,mean=sums.200.harvest$mean,high=sums.200.harvest$high,
                  twentyfive=sums.200.harvest$twentyfive,median=sums.200.harvest$median,seventyfive=sums.200.harvest$seventyfive,
                  Harvest="200",Group="Normal Flow",Kokanee="Mean",CaRegs="Catch Release") 

sim.res30<-tibble(yrs,low=sums.100.harvest$low,mean=sums.100.harvest$mean,high=sums.100.harvest$high,
                  twentyfive=sums.100.harvest$twentyfive,median=sums.100.harvest$median,seventyfive=sums.100.harvest$seventyfive,
                  Harvest="100",Group="Normal Flow",Kokanee="Mean",CaRegs="Catch Release") 

sim.res31<-tibble(yrs,low=sums.600.harvest$low,mean=sums.600.harvest$mean,high=sums.600.harvest$high,
                  twentyfive=sums.600.harvest$twentyfive,median=sums.600.harvest$median,seventyfive=sums.600.harvest$seventyfive,
                  Harvest="600",Group="Normal Flow",Kokanee="Mean",CaRegs="Catch Release") 

#sim.res<-full_join(sim.res28,sim.res29)
sim.res<-full_join(sim.res,sim.res28)
sim.res<-full_join(sim.res,sim.res29)
sim.res<-full_join(sim.res,sim.res30)
sim.res<-full_join(sim.res,sim.res31)

#cbPalette <- c("#009E73", "#56B4E9", "#D55E00","grey") # setting color palette for ggplot

ggplot(data=sim.res, aes(x=yrs,y=median,group=Harvest,fill=Harvest,color=Harvest))+facet_wrap(~Kokanee)+
  geom_ribbon(aes(ymin=twentyfive,ymax=seventyfive,group=Harvest,fill=Harvest,color=Harvest),alpha=0.2,color=NA)+
  geom_line(size=1.5)+
  #scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()






# another simulation 
# co-variates = 400 harvest, max kokanee availability, canadian harvest not allowed (ie catch and release only)
results <- as.data.frame(c())
iterations<-100 

### takes some time to run all of these simulations 
for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 # DATA
                 bt.harv = rep(400,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(max(predict$running_biomass_scale),54), # max kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54), # mean flow
                 ca.reg.changes = rep(0,54), # catch-and-release in canada
                 mean.fish.per.redd = 0.5,
                 
                 # PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant 
  
  # parameters monitored 
  params <- c("Nadults","yr1subs","lambdaA") # add lambda to get median lambda over last 20 yrs to see if population is declining 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
} 

# Results
#results.sa<-results[,55:104] # new sub-adults - starts at year 5 (1984)
#results.lambdaA<-results[,105:157] # lambda
results<-results[,1:54] # adults

results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
head(results.tibs)
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.400.harvest.maxkokanee <- results.tibs %>% 
  rowwise() %>% # sums across all iterations/rows
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>%
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))%>%
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))

# view to make sure everything is working
ggplot(data=sums.400.harvest.maxkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()




#############################################################################################
# 200 harvest, max kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(200,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(max(predict$running_biomass_scale),54), # max kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.200.harvest.maxkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))%>% # 
  mutate(high = quantile(c_across(1:100),probs=0.975))

ggplot(data=sums.200.harvest.maxkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()



#############################################################################################
# 100 harvest, max kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(100,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(max(predict$running_biomass_scale),54), # max kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.100.harvest.maxkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))


ggplot(data=sums.100.harvest.maxkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()




#############################################################################################
# 600 harvest, max kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(600,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(max(predict$running_biomass_scale),54), # max kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.600.harvest.maxkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))


ggplot(data=sums.600.harvest.maxkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()





# combining simulations to plot together
#yrs=seq(2024,2077,1)
sim.res32<-tibble(yrs,low=sums.400.harvest.maxkokanee$low,mean=sums.400.harvest.maxkokanee$mean,high=sums.400.harvest.maxkokanee$high,
                  twentyfive=sums.400.harvest.maxkokanee$twentyfive,median=sums.400.harvest.maxkokanee$median,seventyfive=sums.400.harvest.maxkokanee$seventyfive,
                  Harvest="400",Group="Normal Flow",Kokanee="Max",CaRegs="Catch Release") 

sim.res33<-tibble(yrs,low=sums.200.harvest.maxkokanee$low,mean=sums.200.harvest.maxkokanee$mean,high=sums.200.harvest.maxkokanee$high,
                  twentyfive=sums.200.harvest.maxkokanee$twentyfive,median=sums.200.harvest.maxkokanee$median,seventyfive=sums.200.harvest.maxkokanee$seventyfive,
                  Harvest="200",Group="Normal Flow",Kokanee="Max",CaRegs="Catch Release") 

sim.res34<-tibble(yrs,low=sums.100.harvest.maxkokanee$low,mean=sums.100.harvest.maxkokanee$mean,high=sums.100.harvest.maxkokanee$high,
                  twentyfive=sums.100.harvest.maxkokanee$twentyfive,median=sums.100.harvest.maxkokanee$median,seventyfive=sums.100.harvest.maxkokanee$seventyfive,
                  Harvest="100",Group="Normal Flow",Kokanee="Max",CaRegs="Catch Release") 

sim.res35<-tibble(yrs,low=sums.600.harvest.maxkokanee$low,mean=sums.600.harvest.maxkokanee$mean,high=sums.600.harvest.maxkokanee$high,
                  twentyfive=sums.600.harvest.maxkokanee$twentyfive,median=sums.600.harvest.maxkokanee$median,seventyfive=sums.600.harvest.maxkokanee$seventyfive,
                  Harvest="600",Group="Normal Flow",Kokanee="Max",CaRegs="Catch Release") 

sim.res<-full_join(sim.res,sim.res32)
sim.res<-full_join(sim.res,sim.res33)
sim.res<-full_join(sim.res,sim.res34)
sim.res<-full_join(sim.res,sim.res35)

#cbPalette <- c("#009E73", "#56B4E9", "#D55E00","grey") # setting color palette for ggplot

ggplot(data=sim.res, aes(x=yrs,y=median,group=Harvest,fill=Harvest,color=Harvest))+
  geom_ribbon(aes(ymin=twentyfive,ymax=seventyfive,group=Harvest,fill=Harvest,color=Harvest),alpha=0.2,color=NA)+facet_wrap(~Kokanee)+
  geom_line(size=1.5)+
  #scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()





# another simulation 
# co-variates = 400 harvest, min kokanee availability, canadian harvest not allowed (ie catch and release only)
results <- as.data.frame(c())
iterations<-100 

### takes some time to run all of these simulations 
for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 # DATA
                 bt.harv = rep(400,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(min(predict$running_biomass_scale),54), # min kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54), # mean flow
                 ca.reg.changes = rep(0,54), # catch-and-release in canada
                 mean.fish.per.redd = 0.5,
                 
                 # PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant 
  
  # parameters monitored 
  params <- c("Nadults","yr1subs","lambdaA") # add lambda to get median lambda over last 20 yrs to see if population is declining 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
} 

# Results
results<-results[,1:54] # adults

results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
head(results.tibs)
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.400.harvest.minkokanee <- results.tibs %>% 
  rowwise() %>% # sums across all iterations/rows
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>%
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))%>%
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))

# view to make sure everything is working
ggplot(data=sums.400.harvest.minkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()




#############################################################################################
# 200 harvest, min kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(200,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(min(predict$running_biomass_scale),54), # min kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)

sums.200.harvest.minkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = mean(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))%>% # 
  mutate(high = quantile(c_across(1:100),probs=0.975))

ggplot(data=sums.200.harvest.minkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()



#############################################################################################
# 100 harvest, min kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(100,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(min(predict$running_biomass_scale),54), # min kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.100.harvest.minkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))


ggplot(data=sums.100.harvest.minkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()




#############################################################################################
# 600 harvest, min kokanee availability, canadian harvest not allowed
results <- as.data.frame(c())
iterations<-100

for(i in 1:iterations){
  sim.data<-list(nYears= 54,
                 
                 #DATA
                 bt.harv = rep(600,53), #  fish harvested each year
                 weight = rep(weight.yr$wt[48],54), # current bt weight
                 #kok = rep(predict$running_biomass_scale[44],54), # current kok availability 
                 kok = rep(min(predict$running_biomass_scale),54), # min kokanee
                 low.flow.recruitment = rep(mean(bull.river.low$scale_summerLowflow[6:45]),54),
                 ca.reg.changes = rep(0,54), 
                 mean.fish.per.redd = 0.5,
                 
                 #PARAMETER ESTIMATES - NEED TO SAMPLE
                 alpha_sr = sample(full.mod2$BUGSoutput$sims.list$alpha_sr,size=1),
                 beta_sr = sample(full.mod2$BUGSoutput$sims.list$beta_sr,size=1),
                 
                 tauR_process=(sample(full.mod2$BUGSoutput$sims.list$sigmaR_process,size=1)^-2),
                 tau.S.a = (sample(full.mod2$BUGSoutput$sims.list$sd.S.a,size=1)^-2),
                 
                 S.sa = sample(full.mod2$BUGSoutput$sims.list$S.sa,size=1),
                 a.int = sample(full.mod2$BUGSoutput$sims.list$a.int,size=1),
                 beta.kok = sample(full.mod2$BUGSoutput$sims.list$beta.kok,size=1),
                 beta.dens = sample(full.mod2$BUGSoutput$sims.list$beta.dens,size=1),
                 beta.low.flow = sample(full.mod2$BUGSoutput$sims.list$beta.low.flow,size=1),
                 beta.reg.change = sample(full.mod2$BUGSoutput$sims.list$beta.reg.change,size=1),
                 transition.prob = full.mod2$BUGSoutput$mean$transition.prob) # holding constant
  
  # parameters monitored 
  params <- c("Nadults") 
  
  # run jags
  out <- run.jags(txtstring2, data = sim.data, monitor=params,
                  sample=1, n.chains=1,, summarise=FALSE)
  Simulated <- coda::as.mcmc(out)
  results<-rbind(results,Simulated)
}

# Results
results.tibs<-t(results) # transpose
results.tibs<-data.frame(results.tibs)
x <- seq(1,iterations,1) # number of iterations
colnames(results.tibs) <- x
results.tibs<-results.tibs%>%
  mutate(yrs=seq(2024,2077,1))
results.tibs <- pmax(results.tibs,0) # turns negative numbers to NA
results.tibs[is.na(results.tibs)] <- 0 # turns na to 0 (can't have negative population size)


sums.600.harvest.minkokanee <- results.tibs %>% 
  rowwise() %>% 
  mutate(mean = mean(c_across(1:100)))%>%
  mutate(median = median(c_across(1:100)))%>%
  mutate(low = quantile(c_across(1:100),probs=0.025))%>% # 95% credible interval
  mutate(high = quantile(c_across(1:100),probs=0.975))%>%
  mutate(twentyfive = quantile(c_across(1:100),probs=0.25))%>% # 
  mutate(seventyfive = quantile(c_across(1:100),probs=0.75))


ggplot(data=sums.600.harvest.minkokanee,aes(x=yrs,y=mean))+geom_line()+
  geom_ribbon(aes(ymin=low,ymax=high),alpha=0.1)+labs(x="Year",y="Adult Population Size")+
  theme_minimal()





# combining simulations to plot together
#yrs=seq(2024,2077,1)
sim.res36<-tibble(yrs,low=sums.400.harvest.minkokanee$low,mean=sums.400.harvest.minkokanee$mean,high=sums.400.harvest.minkokanee$high,
                  twentyfive=sums.400.harvest.minkokanee$twentyfive,median=sums.400.harvest.minkokanee$median,seventyfive=sums.400.harvest.minkokanee$seventyfive,
                  Harvest="400",Group="Normal Flow",Kokanee="Min",CaRegs="Catch Release") 

sim.res37<-tibble(yrs,low=sums.200.harvest.minkokanee$low,mean=sums.200.harvest.minkokanee$mean,high=sums.200.harvest.minkokanee$high,
                  twentyfive=sums.200.harvest.minkokanee$twentyfive,median=sums.200.harvest.minkokanee$median,seventyfive=sums.200.harvest.minkokanee$seventyfive,
                  Harvest="200",Group="Normal Flow",Kokanee="Min",CaRegs="Catch Release") 

sim.res38<-tibble(yrs,low=sums.100.harvest.minkokanee$low,mean=sums.100.harvest.minkokanee$mean,high=sums.100.harvest.minkokanee$high,
                  twentyfive=sums.100.harvest.minkokanee$twentyfive,median=sums.100.harvest.minkokanee$median,seventyfive=sums.100.harvest.minkokanee$seventyfive,
                  Harvest="100",Group="Normal Flow",Kokanee="Min",CaRegs="Catch Release") 

sim.res39<-tibble(yrs,low=sums.600.harvest.minkokanee$low,mean=sums.600.harvest.minkokanee$mean,high=sums.600.harvest.minkokanee$high,
                  twentyfive=sums.600.harvest.minkokanee$twentyfive,median=sums.600.harvest.minkokanee$median,seventyfive=sums.600.harvest.minkokanee$seventyfive,
                  Harvest="600",Group="Normal Flow",Kokanee="Min",CaRegs="Catch Release") 

sim.res<-full_join(sim.res,sim.res36)
sim.res<-full_join(sim.res,sim.res37)
sim.res<-full_join(sim.res,sim.res38)
sim.res<-full_join(sim.res,sim.res39)

#cbPalette <- c("#009E73", "#56B4E9", "#D55E00","grey") # setting color palette for ggplot

ggplot(data=sim.res, aes(x=yrs,y=median,group=Harvest,fill=Harvest,color=Harvest))+
  geom_ribbon(aes(ymin=twentyfive,ymax=seventyfive,group=Harvest,fill=Harvest,color=Harvest),alpha=0.2,color=NA)+facet_wrap(~Kokanee)+
  geom_line(size=1.5)+
  #scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()










### final plots, start here ###

#### start here #####
# save all simulations
#saveRDS(sim.res, file="R/sim.res") 
# read in datafile
sim.res <-readRDS(file="R/sim.res")



# mean from 2050 - 2077
mean(sim.res$mean[81:108]) # harvest =300, mean kokanee
mean(sim.res$low[81:108])
mean(sim.res$high[81:108])

mean(sim.res$mean[405:432]) # harvest =300, min kokanee
mean(sim.res$low[405:432])
mean(sim.res$high[405:432])



#cbPalette <- c( "#56B4E9","#999999" ,"#D55E00") #"#009E73",

# change order of factors
sim.res$Kokanee <- factor(sim.res$Kokanee, levels=c("Min", "Mean", "Max"))
sim.res$Harvest <- as.numeric(sim.res$Harvest)

sim.res%>%filter(Group=="Normal Flow") %>%

  ggplot(aes(x=yrs,y=mean,group=Harvest,fill=Harvest,color=Harvest))+
  #geom_ribbon(aes(ymin=low,ymax=high,group=Harvest,fill=Harvest,color=Harvest),alpha=0.4)+ # ,color=NA
  #facet_wrap(~Kokanee, labeller = labeller(Kokanee = 
  #                                          c("Min" = "Min Kokanee",
  #                                           "Current" = "Current Kokanee",
  #                                          "Max" = "Max Kokanee")))+
  facet_grid(CaRegs~Kokanee)+
  geom_line(size=1)+
  #scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()+
  theme(
    #plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
    #strip.background = element_blank(),
    #strip.text.x = element_text(size=0),
    #strip.text.y = element_blank(),
    legend.position = "bottom")

l<-sim.res%>%filter(Group=="Normal Flow") %>%
  ggplot( aes(x=yrs,y=mean,group=Harvest,fill=Harvest,color=Harvest))+
  geom_ribbon(aes(ymin=low,ymax=high,group=Harvest,fill=Harvest,color=Harvest),alpha=0.7)+ #,color=NA
  #facet_wrap(~Kokanee, labeller = labeller(Kokanee = 
  #                                          c("Min" = "Min Kokanee",
  #                                           "Current" = "Current Kokanee",
  #                                          "Max" = "Max Kokanee")))+
  facet_grid(CaRegs~Kokanee,
             labeller = labeller(Kokanee = c("Min" = "Min Kokanee",
                                             "Mean" = "Mean Kokanee",
                                             "Max" = "Max Kokanee")))+
  #geom_line(size=1.5)+
  #scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()+
  theme(
    plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
    strip.background = element_blank(),
    #strip.text.x = element_text(size=0),
    #strip.text.y = element_blank(),
    legend.title = element_text( size=10), legend.text=element_text(size=10),
    legend.key.size = unit(0.3, "cm"),legend.box.margin=margin(0,0,0,0),
    legend.margin=margin(0,0,0,0),
    legend.position = "bottom")
l
#ggsave("Results/summary/simulation.harvest.kokanee.CanadianRegs.png",width = 6, height =4)






p<-sim.res%>%filter(Group=="Normal Flow") %>%
  ggplot( aes(x=yrs,y=mean,group=Harvest,fill=Harvest,color=Harvest))+
  geom_ribbon(aes(ymin=low,ymax=high,group=Harvest,fill=Harvest,color=Harvest),alpha=0.7)+ #,color=NA
  facet_grid(CaRegs~Kokanee)+
  #facet_grid(CaRegs~Kokanee,
  #          labeller = labeller(Kokanee = c("Min" = "Min Kokanee",
  #                                         "Current" = "Current Kokanee",
  #                                        "Max" = "Max Kokanee")))+
  #geom_line(size=1.5)+
  #scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(2030,2050,2070),expand=c(0,0))+
  labs(x="Year",y="Adult Population Size")+theme_bw()+
  theme(
    plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
    strip.background = element_blank(),
    #strip.text.x = element_text(size=0),
    #strip.text.y = element_blank(),
    legend.title = element_text( size=10), legend.text=element_text(size=10),
    legend.key.size = unit(0.3, "cm"),legend.box.margin=margin(0,0,0,0),
    legend.margin=margin(0,0,0,0),
    legend.position = "bottom")+
    scale_color_continuous(breaks=c(0,300,600))+scale_fill_continuous(breaks=c(0,300,600))
#p <- p + guides(shape = guide_legend(override.aes = list(size = 0.2)))
tag_facet(p, open = "", close = "")
#ggsave("Results/summary/simulation.harvest.kokanee.CanadianRegs2.png",width = 6, height =4)





#head(sim.res)
# biomass=wt*Nadult*0.5

### MSY 
spawning.biomass.msy<-4347.2520 # lbs 
current.wt<-weight.yr$wt[48] # current adult bull trout weight = 5.3 lbs

# how many adult fish at maximum sustainable yield?? 
msy<-(spawning.biomass.msy*2)/current.wt # MSY 
msy
# 1636 fish

# current adult bull trout abundance
current.Na<-full.mod2$BUGSoutput$mean$Nadults[44] # 2908

sim.res<-sim.res%>%
  mutate(percentChangeinN = (median-current.Na)/current.Na) # using median instead of mean
#mutate(percentChangeinN.low = (low-current.Na)/current.Na)%>%
#mutate(percentChangeinN.high = (high-current.Na)/current.Na)
head(sim.res)

sim.res.sums<-sim.res%>%filter(yrs>2030)%>%group_by(Harvest,Group,Kokanee,CaRegs)%>%
  summarise(N = mean(mean),
            Change=mean(percentChangeinN),
            Nlow = mean(low),
            N25 =  mean(twentyfive),
            Nmedian = mean(median),
            N75 =  mean(seventyfive),
            Nhigh = mean(high))
head(sim.res.sums)


sim.res.sums %>% filter(CaRegs=="Catch Release")%>%filter(Kokanee=="Mean")%>%filter(Group=="Normal Flow")%>%
  ggplot(aes(x=Harvest,y=Change))+geom_point()
#geom_pointrange(aes(x=Harvest,ymin=Changelow,ymax=Changehigh))
# y limits from -1 to 0.5

sim.res.sums %>% filter(CaRegs=="Catch Release")%>%filter(Kokanee=="Mean")%>%filter(Group=="Normal Flow")%>%
  ggplot(aes(x=Harvest,y=Nmedian))+geom_point(size=3)+
  geom_pointrange(aes(x=Harvest,ymin=Nlow,ymax=Nhigh))+
  geom_pointrange(aes(x=Harvest,ymin=N25,ymax=N75),linewidth=2)
# y limits of 0 to 4500

sim.res.sums %>% filter(CaRegs=="Catch Release")%>%filter(Kokanee=="Min")%>%filter(Group=="Normal Flow")%>%
  ggplot(aes(x=Harvest,y=Nmedian))+geom_point(size=3)+
  geom_pointrange(aes(x=Harvest,ymin=Nlow,ymax=Nhigh))+
  geom_pointrange(aes(x=Harvest,ymin=N25,ymax=N75),linewidth=2)



sim.res.sums %>% filter(Group=="Normal Flow")%>% # filter(CaRegs=="Catch Release")%>%
  ggplot(aes(x=Harvest))+
  geom_hline(yintercept = msy,linetype=2)+
  #annotate("text", x = 1.45, y = msy-130, label = "Max Recruit")+
  geom_pointrange(aes(x=Harvest,y=Nmedian,ymin=Nlow,ymax=Nhigh,color=Kokanee),linewidth=0.25)+
  geom_pointrange(aes(x=Harvest,y=Nmedian,ymin=N25,ymax=N75,color=Kokanee),linewidth=1)+
  geom_point(aes(y=N,color=Kokanee,group=Kokanee))+
  #facet_wrap(~Kokanee)+
  facet_grid(CaRegs~Kokanee,
                       labeller = labeller(Kokanee = c("Min" = "Min Kokanee",
                                                      "Mean" = "Mean Kokanee",
                                                     "Max" = "Max Kokanee")))+
  geom_point(aes(y=(Change+1)*2982.316,color=Kokanee))+
  scale_y_continuous(sec.axis = sec_axis(~./2982.316-1, name="Percent Change")) +
  scale_color_grey(start=0.2,end=0.8)+
  theme_minimal()+
  #theme(legend.position = c(0.2,0.8),strip.text.x = element_text(size=0),
  #                     plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text.x = element_blank(),
        legend.position="bottom")+
  labs(y="Adult Abundance")


### final final plot ####
sim.res.sums %>% filter(Group=="Normal Flow")%>%filter(CaRegs=="Catch Release")%>%
  ggplot(aes(x=Harvest))+
  geom_hline(yintercept = msy,linetype=2)+
  #annotate("text", x = 1.45, y = msy-130, label = "Max Recruit")+
  geom_pointrange(aes(x=Harvest,y=Nmedian,ymin=Nlow,ymax=Nhigh,color=Kokanee),linewidth=0.5)+
  geom_pointrange(aes(x=Harvest,y=Nmedian,ymin=N25,ymax=N75,color=Kokanee),linewidth=1.5)+
  geom_point(aes(y=Nmedian,color=Kokanee,group=Kokanee),size=2.5)+
  facet_wrap(~Kokanee)+
  #geom_point(aes(y=(Change+1)*2982.316,color=Kokanee))+ # y = N/(change+1)
  scale_y_continuous(sec.axis = sec_axis(~./2982.316-1, name="Percent Change")) +
  scale_color_grey(start=0.2,end=0.7)+
  theme_minimal()+
  #theme(legend.position = c(0.2,0.8),strip.text.x = element_text(size=0),
  #                     plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'),
        strip.text.x = element_blank(),
        legend.position="bottom")+
  labs(y="Adult Abundance")
ggsave("Results/summary/PubFigs/Figure 9.png",width = 15, height =8.5, units="cm", dpi=500)



# percent decrease from current to 300 harvest
((2789-2908)/2789)*100 # 4% decrease - mean kokanee
((1941-2908)/1941)*100 # 50% decrease - minimum kokanee
((3363-2908)/3363)*100 # 14% increase - maximum kokanee


### get manuscript numbers
sim.res.sums%>%filter(Harvest==300)

sim.res.sums%>%filter(Harvest==400)




