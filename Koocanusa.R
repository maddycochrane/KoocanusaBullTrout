


##### Bull Trout and Kokanee Data from Koocanusa System 1975 - 2022 ######
### Visually exploring data and trends 

# Maddy Cochrane
# 7/7/23

## 
setwd("C:\\Users\\maddy\\OneDrive - Montana State University\\Koocanusa")

library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)

# read in bull trout catch per unit effort data
bull.dat<-read.csv("Data/BullTroutCPUE.csv",
                  header=TRUE) 

str(bull.dat) # view structure 

bull.dat$Date <- # change date from character to posixct
  as.POSIXct(bull.dat$Date, tz = "", format="%m/%d/%Y") 

bull.dat<-bull.dat%>%
  rename(BullCPUE = CPUE) # rename column
head(bull.dat)

# quick view of bull trout data
ggplot(data=bull.dat)+
  geom_point(aes(x=Year,y=BullCPUE))+
  geom_line(aes(x=Year,y=BullCPUE))

# read in kokanee catch per unit effort data
kok.dat<-read.csv("Data/KokaneeLKCPUE.csv",header=TRUE)
str(kok.dat)

# calculate cpue for 2 sites and combine for total across both
kok.dat1<-kok.dat%>%
  mutate(KokCPUE_Can = KokaneeCount_Canada/NetCount_Canada)%>%
  mutate(KokCPUE_Rex = KokaneeCount_Rexford/NetCount_Rexford)%>%
  mutate(Total_KokCount = KokaneeCount_Canada + KokaneeCount_Rexford)%>%
  mutate(Total_NetCount = NetCount_Canada + NetCount_Rexford)%>%
  mutate(KokCPUE_Total = case_when(is.na(KokCPUE_Rex) ~ NA, # dealing with missing survey years
                                   is.na(KokCPUE_Can) & !is.na(KokCPUE_Rex) ~ KokCPUE_Rex,
                                   !is.na(KokCPUE_Can) & !is.na(KokCPUE_Rex) ~ Total_KokCount / Total_NetCount))

# go from wide to long
long.kok<-kok.dat1%>%select(Year,KokCPUE_Can, KokCPUE_Rex, KokCPUE_Total)
long.kok<-long.kok%>%
  gather(key=Site, value= CPUE, KokCPUE_Can:KokCPUE_Total, factor_key=TRUE)

# quick view of kokanee data
ggplot(data=long.kok)+
  geom_point(aes(x=Year,y=CPUE,color=Site))+
  geom_line(aes(x=Year,y=CPUE,color=Site))



############# 
### combine bull trout and kokanee cpue
bull.dat1<- bull.dat %>%
  rename(CPUE = BullCPUE)%>%
  mutate(Site = "BullCPUE")%>%
  dplyr::select(Year, Site,CPUE)

comb.dat<-full_join(bull.dat1,long.kok) # keeps all observations 

head(comb.dat)
comb.dat<-comb.dat%>%
  mutate(Sp = case_when(Site == "BullCPUE" ~ "Adult Bull",
                        Site != "BullCPUE" ~ "Kokanee"))

# manual color palette
cbp2 <- c("dodgerblue4", "grey65", "grey80", "#000000")


# plot all data together 
ggplot(data=comb.dat)+
  geom_point(data=comb.dat, aes(x=Year,y=CPUE,color=Site,shape=Site),size=2)+
  scale_shape_manual(values=c(15, 1, 2, 16))+
  geom_line(data=comb.dat, aes(x=Year,y=CPUE,color=Site))+
  labs(y="CPUE (ind/net)")+
  scale_colour_manual(values=cbp2)+
  theme_classic()+theme(legend.title = element_blank())

# SCALE Y 
ggplot(data=comb.dat)+
  geom_point(data=comb.dat, aes(x=Year,y=scale(CPUE),color=Site,shape=Site),size=2)+
  scale_shape_manual(values=c(15, 1, 2, 16))+
  geom_line(data=comb.dat, aes(x=Year,y=scale(CPUE),color=Site))+
  labs(y="scaled CPUE")+
  scale_colour_manual(values=cbp2)+
  theme_classic()+theme(legend.title = element_blank())



# 2 plots vertically
p1<-ggplot(data=comb.dat)+
  geom_point(data=comb.dat, aes(x=Year,y=CPUE,color=Site,shape=Site),size=2)+
  scale_shape_manual(values=c(15, 1, 2, 16))+
  geom_line(data=comb.dat, aes(x=Year,y=CPUE,color=Site))+
  labs(y="CPUE (ind/net)")+
  scale_colour_manual(values=cbp2)+
  theme_classic()+theme(legend.title = element_blank())+
  facet_wrap(~Sp,ncol=1, scales = "free")
p1

# save figure
#tiff("Results/bull.kok.cpue.tiff",width = 6, height = 6, units = 'in', res = 300)
#p1
#dev.off()




### bring in juvenile bull trout data 
# read in bull trout catch per unit effort data
juv.bull.dat<-read.csv("Data/JuvenileBullTroutpopulationestimates.csv",
                   header=TRUE) 

str(juv.bull.dat) # view structure 
juv.bull.dat$Fishper100m2<-as.integer(juv.bull.dat$Fishper100m2)

# quick view of juvenile bull trout data
p2<-ggplot(data=juv.bull.dat)+
  geom_point(aes(x=Year,y=PopEst,))+
  geom_line(aes(x=Year,y=PopEst))+
  geom_errorbar(aes(x=Year, ymin=PopEst-X95CI, ymax=PopEst+X95CI), 
                alpha=0.7)+
  labs(y="Juv. Bull Trout Population Estimate")+
  scale_color_grey()+theme_classic()+facet_wrap(~Site)
p2

# save figure
#tiff("Results/juv.bull.pop.tiff",width = 6, height = 6, units = 'in', res = 300)
#p2
#dev.off()


# only site above libby dam
grave<-juv.bull.dat %>% filter(Site == "GraveCreek")%>%
  ggplot()+
  geom_point(aes(x=Year,y=PopEst,))+
  geom_line(aes(x=Year,y=PopEst),alpha=0.2)+
  geom_errorbar(aes(x=Year, ymin=PopEst-X95CI, ymax=PopEst+X95CI), 
                alpha=0.7)+
  labs(y="Juv. Bull Trout Population Estimate")+
  scale_color_grey()+theme_classic()
grave

# save figure
#tiff("Results/graveCreek.juv.bull.pop.tiff",width = 6, height = 6, units = 'in', res = 300)
#grave
#dev.off()




#### read in kokanee spawning enumeration data 
kok.spawn.counts <- read.csv("Data/KoocanusaKokaneeCounts.csv",header=TRUE)
str(kok.spawn.counts)

# go from wide to long 
long.kok.spawn<-kok.spawn.counts%>%
  gather(key=Site, value= Count, SandCreek:AvgPerSurvey, factor_key=TRUE)
head(long.kok.spawn)

# quick view of kokanee spawning data
options(scipen = 999) # don't want scientific notation
ggplot(data=long.kok.spawn)+
  geom_point(aes(x=Year,y=Count,color=Site))+
  geom_line(aes(x=Year,y=Count,color=Site))+
  labs(y="Kokanee Escapement")+
  #scale_color_grey()+
  theme_classic()

p3<-long.kok.spawn %>% filter(Site != "Total")%>%
ggplot()+
  geom_point(aes(x=Year,y=Count))+
  geom_line(aes(x=Year,y=Count),alpha=0.4)+
  labs(y="Kokanee Escapement")+
  scale_color_grey()+theme_classic()+facet_wrap(~Site,nrow=2, scales = "free_y")+
  scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020), labels=c('00', '05','10', '15', '20'))
p3


# save figure
#tiff("Results/kok.spawn.tiff",width = 10, height = 6, units = 'in', res = 300)
#p3
#dev.off()




#### read in bull trout redd data
bt.redd.counts <- read.csv("Data/BullTroutRedds.csv",header=TRUE)
str(bt.redd.counts)

# quick view of data

bt1<-bt.redd.counts %>%
  filter(CoreArea == "Lake Koocanusa")%>% # not including Kootenai River
  filter(Waterbody != "Wigwam River Canada")%>%
  filter(LocalPopulation != "Blackfoot Creek")%>% # data from 17 + has 3 entries
ggplot()+
  geom_point(aes(x=Year,y=Redds,color=LocalPopulation,group=Waterbody))+
  geom_line(aes(x=Year,y=Redds,color=LocalPopulation,group=Waterbody),alpha=0.3)+
  labs(y="Bull Trout Redd Count")+
  #scale_color_grey()+
  theme_classic()+ 
  scale_y_continuous(expand = c(0, 0)) # forcing y axis to start at 0
bt1

# view with facet wrap
redd.fig<-bt.redd.counts %>%
  filter(CoreArea == "Lake Koocanusa")%>% # not including Kootenai River
  #filter(Waterbody != "Wigwam River Canada")%>%
  ggplot()+
  geom_point(aes(x=Year,y=Redds,color=Waterbody))+
  geom_line(aes(x=Year,y=Redds,color=Waterbody),alpha=0.5)+
  labs(y="Bull Trout Redd Count")+
  #scale_color_grey()+
  theme_classic()+ 
  scale_y_continuous(expand = c(0, 0))+ # forcing y axis to start at 0
  facet_wrap(~LocalPopulation, scales = "free_y")
redd.fig

bt2<-bt.redd.counts %>%
  filter(CoreArea == "Lake Koocanusa")%>% # not including Kootenai River
  filter(Waterbody == "Wigwam River Canada")%>%
  ggplot()+
  geom_point(aes(x=Year,y=Redds,color=Waterbody,group=Waterbody))+
  geom_line(aes(x=Year,y=Redds,color=Waterbody,group=Waterbody),alpha=0.3)+
  scale_color_manual(values=c("#619CFF"))+
  labs(y="Bull Trout Redd Count")+
  #scale_color_grey()+
  theme_classic()+ 
  scale_y_continuous(expand = c(0, 0)) # forcing y axis to start at 0
bt2



duplot<-cowplot::plot_grid(bt1 +
                             theme(axis.title.x = element_blank()),
                           
                           bt2 +
                             theme(),
                           
                           nrow = 2, labels = "auto", align = "v")
duplot

# save figure
#tiff("Results/bt.redd.counts.tiff",width = 10, height = 6, units = 'in', res = 300)
#redd.fig
#dev.off()





#########################################################
############ fitting a stock recruitment curve #########

# x = stock level / spawners 
# y = recruitment level / juveniles 

# juvenile data are age 1+ but not migratory adults > 75 mm
# thus we need to lag the redd counts by 2 yrs 


#### grave creek 

# juvenile bull trout data 
head(juv.bull.dat)
grave.juv.bull.dat <- juv.bull.dat %>% filter(Site=="GraveCreek") # grave creek only 
# only the mainstem of grave creek is included

# bt redd data 
head(bt.redd.counts)
grave.bt.redd.dat<- bt.redd.counts %>% filter(LocalPopulation == "Grave Creek")
# includes the mainstem of grave, clarence creek, and blue sky creek
head(grave.bt.redd.dat)

# combining redd counts across grave creek tribs into one total value 
grave.bt.redd.tot <- grave.bt.redd.dat %>% 
  group_by(LocalPopulation,Year)%>%
  summarise(sum(Redds))%>%
  rename(TotalRedds = `sum(Redds)`)%>%
  arrange(Year)%>%
  mutate(lagRedds = lag(TotalRedds, n=2))%>%
  ungroup()%>%
  dplyr::select(Year,lagRedds)
head(grave.bt.redd.tot)

# now just want grave mainstem data
graveonly.bt.redd <- grave.bt.redd.dat%>%
  filter(Waterbody == "Grave Creek")%>%
  arrange(Year)%>%
  mutate(lagRedds_Graveonly = lag(Redds, n=2))%>%
  dplyr::select(Year,lagRedds_Graveonly)
head(graveonly.bt.redd)

# combine grave only to grave and both tribs 
grave.bt.redd.tot<-left_join(grave.bt.redd.tot, graveonly.bt.redd,by="Year")

# join juv and redd data 
head(grave.bt.redd.tot)
head(grave.juv.bull.dat)

grave.sr<-left_join(grave.juv.bull.dat, grave.bt.redd.tot, by="Year")
head(grave.sr)

# need to calculate density for redd counts
grave.sr <- grave.sr %>%
  mutate(reddDensity=lagRedds/37.787397)%>% # value is river length surveyed 
  mutate(reddDensity.graveonly = lagRedds_Graveonly/25.427635)%>%
  mutate(Fishper1000m = 10*sqrt(Fishper100m2)) # want everything to be in km
head(grave.sr)

# plotting stock-recruit function 
grave.sr.plot<-ggplot(data=grave.sr,aes(x=lagRedds, y=log(PopEst/lagRedds )))+
  geom_point()+
  labs(x="Redd Count (2 yr lag)", y="log(Juvenile Population Est. / Redd Count)",
       title="Grave Creek & tributaries")+
  geom_smooth(data=grave.sr, aes(x=lagRedds, y=log(PopEst/lagRedds )),
              method = lm, formula = y ~ x)+
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq","R2","p")), formula = y ~ x)
grave.sr.plot
# save image
#ggsave("Results/stockRecruitFnct/grave.sr.fnct.png",width = 4, height =4)


# excluding tributaries 
grave.only.sr.plot<-ggplot(data=grave.sr,aes(x=lagRedds_Graveonly , 
                                             y=log(PopEst/lagRedds_Graveonly  )))+
  geom_point()+
  labs(x="Redd Count (2 yr lag)", y="",
       title="Grave Creek only")+
  geom_smooth(data=grave.sr, aes(x=lagRedds_Graveonly , y=log(PopEst/lagRedds_Graveonly  )),
              method = lm, formula = y ~ x)+
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq","R2","p")), formula = y ~ x)
grave.only.sr.plot

# save images together
pp<-grid.arrange(grave.sr.plot,grave.only.sr.plot,ncol=2)  

#save_plot("Results/stockRecruitFnct/sr.function.png", pp, ncol = 2)





# plotting stock-recruit function but by DENSITY and not count
grave.sr.plot.bydensity<-ggplot(data=grave.sr,aes(x=reddDensity, 
                                                  y=log(Fishper1000m/reddDensity )))+
  geom_point()+
  labs(y="log(Juvenile Density / Redd Density)",
       title="Grave Creek & tributaries")+
  xlab(bquote('Redd Density '('n km' ^-1))) + 
  geom_smooth(data=grave.sr, aes(x=reddDensity, y=log(Fishper1000m/reddDensity )),
              method = lm, formula = y ~ x)+
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq","R2","p")), formula = y ~ x)
grave.sr.plot.bydensity


# excluding tributaries 
head(grave.sr)
grave.only.sr.bydensity.plot<-ggplot(data=grave.sr,aes(x=reddDensity.graveonly , 
                                             y=log(Fishper1000m/reddDensity.graveonly  )))+
  geom_point()+
  labs(y="",
       title="Grave Creek only")+
  xlab(bquote('Redd Density '('n km' ^-1))) + 
  geom_smooth(data=grave.sr, aes(x=reddDensity.graveonly , 
                                 y=log(Fishper1000m/reddDensity.graveonly  )),
              method = lm, formula = y ~ x)+
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq","R2","p")), formula = y ~ x)
grave.only.sr.bydensity.plot

# save images together
qq<-grid.arrange(grave.sr.plot.bydensity,grave.only.sr.bydensity.plot,ncol=2)  

save_plot("Results/stockRecruitFnct/sr.function.by.density.png", qq, ncol = 2)




# linearized Ricker model 
mod<- lm(log(PopEst/lagRedds) ~ lagRedds, data = grave.sr) 
summary(mod)
# beta = -0.009735 , p < 0.001, F stat = 55.33, df = 21

# grave only model 
mod.graveonly<- lm(log(PopEst/lagRedds_Graveonly) ~ 
                     lagRedds_Graveonly, data = grave.sr) 
summary(mod.graveonly)
# beta = -0.012702 , p < 0.001, F stat = 58.59, df = 21


# view log(R/S) over time 
logRS<-ggplot(data=grave.sr,aes(x=Year, y=log(PopEst/lagRedds )))+
         geom_point()+
  labs(y="log(Juvenile Population Est. / Redd Count)",
       title="Grave and tributaries")

logRS2<-ggplot(data=grave.sr,aes(x=Year, y=log(PopEst/lagRedds_Graveonly )))+
  geom_point()+
  labs(y="log(Juvenile Population Est. / Redd Count)",
       title="Grave only")

logRS_plots <- cowplot::plot_grid(logRS +
                                   theme(),
                                          
                                  logRS2 +
                                       theme(axis.title.y = element_blank()),
                                          
                                          nrow = 1, labels = "auto", align = "v")
logRS_plots

# save image
ggsave("Results/stockRecruitFnct/sr.time.png",width = 8, height =4)




# plot raw R vs S

# predict by count
redd.pred <- tibble(seq(0,250,1))
redd.pred <- redd.pred %>%
  rename(lagRedds = `seq(0, 250, 1)`)%>%
  mutate(lagRedds_Graveonly = lagRedds)

redd.pred$pred<-predict(mod, newdata = redd.pred)
redd.pred$pred.graveonly<-predict(mod.graveonly, newdata = redd.pred)
head(redd.pred)

redd.pred <- redd.pred %>%
  mutate(recruit = lagRedds * (exp(pred)))%>% # transform 
  mutate(recruit.graveonly = lagRedds * (exp(pred.graveonly)))


require(gridExtra)
p1<-ggplot(data=grave.sr)+
  geom_point(aes(x=lagRedds, y=PopEst ))+
  
  labs(y="Juvenile Population Estimate", x="Redd Count (2 yr lag)",
       title="Grave Creek and tributaries")+
  
  geom_line(data=redd.pred, aes(x=lagRedds, y =recruit ), size=1.25)
p1

p2<-ggplot(data=grave.sr)+
  geom_point(aes(x=lagRedds_Graveonly, y=PopEst ),color="blue")+
  
  labs(y="", x="Redd Count (2 yr lag)",
       title="Grave Creek only")+
  
  geom_line(data=redd.pred, aes(x=lagRedds, y =recruit.graveonly ),
            color="blue", linetype= 2, size=1.25) +
  scale_x_continuous(limits=c(0,200))
p3<-grid.arrange(p1,p2,ncol=2)  

#save_plot("Results/stockRecruitFnct/ricker.grave.png", p3, ncol = 2)



# predict by DENSITY 
head(grave.sr)

redd.pred1 <- tibble(seq(0,7,0.1))
redd.pred1 <- redd.pred1 %>%
  rename(reddDensity = `seq(0, 7, 0.1)`)%>%
  mutate(reddDensity.graveonly = reddDensity)
head(redd.pred1)



# linearized Ricker model for density
mod1<- lm(log(Fishper1000m/reddDensity) ~ reddDensity, data = grave.sr) 
summary(mod1)
# beta = -0.3524 , p < 0.001, F stat = 97.43, df = 21

# grave only model 
mod1.graveonly<- lm(log(Fishper1000m/reddDensity.graveonly) ~ 
                      reddDensity.graveonly, data = grave.sr) 
summary(mod1.graveonly)
# beta = -0.30538 , p < 0.001, F stat = 106.1, df = 21


redd.pred1$pred1<-predict(mod1, newdata = redd.pred1)
redd.pred1$pred1.graveonly<-predict(mod1.graveonly, newdata = redd.pred1)
head(redd.pred1)


redd.pred1 <- redd.pred1 %>%
  mutate(recruit = reddDensity  * (exp(pred1)))%>% # transform 
  mutate(recruit.graveonly = reddDensity.graveonly  * (exp(pred1.graveonly)))
head(redd.pred1)


#### add replacement rate lines (2:1, 2.5:1, 3:1)
stock<-tibble(seq(0,7,1))
  stock <- stock %>%
    rename(reddDensity=`seq(0, 7, 1)`)%>%
    mutate(recruit2 = reddDensity*2)%>%
    mutate(recruit2.5 = reddDensity*2.5)%>%
    mutate(recruit3 = reddDensity*3)


head(grave.sr)
q1<-ggplot(data=grave.sr,aes(x=reddDensity , y=Fishper1000m ))+
  geom_point()+
  labs(title="Grave Creek and tributaries")+
  xlab(bquote('Redd Density '('n km' ^-1))) + 
  ylab(bquote('Juvenile Density  '('n km' ^-1))) + 
  geom_line(data=redd.pred1, aes(x=reddDensity, y =recruit ), size=1.25)+
  
  geom_line(data=stock, aes(x=reddDensity, y =recruit2),linetype="dashed",
            size=0.5)+
  geom_line(data=stock, aes(x=reddDensity, y =recruit2.5),linetype="dashed",
            size=0.5)+
geom_line(data=stock, aes(x=reddDensity, y =recruit3),linetype="dashed",
          size=0.5)+ 
  annotate("text", x = 3.7, y = 14, label = "3:1",size=3.5)+
  annotate("text", x = 5, y = 12.5, label = "2.5:1",size=3.5)+
  annotate("text", x = 6.5, y = 10.5, label = "2:1",size=3.5)

q1


q2<-ggplot(data=grave.sr,aes(x=reddDensity.graveonly , y=Fishper1000m ))+
  geom_point(color="blue")+
  labs(title="Grave Creek only")+
  xlab(bquote('Redd Density '('n km' ^-1))) + 
  ylab("") + 
  geom_line(data=redd.pred1, aes(x=reddDensity.graveonly, y =recruit ), size=1.25,
            color="blue")+
  
  geom_line(data=stock, aes(x=reddDensity, y =recruit2),linetype="dashed",
            size=0.5)+
  geom_line(data=stock, aes(x=reddDensity, y =recruit2.5),linetype="dashed",
            size=0.5)+
  geom_line(data=stock, aes(x=reddDensity, y =recruit3),linetype="dashed",
            size=0.5)+ 
  annotate("text", x = 3.7, y = 14, label = "3:1",size=3.5)+
  annotate("text", x = 5, y = 12.5, label = "2.5:1",size=3.5)+
  annotate("text", x = 6.5, y = 10.5, label = "2:1",size=3.5)
q2

qw2<-grid.arrange(q1,q2,ncol=2)  

save_plot("Results/stockRecruitFnct/ricker.grave.density.wreplacement.png", qw2, ncol = 2)







# residual plot
year<-c(seq(1997,2010,1),seq(2012,2015,1),seq(2018,2022,1))
# year is from the juvenile population estimate 
ricker.residuals<-mod$resid
ricker.residuals.graveonly<-mod.graveonly$resid

ricker.residuals.density<-mod1$resid
ricker.residuals.density.graveonly<-mod1.graveonly$resid

resids<-tibble(year,ricker.residuals,ricker.residuals.graveonly,
               ricker.residuals.density, ricker.residuals.density.graveonly)

r1<-ggplot(resids, aes (x = year, y=ricker.residuals))+
  geom_point()+
  labs(x="Year", y="Resdiduals from Count", 
       title="Grave Creek and tributaries")+
  geom_hline(yintercept=0,linetype="dashed")
r1
r2<-ggplot(resids, aes (x = year, y=ricker.residuals.graveonly))+
  geom_point()+
  labs(x="Year", y="", 
       title="Grave Creek only")+
  geom_hline(yintercept=0,linetype="dashed")
r2


r3<-ggplot(resids, aes (x = year, y=ricker.residuals.density))+
  geom_point()+
  labs(x="Year", y="Resdiduals from Density ", title="Grave Creek and tributaries")+
  geom_hline(yintercept=0,linetype="dashed")
r3

r4<-ggplot(resids, aes (x = year, y=ricker.residuals.density.graveonly))+
  geom_point()+
  labs(x="Year", y="", title="Grave Creek only")+
  geom_hline(yintercept=0,linetype="dashed")
r4


rr<-grid.arrange(r1,r2,r3,r4,ncol=2)  

save_plot("Results/stockRecruitFnct/residuals.png", rr, ncol = 2,nrow = 2)















