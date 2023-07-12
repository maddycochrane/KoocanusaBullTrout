


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

# bt redd data 
head(bt.redd.counts)
grave.bt.redd.dat<- bt.redd.counts %>% filter(LocalPopulation == "Grave Creek")

# combining redd counts across grave creek tribs into one total value 
grave.bt.redd.tot <- grave.bt.redd.dat %>% 
  group_by(LocalPopulation,Year)%>%
  summarise(sum(Redds))%>%
  rename(TotalRedds = `sum(Redds)`)%>%
  arrange(Year)%>%
  mutate(lagRedds = lag(TotalRedds, n=2))%>%
  dplyr::select(Year,lagRedds)

# join juv and redd data 
head(grave.bt.redd.tot)
head(grave.juv.bull.dat)

grave.sr<-left_join(grave.juv.bull.dat, grave.bt.redd.tot, by="Year")
head(grave.sr)


# plotting stock-recruit function 
grave.sr<-ggplot(data=grave.sr,aes(x=lagRedds, y=log(PopEst/lagRedds )))+
  geom_point()+
  labs(x="Redd Count (2 yr lag)", y="log(Juvenile Population Est. / Redd Count)",
       title="Grave Creek Stock-Recruit Function")+
  geom_smooth(data=grave.sr, aes(x=lagRedds, y=log(PopEst/lagRedds )),
              method = lm, formula = y ~ x)+
  ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq","R2","p")), formula = y ~ x)
grave.sr
# save image
#ggsave("Results/grave.sr.fnct.png",width = 4, height =4)

# linearized Ricker model 
mod<- lm(log(PopEst/lagRedds) ~ lagRedds, data = grave.sr) 
summary(mod)
# beta = -0.009735 , p < 0.001, F stat = 55.33



# view log(R/S) over time 
ggplot(data=grave.sr,aes(x=Year, y=log(PopEst/lagRedds )))+
         geom_point()+
  labs(y="log(Juvenile Population Est. / Redd Count)")
# save image
ggsave("Results/sr.time.png",width = 4, height =4)








