

### bull trout integrated population model ###

# 5.29.25
# madaline cochrane


library(R2jags)
library(dplyr)
library(ggplot2)
library(mcmcplots)
library(scales)
library(tidyhydat)
library(bayestestR)

# set working directory



####### summarize bull trout catch per unit effort data ####
# read in bull trout gill-netcatch data by year
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
head(adult.bull.net.sum)

#ggplot(data=bull.raw.net.dat, aes(x=Length_mm))+geom_histogram()+facet_wrap(~YYYY)
# what percentage of fish are > 75cm after restriction in 2015 goes into effect
bull.raw.net.dat%>%filter(YYYY>2014)%>%
  mutate(big_fish = if_else(Length_mm>=750,1,0))%>%group_by(big_fish)%>%
  count()
16/542 # only 3% of fish are larger or equal to 75cm

#########################################
## look at condition across years
adult.bull.net.sum2<-bull.raw.net.dat%>%
  mutate(Stage = if_else(Length_mm >= 500,"A","SA"))%>%  # adults are >= 500 mm
  filter(Stage == "A")%>%
  group_by(Section, YYYY,NetNumb)

ggplot(data=adult.bull.net.sum2, aes(x=log(Length_mm),y=log(Weight_g)))+geom_point()
bod.cond.bt <- lm(log(Weight_g) ~ log(Length_mm), data=adult.bull.net.sum2)
summary(bod.cond.bt)
# as size increases so does body condition 



condition<-adult.bull.net.sum2 %>% group_by(YYYY) %>% tidyr::drop_na(ConditionFultonsK) %>% filter(ConditionFiltered == TRUE)%>%
  filter(ConditionFultonsK >0)%>% filter(YYYY >1979) %>% dplyr::select(YYYY, Length_mm,Weight_g,Date,ConditionFultonsK,
                                                                       ConditionFiltered)%>%
  mutate(FultonK = (10^5*Weight_g) / (Length_mm^3))%>%
  mutate(lg.avg.weight = -11.72463 + 3.03431*log(Length_mm))%>%
  mutate(resid.body.condition = log(Weight_g)-lg.avg.weight) # positive is more weight than average


condition$ConditionFultonsK = as.numeric(condition$ConditionFultonsK)

condition.sums <- condition %>% 
  summarise(meanCondition = mean(resid.body.condition,na.rm =T), sdCondition = sd(resid.body.condition,na.rm =T))

ggplot(data=condition.sums, aes(x=YYYY, y=meanCondition))+geom_point()+
  geom_errorbar(aes(ymin=meanCondition-sdCondition, ymax=meanCondition+sdCondition),alpha=0.2)+
  geom_smooth()+
  #geom_line()+
  theme_minimal()+
  labs(x="Year",y="Residual body condition")+
  theme(axis.text.x=element_text(angle=45,hjust=1),#text = element_text(size = 20),
        plot.background = element_rect(fill = 'white', colour = 'white'))

condition <- condition %>% summarise(meanCondition = mean(FultonK,na.rm =T), sdCondition = sd(FultonK,na.rm =T))

ggplot(data=condition, aes(x=YYYY, y=meanCondition))+geom_point()+
  geom_errorbar(aes(ymin=meanCondition-sdCondition, ymax=meanCondition+sdCondition),alpha=0.2)+
  geom_smooth()+
  #geom_line()+
  theme_minimal()+
  labs(x="Year",y="Fulton's condition factor")+
  theme(axis.text.x=element_text(angle=45,hjust=1),#text = element_text(size = 20),
        plot.background = element_rect(fill = 'white', colour = 'white'))
#ggsave("Results/summary/bt-condition.png",width = 4, height = 4)




# read in spring net soak times
set.times<-read.csv("Data/updated.set.times.spring.nets.csv",header=TRUE) 
set.times <- set.times %>% 
  dplyr::select(YYYY,Date,Section,LOC_Name,TimeSet,TimePulled,HoursSet, NetNumb)
head(set.times)

# join adult catch and soak time datasets together to calculate CPUE
adult.bull.net.sums<-left_join(set.times,adult.bull.net.sum,by=c("Section","YYYY","NetNumb"))

# Replace NAs with 0s
adult.bull.net.sums$n[is.na(adult.bull.net.sums$n)] <- 0

adult.bull.net.sums <- adult.bull.net.sums %>% 
  mutate(cpue = round(n/HoursSet*16,0))%>% # standardizing CPUE by 16 hour soak time
  tibble::add_row(YYYY = 1983)%>% # no surveys this year
  dplyr::select(YYYY,cpue,NetNumb)%>%ungroup()%>%arrange(YYYY,NetNumb)

# add sample number by year as a new column
adult.bull.net.sums$sample <- data.table::rowid(adult.bull.net.sums$YYYY)
head(adult.bull.net.sums)

adult.bull.net.sums <- adult.bull.net.sums %>% dplyr::select(-NetNumb) # remove NetNumb column
adult.bull.net.final2<-tidyr::spread(adult.bull.net.sums,key=YYYY, value=cpue) # go from long to wide, so each year is a column
head(adult.bull.net.final2)

# nets per year summary
print(adult.bull.net.sums %>% group_by(YYYY)%>%summarise(n=n()),n=100)

adult<-adult.bull.net.sums%>%filter(YYYY>1979)%>%
  ggplot()+
  geom_boxplot(aes(x=YYYY,y=cpue,group=YYYY))+
  scale_x_continuous(limits=c(1979,2024), minor_breaks = seq(1980,2022,2),
                     breaks=seq(1980,2022,4),expand=c(0,0))+
  scale_y_continuous(limits=c(0,9),breaks=seq(0,10,2),expand=c(0.01,0))+
  labs(x="Year",y="Adult CPUE")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text = element_text(size = 20),
        plot.background = element_rect(fill = 'white', colour = 'white'))
adult
# save adult cpue by year figure
#ggsave("Results/summary/adult.bt.cpue.png",width = 8, height = 6) 

# view adult cpue distribution by year (take home = many zeros)
ggplot(data=adult.bull.net.sums, aes(x=cpue))+geom_histogram()+facet_wrap(~YYYY)


# want to summarize sub-adults by net and year
subadult.bull.net.sum<-bull.raw.net.dat%>%
  mutate(Stage = if_else(Length_mm >= 500,"A","SA"))%>% # sub-adults are < 500 mm
  filter(Stage == "SA")%>%
  group_by(Section, YYYY,NetNumb) %>%
  summarise(n = n())

# join sub-adult catch and soak time datasets together to calculate CPUE
subadult.bull.net.sums<-left_join(set.times,subadult.bull.net.sum,by=c("Section","YYYY","NetNumb"))

# Replace NAs with 0s
subadult.bull.net.sums$n[is.na(subadult.bull.net.sums$n)] <- 0

subadult.bull.net.sums <- subadult.bull.net.sums %>% 
  mutate(cpue = round(n/HoursSet*16,0))%>% # standardizing CPUE by 16 hour soak time
  tibble::add_row(YYYY = 1983)%>% # no surveys this year
  dplyr::select(YYYY,cpue,NetNumb)%>%ungroup()%>%arrange(YYYY,NetNumb)
tail(subadult.bull.net.sums)

# view sub-adult cpue by year
sub.cpue<-ggplot(data=subadult.bull.net.sums, aes(x=YYYY,y=cpue,group=YYYY))+
  geom_boxplot()+
  scale_x_continuous(limits=c(1979,2024), minor_breaks = seq(1980,2022,2),
                     breaks=seq(1980,2022,4),expand=c(0,0))+
  scale_y_continuous(limits=c(0,13),breaks=seq(0,12,2),expand=c(0,0))+
  labs(x="Year",y="Sub-Adult CPUE")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1),text = element_text(size = 20),
              plot.background = element_rect(fill = 'white', colour = 'white'))
sub.cpue
# save adult cpue by year figure
#ggsave("Results/summary/adult.bt.cpue.png",width = 8, height = 6) 


### stacked CPUE plots
duplot<-cowplot::plot_grid(adult +theme(axis.title.x = element_blank(),
                                        axis.text.x = element_blank()),
                           sub.cpue +theme(), 
                           nrow = 2, ncol =1, labels = "auto", align = "v")
duplot
#ggsave("Results/summary/bt.cpue.png",width = 10, height =10,units=c("in"),dpi=300)




# add capture number as id column
subadult.bull.net.sums$sample <- data.table::rowid(subadult.bull.net.sums$YYYY)
head(subadult.bull.net.sums)

subadult.bull.net.sums <- subadult.bull.net.sums %>% dplyr::select(-NetNumb)
subadult.bull.net.final2<-tidyr::spread(subadult.bull.net.sums,key=YYYY, value=cpue) # go from long to wide, so each year is a column
head(subadult.bull.net.final2)

# view sub-adult cpue distribution by year (take home = many zeros)
#ggplot(data=subadult.bull.net.sums, aes(x=cpue))+geom_histogram()+facet_wrap(~YYYY)



########################################################
# weight estimates by year from adult bull trout cpue data
weight.yr <- bull.raw.net.dat %>% 
  group_by(YYYY)%>%
  mutate(Stage = if_else(Length_mm >= 500,"A","SA"))%>% 
  filter(Stage == "A")%>%
  mutate(wt_lbs = Weight_g*0.00220462)%>% # converging to lbs
  mutate(wt_kg = Weight_g*0.001)%>% # converging to lbs
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

# view adult bull trout weight by year
ggplot(data=weight.yr,aes(x=Season,y=wt))+geom_point()+geom_line()


# average weights
bull.raw.net.dat %>% 
  mutate(Stage = if_else(Length_mm >= 500,"A","SA"))%>% 
  #filter(Stage == "A")%>%
  group_by(Stage)%>%
  mutate(wt_kg = Weight_g*0.001)%>%
  summarise(mean = mean(wt_kg,na.rm=T),sd = sd(wt_kg,na.rm=T))




#########################################
## read in bull trout redd data
bt.redds<-read.csv("Data/BullTroutRedds.csv",header=TRUE) 

bt.redds <- bt.redds %>%
  filter(CoreArea== "Lake Koocanusa")%>% 
  filter(LocalPopulation != "White River")%>%
  filter(LocalPopulation !=  "Blackfoot Creek")%>%
  mutate(PatchWB = paste(BLTPatchID,gsub(" ",'',Waterbody),sep='.')) %>%
  add_row(Year=2016,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  add_row(Year=2016,CoreArea="Lake Koocanusa", 
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  add_row(Year=1995,CoreArea="Lake Koocanusa", 
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  add_row(Year=2016,CoreArea="Lake Koocanusa", 
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=2016,CoreArea="Lake Koocanusa", 
          LocalPopulation = "Wigwam River",Waterbody="Wigwam River US",Redds=NA)%>%
  add_row(Year=2019,CoreArea="Lake Koocanusa", 
          LocalPopulation = "Wigwam River",Waterbody="Wigwam River US",Redds=NA)%>% 
  add_row(Year=2020,CoreArea="Lake Koocanusa", 
          LocalPopulation = "Wigwam River",Waterbody="Wigwam River US",Redds=NA)%>%
  
  add_row(Year=1985,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  add_row(Year=1986,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  add_row(Year=1987,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  add_row(Year=1988,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  add_row(Year=1989,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  add_row(Year=1990,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  add_row(Year=1991,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  add_row(Year=1992,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Blue Sky Creek",Redds=NA)%>%
  
  add_row(Year=1986,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  add_row(Year=1987,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  add_row(Year=1988,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  add_row(Year=1989,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  add_row(Year=1990,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  add_row(Year=1991,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  add_row(Year=1992,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Clarence Creek",Redds=NA)%>%
  
  add_row(Year=1986,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1987,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1988,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1989,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1990,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1991,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1992,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1993,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1994,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1976,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1977,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1978,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1979,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1980,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1981,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  add_row(Year=1982,CoreArea="Lake Koocanusa", # no surveys these years
          LocalPopulation = "Grave Creek",Waterbody="Grave Creek",Redds=NA)%>%
  arrange(Waterbody,Year)%>%
  ungroup()%>%
  group_by(Waterbody)

#Interpolate missing years (average of previous and next year)
bt.redds$Redds[which(is.na(bt.redds$Redds))] <- lapply(which(is.na(bt.redds$Redds)),FUN=function(x){
  mean(c(bt.redds$Redds[x-1],bt.redds$Redds[x+1]))
}) %>% unlist()

bt.redds1<-bt.redds%>%dplyr::select(Year,LocalPopulation,Waterbody,Redds) # get rid of unnecessary columns
bt.redd.wide <- tidyr::spread(data=bt.redds1, key=Year, value=Redds) # from long to wide (year as columns)
bt.redd.wide$Site <- seq(1,7) # adding site number in addition to name
#head(bt.redd.wide)



# combining bt redd counts into 2 streams instead of 7 streams (by population)
# grave creek and tributaries includes blue sky and clarence creeks
# wigwam counts from usa and canada are combined
bt.redds.combined<-bt.redds %>% 
  ungroup()%>%
  tidyr::drop_na(Redds)%>% # interpolating missing redds for combo streams using most recent value
  add_row(LocalPopulation = "Grave Creek", Waterbody="Blue Sky Creek",Year=1985,Redds = 1)%>%
  add_row(LocalPopulation = "Grave Creek", Waterbody="Grave Creek",Year=1993,Redds = 15)%>%
  add_row(LocalPopulation = "Grave Creek", Waterbody="Grave Creek",Year=1994,Redds = 15)%>%
  add_row(LocalPopulation = "Wigwam River", Waterbody="Wigwam River US",Year=1995,Redds = 12)%>%
  add_row(LocalPopulation = "Wigwam River", Waterbody="Wigwam River US",Year=2019,Redds = 4)%>% # avg btw 2,6
  add_row(LocalPopulation = "Wigwam River", Waterbody="Wigwam River US",Year=2020,Redds = 4)%>% # avg btw 2,6
  group_by(LocalPopulation,Year)%>%
  arrange(Waterbody,Year)%>%
  summarise(Redds = sum(Redds))%>%
  ungroup()%>%
  add_row(LocalPopulation="Grave Creek",Year=1976,Redds=NA)%>% # no surveys these years
  add_row(LocalPopulation="Grave Creek",Year=1977,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1978,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1979,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1980,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1981,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1982,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1986,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1987,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1988,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1989,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1990,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1991,Redds=NA)%>%
  add_row(LocalPopulation="Grave Creek",Year=1992,Redds=NA)%>%
  arrange(LocalPopulation,Year)

bt.redds.combined.wide<-tidyr::spread(data=bt.redds.combined, key=Year, value=Redds) # go from long to wide
# just using wigwam and grave for now 
bt.redds.combined.wide <- bt.redds.combined.wide[-2,] # removing skook
bt.redds.combined.wide <- bt.redds.combined.wide[-3,]  # removing wildhorse
head(bt.redds.combined.wide) # years are columns

# summary figure of wigwam and grave creek redd counts over time
bt.redds.combined %>% filter(LocalPopulation == "Grave Creek" | LocalPopulation == "Wigwam River")%>%
  filter(Year > 1982)%>%
  ggplot(aes(x=Year,y=Redds,color=LocalPopulation))+geom_point()+geom_line()+
  theme_bw()+labs(y="Redd Count")+
  scale_y_continuous(expand=c(0,0),limits=c(0,2400))+
  scale_x_continuous(expand=c(0.01,0),breaks=c(1985,1990,1995,2000,2005,2010,2015,2020))+theme_minimal()+
  theme(legend.title=element_blank(),
        legend.position = c(0.2,0.8),
        plot.background = element_rect(fill = 'white', colour = 'white'))
#ggsave("Results/summary/bt.redd.data.png",width = 6, height = 4)






#############################################################
# reading in bull trout harvest data
bt.harvest<-read.csv("Data/KoocanusaBTHarvest.csv",header=TRUE) 
head(bt.harvest)

bt.harvest <- bt.harvest %>% # adding in years when no harvest took place to match cpue and redd datasets
  add_row(Season=1975,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1976,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1977,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1978,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1979,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1980,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1981,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1982,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1983,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1984,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1985,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1986,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1987,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1988,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1989,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1990,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1991,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1992,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1993,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1994,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1995,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  add_row(Season=1996,BullTroutHarvested=0,Lower=0,Upper=0)%>%
  arrange(Season)

# view harvest data
ggplot(data=bt.harvest, aes(x=Season, y=BullTroutHarvested))+geom_point()+geom_line()+
  #geom_ribbon(aes(ymin=Lower,ymax=Upper),alpha=0.1)+
  scale_y_continuous(limits=c(0,700),expand=c(0,0))+labs(y="Harvest (n)")+theme_bw()




########## co-variate data ################
#### download bull river, british columbia, Canada discharge dataset 

# download_hydat() # might need to update HYDAT 
hy_daily_flows(station_number = "08NG002") # bull river near wardner

bull.river.q<-hy_daily_flows(station_number = "08NG002", 
                             start_date = "1976-01-01", 
                             end_date = "2022-12-31") # only archived to end of 2022 (as of this analysis 7.3.24)
# value = daily mean discharge in cms
tail(bull.river.q)

bull.river.q$Date <- # change date from character to posixct
  as.POSIXct(bull.river.q$Date, tz = "", format="%Y-%m-%d") 

# exported 2023 (un-archived) bull river data (real-time) to complete dataset
bull.river.q.23<-read.csv("Data/discharge/bullriver_2023.csv",header=TRUE) 
head(bull.river.q.23) 
tail(bull.river.q.23)

# Value is discharge in cms
bull.river.q.23$Date <- # change date from character to posixct
  as.POSIXct(bull.river.q.23$Date, tz = "", format="%m/%d/%Y %H:%M") 
str(bull.river.q.23) # check structure

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

# view hydrograph for bull river
ggplot(data=bull.river.q2, aes(x=Date,y=Value))+geom_line()+labs(y="Discharge (cms)")



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
  filter(month ==6 | month == 7 | month == 8 | month == 9)%>% # only want summer months
  group_by(Year)%>%
  summarise(min(rollingQ2))%>%
  dplyr::rename(summerLowflow = `min(rollingQ2)`)%>%
  mutate(scale_summerLowflow = scale(summerLowflow)) # scaling discharge to make comparable to other covariates
tail(bull.river.low)

# 3-yr running average of 14-day low flow 
rolling.low.bull<- zoo::rollapplyr(bull.river.low$summerLowflow,3, mean, align='right') 
rolling.low.bull<-tibble(rolling.low.bull[1:45],Year=seq(1979,2023))

bull.river.low<-left_join(bull.river.low,rolling.low.bull)
bull.river.low <- bull.river.low %>%
  rename(rolling.low = `rolling.low.bull[1:45]`)%>%
  mutate(scale_rolling_low = scale(rolling.low))%>%
  add_row(Year = 2024, scale_summerLowflow=0)
head(bull.river.low)

# view bull trout summer low flow calculations
ggplot(data=bull.river.low, aes(x=Year, y =summerLowflow))+geom_point()+geom_line()+
  theme_bw()+
  labs(x="Year", y="Lowest 14-day Summer Discharge (CMS)",title="Bull River")
#ggsave("Results/summary/low.flow.png",width = 6, height = 4)



### calculate peak winter discharge for bull river
bull.q.winter.peak <- bull.river.q2 %>%
  mutate(month = lubridate::month(Date))%>%
  filter(month > 9 | month < 5)%>% # only october through april
  mutate(Year = lubridate::year(Date))%>%
  ungroup()%>%
  mutate(Survey_year = ifelse(month < 4,Year-1,Year))%>% # survey year defined by the first year (start of winter)
  group_by(Survey_year)%>%
  summarise(max(Value))%>%
  rename(maxWinterq = `max(Value)`)%>%
  rename(Year = Survey_year)%>%
  mutate(scale_maxWinterq = scale(maxWinterq)) # scaling max discharge

# view peak winter discharge for bull river over time
ggplot(data=bull.q.winter.peak, aes(x=Year, y =maxWinterq))+geom_point()+geom_line()+
  theme_bw()+
  labs(x="Year", y="Peak Winter Discharge (CMS)",title="Bull River")
#ggsave("Results/summary/peak.flow.png",width = 6, height = 4)




############################ read in kokanee estimates from separate r script / analysis ###################
predict <-readRDS(file="R/kokanee_biomass")
head(predict)
# for 3-yr running average of kokanee we need to give the first 2 yrs (1980, 81) the same values as 1983
predict$running_accessible_biomass[1]<-predict$running_accessible_biomass[3] # giving first and second values the same as third
predict$running_accessible_biomass[2]<-predict$running_accessible_biomass[3]

predict <- predict %>%
  mutate(running_biomass_scale = scale(running_accessible_biomass)) # scaling kokanee biomass to make comptabile with discharge, pool elevation data


################# stage-structured ipm model for bull trout  #####################

# group all data sets together into a list to pass to jags
jags.data<-list(y.adult=adult.bull.net.final2[,6:49], # adult bt cpue, 1980 - 2023
                y.subadult=subadult.bull.net.final2[,6:49], # sub-adult bt cpue, 1980 - 2023
                nYears=ncol(adult.bull.net.final2[,6:49]), # 44  years
                nReps = nrow(adult.bull.net.final2[,6:49]), # 45 max nets/yr
                y.redd = bt.redds.combined.wide[,6:49], # redd counts for wigwam and grave + tribs, 1980 - 2023
                nSites = nrow(bt.redds.combined.wide[,6:49]), # 2 sites
                bt.harv = bt.harvest$BullTroutHarvested[6:48], # lake kook adult bt harvest, 1980-2022 (season ends in march of nxt yr), so 2004 harvest yr will affect survival from 2004-2005
                weight =weight.yr$wt[5:44], # adult bull trout weight, 1980 - 2019
                
                kok = predict$running_biomass_scale[1:43], # available kokanee biomass (sub-adult and adult), 1980 kok biomass will affect survival from 1980 - 1981
                
                low.flow.recruitment = bull.river.low$scale_summerLowflow[6:45], # low summer flow data from 1981 - 2020 for bull river, accounting for 3 yr lag to NSub estimate
                #low.flow.recruitment2 = bull.river.low$scale_summerLowflow[7:46], # low summer flow data from 1982 - 2021 
                peak.flow = bull.q.winter.peak$scale_maxWinterq[5:44], # max winter discharge from 1980 - 2019 for bull river
                
                #pool.elevation.recruit = c(pool.elev.sum$scale_avgPoolElevation[5:44]), # pool elevation data for 1983 - 2022
                ca.reg.changes = c(rep(1,16),rep(0,27))) # 1980 - 1995 no sig. fishing regs, 1996 - 2022 catch and release only 


# jags model 
cat(file="jags/s-r-new2.txt", "
# likelihood
model{ 

# Process Likelihood
  for (t in 1:(nYears-4)) { # 1984 - 2023 (40 observtions) as recruitment is lagged by 4 yrs

  # covariate effects on alpha parameter (density-independent) in Ricker model
    alpha[t] <- exp(alpha_sr + beta.low.flow*low.flow.recruitment[t] + beta.peak.winter*peak.flow[t] ) 
    
  # Ricker model for sub-adult BT recruitment
    log_predR[t+4] <- log(biomass[t]) + log(alpha[t]) - beta_sr*biomass[t] 
    
    biomass[t]<-weight[t]*Nadults[t]*mean.fish.per.redd  
    
    logR_Obs[t+4] ~ dnorm(log_predR[t+4], tauR_process) # recruitment includes process error/unaccounted for variation
    yr1subs[t+4] <- exp(logR_Obs[t+4]) # 4-yr lag and log transformation
    Nsubadults[t+4] <- yr1subs[t+4] + (1-transition.prob)*S.sa*Nsubadults[t+3] 
  }
  
  for(t in 1:(nYears-1)){
  # covariate effects on adult BT survival
    S.a.step1[t] <- a.int + beta.dens*(NumbLake[t]/10000)+ 
                  beta.kok*kok[t] + beta.reg.change*ca.reg.changes[t] # scaling raw abundance
                  
    mean.S.a[t] <- 1/(1+exp(-S.a.step1[t])) #  logit transformation
    S.a[t] ~ dnorm(mean.S.a[t], tau.S.a)T(0,1) # yr-specific survival includes process error 
    
    NA[t]<-max(Nadults[t],bt.harv[t]) # to keep Nadults above 0 (simulation problem when N can't go below zero)
    Nadults[t+1] <- (NA[t]-bt.harv[t])*S.a[t] + transition.prob*S.sa*Nsubadults[t]
    
  # Derived Parameters
    lambdaA[t] <-Nadults[t+1]/Nadults[t] # change in adult population size over time
    lambda[t] <-Nsubadults[t+1]/Nsubadults[t] # change in sub-adult population size over time
    NumbLake[t]<-Nsubadults[t]+Nadults[t] # using total bt population to understand density dependence
  }
  
  
# Observation Likelihood
  for (t in 1:nYears){
    for(i in 1:nReps){ 
# CPUE 
      y.adult[i,t] ~ dnegbin(p.a[t], r.a) # includes observation error
      y.subadult[i,t] ~ dnegbin(p.sa[t], r.sa) # includes observation error
    }
    p.a[t] <- r.a/(r.a+(Nadults[t]*q.lake))
    p.sa[t] <- r.sa/(r.sa+(Nsubadults[t]*q.lake.sa))
    
    for(j in 1:nSites){
# Redd counts
      y.redd[j,t] ~ dnegbin(p.a.redd[j,t],r.a.redd)
      p.a.redd[j,t] <- r.a.redd/(r.a.redd+Nspawners[t]*q[j])
    }
    Nspawners[t] <- Nadults[t]*mean.fish.per.redd
    }
    
  
  # PRIORS
  mean.fish.per.redd <- 0.5 # 2 fish/redd
  
  # intrinsic productivity (recruitment prior)
  alpha_sr ~ dunif(-10,2.5)
  
  # strength of dens depend (recruitment prior)
  beta_inv ~ dunif(0,5000)
  beta_sr <- 1/beta_inv 
  
  # priors for covariate effects on sub-adult recruitment/alpha parameter
  beta.low.flow ~ dnorm(0,1)
  beta.peak.winter ~ dnorm(0,1)
  beta.pool.recruitment ~ dnorm(0,1)
  ind ~ dbeta(1,3) # weighted towards age-1 effect more than age-2 for this indicator variable
  
  # unk process variance for recruits model 
  sigmaR_process ~ dnorm(0, 1/3^2)T(0,) # sd = 3
  tauR_process <- pow(sigmaR_process, -2)
  
  # catchability coefficients (scaling cpue to true abundance)
  q.lake ~ dnorm(0.001,100)T(0,)
  q.lake.sa ~ dnorm(0.001,100)T(0,)
  q ~ ddirch(c(0.08,0.92)) 
  
  # over-disperstion parameter for neg binomial distribution
  r.a ~ dgamma(5,1)
  r.sa ~ dgamma(5,1)
  r.a.redd ~ dgamma(1,1)
  
  # Adult bt survival priors
  sd.S.a ~ dnorm(0.1,2500)T(0,) # sd = 0.02
  tau.S.a <- pow(sd.S.a, -2)
  a.int ~ dnorm(2.2,1/0.5^2)
  beta.kok ~ dnorm(0.5,4) # effect of kokanee abundance on adult bt survival
  beta.dens ~ dnorm(-1,4)T(,0) # dd effect (forcing negative)
  beta.reg.change ~ dnorm(-0.5,1) # when canadian harvest open, it will negatively affect survival

  # Transition probability (from sub-adult to adult bt) prior
  transition.prob~ dnorm(0.5,400)T(0,1) # sd = 0.05
  S.sa ~ dnorm(0.6,400)T(0,1) # sd = 0.05
  
  Nadults[1] ~ dunif(500,1000) 
  
  # initial sub-adult population size estimates when recruitment model can't be fit due to 4-yr lag
  Nsubadults[1] ~ dunif(500,2500)
  Nsubadults[2] ~ dunif(500,2500)
  Nsubadults[3] ~ dunif(500,2500)
  Nsubadults[4] ~ dunif(500,2500)
} 
")

# parameters monitored 
params <- c("Nadults","r.a","Nsubadults", "lambda","q", "r.sa","r.a.redd",
            "q.lake.sa","q.lake","sigmaR_process", "alpha_sr","beta_sr", 
            "lambdaA","S.a", "sd.S.a","transition.prob","a.int",
            "beta.kok","beta.dens","sd.lambda","mean.lambda", "transition.prob",
            "alpha","beta.low.flow","beta.peak.winter",
            "beta.reg.change","S.sa","yr1subs") 

# initial parameter draws
inits <- function(){
  list(
    alpha_sr = runif(1,0.2,0.5),
    beta_inv = runif(1,4500,4600),
    sigmaR_process = runif(1,0.4,0.5),
    q.lake = runif(1,0.0001,0.0009),
    q.lake.sa = runif(1,0.001,0.0012),
    q = c(0.08,0.92),
    r.a = runif(1,2,5),
    r.sa = runif(1,2,4),
    r.a.redd = runif(1,7.1,7.5),
    'Nadults[1]' = runif(1,500,1000),
    sd.S.a = runif(1,0.08,0.12),
    a.int = runif(1,2,2.5),
    S.sa = runif(1,0.55,0.6),
    beta.kok = runif(1,0.5,2),
    beta.dens = runif(1,-1,-0.5),
    beta.reg.change = runif(1,-2,-1),
    transition.prob = runif(0,0.48,0.51),
    'Nsubadults[1]' = runif(1,500,1500),
    'Nsubadults[2]' = runif(1,1000,2000),
    'Nsubadults[3]' = runif(1,1000,2000),
    'Nsubadults[4]' = runif(1,1000,1500),
    #beta.low.flow.redd = runif(1,-0.5,0.5),
    beta.low.flow = runif(1,0,0.5),
    beta.peak.winter = runif(1,-0.5,0.5))
}

# run model - this is going to take awhile 
full.mod2 <- jags.parallel(data = jags.data, model.file = "jags/s-r-new2.txt", parameters.to.save = params,
              inits=inits, n.chains = 3, n.iter = 25000, n.burnin = 23000, n.thin=5)

# view all results
options(max.print=1000000)
print(full.mod2,digits = 3)
mcmcplot(full.mod2) # check for convergence, rhat values, etc


full.mod2$BUGSoutput$mean$beta_sr
quantile(full.mod2$BUGSoutput$sims.list$beta_sr,0.025)
quantile(full.mod2$BUGSoutput$sims.list$beta_sr,0.975)

full.mod2$BUGSoutput$mean$q.lake
quantile(full.mod2$BUGSoutput$sims.list$q.lake,0.025)
quantile(full.mod2$BUGSoutput$sims.list$q.lake,0.975)

full.mod2$BUGSoutput$mean$q.lake.sa
quantile(full.mod2$BUGSoutput$sims.list$q.lake.sa,0.025)
quantile(full.mod2$BUGSoutput$sims.list$q.lake.sa,0.975)

1/(1+exp(-full.mod2$BUGSoutput$mean$a.int))
1/(1+exp(-quantile(full.mod2$BUGSoutput$sims.list$a.int,0.025)))
1/(1+exp(-quantile(full.mod2$BUGSoutput$sims.list$a.int,0.975)))

# save jag model results
#saveRDS(full.mod2, file="R/full.mod.new.ipm") # new one to use
#full.mod2 <-readRDS(file="R/full.mod.new.ipm")



Nadult.est<-c(full.mod2$BUGSoutput$mean$Nadults) # adult bt abundance estimates
Nsubadult.est<-c(full.mod2$BUGSoutput$mean$Nsubadults[1:44]) # sub-adult bt abundance estimates
Nlake.est<-Nadult.est+Nsubadult.est # combine adult and sub-adult abundances
meanNlake<-mean(Nlake.est) 



#######################################
# looking at stock-recruit replacement 
Nreplace <- Nadult.est[7:44] # lagging adults by 6 years
Nadult <- Nadult.est[1:38]
line<- seq(0,5000,length.out=38)

res2<-tibble(Nadult,Nreplace,yrs=seq(1986,2023),line)
res2 <- res2 %>% # only want last two digits of years
  mutate(id=substr(yrs,nchar(yrs)-1,nchar(yrs)))

# plotting
ggplot(data=res2, aes(x=Nadult,y=Nreplace))+
  geom_text(aes(label=id))+
  geom_line(aes(x=line,y=line),color="grey")+
  labs(x="Adult stock (n)", y="Adults lagged 6-yrs (n)")+theme_minimal()+
  annotate("text", x = 4800, y = 4400, label = "1:1",color="darkgrey")+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))
#plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
#    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
#   axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
# yr represents the replacement year not the initial adult stock size 
#ggsave("Results/summary/PubFigs/Figure 2.png",width = 12, height =8.5, units="cm", dpi=500)


# need to lag sub-adult estimate for comparison in figure by 4 yrs
Nsubadult.est<-c(full.mod2$BUGSoutput$mean$Nsubadults)
Nsubadult.est2 <- Nsubadult.est[5:44] 
Nadult.est2 <- Nadult.est[1:40]
wt<-weight.yr$wt[5:44] # 1980-2019
a <-c(exp(full.mod2$BUGSoutput$mean$alpha_sr))
b<-c(full.mod2$BUGSoutput$mean$beta_sr)
yr1subs<-c(full.mod2$BUGSoutput$mean$yr1subs[5:44] )
#yr1subs<-Nsubadult.est[5:44] - (1-transition)*S.a*Nsubadult.est[4:43]


res<-tibble(Nadult.est2,wt,Nsubadult.est2,yrs=seq(1984,2023),yr1subs) # yrs are for recruit yr
res <- res %>% # only want last two digits of years
  mutate(id=substr(yrs,nchar(yrs)-1,nchar(yrs)))

res <- res %>%
  mutate(biomass=wt*Nadult.est2*0.5)%>%
  mutate(model_recruitment1 = log(biomass) + log(a) - b*biomass)%>%
  
  # Nsubadults[t+4] <- yr1subs[t+4] + (1-transition.prob)*S.a[t+3]*Nsubadults[t+3] 
  mutate(model_recruitment = exp(model_recruitment1))%>%
  mutate(residual = yr1subs - model_recruitment)

# how are residuals from stock-recruitment changing over time
ggplot(data=res,aes(x=yrs,y=residual))+
  geom_point()+
  #geom_text(aes(label=id))+
  #stat_smooth()+ # method="lm"
  scale_y_continuous(limits=c(-2500,2700))+
  labs(x="Year",y="Residual")+theme_minimal()+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))
# negative values = Nestimate is lower than expected from best fit line
#ggsave("Results/summary/PubFigs/Figure S2.png",width = 4, height = 4)


#ggplot(data=res, aes(x=biomass,y=yr1subs))+geom_point()
#ggplot(data=res, aes(x=biomass,y=yr1subs))+ geom_text(aes(label=id))





######################################
### visualize results of stock-recruitment curve
transition<-c(full.mod2$BUGSoutput$mean$transition.prob)
S.sa<-c(full.mod2$BUGSoutput$mean$S.sa)
transition # 0.47
# Nsubadults[t+4] <- exp(logR_Obs[t+4]) + (1-transition.prob)*Nsubadults[t+3]
# Nsubadults[t+4] <- yr1subs[t+4] + (1-transition.prob)*S.a[t+3]*Nsubadults[t+3] 

Nadult.est2 <- Nadult.est[1:40]
wt<-weight.yr$wt[5:44] # 1980-2019
Nsubadult.est2 <- Nsubadult.est[5:44] 
yr1subs.est2 <- yr1subs[5:44] 

#surv<-c(full.mod2$BUGSoutput$mean$S.a)

res<-tibble(Nadult.est2,wt,Nsubadult.est2,yrs=seq(1984,2023), yr1subs.est2 ) # yrs are for recruit yr  
head(res)
res <- res %>% # only want last two digits of years
  mutate(id=substr(yrs,nchar(yrs)-1,nchar(yrs)))%>%
  mutate(biomass=wt*Nadult.est2*0.5)
#mutate(prvYrrecruits = (lag(Nsubadult.est2)*S.a )*(1-transition))%>%
#mutate(yr1recruits = Nsubadult.est2 - prvYrrecruits)
mean.biomass<-c(mean(res$biomass))

head(res)

sr.sim.dat <- as.matrix(data.frame(constant = 1,
                                   biomass = seq(from = 0.001, to = max(res$biomass), 
                                                 length.out=20)))


# posterior draws
posterior.draws<-as.matrix(as.mcmc(full.mod2))
posterior.draws<-as_tibble(posterior.draws)
posterior.draws3<-posterior.draws %>%
  dplyr::select(alpha_sr,beta_sr)%>% # make sure this is in the same order as sim.dat # 
  mutate(lg_alpha = log(exp(alpha_sr)))%>%
  mutate(beta = -beta_sr)%>%
  dplyr::select(lg_alpha,beta)
#head(posterior.draws3)
posterior.draws3<-as.matrix(posterior.draws3)

Xb <- t(sr.sim.dat %*% t(posterior.draws3)) # multiplying all posterior draws by each row
#head(Xb)
sr.pp <- as.data.frame(Xb)

names(sr.pp) <- sr.sim.dat[, "biomass"]
#head(sr.pp)
## these are draws for alpha - beta (still need to add log(biomass))

sr.pp.long <- tidyr::gather(sr.pp, key = "biomass")
#head(sr.pp.long)
sr.pp.long$biomass<-as.numeric(sr.pp.long$biomass)
sr.pp.long <- sr.pp.long %>% 
  mutate(value = value + log(biomass))

sr.pp.sum <- summarize(group_by(sr.pp.long, biomass), # summarize results
                       median_pp = median(value), 
                       mean_pp = mean(value),
                       lower_pp = quantile(value, probs = c(0.025)), 
                       upper_pp = quantile(value, probs = c(0.975)),
                       low10_pp = quantile(value, probs = c(0.1)),
                       high90_pp = quantile(value, probs = c(0.90)),
                       low25_pp = quantile(value, probs = c(0.25)),
                       high75_pp = quantile(value, probs = c(0.75)))

sr.pp.sum$biomass <- as.numeric(sr.pp.sum$biomass)

# plotting
ggplot(data = sr.pp.sum, aes(x = biomass*0.453592, y = exp(mean_pp))) +
  geom_ribbon(aes(ymin=exp(lower_pp),ymax=exp(upper_pp)),alpha=0.1)+
  #geom_point(size=3)+
  geom_line(size=1.5,col="grey") + 
  geom_point(data=res,aes(x=biomass*0.453592, y=yr1subs.est2))+
  scale_x_continuous(expand=c(0,80))+scale_y_continuous(expand=c(0,50))+
  theme_minimal()+
  labs(x="Spawning Biomass (kg)", y="New Sub-Adult Recruits")+
  theme(plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
#ggsave("Results/summary/PubFigs/Figure 3.png",width = 12, height =8.5, units="cm", dpi=500)




# plotting with years instead of points
ggplot(data = sr.pp.sum, aes(x = biomass, y = exp(mean_pp))) +
  geom_ribbon(aes(ymin=exp(lower_pp),ymax=exp(upper_pp)),alpha=0.1)+
  geom_line(size=1.5) + 
  #scale_x_continuous(expand=c(0,80))+scale_y_continuous(expand=c(0,50))+
  theme_minimal()+
  geom_text(data=res,aes(x=biomass, y=yr1subs.est2,label=id))+
  labs(x="Spawning Biomass (lbs)", y="New Sub-Adult Recruits")+
  theme( #plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
    plot.background = element_rect(fill = 'white', colour = 'white'),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
#ggsave("Results/summary/s-r-curve-w-dates.png",width = 6, height = 4)



# effect size 

# maximize recruitment 
exp(max(sr.pp.sum$mean_pp)) # 1696 new sub-adult individuals
# biomass = 4312 lb for max (from sr.pp.sum)
4347.2520*0.453592 # 1972 kg (spawning biomass) for maximum sustainable yield
tail(weight.yr) # current weight = 5.31 lb
4347.2520*2/5.31 # 1637 = number of adults that maximizes recruitment with current bull trout weight


# highest biomass = 14105 lb
13766.2958*0.453592 # 6244 kg
exp(6.497433) # new individuals at max biomass = 663

# number decrease
1696.343- 663.4364 # 1033
# percent decrease
((663.4364-1696.343)/1696.343)*100 # 61% decrease






full.mod2$BUGSoutput$mean$beta_sr # 0.00022
quantile(full.mod2$BUGSoutput$sims.list$beta_sr,0.025) # 0.00020
quantile(full.mod2$BUGSoutput$sims.list$beta_sr,0.975) # 0.00026

full.mod2$BUGSoutput$mean$q.lake # 0.00055
quantile(full.mod2$BUGSoutput$sims.list$q.lake,0.025) # 0.00047
quantile(full.mod2$BUGSoutput$sims.list$q.lake,0.975) # 0.00063

full.mod2$BUGSoutput$mean$q.lake.sa # 0.00115
quantile(full.mod2$BUGSoutput$sims.list$q.lake.sa,0.025) # 0.00086
quantile(full.mod2$BUGSoutput$sims.list$q.lake.sa,0.975) # 0.00154

1/(1+exp(-full.mod2$BUGSoutput$mean$a.int)) # 0.92
1/(1+exp(-quantile(full.mod2$BUGSoutput$sims.list$a.int,0.025))) # 0.85
1/(1+exp(-quantile(full.mod2$BUGSoutput$sims.list$a.int,0.975))) # 0.96



## residuals
data.grave.redd <- bt.redds.combined %>% filter(LocalPopulation == "Grave Creek")%>%
  filter(Year>1989)%>%rename(Redds.grave = Redds)%>%dplyr::select(-LocalPopulation)
data.wigwam.redd <- bt.redds.combined %>% filter(LocalPopulation == "Wigwam River")%>%
  filter(Year>1989)%>%rename(Redds.wigwam = Redds)%>%dplyr::select(-LocalPopulation)
head(data.wigwam.redd)
p.grave <- full.mod2$BUGSoutput$mean$q[1]
p.wig<-full.mod2$BUGSoutput$mean$q[2]
q.lake<-full.mod2$BUGSoutput$mean$q.lake
model.NA = full.mod2$BUGSoutput$mean$Nadults
model.NA.low =full.mod2$BUGSoutput$summary[1:44,3]
model.NA.high =full.mod2$BUGSoutput$summary[1:44,7]
median.cpue<-adult.bull.net.sums %>% group_by(YYYY)%>%filter(YYYY>1979)%>%rename(Year=YYYY)%>%
  summarise(median.cpue = median(cpue))

resid.df <- tibble( model.NA=model.NA,Year=seq(1980,2023,1),model.NA.low=model.NA.low,model.NA.high=model.NA.high)
resid.df<-left_join(resid.df,data.grave.redd)
resid.df<-left_join(resid.df,data.wigwam.redd)
resid.df<-left_join(resid.df,median.cpue,by="Year")
head(resid.df)
resid.df <- resid.df %>% 
  mutate(model.grave = (model.NA/2)*p.grave)%>%
  mutate(model.wigwam = (model.NA/2)*p.wig)%>%
  #mutate(data.cpue = median.cpue/q.lake)%>%
  mutate(resid.grave = Redds.grave-model.grave)%>%
  mutate(resid.wig = Redds.wigwam-model.wigwam)%>%
  mutate(resid.cpue = median.cpue-(model.NA*c(q.lake)))
head(resid.df)

ggplot(data=resid.df, aes(x=Year,y=resid.wig))+geom_point()+#geom_line()+
  geom_point(aes(y=resid.grave),color="blue")+
  #geom_line(aes(y=resid.grave),color="blue")+
  #geom_point(aes(y=resid.cpue),color="red")+geom_line(aes(y=resid.cpue),color="red")+
  scale_x_continuous(limits=c(1994,2023))+
  scale_y_continuous(limits=c(-1700,1700))+
  theme_minimal()+labs(y="Difference between data and model")





ggplot(data=resid.df, aes(x=Year,y=Redds.grave+Redds.wigwam))+geom_point()+geom_line()+
  geom_line(aes(x=Year,y=model.NA/2),color="grey")+
  geom_point(aes(x=Year,y=model.NA/2),color="grey")+
  geom_ribbon(aes(x=Year,ymin=model.NA.low/2, ymax=model.NA.high/2),alpha=0.1)+
  theme_minimal()+labs(y="Redds")+
  annotate("text", x = 1985, y = 520, label = "Model",color="grey")+
  annotate("text", x = 2006, y = 710, label = "Data")
#ggsave("Results/summary/PubFigs/Figure S4.png",width = 12, height =8.5, units="cm", dpi=500)


adult.bull.net.sums%>%filter(YYYY>1979)%>%
  ggplot()+
  geom_boxplot(aes(x=YYYY,y=cpue,group=YYYY))+
  scale_x_continuous(limits=c(1979,2024), minor_breaks = seq(1980,2022,2),
                     breaks=seq(1980,2022,4),expand=c(0,0))+
  scale_y_continuous(limits=c(0,9),breaks=seq(0,10,2),expand=c(0.01,0))+
  labs(x="Year",y="Adult CPUE")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1), #text = element_text(size = 20),
        plot.background = element_rect(fill = 'white', colour = 'white'))+
  geom_point(data=resid.df, aes(x=Year,y=model.NA*c(q.lake)),color="grey")+
  geom_line(data=resid.df, aes(x=Year,y=model.NA*c(q.lake)),color="grey")+
  geom_ribbon(data=resid.df, aes(x=Year,ymin=model.NA.low*c(q.lake), ymax=model.NA.high*c(q.lake)),alpha=0.1)
#ggsave("Results/summary/PubFigs/Figure S5.png",width = 12, height =8.5, units="cm", dpi=500)


ggplot(data=resid.df, aes(x=Year,y=resid.cpue))+geom_point()+geom_line()


# qq plot equivalent
grave.fig<-ggplot(data=resid.df, aes(x=Redds.grave,y=model.grave))+geom_point()+
  scale_y_continuous(limits=c(0,250))+labs(x="Redds in Grave",y="Modeled Redds")+
  geom_abline (slope=1, linetype = "dashed", color="Red")+theme_minimal()
wig.fig<-ggplot(data=resid.df, aes(x=Redds.wigwam,y=model.wigwam))+geom_point()+
  scale_y_continuous(limits=c(0,2200))+labs(x="Redds in Wigwam",y="Modeled Redds")+
  geom_abline (slope=1, linetype = "dashed", color="Red")+theme_minimal()


#ggplot(data=resid.df, aes(x=data.cpue,y=model.NA))+geom_point()
cpue.fig<-ggplot(data=resid.df, aes(x=median.cpue,y=model.NA*c(q.lake)))+geom_point()+
  scale_y_continuous(limits=c(0,3.5))+labs(x="Median CPUE",y="Modeled CPUE")+
  geom_abline (slope=1, linetype = "dashed", color="Red")+theme_minimal()



### stacked CPUE plots
duplot<-cowplot::plot_grid(grave.fig +theme(),
                           wig.fig +theme(), 
                           
                           cpue.fig,
                           nrow = 1, ncol =3, labels = "auto", align = "v")
duplot
#ggsave("Results/summary/bt.resids.png",width = 14, height =6,units=c("in"),dpi=500)




# r(1-p)/p^2 = expected variance in negative binomial 
r.a<-full.mod2$BUGSoutput$mean$r.a
obs.var<-adult.bull.net.sums %>% group_by(YYYY)%>%filter(YYYY>1979)%>%
  summarise(obs.var = var(cpue))
head(obs.var)
exp.var <- tibble(model.NA = model.NA,q.lake=q.lake,r.a=r.a,obs.var=obs.var$obs.var)
exp.var <-exp.var %>%
  mutate(p.a =  r.a/(r.a+(model.NA*q.lake)))%>%
  mutate(exp.var = r.a*(1-p.a)/(p.a^2))
head(exp.var)

## overdispersion parameter - good fit??
png(filename="Results/summary/observedVSexpectedVARIANCE.png")
plot(exp.var$obs.var,exp.var$exp.var,xlab="Observed variance",ylab="Expected variance",main="CPUE",
     ylim=c(0,5))
dev.off()


png(filename="Results/summary/dispersion.posteriors.png",width = 880, height = 480, units = "px")
par(mfrow=c(1,2))
hist(full.mod2$BUGSoutput$sims.list$r.a, #xlim=c(1,10),
     main="CPUE",xlab="Parameter estimate") # posterior distribution of overdispersion parameter (>1 = overdispersion)
hist(full.mod2$BUGSoutput$sims.list$r.a.redd,
     main="Redd",xlab="Parameter estimate")
dev.off()

r.a.redd<-full.mod2$BUGSoutput$mean$r.a.redd

exp.var.redd <- tibble(model.NA = model.NA,q.wig = p.wig, q.grave = p.grave, r.a=r.a.redd)
head(exp.var.redd)
exp.var.redd <-exp.var.redd %>%
  mutate(p.a =  r.a/(r.a+(model.NA/2*q.wig)))%>%
  mutate(exp.var = r.a*(1-p.a)/(p.a^2))
mean(exp.var.redd$exp.var) # 203,050
var(resid.df$Redds.wigwam,na.rm=TRUE) #  26277

exp.var.redd <-exp.var.redd %>%
  mutate(p.a =  r.a/(r.a+(model.NA/2*q.grave)))%>%
  mutate(exp.var = r.a*(1-p.a)/(p.a^2))
mean(exp.var.redd$exp.var) # 1686
var(resid.df$Redds.grave,na.rm=TRUE)  # 3631





### probability of distribution for beta values
median(full.mod2$BUGSoutput$sims.list$beta.kok) # positive
p_direction(x=full.mod2$BUGSoutput$sims.list$beta.kok, method="direct",null=0)
# 99% proportion of the posterior distribution that is of the median


median(full.mod2$BUGSoutput$sims.list$beta.dens) # negative
p_direction(x=full.mod2$BUGSoutput$sims.list$beta.dens, method="direct",null=0)
# 100% proportion of the posterior distribution that is of the median’s sign

median(full.mod2$BUGSoutput$sims.list$beta.low.flow) # positive 
p_direction(x=full.mod2$BUGSoutput$sims.list$beta.low.flow, method="direct",null=0)
# 94% proportion of the posterior distribution that is of the median’s sign


median(full.mod2$BUGSoutput$sims.list$beta.peak.winter) # negative (barely) 
p_direction(x=full.mod2$BUGSoutput$sims.list$beta.peak.winter, method="direct",null=0)
# 62% proportion of the posterior distribution that is of the median’s sign

median(full.mod2$BUGSoutput$sims.list$beta.reg.change) # negative
p_direction(x=full.mod2$BUGSoutput$sims.list$beta.reg.change, method="direct",null=0)
# 100% proportion of the posterior distribution that is of the median’s sign




############# visualize co-variate relationships on adult survival 

Nadult.est<-c(full.mod2$BUGSoutput$mean$Nadults) # adult bt abundance estimates
Nsubadult.est<-c(full.mod2$BUGSoutput$mean$Nsubadults[1:44]) # sub-adult bt abundance estimates
Nlake.est<-Nadult.est+Nsubadult.est # combine adult and sub-adult abundances
meanNlake<-mean(Nlake.est) 

#a.int<-c(full.mod2$BUGSoutput$mean$a.int) # maximum adult survival rate
#beta.dens<-full.mod2$BUGSoutput$mean$beta.dens # effect of bull trout density on adult survival
#beta.kok<-c(full.mod2$BUGSoutput$mean$beta.kok) # kokanee biomass effect on adult survival
#beta.reg.change<-c(full.mod2$BUGSoutput$mean$beta.reg.change) # effect of canadian fishing regulation change on adult survival




##### varying available kokanee biomass
kok.sim.dat <- as.matrix(data.frame(constant = 1,
                                    NumbLake = meanNlake/10000,   
                                    kok = seq(from = min(predict$running_biomass_scale,na.rm=T), 
                                              to = max(predict$running_biomass_scale,na.rm=T), 
                                              length.out=20),
                                    ca.reg.changes=1)) # open / NOT catch and release (second half of study)
#head(kok.sim.dat)

# posterior draws
posterior.draws<-as.matrix(as.mcmc(full.mod2))
posterior.draws<-as_tibble(posterior.draws)

#head(posterior.draws)
posterior.draws<-posterior.draws %>%
  dplyr::select(a.int,beta.dens,beta.kok,beta.reg.change) # 
posterior.draws<-as.matrix(posterior.draws)


Xb <- t(kok.sim.dat %*% t(posterior.draws)) # multiplying all posterior draws by possible co-variate value
#head(Xb)
Xc<-plogis(Xb)# still need to transform
#head(Xc)
kok.pp <- as.data.frame(Xc)

names(kok.pp) <- kok.sim.dat[, "kok"] # 
#head(kok.pp)

kok.pp.long <- tidyr::gather(kok.pp, key = "kok") # from wide to long
#head(kok.pp.long)
kok.pp.sum2 <- summarize(group_by(kok.pp.long, kok), # summarize results 
                         mean_pp = mean(value), 
                         median_pp = median(value), 
                         lower_pp = quantile(value, probs = c(0.025)), 
                         upper_pp = quantile(value, probs = c(0.975)),
                         low10_pp = quantile(value, probs = c(0.1)),
                         high90_pp = quantile(value, probs = c(0.90)),
                         low25_pp = quantile(value, probs = c(0.25)),
                         high75_pp = quantile(value, probs = c(0.75)))
#head(kok.pp.sum2)

kok.pp.sum2$kok <- as.numeric(kok.pp.sum2$kok)

# view results with untransformed x-axis
ggplot(data = kok.pp.sum2, aes(x = kok, y = mean_pp)) +
  geom_ribbon(aes(ymin=lower_pp,ymax=upper_pp),alpha=0.1)+
  geom_line() + 
  theme_bw()+ xlab("Kokanee")+ ylab("Survival")

# convert scaled kokanee back to actual biomass 
# scale subtracts the mean and divides by the sd
# so multifply by sd then add mean to transform back
m<-mean(predict$running_accessible_biomass,na.rm=T)
s<-sd(predict$running_accessible_biomass,na.rm=T)

#head(kok.pp.sum2)
predicted_prob_sum2 <- kok.pp.sum2 %>%
  mutate(step1 = kok*s)%>%
  mutate(step2 = step1+m)%>%
  rename(AvailKok = step2)%>%
  dplyr::select(-step1)%>%
  mutate("Release" = "N")
head(predicted_prob_sum2)


##  plot to view effect of kokanee availability of adult bull trout survival
ggplot(data = predicted_prob_sum2, aes(x =  AvailKok, y = mean_pp)) +
  geom_ribbon(aes(ymin=lower_pp,ymax=upper_pp),alpha=0.1)+
  geom_smooth(method = 'loess', formula = 'y ~ x',color="black",size=1,se=FALSE) + 
  theme_bw()+ 
  #geom_point(data=dat,aes(x=kok,y=survival))+
  labs(x="Available Kokanee Biomass (lbs)",y="Adult Bull Trout Survival") +
  scale_x_continuous(expand=c(0,0),labels = comma)



### plot the effect of kokanee on survival across different fishing regulations
kok.sim.dat <- as.matrix(data.frame(constant = 1,
                                    NumbLake = meanNlake/10000,   
                                    kok = seq(from = min(predict$running_biomass_scale,na.rm=T), 
                                              to = max(predict$running_biomass_scale,na.rm=T), 
                                              length.out=20),
                                    ca.reg.changes=0)) #  significant fishing regs
#head(kok.sim.dat)

# posterior draws
posterior.draws<-as.matrix(as.mcmc(full.mod2))
posterior.draws<-as_tibble(posterior.draws)

#head(posterior.draws)
posterior.draws<-posterior.draws %>%
  dplyr::select(a.int,beta.dens,beta.kok,beta.reg.change) # 
#head(posterior.draws)
posterior.draws<-as.matrix(posterior.draws)


Xb <- t(kok.sim.dat %*% t(posterior.draws)) # multiplying all posterior draws by covariate value
#head(Xb)
Xc<-plogis(Xb)# still need to transform
#head(Xc)
kok.pp <- as.data.frame(Xc)

names(kok.pp) <- kok.sim.dat[, "kok"] # 
#head(kok.pp)

kok.pp.long <- tidyr::gather(kok.pp, key = "kok")
#head(kok.pp.long)
kok.pp.sum <- summarize(group_by(kok.pp.long, kok), 
                        mean_pp = mean(value), 
                        median_pp = median(value), 
                        lower_pp = quantile(value, probs = c(0.025)), 
                        upper_pp = quantile(value, probs = c(0.975)),
                        low10_pp = quantile(value, probs = c(0.1)),
                        high90_pp = quantile(value, probs = c(0.90)),
                        low25_pp = quantile(value, probs = c(0.25)),
                        high75_pp = quantile(value, probs = c(0.75)))

kok.pp.sum$kok <- as.numeric(kok.pp.sum$kok)

ggplot(data = kok.pp.sum, aes(x = kok, y = mean_pp)) +
  geom_ribbon(aes(ymin=lower_pp,ymax=upper_pp),alpha=0.1)+
  geom_line() + 
  theme_bw()+ xlab("Kokanee")+ ylab("Survival")



# convert scaled kokanee to actual biomass
# scale subtracts the mean and divides by the sd
# so multifply by sd then add mean to transform back
m<-mean(predict$running_accessible_biomass,na.rm=T)
s<-sd(predict$running_accessible_biomass,na.rm=T)

head(kok.pp.sum)
predicted_prob_sum <- kok.pp.sum %>%
  mutate(step1 = kok*s)%>%
  mutate(step2 = step1+m)%>%
  rename(AvailKok = step2)%>%
  dplyr::select(-step1)%>%
  mutate("Release" = "Y")

predicted_prob_sum3<-full_join(predicted_prob_sum,predicted_prob_sum2)
predicted_prob_sum3$Release <- as.factor(predicted_prob_sum3$Release)

# plotting
ggplot(data = predicted_prob_sum3, aes(x =  AvailKok, y = mean_pp,group=`Release`,
                                       color=`Release`,linetype=`Release`)) +
  geom_ribbon(aes(ymin=lower_pp,ymax=upper_pp),alpha=0.1)+
  geom_line(size=1.5)+
  theme_bw()+ 
  theme(legend.position = c(0.8,0.2))+
  labs(x="Kokanee Biomass Available (lbs)",y="Adult Bull Trout Survival") +
  scale_x_continuous(expand=c(0,0),labels = comma)+ scale_color_manual(values=c("#999999", "black"))

predicted_prob_sum3 %>% filter(Release=="Y")%>%
  ggplot(aes(x =  AvailKok*0.453592, y = mean_pp)) +
  geom_ribbon(aes(ymin=lower_pp,ymax=upper_pp),alpha=0.1)+
  geom_smooth(method = 'loess', formula = 'y ~ x',color="black",size=1,se=FALSE) + 
  theme_minimal()+ 
  #geom_point(data=dat,aes(x=kok,y=survival))+
  labs(x="Available Kokanee Biomass (kg)",y="Survival") +
  scale_x_continuous(expand=c(0,0),labels = comma)+
  theme(plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
#ggsave("Results/summary/PubFigs/Figure 6.png",width = 8.5, height =8.5, units="cm", dpi=500)


### want to get median increase in adult survival between two regulatory periods
#head(predicted_prob_sum3)
predicted_prob_sum3%>%filter(kok>0 & kok <0.2)
((0.889-0.518)/0.518)*100 # 72% increase 


# percent increase in survival
# only using later half of study, e.g. regulations = yes
predicted_prob_sum3<-predicted_prob_sum3%>%filter(Release=="Y")%>%arrange(kok)
predicted_prob_sum3$median_pp[1] # 75% survival - lowest kokanee biomass
predicted_prob_sum3$median_pp[20] # 97% survival - highest kokanee biomass
((predicted_prob_sum3$median_pp[20]-predicted_prob_sum3$median_pp[1])/predicted_prob_sum3$median_pp[1])*100
# 31% increase 







######################################
### effect of variable bull trout density on adult survival
dd.sim.dat <- as.matrix(data.frame(constant = 1,   
                                   kok = 0, 
                                   NumbLake = seq(from = 0, to = max(Nlake.est)/10000, 
                                                  length.out=20),
                                   ca.reg.changes=0)) # catch and release


# posterior draws
posterior.draws<-as.matrix(as.mcmc(full.mod2))
posterior.draws<-as_tibble(posterior.draws)
posterior.draws2<-posterior.draws %>%
  dplyr::select(a.int,beta.kok,beta.dens,beta.reg.change) # make sure this is in the same order as sim.dat # 
posterior.draws2<-as.matrix(posterior.draws2)


Xb <- t(dd.sim.dat %*% t(posterior.draws2)) # multiplying all posterior draws by each dd row
Xc<-plogis(Xb)# still need to use logit transform
dd.pp <- as.data.frame(Xc)

names(dd.pp) <- dd.sim.dat[, "NumbLake"]

dd.pp.long <- tidyr::gather(dd.pp, key = "NumbLake")
#head(dd.pp.long)
dd.pp.sum <- summarize(group_by(dd.pp.long, NumbLake), 
                       median_pp = median(value), 
                       mean_pp = mean(value),
                       lower_pp = quantile(value, probs = c(0.025)), 
                       upper_pp = quantile(value, probs = c(0.975)),
                       low10_pp = quantile(value, probs = c(0.1)),
                       high90_pp = quantile(value, probs = c(0.90)),
                       low25_pp = quantile(value, probs = c(0.25)),
                       high75_pp = quantile(value, probs = c(0.75)))
#head(dd.pp.sum)

dd.pp.sum$NumbLake <- as.numeric(dd.pp.sum$NumbLake)

# plotting
ggplot(data = dd.pp.sum, aes(x = NumbLake*10000, y = mean_pp)) +
  geom_ribbon(aes(ymin=lower_pp,ymax=upper_pp),alpha=0.1)+
  geom_line(size=1.5) + 
  scale_x_continuous(expand=c(0,0))+
  theme_minimal()+ xlab("Bull Trout in Reservoir (n)")+ ylab("Survival")+
  theme(plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
#ggsave("Results/summary/PubFigs/Figure 7.png",width = 4, height = 4)

max(Nlake.est) #  7051
min(Nlake.est) # 1212
dd.pp.sum$median_pp[2] # 91% survival - uses 2 b/c population was never 0
dd.pp.sum$median_pp[20] # 86% survival
((dd.pp.sum$median_pp[2]-dd.pp.sum$median_pp[20])/dd.pp.sum$median_pp[20])*100
# 7% decrease



#### visualize how survival rate changed over time

## survival rate over time
survival<-full.mod2$BUGSoutput$mean$S.a
surv<-tibble(survival,yrs=seq(1981,2023))

ggplot(data=surv,aes(x=yrs,y=survival))+geom_point()+geom_line()+
  labs(x="Year",y="Adult Bull Trout Survival")+
  geom_hline(yintercept=mean(full.mod2$BUGSoutput$mean$S.a),linetype="dashed") # mean survival
#ggsave("Results/summary/bt.survival.png",width = 4, height = 4)


# avg survival from 1980 - 1995
mean(surv$survival[1:16]) # 44
sd(surv$survival[1:16]) # 0.10

# avg survival from 1998 - 2004
mean(surv$survival[19:25]) # 92
sd(surv$survival[19:25]) # 0.01

# overall mean
mean(survival) # 0.71
sd(survival) # 0.22

min(survival) # 0.30
max(survival) # 0.94

# percent increase from pre-reg to regulatory period 
(( mean(surv$survival[17:43])-mean(surv$survival[1:16]))/mean(surv$survival[1:16]))*100
# 100% increase 



######################################
### relationship between low flow and sub-adult recruitment  ####
low.flow.recruitment2 = bull.river.low$scale_summerLowflow[6:45]

sr.sim.dat <- as.matrix(data.frame(constant = 1,
                                   low.flow.recruitment = seq(from = min(low.flow.recruitment2), to = max(low.flow.recruitment2), 
                                                              length.out=20),
                                   peak.flow = 0)) # set at averages pool.elevation.recruit = 0

# posterior draws
posterior.draws<-as.matrix(as.mcmc(full.mod2))
posterior.draws<-as_tibble(posterior.draws)
posterior.draws4<-posterior.draws %>%
  dplyr::select(alpha_sr,beta.low.flow, beta.peak.winter) # make sure this is in the same order as sim.dat #  ,beta.pool.recruitment
#head(posterior.draws4)
posterior.draws4<-as.matrix(posterior.draws4)

Xb <- t(sr.sim.dat %*% t(posterior.draws4)) # multiplying all posterior draws by each row
#head(Xb)
Xc<-exp(Xb)# transform
#head(Xc)
sr.pp <- as.data.frame(Xc)

names(sr.pp) <- sr.sim.dat[, "low.flow.recruitment"]
#head(sr.pp)
## these are draws for alpha - beta (still need to add log(biomass))

sr.pp.long <- tidyr::gather(sr.pp, key = "low.flow.recruitment")
sr.pp.long$low.flow.recruitment<-as.numeric(sr.pp.long$low.flow.recruitment)
sr.pp.long <- sr.pp.long %>% 
  mutate(value =  exp( log(mean.biomass) + log(value) - full.mod2$BUGSoutput$mean$beta_sr*mean.biomass) ) # including mean beta estimate

sr.pp.sum <- summarize(group_by(sr.pp.long, low.flow.recruitment), 
                       median_pp = median(value), 
                       mean_pp = mean(value),
                       lower_pp = quantile(value, probs = c(0.025)), 
                       upper_pp = quantile(value, probs = c(0.975)),
                       low10_pp = quantile(value, probs = c(0.1)),
                       high90_pp = quantile(value, probs = c(0.90)),
                       low25_pp = quantile(value, probs = c(0.25)),
                       high75_pp = quantile(value, probs = c(0.75)))

sr.pp.sum$low.flow.recruitment <- as.numeric(sr.pp.sum$low.flow.recruitment)

# converting back to actual discharge values
m<-mean(bull.river.low$summerLowflow,na.rm=T)
s<-sd(bull.river.low$summerLowflow,na.rm=T)

#head(sr.pp.sum)
sr.pp.sum <- sr.pp.sum %>%
  mutate(step1 = low.flow.recruitment*s)%>%
  mutate(step2 = step1+m)%>%
  rename(LowFlow = step2)%>%
  dplyr::select(-step1)
#head(sr.pp.sum) # discharge in cubic meters per second

# plotting
ggplot(data = sr.pp.sum, aes(x = LowFlow, y = mean_pp)) +
  geom_ribbon(aes(ymin=lower_pp,ymax=upper_pp),alpha=0.1)+
  geom_line(size=1.5) + 
  scale_x_continuous(expand=c(0,0),breaks=c(8,10,12,14,16))+
  theme_minimal()+
  xlab(bquote('Low Flow ('~m^3~s^-1*' )'))+
  labs(y="Sub-Adult Recruits")+
  theme(plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
#ggsave("Results/summary/PubFigs/Figure 4.png",width = 8.5, height =8.5, units="cm", dpi=500)


# effect size 
sr.pp.sum$median_pp[1]-sr.pp.sum$median_pp[20] # - 1085
sr.pp.sum$lower_pp[1]-sr.pp.sum$lower_pp[20] # - 544
sr.pp.sum$upper_pp[1]-sr.pp.sum$upper_pp[20] # - 1961

((sr.pp.sum$median_pp[20]-sr.pp.sum$median_pp[1])/sr.pp.sum$median_pp[1])*100 #≈ 87% increase
((sr.pp.sum$lower_pp[20]-sr.pp.sum$lower_pp[1])/sr.pp.sum$lower_pp[1])*100 #≈ 77% increase
((sr.pp.sum$upper_pp[20]-sr.pp.sum$upper_pp[1])/sr.pp.sum$upper_pp[1])*100 #≈ 92% increase






##### use bayesplot to view all beta parameter estimates
pe <- as.array(full.mod2$BUGSoutput$sims.array)
#dimnames(pe)

bayesplot::mcmc_intervals(pe, pars= c("beta.dens",
                                      "beta.kok",
                                      "beta.low.flow",
                                      "beta.peak.winter"),
                          prob = 0.75, 
                          prob_outer = 0.9) + 
  scale_y_discrete(labels = c("beta.dens" = expression(paste(beta[S-density])), 
                              "beta.kok" = expression(paste(beta[S-kok])),
                              "beta.low.flow" = expression(paste(beta[R-lowQ])),
                              "beta.peak.winter" = expression(paste(beta[R-peakQ])))) +
  labs( x = "Credible Interval", y = "Parameter") +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype="dashed") +
  theme_bw() + 
  theme(axis.text=element_text(size=12), # make text larger in plots
        axis.title=element_text(size=11,face="bold")) 
# R stands for recruitment, S stands for survival
#ggsave("Results/summary/bt.beta.significance.png",width = 6, height = 6)






####### lambda over time
lambdaSA<-full.mod2$BUGSoutput$mean$lambda
lambdaA<-full.mod2$BUGSoutput$mean$lambdaA
yr <- seq(1981,2023,1)

l.tib<-tibble(yr,lambdaA,lambdaSA,Nlake=as.integer(Nlake.est[1:43]),
              Nadult = as.integer(Nadult.est[1:43]))

ggplot(data=l.tib, aes(x=yr, y=lambdaA))+geom_point()+geom_line()+ # black = adult
  geom_point(aes(x=yr, y=lambdaSA),color="grey")+ # blue = subadult
  geom_line(aes(x=yr, y=lambdaSA),color="grey")+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme_classic()+labs(x="Year",y="Lambda")+
  annotate("text", x = 2000, y = 2.5, label = "Sub-Adult",color="grey")+
  annotate("text", x = 1987, y = 0.6, label = "Adult")
#ggsave("Results/summary/bt.lambda.png",width = 4, height = 4)

# average lambda over time (with sd)
mean(full.mod2$BUGSoutput$mean$lambda) # 1.12
sd(full.mod2$BUGSoutput$mean$lambda) # 0.41

mean(full.mod2$BUGSoutput$mean$lambdaA) # 1.06
sd(full.mod2$BUGSoutput$mean$lambdaA)  # 0.19





######### population size estimates ###################
qa<-c(full.mod2$BUGSoutput$mean$q.lake) # catchability coefficients
qsa<-c(full.mod2$BUGSoutput$mean$q.lake.sa)

grave<-full.mod2$BUGSoutput$mean$q[1] # proportion of the population from grave creek
grave.low<-quantile(full.mod2$BUGSoutput$sims.list$q[1],0.025)
grave.high<-quantile(full.mod2$BUGSoutput$sims.list$q[1],0.975)

wig<-full.mod2$BUGSoutput$mean$q[2]# proportion of the population from wigwam river

yr <- seq(1980,2023,1)

#sdNadults<-full.mod2$BUGSoutput$sd$Nadults
#sdNsubadults<-full.mod2$BUGSoutput$sd$Nsubadults
Nadult.low<-full.mod2$BUGSoutput$summary[1:44,3] # 2.5%
Nadult.high<-full.mod2$BUGSoutput$summary[1:44,7] # 97.5%
Nsubadult.low<-full.mod2$BUGSoutput$summary[45:88,3] # 2.5%
Nsubadult.high<-full.mod2$BUGSoutput$summary[45:88,7] # 97.5%

tot.wt<-weight.yr$wt[5:48] # 1980-2023

out.result<-tibble(yr,Nsubadult.est=Nsubadult.est[1:44],Nadult.est=Nadult.est[1:44],
                   tot.wt, Nadult.low, Nadult.high, # sdNadults,sdNsubadults=sdNsubadults[1:44],
                   Nsubadult.low, Nsubadult.high)
out.result <- out.result %>% 
  #mutate(Nadult.low = Nadult.est-sdNadults)%>%
  #mutate(Nadult.high = Nadult.est+sdNadults)%>%
  #mutate(Nsubadult.low = Nsubadult.est-sdNsubadults)%>%
  #mutate(Nsubadult.high = Nsubadult.est+sdNsubadults)%>%
  mutate(Total.low = Nadult.low+Nsubadult.low)%>%
  mutate(Total.high = Nadult.high+Nsubadult.high)%>%
  mutate(Grave_A = Nadult.est*grave*0.5)%>%
  mutate(Grave_Alow = Nadult.low*grave.low*0.5)%>%
  mutate(Grave_Ahigh = Nadult.high*grave.high*0.5)%>%
  
  mutate(Wigwam_A = Nadult.est*wig*0.5)%>%
  mutate(SA_reduc = Nsubadult.est*qsa)%>%
  mutate(A_reduc = Nadult.est*qa)%>%
  mutate(adult_biomass = Nadult.est*tot.wt)



btpop<-ggplot(data=out.result, aes(x=yr, y=Nsubadult.est+Nadult.est))+
  geom_point()+geom_line()+
  geom_point(aes(x=yr, y=Nsubadult.est),color="darkgrey")+
  geom_line(aes(x=yr, y=Nsubadult.est),color="darkgrey")+
  geom_point(aes(x=yr, y=Nadult.est),color="blue")+
  geom_line(aes(x=yr, y=Nadult.est),color="blue",linetype="dashed") +
  geom_ribbon(aes(ymin=Nsubadult.low, ymax=Nsubadult.high),alpha=0.1)+
  geom_ribbon(aes(ymin=Nadult.low, ymax=Nadult.high),alpha=0.1,)+
  geom_ribbon(aes(ymin=Total.low, ymax=Total.high),alpha=0.1,)+
  #geom_point(data=bt.harvest[5:48,], aes(x=Season,y=BullTroutHarvested),
  # color="orange")+
  #geom_line(data=bt.harvest[5:48,],aes(x=Season,y=BullTroutHarvested),color="orange")+
  labs(x="Year",y="Bull Trout Population Size")+
  theme_classic()+ 
  #annotate("text", x = 2018, y = 600, label = "Harvest",color="orange")+
  annotate("text", x = 2017, y = 400, label = "Sub-Adult",color="darkgrey")+
  annotate("text", x = 1999, y = 500, label = "Adult",color="blue")+
  annotate("text", x = 1987, y = 3800, label = "Total",color="black")+
  scale_x_continuous(expand = c(0, 0),breaks = c(1980,1985,1990,1995,2000,2005,2010,2015,2020))
btpop


# plotting adult biomass over time
ggplot(data=out.result, aes(x=yr,y=adult_biomass))+
  geom_point()+geom_line()


# adult population size only
ggplot(data=out.result,aes(x=yr, y=Nadult.est))+
  geom_point(color="blue")+
  geom_line(color="blue") +
  geom_ribbon(aes(ymin=Nadult.low, ymax=Nadult.high),alpha=0.1)+
  
  labs(x="Year",y="Bull Trout Population")+theme_classic()+ 
  scale_x_continuous(expand = c(0, 0))


# make better plot, go from wide to long for legends
out.result.long<-  out.result%>% dplyr::select(yr,Nsubadult.est,Nadult.est) %>%
  rename("Sub-Adult" = "Nsubadult.est")%>%
  rename("Adult" = "Nadult.est")%>%
  tidyr::gather(key=population, value=N,-yr )

out.result.long.low<-  out.result%>% dplyr::select(yr,Nsubadult.low,Nadult.low) %>%
  rename("Sub-Adult" = "Nsubadult.low")%>%
  rename("Adult" = "Nadult.low")%>%
  tidyr::gather(key=population, value=N.low,-yr )

out.result.long.high<-  out.result%>% dplyr::select(yr,Nsubadult.high,Nadult.high) %>%
  rename("Sub-Adult" = "Nsubadult.high")%>%
  rename("Adult" = "Nadult.high")%>%
  tidyr::gather(key=population, value=N.high,-yr )

out.result.long<-left_join(out.result.long,out.result.long.low, by=c("yr","population"))
out.result.long<-left_join(out.result.long,out.result.long.high, by=c("yr","population"))


bt<-ggplot(data=out.result.long, aes(x=yr, y=N,group=population,color=population,linetype=population))+
  geom_point(size=1.5)+
  geom_line(size=1.1)+
  geom_ribbon(aes(ymin=N.low,ymax=N.high),alpha=0.1)+ # ,color=NA
  scale_x_continuous(expand = c(0, 0))+ # ,breaks = c(1980,1985,1990,1995,2000,2005,2010,2015,2020)
  scale_y_continuous(expand=c(0.06,0))+
  theme_minimal()+
  scale_color_grey()+
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8))+
  labs(x="Year",y="Bull Trout (n)")+
  theme(plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
bt
#ggsave("Results/summary/PubFigs/Figure 1.png",width = 12, height =8.5, units="cm", dpi=500)


# population size in 1980s
#head(out.result)
mean(out.result$Nadult.est[1:10]) # 532
max(out.result$Nadult.est)# 4902
min(out.result$Nadult.est) # 371

mean(out.result$Nadult.est[32:44]) # 2918
min(out.result$Nadult.est[32:44])# 2639
max(out.result$Nadult.est[32:44])# 3230



####################### compare known harvest on lake koocanusa as proportion of mortality
bt.harv = bt.harvest$BullTroutHarvested[6:48]
Harvest<-bt.harv/Nadult.est[2:44]
Natural<-1-survival
Natural<-1-surv$survival
mort.og<-tibble(Natural, Harvest, total=Natural+Harvest, 
                proportion.harv= Harvest/total, yrs=seq(1980,2022))

mort<-tibble(mort=Natural, yrs=seq(1980,2022),type="Natural")
mort2<-tibble(mort=Harvest, yrs=seq(1980,2022),type="Known Harvest")

mort2$mort[25] # 2004 year when catch was high
mort2$mort[38]# 2017 when catch was low 
mort2$mort[25]/(mort$mort[25]+mort2$mort[25]) # 65% of mortality was due to known harvest
mort2$mort[38]/(mort$mort[38]+mort2$mort[38]) # only 15% of mortality was due to known harvest

mortfull<-full_join(mort,mort2)
#head(mortfull)

ggplot(data=mort.og, aes(x=yrs))+
  geom_area(aes(y=Natural+Harvest, fill="Total"))+ # maybe switch to total
  geom_area(aes(y=Harvest,fill="Harvest"))+
  scale_y_continuous(expand=c(0,0),minor_breaks = seq(0,0.6,0.1))+
  scale_x_continuous(expand=c(0,0),breaks=c(1980,1990,2000,2010,2020))+
  scale_fill_grey()+
  theme_minimal()+
  geom_vline(xintercept = 1996,linetype=2)+
  labs(x="Year",y="Mortality Rate")+  
  theme(legend.title = element_blank(),
        legend.position = c(0.8,0.8),
        #panel.grid.minor = element_line(size=0.5),
        #plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
        plot.background = element_rect(fill = 'white', colour = 'white'))
#axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
#axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
ggsave("Results/summary/PubFigs/Figure 5.png",width = 8, height = 6,units="cm", dpi=500) 

mort.og2 <- mort.og %>%
  mutate(Harvest = if_else(Harvest >0,Harvest,NA))%>%
  mutate(Natural  = if_else(Harvest>0,Natural,NA))
mort.og2$Harvest[33]<-0
mort.og2$Harvest[34]<-0
mort.og2$Harvest[35]<-0
mort.og2$Harvest[36]<-0
mort.og2$Natural[33]<-mort.og2$total[33]
mort.og2$Natural[34]<-mort.og2$total[34]
mort.og2$Natural[35]<-mort.og2$total[35]
mort.og2$Natural[36]<-mort.og2$total[36]

# different visualization of previous
ggplot(data=mort.og2, aes(x=yrs))+
  geom_area(aes(y=Natural+Harvest, fill="Natural"))+
  geom_area(aes(y=Harvest,fill="Harvest"))+
  geom_line(aes(y=total))+
  #geom_point(aes(y=total))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0),breaks=c(1980,1985,1990,1995,2000,2005,2010,2015,2020))+
  scale_fill_grey()+
  theme_bw()+
  labs(x="Year",y="Mortality Rate")+  
  theme(legend.title = element_blank(),
        legend.position = c(0.8,0.6))+
  annotate("text", x = 1990, y = 0.6, label = "Total",color="black")
#ggsave("Results/summary/mortality.rates.png",width = 6, height = 4) 







#### what unexplained mortality is left after dd, kokanee and regulation change
beta.dens<-full.mod2$BUGSoutput$mean$beta.dens
beta.kok<-full.mod2$BUGSoutput$mean$beta.kok
beta.reg<-full.mod2$BUGSoutput$mean$beta.reg.change
Nsubadult.est
Nadult.est
NumbLake <- Nsubadult.est+Nadult.est
a.int<-full.mod2$BUGSoutput$mean$a.int
kok = predict$running_biomass_scale[1:43] # 1980-2022
ca.reg.changes = c(rep(1,16),rep(0,27))

S.a.step3 <- a.int +  beta.kok*kok + beta.dens*(NumbLake[1:43]/10000) + beta.reg*ca.reg.changes
mean.S.a <- 1/(1+exp(-S.a.step3)) # same as logit

a.int.low<-quantile(full.mod2$BUGSoutput$sims.list$a.int,0.025)
beta.kok.low<-quantile(full.mod2$BUGSoutput$sims.list$beta.kok,0.025)
beta.dens.low<-quantile(full.mod2$BUGSoutput$sims.list$beta.dens,0.025)
beta.reg.low<-quantile(full.mod2$BUGSoutput$sims.list$beta.reg.change,0.025)

S.a.step3.low <- a.int + beta.kok.low*kok + beta.dens.low*(NumbLake[1:43]/10000) + beta.reg.low*ca.reg.changes
mean.S.a.low <- 1/(1+exp(-S.a.step3.low)) 

a.int.high<-quantile(full.mod2$BUGSoutput$sims.list$a.int,0.975)
beta.kok.high<-quantile(full.mod2$BUGSoutput$sims.list$beta.kok,0.975)
beta.dens.high<-quantile(full.mod2$BUGSoutput$sims.list$beta.dens,0.975)
beta.reg.high<-quantile(full.mod2$BUGSoutput$sims.list$beta.reg.change,0.975)

S.a.step3.high <- a.int +  beta.kok.high*kok + beta.dens.high*(NumbLake[1:43]/10000) + beta.reg.high*ca.reg.changes
mean.S.a.high <- 1/(1+exp(-S.a.step3.high)) 

S.a<-full.mod2$BUGSoutput$mean$S.a # total
S.a.low<-full.mod2$BUGSoutput$summary[89:131,3] # 2.5%
S.a.high<-full.mod2$BUGSoutput$summary[89:131,7] # 97.5%

# diff btw known and unknown cause of variation in survival
diff <- S.a-mean.S.a
diff.low <- S.a.low-mean.S.a.high
diff.high <- S.a.high-mean.S.a.low

diffs<-tibble(yrs=seq(1981,2023),mean.S.a.low, mean.S.a, mean.S.a.high, S.a.high, S.a, S.a.low,
              diff.low,diff,diff.high)

ggplot(diffs, aes(x=yrs,y=diff))+geom_point()+geom_line()+
  geom_ribbon(aes(ymin=diff.low,ymax=diff.high),alpha=0.1)+
  geom_hline(aes(yintercept=0))+
  labs(x="Year",y="Difference in Survival Rate")+theme_bw()
# black line negative = unknown variation reduced survival rate
# black line positive  = unknown variation increased survival rate
#ggsave("Results/summary/unexplained.survival.png",width = 6, height = 4)



## summary stats on kokanee biomass available 
max(predict$running_accessible_biomass)*0.453592 # 763545 kg
mean(predict$running_accessible_biomass)*0.453592 # 386693
min(predict$running_accessible_biomass)*0.453592 # 12512

mean(predict$running_accessible_biomass_low)*0.453592 # 221544
mean(predict$running_accessible_biomass_high)*0.453592 # 902550


############## summary stats from raw data ###
# bull trout cpue
head(adult.bull.net.sums)
cpue.sums <- adult.bull.net.sums %>% filter(YYYY>1979) %>% group_by(YYYY)%>%
  tidyr::drop_na()%>%summarise(median=median(cpue),max=max(cpue))
head(cpue.sums)
median(cpue.sums$median) # 1

head(subadult.bull.net.sums)
sa.cpue.sums <- subadult.bull.net.sums %>% filter(YYYY>1979) %>% group_by(YYYY)%>%tidyr::drop_na()%>%
  summarise(median=median(cpue),max=max(cpue))
median(cpue.sums$median) # 1

### redd counts
redds.sums<-bt.redds.combined %>% group_by(LocalPopulation) %>% filter(Year>1979) %>% tidyr::drop_na()%>%
  summarise(mean = mean(Redds), sd=sd(Redds))
# view df for sums

# known harvest data
bt.harvest %>% filter(BullTroutHarvested>0)%>%summarise(mean = mean(BullTroutHarvested), sd = sd(BullTroutHarvested))

# length and weight data 
bull.raw.net.dat %>%
  mutate(Stage = if_else(Length_mm >= 500,"A","SA"))%>% 
  group_by(Stage)%>%
  dplyr::select(YYYY,Stage, Length_mm, Weight_lbs)%>%
  summarise(mean.wt = mean(Weight_lbs, na.rm=T), sd.wt = sd(Weight_lbs, na.rm=T),
            mean.length = mean(Length_mm, na.rm=T), sd.length = sd(Length_mm, na.rm=T))

# co-variates for discharge
bull.q.winter.peak %>% filter(Year > 1979) %>% summarise(mean = mean(maxWinterq), sd = sd(maxWinterq))
bull.river.low %>% filter(Year > 1979) %>% tidyr::drop_na(summerLowflow)%>% 
  summarise(mean = mean(summerLowflow), sd = sd(summerLowflow))



