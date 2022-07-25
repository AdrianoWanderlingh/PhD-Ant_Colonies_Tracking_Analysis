####################################################################################
#### THIS SCRIPT CONTAINS:
####
####################################################################################
library(FortMyrmidon)
library(ggplot2)
library(lubridate)
library(plotrix)
library(scales)
library(car)
library(lme4)
library(Hmisc)
library(viridis)
library(stringr)
library(dplyr)

WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data"

SpaceUse <- read.table(paste(WORKDIR,"/Data/AntTasks_SpaceUse.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

## Colony sizes info
SpaceUse$treat <- str_sub(SpaceUse$colony,-2,-1)
Mean_ants_exp <- aggregate(colony_size ~ treat, FUN=mean, na.action=na.omit, SpaceUse)
SD_ants_exp <- aggregate(colony_size ~ treat, FUN=sd, na.action=na.omit, SpaceUse); colnames(SD_ants_exp) [match("colony_size",colnames(SD_ants_exp))] <- "SD_received"
Colony.size    <-  plyr::join(x=Mean_ants_exp, y=SD_ants_exp, type = "full", match = "all")
data.frame(treat=Colony.size$treat, Colony_size=sprintf("%.2f \U00B1 %.2f",Colony.size$colony_size,Colony.size$SD_received))


#remove queens
SpaceUse <- SpaceUse[which(SpaceUse$IsQueen=="no"),]
SpaceUse <- SpaceUse[which(!is.na(SpaceUse$Exposed)),]

N_exp <- as.data.frame(table(SpaceUse$Exposed))

#### PRODUCE BOXPLOT OF USED SPACE FOR NURSES DELTA TIME INSIDE
ggplot(SpaceUse, aes(x=treatment, y=delta_time_inside,color=AntTask))+
 # geom_jitter(aes(fill = AntTask))+
  #geom_point(position=position_jitterdodge())+
  geom_boxplot()+
  facet_wrap(~Exposed) +
  labs(title= "SpaceUse",
       subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],N_exp[1,1],":",N_exp[2,2]))


