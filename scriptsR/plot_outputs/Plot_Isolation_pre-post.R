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

SpaceUse <- read.table(paste(DATADIR,"/AntTasks_SpaceUse_july2022.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

SpaceUse$REP_treat <- SpaceUse$colony


## Colony sizes info
SpaceUse$treat <- str_sub(SpaceUse$colony,-2,-1)
Mean_ants_exp <- aggregate(colony_size ~ treat, FUN=mean, na.action=na.omit, SpaceUse)
SD_ants_exp <- aggregate(colony_size ~ treat, FUN=sd, na.action=na.omit, SpaceUse); colnames(SD_ants_exp) [match("colony_size",colnames(SD_ants_exp))] <- "SD_received"
Colony.size    <-  plyr::join(x=Mean_ants_exp, y=SD_ants_exp, type = "full", match = "all")
data.frame(treat=Colony.size$treat, Colony_size=sprintf("%.2f \U00B1 %.2f",Colony.size$colony_size,Colony.size$SD_received))


#select some SpaceUse cols
SpaceUse <- SpaceUse[,c("REP_treat","antID","AntTask","delta_time_inside")]

### Get info from metadata
SpaceMeta <- list(SpaceUse,metadata)
SpaceMeta <- Reduce(function(x, y) merge(x, y, all=TRUE), SpaceMeta)
#clean
SpaceMeta <- SpaceMeta[which(!is.na(SpaceMeta$antID)),]

N_exp <- as.data.frame(table(SpaceMeta$Exposed))

#LEVELS RENAMING
# Rename by name
SpaceMeta$treatment <- as.factor(SpaceMeta$treatment)
levels(SpaceMeta$treatment)[levels(SpaceMeta$treatment)=="BS"] <- "Big Sham"
levels(SpaceMeta$treatment)[levels(SpaceMeta$treatment)=="BP"] <- "Big Pathogen"
levels(SpaceMeta$treatment)[levels(SpaceMeta$treatment)=="SS"] <- "Small Sham"
levels(SpaceMeta$treatment)[levels(SpaceMeta$treatment)=="SP"] <- "Small Pathogen"

# Rename levels
N_exp$Var1 <- as.factor(N_exp$Var1)
levels(N_exp$Var1)[levels(N_exp$Var1)==TRUE] <- "exposed"
levels(N_exp$Var1)[levels(N_exp$Var1)==FALSE] <- "non-exposed"


#CHECK IF THERE ARE ANY EXPOSED NURSES LABELED AS FORAGERS
SpaceMeta[which(SpaceMeta$AntTask=="forager" & SpaceMeta$Exposed==TRUE),] #only 3 individuals, no issues in code!
#change task for the sake of plotting
SpaceMeta[which(SpaceMeta$AntTask=="forager" & SpaceMeta$Exposed==TRUE),"AntTask"] <- "nurse"

#aggregg for boxplot - ants means, by colony
#mean by colony, AntTask, exposed, treatment
Mean_SpaceMeta <- aggregate(delta_time_inside ~ REP_treat + AntTask + Exposed + treatment, FUN=mean, na.rm=T, na.action=na.pass, SpaceMeta)
SD_SpaceMeta <- aggregate(delta_time_inside ~ REP_treat + AntTask + Exposed + treatment, FUN=sd, na.rm=T, na.action=na.pass, SpaceMeta); colnames(SD_SpaceMeta) [match("delta_time_inside",colnames(SD_SpaceMeta))] <- "SD_delta_time_inside"

#merge dfs
Mean_SpaceMeta <- plyr::join (x = Mean_SpaceMeta , y=SD_SpaceMeta, type = "right", match = "all")
# Rename levels
Mean_SpaceMeta$Exposed <- as.factor(Mean_SpaceMeta$Exposed)
levels(Mean_SpaceMeta$Exposed)[levels(Mean_SpaceMeta$Exposed)==TRUE] <- "exposed"
levels(Mean_SpaceMeta$Exposed)[levels(Mean_SpaceMeta$Exposed)==FALSE] <- "non-exposed"


### for barplot - colony means
Mean_SpaceMeta_Rep <- aggregate(delta_time_inside ~ AntTask + Exposed + treatment, FUN=mean, na.rm=T, na.action=na.pass, Mean_SpaceMeta)
SE_SpaceMeta_Rep   <- aggregate(delta_time_inside ~ AntTask + Exposed + treatment, FUN=std.error, na.rm=T, na.action=na.pass, Mean_SpaceMeta); colnames(SE_SpaceMeta_Rep) [match("delta_time_inside",colnames(SE_SpaceMeta_Rep))] <- "SE_delta_time_inside"
#merge dfs
Mean_SpaceMeta_Rep <- plyr::join (x = Mean_SpaceMeta_Rep , y=SE_SpaceMeta_Rep, type = "right", match = "all")

##PLOTTING
STYLE_NOVIR <- list(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
                    theme_bw(),
                    scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)


#### BOXPLOT OF USED SPACE FOR DELTA TIME INSIDE
ggplot(Mean_SpaceMeta, aes(x=treatment, y=delta_time_inside,color=AntTask))+
  # geom_errorbar( aes(x=treatment,ymin=delta_time_inside-SD_delta_time_inside, ymax=delta_time_inside+SD_delta_time_inside),position=position_dodge2(width=0.8, preserve = "single"))+
  # geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
 # geom_jitter(aes(fill = AntTask))+
  geom_boxplot(position = position_dodge(width = 0.8, preserve = "single"))+
  #geom_point(size=0.8,position=position_dodge2(width = 0.8, preserve = "single"))+
  facet_wrap(~Exposed) +
  STYLE_NOVIR +
  labs(title= "Space Use, ant means by colony",
       subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))




#### BARPLOT OF USED SPACE FOR DELTA TIME INSIDE
ggplot(Mean_SpaceMeta_Rep, aes(x=treatment, y=delta_time_inside,fill=AntTask))+
  geom_errorbar( aes(x=treatment,ymin=delta_time_inside-SE_delta_time_inside, ymax=delta_time_inside+SE_delta_time_inside),position=position_dodge2(width=0.8, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  facet_wrap(~Exposed) +
  STYLE_NOVIR +
  labs(title= "Space Use, colony means w/ SE",
       subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))

