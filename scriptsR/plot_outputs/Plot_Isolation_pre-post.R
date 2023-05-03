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
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data"

#SpaceUse <- read.table(paste(DATADIR,"/AntTasks_SpaceUse_july2022.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021_2023-02-27.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

#SpaceUse$REP_treat <- SpaceUse$colony
#SpaceUse$delta_time_inside <- NULL
#SpaceUse$delta_time_outside <-SpaceUse$prop_inside_24hPRE -  SpaceUse$prop_inside_24hPOST


## Colony sizes info
SpaceUse$size_treat <- str_sub(SpaceUse$colony,-2,-1)
Mean_ants_exp <- aggregate(colony_size ~ size_treat, FUN=mean, na.action=na.omit, SpaceUse)
SD_ants_exp <- aggregate(colony_size ~ size_treat, FUN=sd, na.action=na.omit, SpaceUse); colnames(SD_ants_exp) [match("colony_size",colnames(SD_ants_exp))] <- "SD_received"
Colony.size    <-  plyr::join(x=Mean_ants_exp, y=SD_ants_exp, type = "full", match = "all")
data.frame(size_treat=Colony.size$size_treat, Colony_size=sprintf("%.2f \U00B1 %.2f",Colony.size$colony_size,Colony.size$SD_received))


#select some SpaceUse cols
SpaceUse <- SpaceUse[,c("REP_treat","antID","delta_time_outside")] #remove the ant_task as it is old comapred to the one in metadata

### Get info from metadata
SpaceMeta <- list(SpaceUse,metadata)
SpaceMeta <- Reduce(function(x, y) merge(x, y, all=TRUE), SpaceMeta)
#clean
SpaceMeta <- SpaceMeta[which(!is.na(SpaceMeta$antID)),]

N_exp <- as.data.frame(table(SpaceMeta$Exposed))

#LEVELS RENAMING
# Rename by name
SpaceMeta$size_treat <- as.factor(SpaceMeta$size_treat)
levels(SpaceMeta$size_treat)[levels(SpaceMeta$size_treat)=="BS"] <- "Big Sham"
levels(SpaceMeta$size_treat)[levels(SpaceMeta$size_treat)=="BP"] <- "Big Pathogen"
levels(SpaceMeta$size_treat)[levels(SpaceMeta$size_treat)=="SS"] <- "Small Sham"
levels(SpaceMeta$size_treat)[levels(SpaceMeta$size_treat)=="SP"] <- "Small Pathogen"

#remove queens and dead ants
SpaceMeta <- SpaceMeta[which(SpaceMeta$AntTask!= "queen"),]
SpaceMeta <- SpaceMeta[which(SpaceMeta$IsAlive== TRUE),]


# Rename levels
N_exp$Var1 <- as.factor(N_exp$Var1)
levels(N_exp$Var1)[levels(N_exp$Var1)==TRUE] <- "treated"
levels(N_exp$Var1)[levels(N_exp$Var1)==FALSE] <- "untreated"


#CHECK IF THERE ARE ANY EXPOSED NURSES LABELED AS FORAGERS
SpaceMeta[which(SpaceMeta$AntTask=="forager" & SpaceMeta$Exposed==TRUE),] #only 3 individuals, no issues in code!
#change task for the sake of plotting
SpaceMeta[which(SpaceMeta$AntTask=="forager" & SpaceMeta$Exposed==TRUE),"AntTask"] <- "nurse"

#aggregg for boxplot - ants means, by colony
#mean by colony, AntTask, exposed, size_treat
Mean_SpaceMeta <- aggregate(delta_time_outside ~ REP_treat + AntTask + Exposed + size_treat, FUN=mean, na.rm=T, na.action=na.pass, SpaceMeta)
SD_SpaceMeta <- aggregate(delta_time_outside ~ REP_treat + AntTask + Exposed + size_treat, FUN=sd, na.rm=T, na.action=na.pass, SpaceMeta); colnames(SD_SpaceMeta) [match("delta_time_outside",colnames(SD_SpaceMeta))] <- "SD_delta_time_outside"

#merge dfs
Mean_SpaceMeta <- plyr::join (x = Mean_SpaceMeta , y=SD_SpaceMeta, type = "right", match = "all")
# Rename levels
Mean_SpaceMeta$Exposed <- as.factor(Mean_SpaceMeta$Exposed)
levels(Mean_SpaceMeta$Exposed)[levels(Mean_SpaceMeta$Exposed)==TRUE] <- "treated"
levels(Mean_SpaceMeta$Exposed)[levels(Mean_SpaceMeta$Exposed)==FALSE] <- "untreated"

Mean_SpaceMeta$Status <- paste(Mean_SpaceMeta$Exposed, Mean_SpaceMeta$AntTask)

### for barplot - colony means
Mean_SpaceMeta_Rep <- aggregate(delta_time_outside ~ Status + size_treat, FUN=mean, na.rm=T, na.action=na.pass, Mean_SpaceMeta)
SE_SpaceMeta_Rep   <- aggregate(delta_time_outside ~ Status + size_treat, FUN=std.error, na.rm=T, na.action=na.pass, Mean_SpaceMeta); colnames(SE_SpaceMeta_Rep) [match("delta_time_outside",colnames(SE_SpaceMeta_Rep))] <- "SE_delta_time_outside"
#merge dfs
Mean_SpaceMeta_Rep <- plyr::join (x = Mean_SpaceMeta_Rep , y=SE_SpaceMeta_Rep, type = "right", match = "all")

##PLOTTING
STYLE_NOVIR <- list(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
                    theme_bw(),
                    scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)


#### BOXPLOT OF USED SPACE FOR DELTA TIME INSIDE
ggplot(Mean_SpaceMeta, aes(x=size_treat, y=delta_time_outside,color=Status))+
  # geom_errorbar( aes(x=size_treat,ymin=delta_time_outside-SD_delta_time_outside, ymax=delta_time_outside+SD_delta_time_outside),position=position_dodge2(width=0.8, preserve = "single"))+
  # geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
 # geom_jitter(aes(fill = Status))+
  geom_boxplot(position = position_dodge(width = 0.8, preserve = "single"))+
  #geom_point(size=0.8,position=position_dodge2(width = 0.8, preserve = "single"))+
  #facet_wrap(~Exposed) +
  STYLE_NOVIR +
  labs(title= "Space Use, ant means by colony",
       subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))






# Define custom color palette
custom_palette <- c("treated nurse" = "#333333", "untreated forager" = "#FDE725", "untreated nurse" = "#1F9E89")

# Create the boxplot
ggplot(Mean_SpaceMeta, aes(x=size_treat, y=delta_time_outside, color=Status)) +
  geom_boxplot(position = position_dodge(width = 0.8, preserve = "single")) +
  scale_color_manual(values = custom_palette) + # Apply custom color palette
  STYLE_NOVIR +
  labs(x="")
  #labs(title= "Space Use, ant means by colony",
  #     subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))



#### BARPLOT OF USED SPACE FOR DELTA TIME INSIDE
ggplot(Mean_SpaceMeta_Rep, aes(x=size_treat, y=delta_time_outside,fill=Status))+
  geom_errorbar( aes(x=size_treat,ymin=delta_time_outside-SE_delta_time_outside, ymax=delta_time_outside+SE_delta_time_outside),position=position_dodge2(width=0.8, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  #facet_wrap(~Exposed) +
  scale_fill_manual(values = custom_palette) + # Apply custom color palette
  STYLE_NOVIR +
  labs(x="")
  #+
  #labs(title= "Space Use, colony means w/ SE",
  #     subtitle=paste("N",N_exp[1,1],":",N_exp[1,2],";","N",N_exp[2,1],":",N_exp[2,2]))

