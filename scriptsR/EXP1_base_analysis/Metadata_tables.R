##################################################################
################## METADATA CHECKS AND TABLES ####################
##################################################################

### tables from metadata information
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

USER <- "supercompAdriano"

if (USER=="supercompAdriano") {
  WORKDIR <- "/media/cf19810/DISK4/ADRIANO" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  #SCRIPTDIR <- "/media/cf19810/DISK4/EXP1_base_analysis/EXP1_analysis scripts"
}

###### LOAD METADATA
metadata_present <- read.table(file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt"),header=T,stringsAsFactors = F, sep=",")

##################################################################
################## THRESHOLD FOR ANT_TASK ########################
##################################################################

# calculate AntTasks with different prop_time_outside limits
MAX <- c(0.01, 0.02, 0.03, 0.04)

#create new column per each prop_time_outside limit
for (i in seq_along(MAX)) { 
metadata_present <- metadata_present %>% mutate(!!paste0("AntTask", MAX[i]*100, "perc") := NA)
}

#list of plots
plots <- list()
#loop through max  prop_time_outside limits 
for (i in seq_along(MAX)) {
  #call new column for the loop
  ANT_TASK_PERC <- paste0("AntTask", MAX[i]*100, "perc")
  metadata_present[which(metadata_present$prop_time_outside <= MAX[i]), ANT_TASK_PERC] <- "nurse"
  metadata_present[which(metadata_present$prop_time_outside > MAX[i]), ANT_TASK_PERC] <- "forager"
  metadata_present[which(metadata_present$IsQueen ==T), ANT_TASK_PERC] <- "queen"
#remove dead ants
#metadata_present <- metadata_present[which(metadata_present$IsAlive==TRUE),]

#AntTasks_by_col <- dcast(metadata_present, size_treat + REP_treat ~ AntTask, fun.aggregate = length)
AntTasks_by_col <- metadata_present %>%
  group_by(!!sym(ANT_TASK_PERC), size_treat, REP_treat) %>%
  summarise(occurrences = n())

#remove NAs
#AntTasks_by_col <- AntTasks_by_col[!is.na(AntTasks_by_col[[ANT_TASK_PERC]]),]
#remove Q
#AntTasks_by_col <- AntTasks_by_col[which(AntTasks_by_col[[ANT_TASK_PERC]]!="queen"),]

#proportion of ants by task group
AntTasks_by_col_prop <- AntTasks_by_col %>%
  group_by(REP_treat, size_treat, !!sym(ANT_TASK_PERC)) %>%
  summarise(total_occurrences = sum(occurrences)) %>%
  mutate(proportion = total_occurrences/sum(total_occurrences))

#plotting
AntTasks_by_col_prop_nurses <- AntTasks_by_col_prop[which(AntTasks_by_col_prop[[ANT_TASK_PERC]]=="nurse"),]
p <- ggplot(AntTasks_by_col_prop_nurses, aes(x = size_treat, y = proportion, fill = AntTasks_by_col_prop_nurses[[ANT_TASK_PERC]])) +
  geom_boxplot() +
  geom_jitter()+
  facet_wrap(. ~  AntTasks_by_col_prop_nurses[[ANT_TASK_PERC]]) +
  ylim(0, 1) + # set y-axis limits to the same range
  labs(x = "REP_treat", y = "Proportion", fill = ANT_TASK_PERC) +
  ggtitle( paste0("nurses prop if < ", MAX[i]*100, "%\ntime outside")) +
  theme_bw() +
  theme(legend.position = "none")

plots[[i]] <- p
}

# Combine ggplots into a grid
grid.arrange(grobs = plots, nrow = 1, ncol = 4)

metadata_present <- select(metadata_present, -c("AntTask1perc","AntTask3perc","AntTask4perc"))
colnames(metadata_present)[which(colnames(metadata_present)=="AntTask2perc")] <- "AntTask"

write.table(metadata_present,file=file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt"),append=F,col.names=T,row.names=F,quote=T,sep=",")


##################################################################
##################              ########################
##################################################################



