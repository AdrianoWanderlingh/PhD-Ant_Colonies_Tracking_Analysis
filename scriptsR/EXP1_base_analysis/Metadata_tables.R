### tables from metadata information
library(reshape2)
library(dplyr)
library(ggplot2)

USER <- "supercompAdriano"

if (USER=="supercompAdriano") {
  WORKDIR <- "/media/cf19810/DISK4/ADRIANO" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  #SCRIPTDIR <- "/media/cf19810/DISK4/EXP1_base_analysis/EXP1_analysis scripts"
}

 

metadata_present <- read.table(file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt"),header=T,stringsAsFactors = F, sep=",")

#AntTasks_by_col <- dcast(metadata_present, size_treat + REP_treat ~ AntTask, fun.aggregate = length)

AntTasks_by_col <- metadata_present %>%
  group_by(AntTask, size_treat, REP_treat) %>%
  summarize(occurrences = n())

ggplot(AntTasks_by_col, aes(x=size_treat, y=occurrences)) +
  geom_boxplot() +
  facet_wrap( . ~AntTask)

#remove NAs
AntTasks_by_col <- AntTasks_by_col[!is.na(AntTasks_by_col$AntTask),]
#remove Q
AntTasks_by_col <- AntTasks_by_col[which(AntTasks_by_col$AntTask!="queen"),]

AntTasks_by_col_prop <- AntTasks_by_col %>%
  group_by(REP_treat, size_treat, AntTask) %>%
  summarise(total_occurrences = sum(occurrences)) %>%
  mutate(proportion = total_occurrences/sum(total_occurrences))

ggplot(AntTasks_by_col_prop, aes(x = size_treat, y = proportion, fill = AntTask)) +
  geom_boxplot() +
  geom_jitter()+
  facet_wrap(. ~ AntTask) +
  labs(x = "REP_treat", y = "Proportion", fill = "AntTask") +
  ggtitle("AntTask = nurse if <1% time outside") +
  theme_bw()


