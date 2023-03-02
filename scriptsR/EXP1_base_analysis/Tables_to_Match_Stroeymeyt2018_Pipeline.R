##################################################################
####### TABLES TO MATCH STROEYMEYT ET AL., 2018, SCIENCE #########
##################################################################

### tables from metadata information
library(reshape2)
library(dplyr)

USER <- "supercompAdriano"

if (USER=="supercompAdriano") {
  WORKDIR <- "/media/cf19810/DISK4/ADRIANO" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  SAVEDIR <- "/media/cf19810/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment/original_data"
  #SCRIPTDIR <- "/media/cf19810/DISK4/EXP1_base_analysis/EXP1_analysis scripts"
}

###### LOAD METADATA
metadata_present <- read.table(file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt"),header=T,stringsAsFactors = F, sep=",")
pathogen_load <- read.csv("/media/cf19810/DISK4/EXP1_base_analysis/Personal_Immunity/Pathogen_Quantification_Data/Adriano_qPCR_pathogen_load_MASTER_REPORT.csv")

#remove dead ants
metadata_present <- metadata_present[which(metadata_present$IsAlive==TRUE),]
metadata_present$colony <- NA

##### ADD EXTRA COLS TO METADATA
#conform naming to science2018
# colony code
metadata_present$REP_NUM           <- substring(metadata_present$REP_treat, 2, nchar(metadata_present$REP_treat)-1)
metadata_present$N_CHAR            <-  3-nchar(metadata_present$REP_NUM)
#rep function is unhappy if the reps are 0...
#paste("colony", paste(rep(0,N_CHAR),collapse=""), metadata_present$REP_NUM,sep="")
metadata_present[which(metadata_present$N_CHAR==0),"colony"] <- paste("colony", metadata_present[which(metadata_present$N_CHAR==0),"REP_NUM"],sep="")
metadata_present[which(metadata_present$N_CHAR==1),"colony"] <- paste("colony", 1, metadata_present[which(metadata_present$N_CHAR!=0),"REP_NUM"],sep="")
# colony_status
metadata_present$status_char       <- substr(metadata_present$REP_treat, nchar(metadata_present$REP_treat), nchar(metadata_present$REP_treat))
metadata_present$colony_status     <- ifelse(metadata_present$status_char  == "P", "pathogen", ifelse(metadata_present$status_char  == "S", "control", NA))
# treatment code
#period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
metadata_present$size_char         <- substr(metadata_present$REP_treat, nchar(metadata_present$REP_treat)-1, nchar(metadata_present$REP_treat)-1)
metadata_present$size_status       <- ifelse(metadata_present$size_char  == "S", "small", ifelse(metadata_present$size_char  == "B", "big", NA))
metadata_present$treatment_code    <- paste(metadata_present$colony_status ,metadata_present$size_status,sep="_")

##### ADD EXTRA COLS TO PATHOHEN LOAD DATA
# colony code
pathogen_load$REP_NUM           <- substring(pathogen_load$Colony, 2, nchar(pathogen_load$Colony)-1)
pathogen_load$N_CHAR            <-  3-nchar(pathogen_load$REP_NUM)
#rep function is unhappy if the reps are 0...
#paste("colony", paste(rep(0,N_CHAR),collapse=""), pathogen_load$REP_NUM,sep="")
pathogen_load[which(pathogen_load$N_CHAR==0),"colony"] <- paste("colony", pathogen_load[which(pathogen_load$N_CHAR==0),"REP_NUM"],sep="")
pathogen_load[which(pathogen_load$N_CHAR==1),"colony"] <- paste("colony", 1, pathogen_load[which(pathogen_load$N_CHAR!=0),"REP_NUM"],sep="")
# colony_status
pathogen_load$status_char       <- substr(pathogen_load$Colony, nchar(pathogen_load$Colony), nchar(pathogen_load$Colony))
pathogen_load$colony_status     <- ifelse(pathogen_load$status_char  == "P", "pathogen", ifelse(pathogen_load$status_char  == "S", "control", NA))
# treatment code
#period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
pathogen_load$size_char         <- substr(pathogen_load$Colony, nchar(pathogen_load$Colony)-1, nchar(pathogen_load$Colony)-1)
pathogen_load$size_status       <- ifelse(pathogen_load$size_char  == "S", "small", ifelse(pathogen_load$size_char  == "B", "big", NA))
pathogen_load$treatment_code    <- paste(pathogen_load$colony_status ,pathogen_load$size_status,sep="_")


##################################################################
##################     task_groups.txt    ########################
##################################################################

task_groups <- data.frame(colony = metadata_present$colony,
                          tag =  metadata_present$antID,
                          task_group = "decide_threshold",
                          treatment = metadata_present$treatment_code,
                          REP_treat= metadata_present$REP_treat)
warning("task_group = decide_threshold")
write.table(task_groups, file = file.path(SAVEDIR,"task_groups.txt"), append = F, col.names = T, row.names = F, quote = T, sep = ",")

##################################################################
##############     treated_worker_list.txt    ####################
##################################################################

metadata_treated <- metadata_present[which(metadata_present$Exposed==TRUE),]

treated_worker_list <- data.frame(colony = metadata_treated$colony,
                                  tag =  metadata_treated$antID,
                                  survived_treatment = T,
                                  treatment = metadata_treated$treatment_code,
                                  REP_treat= metadata_treated$REP_treat)

write.table(treated_worker_list, file = file.path(SAVEDIR,"treated_worker_list.txt"), append = F, col.names = T, row.names = F, quote = T, sep = ",")

##################################################################
#######################    info.txt    ###########################
##################################################################

metadata_colony <- metadata_present[,c("colony","treatment_code","box","colony_size","REP_treat")]
metadata_colony <- unique(metadata_colony)

info <- data.frame(colony = metadata_colony$colony,
                  treatment = metadata_colony$treatment_code,
                  box = metadata_colony$box,
                  colony_size = metadata_colony$colony_size,
                  Ynestmin = NA,
                  Ynestmax = NA,
                  nb_larvae = NA,
                  nb_pupae = NA,
                  colony_age = NA,
                  REP_treat= metadata_colony$REP_treat)

write.table(info, file = file.path(SAVEDIR,"info.txt"), append = F, col.names = T, row.names = F, quote = T, sep = ",")

##################################################################
#####################    qPCR_file.txt    ########################
##################################################################

head(pathogen_load[,c("colony","treatment_code","antID","Exposed","IsAlive","MbruDNA","above_detection_threshold", "Colony")])

qPCR_file <- data.frame(colony = pathogen_load$colony,
                        treatment = pathogen_load$treatment_code,
                        tag =  pathogen_load$antID,
                        age = NA,
                        status = pathogen_load$Exposed,
                        alive_at_sampling_time = pathogen_load$IsAlive,
                        measured_load_ng_per_uL = pathogen_load$MbruDNA,
                        above_detection_threshold = pathogen_load$above_detection_threshold,
                        REP_treat= pathogen_load$Colony)


write.table(qPCR_file, file = file.path(SAVEDIR,"qPCR","qPCR_file.txt"), append = F, col.names = T, row.names = F, quote = T, sep = ",")
