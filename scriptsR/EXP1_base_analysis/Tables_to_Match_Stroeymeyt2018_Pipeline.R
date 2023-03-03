##################################################################
####### TABLES TO MATCH STROEYMEYT ET AL., 2018, SCIENCE #########
##################################################################

### tables from metadata information
library(reshape2)
library(dplyr)
library(FortMyrmidon)

#### FUNCTIONS
#list files recursive up to a certain level (level defined by "n" parameter)
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}


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
metadata_present$REP_NUM           <- substring(metadata_present$REP_treat, 2, nchar(metadata_present$REP_treat))
metadata_present$N_CHAR            <-  4-nchar(metadata_present$REP_NUM)
#rep function is unhappy if the reps are 0...
#paste("colony", paste(rep(0,N_CHAR),collapse=""), metadata_present$REP_NUM,sep="")
metadata_present[which(metadata_present$N_CHAR==0),"colony"] <- paste("colony", metadata_present[which(metadata_present$N_CHAR==0),"REP_NUM"],sep="")
metadata_present[which(metadata_present$N_CHAR==1),"colony"] <- paste("colony", 0, metadata_present[which(metadata_present$N_CHAR!=0),"REP_NUM"],sep="")
# colony_status
metadata_present$status_char       <- substr(metadata_present$REP_treat, nchar(metadata_present$REP_treat), nchar(metadata_present$REP_treat))
metadata_present$colony_status     <- ifelse(metadata_present$status_char  == "P", "pathogen", ifelse(metadata_present$status_char  == "S", "control", NA))
# treatment code
#period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
metadata_present$size_char         <- substr(metadata_present$REP_treat, nchar(metadata_present$REP_treat)-1, nchar(metadata_present$REP_treat)-1)
metadata_present$size_status       <- ifelse(metadata_present$size_char  == "S", "small", ifelse(metadata_present$size_char  == "B", "big", NA))
metadata_present$treatment_code    <- paste(metadata_present$colony_status ,metadata_present$size_status,sep=".")

##### ADD EXTRA COLS TO PATHOHEN LOAD DATA
# colony code
pathogen_load$REP_NUM           <- substring(pathogen_load$Colony, 2, nchar(pathogen_load$Colony))
pathogen_load$N_CHAR            <-  4-nchar(pathogen_load$REP_NUM)
#rep function is unhappy if the reps are 0...
#paste("colony", paste(rep(0,N_CHAR),collapse=""), pathogen_load$REP_NUM,sep="")
pathogen_load[which(pathogen_load$N_CHAR==0),"colony"] <- paste("colony", pathogen_load[which(pathogen_load$N_CHAR==0),"REP_NUM"],sep="")
pathogen_load[which(pathogen_load$N_CHAR==1),"colony"] <- paste("colony", 0, pathogen_load[which(pathogen_load$N_CHAR!=0),"REP_NUM"],sep="")
# colony_status
pathogen_load$status_char       <- substr(pathogen_load$Colony, nchar(pathogen_load$Colony), nchar(pathogen_load$Colony))
pathogen_load$colony_status     <- ifelse(pathogen_load$status_char  == "P", "pathogen", ifelse(pathogen_load$status_char  == "S", "control", NA))
# treatment code
#period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
pathogen_load$size_char         <- substr(pathogen_load$Colony, nchar(pathogen_load$Colony)-1, nchar(pathogen_load$Colony)-1)
pathogen_load$size_status       <- ifelse(pathogen_load$size_char  == "S", "small", ifelse(pathogen_load$size_char  == "B", "big", NA))
pathogen_load$treatment_code    <- paste(pathogen_load$colony_status ,pathogen_load$size_status,sep=".")


##################################################################
##################     task_groups.txt    ########################
##################################################################

task_groups <- data.frame(colony = metadata_present$colony,
                          tag =  metadata_present$antID,
                          task_group = metadata_present$AntTask,
                          treatment = metadata_present$treatment_code,
                          REP_treat= metadata_present$REP_treat)
write.table(task_groups, file = file.path(SAVEDIR,"task_groups.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

##################################################################
##############     treated_worker_list.txt    ####################
##################################################################

metadata_treated <- metadata_present[which(metadata_present$Exposed==TRUE),]

treated_worker_list <- data.frame(colony = metadata_treated$colony,
                                  tag =  metadata_treated$antID,
                                  survived_treatment = T,
                                  treatment = metadata_treated$treatment_code,
                                  REP_treat= metadata_treated$REP_treat)

write.table(treated_worker_list, file = file.path(SAVEDIR,"treated_worker_list.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

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

write.table(info, file = file.path(SAVEDIR,"info.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

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


write.table(qPCR_file, file = file.path(SAVEDIR,"qPCR","qPCR_file.txt"), append = F, col.names = T, row.names = F, quote = F, sep = "\t")

##################################################################
#####################    tag_files.txt    ########################
##################################################################

### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]

# replicate folder
for (REP.n in 1:length(files_list)) {
  # REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = "CapsuleDef2018_q.myr")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  #replicate file
  for (REP.FILES in REP.filefolder) {
    #get substring in variable until R9BS_
    # REP.FILES <- REP.filefolder[2]
    REP_treat <- sub("\\_.*", "", basename(REP.FILES))
    
    
    #conform naming to science2018
    # colony code
    REP_NUM           <- substring(REP_treat, 2, nchar(REP_treat))
    colony            <- paste("colony", paste(rep(0,4-nchar(REP_NUM)),collapse=""), REP_NUM,sep="")
    # colony_status
    status_char       <- substr(REP_treat, nchar(REP_treat), nchar(REP_treat))
    colony_status     <- ifelse(status_char == "P", "pathogen", ifelse(status_char == "S", "control", NA))
    # treatment code
    size_char         <- substr(REP_treat, nchar(REP_treat)-1, nchar(REP_treat)-1)
    size_status       <- ifelse(size_char == "S", "small", ifelse(size_char == "B", "big", NA))
    treatment_code    <- paste(colony_status,size_status,sep=".")
    
    print(paste(REP_treat,colony,treatment_code,sep=" | ")) ##}}
    #open experiment
    e <- fmExperimentOpen(REP.FILES)
    exp.Ants  <- e$ants
    exp_end   <- fmQueryGetDataInformations(e)$end
    exp_start <- fmQueryGetDataInformations(e)$start

    tag_stats <- fmQueryComputeTagStatistics(e)
    
    ############# CREATE BASE FILE
    tag_file <- NULL
    
    for (ant in exp.Ants){
      for (id in ant$identifications){
        #if the ant died, skip
        if (capture.output(ant$identifications[[1]]$end)=="+∞") {
        tag_file <- rbind(tag_file,data.frame(tag =  ant$ID,
                               count = tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"count"], #count considers also the acclimation time, therefore it is always higher than the N of frames
                               last_det = 51*60*60 * 8, #fixed= sll alive, so it is the n of frames for 51 hours.   tag_stats[which(tag_stats$tagDecimalValue == id$tagValue),"lastSeen"], 
                               rot = round(deg(id$antAngle),3), # assumed to be the the relative angle of the ant to the tag
                               displacement_distance = 0, #id$antPosition one of the two measures
                               displacement_angle = 0,
                               antenna_reach = NA,
                               trapezoid_length = NA,
                               type = "N", #what does it mean?
                               size = 0 ,
                               headwidth = 0,
                               death = 0,  #all alive, dead not included
                               age = 0,
                               final_status = "alive", #ifelse(metadata_present[which(metadata_present$tagIDdecimal == id$tagValue & metadata_present$REP_treat==REP_treat),"IsAlive"]==T,"alive","dead"),
                               group = metadata_present[which(metadata_present$tagIDdecimal == id$tagValue & metadata_present$REP_treat==REP_treat),"AntTask"],
                               REP_treat = REP_treat,
                               stringsAsFactors = F))
        
        tag_file$rot <- round(tag_file$rot,3)
        #identifEnd <- ifelse(capture.output(print(id$end))=="+∞",NA,ifelse(capture.output(print(id$end))))
        #fmQueryGetDataInformations(e)$end - 51*60*60
        
  
          # individual  <- ant$ID
          # #extract metadata info per key
          #   #for more than one row, always last row will be picked (the relevant one if there is a timed change or the default one if there is not)
          #   for(ROW in 1:2) {
          #     #assign metadata_key value when ID corresponds
          #     tag_file[which(tag_file$tag==ant$ID),"final_status"] <- ant$getValues("IsAlive")["IsAlive","values"]
          #     #if the ant died
          #     if (METADATA_KEY=="IsAlive") {
          #       if (ant$getValues("IsAlive")[ROW,"values"]==FALSE) {
          #         metadata[which(metadata$antID==ant$ID),"surviv_time"] <- ant$getValues("IsAlive")[ROW,"times"]
          #         # if didn't die    
          #       }else if (ant$getValues("IsAlive")[ROW,"values"]==TRUE) {
          #         metadata[which(metadata$antID==ant$ID),"surviv_time"] <- as.POSIXct( exp_end,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
          #       }} # IsAlive check
          #   } # ROW
        }
    }
    
      
    }#ant is dead
 
    # #check if the tag_file file exists
    # if(file.exists(file.path(SAVEDIR,"tag_files",paste(colony,treatment_code,".txt")))){
    #     print(paste0(REP_treat," already present in tag_file, skip"))
    #   } else {

        write.table(tag_file, file = file.path(SAVEDIR,"tag_files",paste0(colony,"_",treatment_code,".txt")), append = F, col.names = T, row.names = F, quote = F, sep = "\t")
        
        
        # }
  } # REP folders
} # REP by REP


