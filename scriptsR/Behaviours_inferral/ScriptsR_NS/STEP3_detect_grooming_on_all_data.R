##########################################################################################
############## BEH MAIN Behaviours Analysis ##############################################
##########################################################################################

#### THIS VERSION IS FORT 0.8.1 COMPATIBLE ####

#this should be the version of the script maintained for long term use.
#For previous versions of this script and the sourced one, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral

#NOTE: at the current time (24 March 2022), the computational part of the script works but the saving of the plots does not (needs to be fixed, possibly according to plots_structure.R)

#clean start
rm(list=ls())
gc()

###parameter to set at start
USER <- "Nathalie"
###To set by user
BEH <- "G"
FRAME_RATE <- 8

if (BEH =="G"){
  interactions_of_interest <- list(c("head","body"))
  ###if interested in more than one type, list them as follows: interactions_of_interest <- c(list(c("head","body")),list(c("head","head")))
}

###############################################################################
###### GLOSSARY ###############################################################
###############################################################################

# MAN         : manual interactions/trajectories deriving from the hand annotated data (from the file annotations)
# AUTO        : automatically extracted interactions/trajectories deriving from fmQueryComputeAntInteractions
# REPLICATE   : nest ID, either "R3SP" or "R9SP" (SP: small pathogen)
# PERIOD      : the treatment period, either "pre" or "post" pathogen exposure
# REP_PER     : each of the 4 blocks analised, crossing the REPLICATE and PERIOD

###############################################################################
###### LOAD LIBRARIES AND FUNCTIONS #####################################
###############################################################################

####### navigate to folder containing myrmidon file
if (USER=="Adriano") {
  WORKDIR <- "/home/cf19810/Documents/Ants_behaviour_analysis"
  DATADIR <- paste(WORKDIR,"Data",sep="/")
  SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")
  EXPDATADIR <- "/media/cf19810/DISK4/EXPERIMENT_DATA"
}
if (USER=="Tom")     {
  WORKDIR <- "/media/tom/MSATA/Dropbox/Ants_behaviour_analysis"
  DATADIR <- paste(WORKDIR,"Data",sep="/")
  SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")
  EXPDATADIR <- ""
}
if (USER=="Nathalie"){
  WORKDIR <- "/media/bzniks/DISK3/Ants_behaviour_analysis"
  DATADIR <- paste(WORKDIR,"Data",sep="/")
  SCRIPTDIR <- "/home/bzniks/Dropbox/SeniorLectureship_Bristol/Students_postdocs/PhD_students/2019 Adriano Wanderlingh/code/PhD-exp1-data-analysis-main/scriptsR/Behaviours_inferral_Nathalie_NEW"
  EXPDATADIR <- "/media/bzniks/DISK3/ADRIANO/EXPERIMENT_DATA"
}

###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"BEH_Extract_movement_variables_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_Auto_Man_agreement_matrix_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_PCA_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_self_defined_functions.R",sep="/"))
source(paste(SCRIPTDIR,"interaction_detection.R",sep="/"))
source(paste(SCRIPTDIR,"trajectory_extraction.R",sep="/"))
suppressMessages(source(paste(SCRIPTDIR,"BEH_libraries.R",sep="/")))
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste(SCRIPTDIR,"add_angles.cpp",sep="/"))
sourceCpp(paste(SCRIPTDIR,"merge_interactions.cpp",sep="/"))

###############################################################################
###### PARAMETERS #############################################################
###############################################################################
#plotting limits used for the coordinates plotting
Xmin <- 2000
Xmax <- 7900
Ymin <- 0000
Ymax <- 5500

#general
N_DECIMALS                  <- 3 ## when assigning time in seconds, the number of decimals to preserve when rounding to match between interactions & collisions
FUZZY_MATCH                 <- TRUE  ## fuzzy matching between data frames in collisions detection

#trajectories cutting gap, relevant for fmQueryComputeAntTrajectories
max_gap                     <- fmHour(24*365)   ## important parameter to set! Set a maximumGap high enough that no cutting will happen. For example, set it at an entire year: fmHour(24*365)

###NATH_FLAG: don't just use max and min, give some extra 
###NATH_FLAG: these parameters should be adjusted for each box as ant length may depend on focus /  camera distance
# AntDistanceSmallerThan      <- 300 #for higher accuracy, recalculate it from: max(interaction_MANUAL$straightline_dist_px,na.rm = T)
# AntDistanceGreaterThan      <- 70 #for higher accuracy, recalculate it from: min(interaction_MANUAL$straightline_dist_px,na.rm = T)
ANT_LENGHT_PX               <- 153 #useful for matcher::meanAntDisplacement mean and median value are similar
# minimumGap                  <- fmSecond(0.5) ## THIS OPTION DOES NOT WORK IN INTERACTIONS SO DISABLED! for a given pair of interacting ants, when interaction is interrupted by more than minimumGap, interaction will check whether ants have moved since - and if so, will create new interaction

DISAGREEMENT_THRESH <- 0.5
###Fixed parameter
set.seed(2)       ###define I(arbitrary) seed so results can be reproduced over multiple runs of the program

###############################################################################
#### READ CHOSEN METHOD #######################################################
###############################################################################
chosen <- read.table(paste(WORKDIR,"/Data/MachineLearning_outcomes/quality_scores_CHOSEN.txt",sep=""),header=T,stringsAsFactors = F)

###############################################################################
###### EXTRACT CHOSEN PARAMETERS FROM CHOSEN ##################################
###############################################################################
# #####Arguments to loop over - to comment out when running the loop
subDir                      <- paste0("Loop_ID_",chosen[,"Loop_ID"])
CAPSULE_FILE                <- chosen[,"CAPSULE_FILE"]
DT_dist_THRESHOLD           <- chosen[,"DT_dist_THRESHOLD"]
MAX_INTERACTION_GAP         <- chosen[,"MAX_INTERACTION_GAP"]
DISAGREEMENT_THRESH         <- chosen[,"DISAGREEMENT_THRESH"]
trim_length_sec             <- chosen[,"trim_length_sec"]
DT_frame_THRESHOLD          <- chosen[,"DT_frame_THRESHOLD"]

###Load BN_list
BN_list <- dget (  file.path(DATADIR, "MachineLearning_outcomes",subDir,"BN_object_list.dat") )
###Load_classifier
classifier        <- list(readRDS (file.path(DATADIR, "MachineLearning_outcomes",subDir,"fits",paste(chosen[,"classifier"],".rds",sep=""))))
names(classifier) <- chosen[,"classifier"]

###############################################################################
###### LIST MYRMIDON FILES ON WHICH TO PERFORM GROOMING DETECTION ##################################
###############################################################################
setwd(EXPDATADIR)
myrmidon_files <- list.files(pattern="_AutoOriented_withMetaData_NS",recursive = T)
myrmidon_files <- myrmidon_files[which(grepl(CAPSULE_FILE,myrmidon_files))]

to_keep <- c(ls(),"to_keep","myrmidon_file")
for (myrmidon_file in myrmidon_files){
  print(myrmidon_file)
  rm(list=ls()[which(!ls()%in%to_keep)])
  gc()
  
  REPLICATE <- unlist(strsplit(myrmidon_files,split="_"))[grepl("SP|SS|BP|BS",unlist(strsplit(myrmidon_files,split="_")))]
  
  e <- fmExperimentOpen(myrmidon_file)
  
  PERIOD    <- "whole_experiment" ###to gain time you might want to loop over pre/post 24 hour periods and thus use different time_start and time_stop and different values for PERIOD
  time_start <- fmTimeCreate(offset=fmQueryGetDataInformations(e)$start) 
  time_stop  <- fmTimeCreate(offset=fmQueryGetDataInformations(e)$end)

  IdentifyFrames      <- fmQueryIdentifyFrames(e,start=time_start, end=time_stop,showProgress = FALSE)
  IF_frames           <- IdentifyFrames$frames
  rm(list=c("IdentifyFrames"))
  gc()
  # Assign a frame to each time since start and use it as baseline for all matching and computation
  IF_frames$frame_num <- as.numeric(seq.int(nrow(IF_frames)))
  
  # assign a time in sec to match annotation time (UNIX hard to use for class mismatch)
  IF_frames$time_sec <- round(as.numeric(IF_frames$time),N_DECIMALS)
  IF_frames$cum_diff <- c(0, with(IF_frames, diff(time)))
  IF_frames$cum_diff <- cumsum(IF_frames$cum_diff)
  
  ###get trajectories using self-written function which automatically performs the extraction on successive chunks of 12 hours to avoid crashes, and merge all trajectories for a single ant into a single trajectory
  positions <- extract_trajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = T,showProgress = F,IF_frames=IF_frames) #set true to obtain the zone of the ant
  
  ## immediately code hard link between  positions$trajectories_summary and  positions$trajectories
  positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
  names(positions$trajectories)            <- positions$trajectories_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
  
  ###add frame number information to trajectories
  positions$trajectories_summary$frame_num <- NA
  #assign starting frame number
  positions$trajectories_summary["frame_num"] <- lapply(positions$trajectories_summary["start"], function(x) IF_frames$frame_num[match(x, IF_frames$time)])
  
  ## Add ant_x names and times to the positions to convert from "FRAME since the start of the experiment", to FRAMES
  for (traj_idx in 1:nrow(positions$trajectories_summary))
  {
    AntID                         <- paste("ant_",positions$trajectories_summary[traj_idx,"antID"],sep="") ###NATH_FLAG: why not loop directly on antID_str which you just created?
    First_Obs_Frame               <- positions$trajectories_summary[traj_idx,"frame_num"] 
    IF_frames$new_zero_diff  <- IF_frames$cum_diff - IF_frames[IF_frames$frame_num==First_Obs_Frame,"cum_diff"] #subtracting the $frames zeroed-time  corresponding to the $start time from the zeroed-time column itself (New-Zeroed-time)
    # print(paste("Adding first obs FRAME", First_Obs_Frame, "to the time-zeroed trajectory of ant", AntID))
    #assign corresponding frame N when the New-Zeroed-time and $time correspond, closest.match 0.05 (well inside the 0.125 frame length in sec)
    positions$trajectories[[traj_idx]]$frame <- match.closest(x = positions$trajectories[[traj_idx]]$time, table = IF_frames$new_zero_diff, tolerance = 0.05)
    IF_frames$new_zero_diff <- NA
  }
  
  capsules  <- e$antShapeTypeNames
  names(capsules) <- as.character( 1:length(capsules))
  
  ALL_CAPS_MATCHERS <- list()
  for (interaction_of_interest in interactions_of_interest){
    caps_matcher <- c()
    for (caps in interaction_of_interest[c(1:2)]){
      caps_matcher <- c(caps_matcher,as.numeric(names(capsules)[[which(capsules==caps)]]) )
    }
    ALL_CAPS_MATCHERS <- c(ALL_CAPS_MATCHERS,list(caps_matcher))
  }
  # head_id <- as.numeric(names(capsules)[[which(capsules=="head")]])
  # body_id <- as.numeric(names(capsules)[[which(capsules=="body")]])
  gc()
  interac_start <- Sys.time()
  interacts_AUTO_REP_PER <- interaction_detection (e=e
                                                   ,start=time_start
                                                   ,end=time_stop
                                                   ,max_time_gap = MAX_INTERACTION_GAP
                                                   ,max_distance_moved = 2*ANT_LENGHT_PX
                                                   ,capsule_matcher=ALL_CAPS_MATCHERS
                                                   ,IF_frames=IF_frames
  )
  interac_stop <- Sys.time()

  ## Add ant_x names and times to the interacts_AUTO_REP_PER to convert from FRAME since the start of the experiment, to FRAMES
  ##creates a ID string for each ant in $interactions
  interacts_AUTO_REP_PER$ant1ID_str            <- paste("ant_",interacts_AUTO_REP_PER$ant1,sep="")
  interacts_AUTO_REP_PER$ant2ID_str            <- paste("ant_",interacts_AUTO_REP_PER$ant2,sep="")
  # Assign interaction pair
  interacts_AUTO_REP_PER$pair <- paste(interacts_AUTO_REP_PER$ant1, interacts_AUTO_REP_PER$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...
  
  ##convert times to different formats
  interacts_AUTO_REP_PER$T_start_UNIX <- as.POSIXct(interacts_AUTO_REP_PER$start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  interacts_AUTO_REP_PER$T_stop_UNIX  <- as.POSIXct(interacts_AUTO_REP_PER$end,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  #assign time in sec to avoid issues on time management and matching
  interacts_AUTO_REP_PER$T_start_sec <- round(as.numeric(interacts_AUTO_REP_PER$T_start_UNIX),N_DECIMALS)
  interacts_AUTO_REP_PER$T_stop_sec <- round(as.numeric(interacts_AUTO_REP_PER$T_stop_UNIX),N_DECIMALS)
  
  #assign start and end frame number
  interacts_AUTO_REP_PER$frame_start <- match.closest(x = interacts_AUTO_REP_PER$T_start_sec, table = IF_frames$time_sec, tolerance = 0.05)
  interacts_AUTO_REP_PER$frame_stop  <- match.closest(x = interacts_AUTO_REP_PER$T_stop_sec,  table = IF_frames$time_sec, tolerance = 0.05)
  
  #calc duration (including start frame)
  interacts_AUTO_REP_PER$duration <- interacts_AUTO_REP_PER$T_stop_sec - interacts_AUTO_REP_PER$T_start_sec + 1/FRAME_RATE
  #hist(interacts_AUTO_REP_PER$interactions$Duration,breaks = 30)
  interacts_AUTO_REP_PER$unique_interaction_id <- 1:nrow(interacts_AUTO_REP_PER)
  
  
  #~~~~~EXTRACT MOVEMERNT VARIABLES FROM DETECTED INTERACTIONS
  print("Extracting movement variables...")
  extraction_start <- Sys.time()
  candidate_groomings <- extract_from_object(interacts_AUTO_REP_PER
                                                  ,IF_frames
                                                  ,positions
                                                  ,BEH
                                                  ,REPLICATE
                                                  ,PERIOD
                                                  ,selected_variables = unlist(lapply(names(BN_list),function(x)unlist(strsplit(x,split="\\."))[1]))
  )[["summary_variables"]]
  extraction_stop <- Sys.time()
  
  #~~~~~PREDICT GROONING USING CLASSIFIED
  print("Predicting grooming...")
  prediction_start <- Sys.time()
  candidate_groomings["predicted_Hit"] <- predict_class  (summary_AUTO =    candidate_groomings
                                                         ,BN_list      =  BN_list
                                                         ,classifier   = classifier
  )
  prediction_stop <- Sys.time()
  
  #~~~~~COPY USEFUL INFO FROM CANDIDATE GROOMINGS TO INTERACTS
  interacts_AUTO_REP_PER <- cbind(interacts_AUTO_REP_PER,candidate_groomings[match(interacts_AUTO_REP_PER$unique_interaction_id,candidate_groomings$unique_interaction_id),c("Act_Name","Rec_Name","predicted_Hit")])                   
  
  ###FINAL GROOMING TABLE, TO SAVE
  inferred_groomings     <- interacts_AUTO_REP_PER[which(interacts_AUTO_REP_PER$predicted_Hit==1),]
  

  print(paste("Interaction detection took",round((as.numeric(interac_stop)-as.numeric(interac_start))/60,digits=2),"minutes."))
  print(paste("Extraction of movement variables took",round((as.numeric(extraction_stop)-as.numeric(extraction_start))/3600,digits=2),"hours."))
  print(paste("Prediction took",round((as.numeric(prediction_stop)-as.numeric(prediction_start)),digits=2),"seconds."))
  
}