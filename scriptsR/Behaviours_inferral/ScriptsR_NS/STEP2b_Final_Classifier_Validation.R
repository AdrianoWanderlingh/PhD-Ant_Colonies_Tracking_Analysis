##########################################################################################
############## Final Classifier Validation  ##############################################
##########################################################################################

#### THIS VERSION IS FORT 0.8.2 COMPATIBLE ####

# Script created by Adriano Wanderlingh, Nathalie Stroeymeyt and Tom Richardson, with contributions by Enrico Gavagnign

#this should be the version of the script maintained for long term use.
#For previous versions of this script and the sourced one, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral


#clean start
rm(list=ls())
gc()

library(multcomp)

###parameter to set at start
USER <- "Nathalie"#"Adriano"
ANNOT_DATASET <- "Vasudha2021" # "Adriano2022" |"Vasudha2021
BEH <- "G"
FRAME_RATE <- 8
measure <- "CSI" ##""Fbeta" "CSI"
###TO DO: also edit the path to BODYLENGTH_FILE when defining folders from line 45 onwards

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
}
if (USER=="Tom")     {
}
if (USER=="Nathalie"){
  DATADIR <- "/media/bzniks/Seagate\ Portable\ Drive/ADRIANO/EXPERIMENT_DATA_EXTRAPOLATED"
  ANNOTATIONDIR <- "/media/bzniks/Seagate\ Portable\ Drive/Ants_behaviour_analysis/Data"
  SCRIPTDIR <- "~/Desktop/ScriptsR"
  MachineLearningOutcome_DIR <- "/media/bzniks/DISK1/MachineLearning_outcomes_NewAnnotations/"
  BODYLENGTH_FILE <- "/media/bzniks/Seagate\ Portable\ Drive/Ants_behaviour_analysis/Data/Mean_ant_length_per_TrackingSystem.txt"
}

###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"BEH_Extract_movement_variables_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_Auto_Man_agreement_matrix_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_PCA_fort081_FUNCTIONS.R",sep="/"))
source(paste(SCRIPTDIR,"BEH_self_defined_functions.R",sep="/"))
source(paste(SCRIPTDIR,"interaction_detection.R",sep="/"))
suppressMessages(source(paste(SCRIPTDIR,"BEH_libraries.R",sep="/")))
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste(SCRIPTDIR,"add_angles.cpp",sep="/"))
sourceCpp(paste(SCRIPTDIR,"merge_interactions.cpp",sep="/"))

#quality scores output file
output_scores <- file.path(ANNOTATIONDIR,paste0("quality_scores_ClassifierValidation_",ANNOT_DATASET,"_",Sys.Date(),".txt"))
output_name <- file.path(ANNOTATIONDIR,paste0("candidate_groomings_",ANNOT_DATASET,"_",Sys.Date(),".txt"))


###############################################################################
###### PARAMETERS #############################################################
###############################################################################
###body length information
all_body_lengths <-read.table(BODYLENGTH_FILE,header=T,stringsAsFactors = F,sep=",")

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
# ANT_LENGHT_PX               <- 153 #useful for matcher::meanAntDisplacement mean and median value are similar
# minimumGap                  <- fmSecond(0.5) ## THIS OPTION DOES NOT WORK IN INTERACTIONS SO DISABLED! for a given pair of interacting ants, when interaction is interrupted by more than minimumGap, interaction will check whether ants have moved since - and if so, will create new interaction

DISAGREEMENT_THRESH <- 0.5
###Fixed parameter
set.seed(2)       ###define I(arbitrary) seed so results can be reproduced over multiple runs of the program

### SELECT ANNOTATION DATASET TO USE
if (ANNOT_DATASET=="Vasudha2021") {
  FOCAL <- F
  ###############################################################################
  ###READ WHOLE ANNOTATIONS DATASET #############################################
  ###############################################################################
  print("Loading manual annotations...")
  #the current annotation file FINAL_script_output_CROSSVAL_25PERC_AND_TROPH.csv underwent cross-validation by Adriano
  annotations_all <- read.csv(paste(ANNOTATIONDIR,"/R3SP_R9SP_All_data_FINAL_script_output_CROSSVAL_25PERC_AND_TROPH.csv",sep = ""), sep = ",")
  
  ##rename some columns
  names(annotations_all)[which(grepl("rep",names(annotations_all)))]    <- "REPLICATE"
  names(annotations_all)[which(grepl("period",names(annotations_all)))] <- "PERIOD"
  
  
  annotations_all$Behaviour     <- as.character(annotations_all$Behaviour)
  annotations_all$Actor         <- as.character(annotations_all$Actor)
  annotations_all$Receiver      <- as.character(annotations_all$Receiver)
  
  ###define new time columns called T_start_UNIX and T_stop_UNIX and delete T_start and T_stop
  annotations_all$T_start_UNIX  <- as.POSIXct(annotations_all$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  annotations_all$T_stop_UNIX   <- as.POSIXct(annotations_all$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  annotations_all$T_start <- NULL
  annotations_all$T_stop <- NULL
  
  ###create object that will contain time limits for annotations_all (contains all annotations for all behaviour so true time limit for period analysed)
  time_window_all        <- merge(aggregate(T_start_UNIX ~ PERIOD + REPLICATE, FUN=min, data=annotations_all),aggregate(T_stop_UNIX  ~ PERIOD + REPLICATE, FUN=max, data=annotations_all))
  names(time_window_all) <- c("PERIOD","REPLICATE","time_start","time_stop")
  
  ###finally, subset annotations_all for the behaviour of interest
  annotations_all <- annotations_all[which(annotations_all$Behaviour==BEH),]
  
}


### SELECT ANNOTATION DATASET TO USE
if (ANNOT_DATASET=="Adriano2022") {
  
  FOCAL <- T
  ###############################################################################
  ###READ NEW ANNOTATIONS DATASET ##############################################
  ###############################################################################
  print("Loading manual annotations...")
  ###create object that will contain time limits for annotations_all (contains all annotations for all behaviour so true time limit for period analysed)
  metadata_info         <- read.csv(paste(ANNOTATIONDIR,"/Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU.csv",sep = ""), sep = ",")
  time_window_all        <- data.frame(
    PERIOD = "post" # "post"
    ,REPLICATE = metadata_info$REP_treat
    ,time_start = metadata_info$ReturnExposed_time
  )
  
  time_window_all$time_start <- as.POSIXct(time_window_all$time_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  ###add two minutes to start
  time_window_all$time_start <- time_window_all$time_start + 2*60
  time_window_all$time_stop  <- time_window_all$time_start + 30*60
  
  
  #the current annotation file FINAL_script_output_CROSSVAL_25PERC_AND_TROPH.csv underwent cross-validation by Adriano
  annotations_all <- read.csv(paste(ANNOTATIONDIR,"/Grooming_Classifier_CrossVal_ANNOTATIONS.csv",sep = ""), sep = ",")
  annotations_all$Behaviour     <- "G"
  
  ##rename some columns
  names(annotations_all)[which(grepl("rep",names(annotations_all)))]    <- "REPLICATE"
  names(annotations_all)[which(grepl("treatment",names(annotations_all)))] <- "PERIOD"
  
  annotations_all$Actor         <- as.character(annotations_all$Actor)
  annotations_all$Receiver      <- as.character(annotations_all$Receiver)
  
  
  ###define new time columns called T_start_UNIX and T_stop_UNIX and delete T_start and T_stop
  annotations_all$T_start_UNIX  <- as.POSIXct(annotations_all$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  annotations_all$T_stop_UNIX   <- as.POSIXct(annotations_all$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  annotations_all$T_start <- NULL
  annotations_all$T_stop <- NULL
  
}

###make sure annotations don't overrun 30 minute window
for (i in 1:nrow(annotations_all)){
  annotations_all[i,"T_stop_UNIX"] <- min(annotations_all[i,"T_stop_UNIX"],time_window_all[which(time_window_all$REPLICATE==annotations_all[i,"REPLICATE"]&time_window_all$PERIOD==annotations_all[i,"PERIOD"]),"time_stop"]-1/FRAME_RATE)
  annotations_all[i,"T_start_UNIX"] <- max(annotations_all[i,"T_start_UNIX"],time_window_all[which(time_window_all$REPLICATE==annotations_all[i,"REPLICATE"]&time_window_all$PERIOD==annotations_all[i,"PERIOD"]),"time_start"]+1/FRAME_RATE)
}



annotations_all$duration      <- as.numeric(annotations_all$T_stop_UNIX - annotations_all$T_start_UNIX)
annotations_all               <- annotations_all[which(annotations_all$duration>0),]

annotations_all$T_start_sec <- round(as.numeric(annotations_all$T_start_UNIX),N_DECIMALS)
annotations_all$T_stop_sec <- round(as.numeric(annotations_all$T_stop_UNIX),N_DECIMALS)


###define name of general output table containing quality scores
chosen_file_name <- file.path(MachineLearningOutcome_DIR, paste("quality_scores_",measure,"_CHOSEN.txt",sep=""))
chosen <- read.table(chosen_file_name,header=T,stringsAsFactors = F)


###############################################################################
###### EXTRACT CHOSEN PARAMETERS FROM CHOSEN ##################################
###############################################################################
# #####Arguments to loop over - to comment out when running the loop
subDir                      <- paste0("Loop_ID_",chosen[,"Loop_ID"])
CAPSULE_FILE                <- chosen[,"CAPSULE_FILE"]
DT_dist_THRESHOLD_BL        <- chosen[,"DT_dist_THRESHOLD_BL"]
MAX_INTERACTION_GAP         <- chosen[,"MAX_INTERACTION_GAP"]
DISAGREEMENT_THRESH         <- chosen[,"DISAGREEMENT_THRESH"]
trim_length_sec             <- chosen[,"trim_length_sec"]
DT_frame_THRESHOLD          <- chosen[,"DT_frame_THRESHOLD"]
beta                        <- chosen[,"beta"]
###Load BN_list
BN_list <- dget (  file.path(MachineLearningOutcome_DIR,subDir,"BN_object_list.dat") )
###Load_classifier
classifier        <- list(readRDS (file.path(MachineLearningOutcome_DIR,subDir,"fits",paste(chosen[,"classifier"],".rds",sep=""))))
names(classifier) <- chosen[,"classifier"]

##################################################################################
##################### EXTRACT VARIABLES FOR MANUAL AND AUTOMATIC INTERACTIONS ####
##################################################################################
###for these particular parameters, evaluate how successful the interaction detetection parameters are at detecting candidate frames
###extract interactions and manual annotations for the whole dataset
all      <- extraction_loop(chunk="all"
                            ,extract_movement_variables=T
                            ,selected_variables = unlist(lapply(names(BN_list),function(x)unlist(strsplit(x,split="\\."))[1]))
                            ,all_body_lengths=all_body_lengths
                            ,focal=FOCAL
)


candidate_groomings <- all[["summary_AUTO"]]

#~~~~~PREDICT GROONING USING CLASSIFIED
print("Predicting grooming...")
prediction_start <- Sys.time()
candidate_groomings["predicted_Hit"] <- predict_class  (summary_AUTO =    candidate_groomings
                                                        ,BN_list      =  BN_list
                                                        ,classifier   = classifier
)
print("Grooming predicted.")

###overall test
aut_man_agreement <- auto_manual_agreement (candidate_groomings[which(candidate_groomings$predicted_Hit==1),]
                                            , all[["summary_MANUAL"]]
                                            , all[["list_IF_Frames"]]
)
`quality_scores_overall` <-round(quality_scores(aut_man_agreement[["true_false_positive_negatives"]],beta),digits=3 )

quality_scores_REP_PER <- NULL
for (REP in unique(candidate_groomings$REPLICATE)) {
  for (PER in unique(candidate_groomings[which(candidate_groomings$REPLICATE==REP),"PERIOD"])){
    
    if (nrow(candidate_groomings[which(candidate_groomings$predicted_Hit==1 & candidate_groomings$REPLICATE==REP & candidate_groomings$PERIOD==PER),])<1) {
      #temporary solution for REPs without candidate grooming
      quality_scores_REP_PER <- data.frame(replicate=REP,period=PER,CSI=NA,Fbeta=NA,precision=NA,sensitivity=NA)
    }else{
      aut_man_agreement_REP_PER <- auto_manual_agreement (candidate_groomings[which(candidate_groomings$predicted_Hit==1 & candidate_groomings$REPLICATE==REP & candidate_groomings$PERIOD==PER),]
                                                          , subset(all[["summary_MANUAL"]],all[["summary_MANUAL"]]$REPLICATE==REP&all[["summary_MANUAL"]]$PERIOD==PER)
                                                          , lapply(all["list_IF_Frames"], function(x) x[grep(paste0("IF_frames_",REP),names(x))])$list_IF_Frames
                                                          #, subset(all["list_IF_Frames"]$list_IF_Frames)
      )
      
      qual_scores <-round(quality_scores(aut_man_agreement_REP_PER[["true_false_positive_negatives"]],beta),digits=3 )
      quality_scores_REP_PER <- rbind(quality_scores_REP_PER,data.frame(replicate=REP,period=PER,as.list(qual_scores)))
    }
    
    
    #add tracking system info, $spaces?
    
    ####  
    
    # if (file.exists(output_scores)){
    #   write.table(quality_scores_REP,file=output_scores,append=T,col.names=F,row.names=F,quote=T)
    # }else{
    #   write.table(quality_scores_REP,file=output_scores,append=F,col.names=T,row.names=F,quote=T)
    # }
    
  }
  
}

if (ANNOT_DATASET=="Adriano2022") {
  quality_scores_REP_PER <- within(  quality_scores_REP_PER,  treatment <- substr(replicate,nchar(replicate)-1,nchar(replicate)))
  quality_scores_REP_PER$treatment <- factor(quality_scores_REP_PER$treatment )
  
  print(aggregate(cbind(precision,sensitivity)~treatment,function(x)cbind(mean(x),std.error(x)),data=quality_scores_REP_PER))
  
  model_sensitivity <- lm(     sensitivity ~ treatment, data=quality_scores_REP_PER)
  shapiro.test(residuals(model_sensitivity))    
  anova(model_sensitivity)
  
  model_precision <- lm(     (precision)^4 ~ treatment, data=quality_scores_REP_PER)
  shapiro.test(residuals(model_precision))  
  anova(model_precision)
  summary(glht(model_precision,linfct = multcomp::mcp(treatment="Tukey")),test=adjusted("BH"))
}

###############################################################################
######        SAVING FINAL GROOMING TABLE        ##############################
###############################################################################
# ### AW
# if (file.exists(output_name)){
#   write.table(inferred_groomings,file=output_name,append=T,col.names=F,row.names=F,quote=T,sep=",")
# }else{
#   write.table(inferred_groomings,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
# }

# write.table(candidate_groomings,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
