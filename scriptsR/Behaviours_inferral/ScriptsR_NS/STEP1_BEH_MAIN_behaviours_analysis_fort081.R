##########################################################################################
############## BEH MAIN Behaviours Analysis ##############################################
##########################################################################################

#### THIS VERSION IS FORT 0.8.1 COMPATIBLE ####

# Script created by Adriano Wanderlingh, Nathalie Stroeymeyt and Tom Richardson, with contributions by Enrico Gavagnign

#this should be the version of the script maintained for long term use.
#For previous versions of this script and the sourced one, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral


#clean start
rm(list=ls())
gc()

###parameter to set at start
USER <- "Adriano"
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
  WORKDIR <- "/media/cf19810/DISK4/Ants_behaviour_analysis"
  DATADIR <- paste(WORKDIR,"Data",sep="/")
  SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")
  SAVEOUTPUT <- "/home/cf19810/Documents"
}
if (USER=="Tom")     {
  WORKDIR <- "/media/tom/MSATA/Dropbox/Ants_behaviour_analysis"
  DATADIR <- paste(WORKDIR,"Data",sep="/")
  SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")
  SAVEOUTPUT <- "/home/bzniks/Documents"
}
if (USER=="Nathalie"){
  WORKDIR <- "/media/bzniks/DISK3/Ants_behaviour_analysis"
  DATADIR <- paste(WORKDIR,"Data",sep="/")
  SCRIPTDIR <- "/home/bzniks/Dropbox/SeniorLectureship_Bristol/Students_postdocs/PhD_students/2019 Adriano Wanderlingh/code/PhD-exp1-data-analysis-main/scriptsR/Behaviours_inferral_Nathalie_NEW"
  SAVEOUTPUT <- "/home/bzniks/Documents"
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
###READ WHOLE ANNOTATIONS DATASET #############################################
###############################################################################
print("Loading manual annotations...")
#the current annotation file FINAL_script_output_CROSSVAL_25PERC_AND_TROPH.csv underwent cross-validation by Adriano
annotations_all <- read.csv(paste(DATADIR,"/R3SP_R9SP_All_data_FINAL_script_output_CROSSVAL_25PERC_AND_TROPH.csv",sep = ""), sep = ",")

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

##############################################################################
######### SPLIT ANNOTATIONS DATASET INTO TEST AND TRAINING  ##################
##############################################################################
### Be careful about this - in the previous you had randomly allocated behaviours to test or training throughout the period
### This means true Hits were wrongly identified as misses when looking at automatic interactions, because the time span of automatic interaction detection was unchanged
### Instead you need to define contiguous periods of time that contain half the events, for each colony/period
print("Splitting manual annotations into test and training chunks...")
###first list nb of events for each behaviour, each replicate and each period
nb_events <- aggregate ( Column ~ Behaviour + PERIOD + REPLICATE, FUN=length, data=annotations_all)
###narrow down to behaviour of interest only
nb_events <- nb_events[which(nb_events$Behaviour==BEH),]

###initialise new test and training objects
annotations_training <- NULL
annotations_test     <- NULL
time_window_training <- NULL
time_window_test     <- NULL

###then loop over nb_events
for (i in 1:nrow(nb_events)){
  ###subset annotations_all to period/replicate of interest
  PERIOD    <- nb_events[i,"PERIOD"]
  REPLICATE <- nb_events[i,"REPLICATE"]
  annotation_subset <- annotations_all[which(annotations_all$Behaviour==BEH & annotations_all$PERIOD==PERIOD & annotations_all$REPLICATE==REPLICATE),]
  
  ### ensure annotation_subset is sorted by time start sec then define a time_limit in between two successive events corresponding to half the events for this period  
  annotation_subset <- annotation_subset[order(annotation_subset$T_start_UNIX),]
  time_limit <- min (c(as.numeric(annotation_subset[floor(nb_events[i,"Column"]/2),"T_stop_UNIX"]),as.numeric(annotation_subset[1+floor(nb_events[i,"Column"]/2),"T_start_UNIX"]) )) + (1/FRAME_RATE) * floor(round(abs(as.numeric(annotation_subset[floor(nb_events[i,"Column"]/2),"T_stop_UNIX"])-as.numeric(annotation_subset[1+floor(nb_events[i,"Column"]/2),"T_start_UNIX"]))/(1/FRAME_RATE))/2)          
  
  ###time_limit will work well if those two events don't overlap
  ###but if they do the corresponding lines will need to be duplicated with one event ending at time limit and one starting just after time_limit
  ###first list events to duplicate (if any)
  to_duplicate  <- annotation_subset[which(  (as.numeric(annotation_subset$T_start_UNIX)<=time_limit & as.numeric(annotation_subset$T_stop_UNIX)>time_limit) ),]
  
  ###perform duplication if necessary
  if (nrow(to_duplicate)>=1){
    to_keep_as_is <- annotation_subset[which(!(as.numeric(annotation_subset$T_start_UNIX)<=time_limit & as.numeric(annotation_subset$T_stop_UNIX)>time_limit) ),]
    to_duplicate_before             <-  to_duplicate
    to_duplicate_before$T_stop_UNIX <- as.POSIXct(time_limit,  origin="1970-01-01", tz="GMT" )
    to_duplicate_before$duration    <- as.numeric(to_duplicate_before$T_stop_UNIX - to_duplicate_before$T_start_UNIX)
    
    to_duplicate_after              <-  to_duplicate
    to_duplicate_after$T_start_UNIX <- as.POSIXct(time_limit,  origin="1970-01-01", tz="GMT" )+1/FRAME_RATE
    to_duplicate_after$duration     <- as.numeric(to_duplicate_after$T_stop_UNIX - to_duplicate_after$T_start_UNIX)
    
    annotation_subset <- rbind(to_keep_as_is,to_duplicate_before,to_duplicate_after)
    annotation_subset <- annotation_subset[order(annotation_subset$T_start_UNIX),]
  }
  
  ### now we split annotation_subset in two chunks before and after time_limit,
  ### randomly allocate each time chunk to test or training dataset,
  ### and store the time windows for those time chunks
  if (runif(1,0,1)<0.5){
    annotations_training <- rbind( annotations_training , annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),]  )
    annotations_test     <- rbind( annotations_test     , annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),])
    time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                     REPLICATE=REPLICATE,
                                                                     time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                     time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
    )
    )
    
    time_window_test     <- rbind( time_window_test     , data.frame(PERIOD=PERIOD, 
                                                                     REPLICATE=REPLICATE,
                                                                     time_start = as.POSIXct(time_limit + 1/FRAME_RATE,origin="1970-01-01", tz="GMT") - 1/FRAME_RATE, 
                                                                     time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
  }else{
    annotations_training <- rbind( annotations_training, annotation_subset[which(annotation_subset$T_start_UNIX>time_limit),]  )
    annotations_test     <- rbind( annotations_test,annotation_subset[which(annotation_subset$T_stop_UNIX<=time_limit),])
    
    time_window_training <- rbind( time_window_training , data.frame(PERIOD=PERIOD, 
                                                                     REPLICATE=REPLICATE,
                                                                     time_start = as.POSIXct(time_limit+1/FRAME_RATE,origin="1970-01-01", tz="GMT")- 1/FRAME_RATE,
                                                                     time_stop  = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_stop"]))
    
    
    time_window_test     <- rbind( time_window_test , data.frame(PERIOD=PERIOD, 
                                                                 REPLICATE=REPLICATE,
                                                                 time_start = time_window_all[which(time_window_all$PERIOD==PERIOD & time_window_all$REPLICATE==REPLICATE),"time_start"],
                                                                 time_stop = as.POSIXct(time_limit,origin="1970-01-01", tz="GMT") + 1/FRAME_RATE
    ))
    
  }
  
}

###finally, recalculate duration for all annotations_all objects and define new time columns express in seconds
for (annotation_object in c("annotations_all","annotations_training","annotations_test")){
  annot <- get(annotation_object)
  #convert Zulu time to GMT
  annot$duration      <- as.numeric(annot$T_stop_UNIX - annot$T_start_UNIX)
  # #transform zulu time in GMT
  # annot$T_start_UNIX <- as.POSIXct(annot$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  # annot$T_stop_UNIX  <- as.POSIXct(annot$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  #assign time in sec to avoid issues on time management and matching
  annot$T_start_sec <- round(as.numeric(annot$T_start_UNIX),N_DECIMALS)
  annot$T_stop_sec <- round(as.numeric(annot$T_stop_UNIX),N_DECIMALS)
  
  assign(annotation_object,annot)
}

#start fresh
Grooming_LDA_output             <- data.frame()
Grooming_LDA_eachRun            <- data.frame()

###initialise general output folder
###remove folder if already exists to amke sure we don't mix things up
if (file.exists(file.path(SAVEOUTPUT, "MachineLearning_outcomes"))){
  unlink(file.path(SAVEOUTPUT, "MachineLearning_outcomes"),recursive=T)
} 
###create folder
dir.create(file.path(SAVEOUTPUT, "MachineLearning_outcomes"),recursive = T)

###define name of general output table containing quality scores
output_name <- file.path(SAVEOUTPUT, "MachineLearning_outcomes","quality_scores.txt")


###############################################################################
###### OUTER PARAMETERS LOOP ##################################################
###############################################################################
#Varying capsule shapes
CAPSULE_FILE_LIST           <- c("CapsuleDef3","CapsuleDef4","CapsuleDef9","CapsuleDef10","CapsuleDef11","CapsuleDef12")

# #####Arguments to loop over - to comment out when running the loop
# CAPSULE_FILE                <- "BODY-HEAD_17March22_auto_oriented_meanlength"
# DT_dist_THRESHOLD           <- 0 # NOT HIGHER THAN 0.5 as it will cut a very large portion of data (most movements are on a very small scale) #tag length is 62 px approx (measured on full size pics in R9SP)
# MAX_INTERACTION_GAP         <- 10 ## in SECONDS, the maximum gap in tracking before cutting the trajectory or interactions in two different object
# DISAGREEMENT_THRESH         <- 0.4#Assign Hit or Miss after matrix difference according to threshold 
# trim_length_sec             <- 1  ####all automated interactions that last trim_length_sec or less will be removed from the automated detection to fit the LDA as they increase noise!
# DT_frame_THRESHOLD          <- 16  #arbitrary value #trajectories jumps/gaps thresholds to avoid getting skewed means summary values (see their use in params extraction scripts: BEH_Traj_from_Man_annotations_fort081.R and BEH_Parameters_Auto_fort081.R)
####
###beta: relates to the calculation of the Fbeta score (https://en.wikipedia.org/wiki/F-score)
###      This is how much we value recall (=sensitivity = TP / (TP+FN)) over precision (=positive predictive value = TP / (TP + FP))
### if beta = 1 the score calculated is a F1 score
### based on lab meeting discussion we felt that we valued precision (minimising FP) over sensitivity to decrease the noise / avoid detecting grooming when there are none
### so we want a beta < 1
### arbitrarily here, it is set at 0.5
# beta <- 0.5 


####initialise Loop_ID
Loop_ID       <- 1
Trunk_Loop_ID  <- 1
for (CAPSULE_FILE in CAPSULE_FILE_LIST) { #list of CAPUSLE FILES TO BE USED  ###NATH_FLAG: CAPSULE_FILE_LIST has not been defined
  #fmQueryComputeAntInteractions matcher for the max time interval in iteraction when the ant pair disengages the interaction. Specific for GROOMING
  #Sequentially vary the interaction gap-filling to check what effect this has on the agreement between the MANUAL & AUTOMATIC interactions
  for (MAX_INTERACTION_GAP in c(15,20,25)) { #IN SECONDS (5,10, 15 never selected)
    trunk_loop_start_time <- Sys.time()
    print(paste("TRUNK LOOP ID:",Trunk_Loop_ID))
    ##################################################################################
    ##################### EXTRACT VARIABLES FOR MANUAL AND AUTOMATIC INTERACTIONS ####
    ##################################################################################
    ###for these particular parameters, evaluate how successful the interaction detetection parameters are at detecting candidate frames
    ###extract interactions and manual annotations for the whole dataset
    print(paste("Evaluating quality of interaction parameters for trunk loop ID ",Trunk_Loop_ID,"...",sep=""))
    all      <- extraction_loop("all",extract_movement_variables=F)
    
    ###Run auto_manual_agreement on all data to get an idea of how well the Loop performs in discovering the candidate interaction frames
    loop_interaction_detection_TPTNFPFN <- auto_manual_agreement (all[["summary_AUTO"]] , all[["summary_MANUAL"]], all[["list_IF_Frames"]] )[["true_false_positive_negatives"]]

    ##THRESHOLD to exclude jitter in the individuals' movement (DISTANCE)
    for (DT_dist_THRESHOLD in c(0,0.2,0.4)) { #NOT HIGHER THAN 0.5 # tag length is 62 px approx (measured on full size pics in R9SP)
      # Assign Hit based on threshold
      # maybe can be put somewhere better as it involves a later stage of the analysis
        for (DT_frame_THRESHOLD in c(32,40)){
          
          ###for these particular parameters, calculate variables of interest for manual annotations and automatic interactions for each of training and test dataset
          print("Extracting movement variables for manually-annotated data and automatic interactions:")
          print("-training data...")
          training <- extraction_loop("training")
          print("-test data...")
          test     <- extraction_loop("test")
          
          ##############################################################################
          ######### MANUAL/AUTO AGREEMENT  ###################################
          ##############################################################################
          ###Run auto_manual_agreement on training data to define Hits and Misses
          training[["summary_AUTO"]] <- auto_manual_agreement (training[["summary_AUTO"]] , training[["summary_MANUAL"]], training[["list_IF_Frames"]] )[["summary_AUTO"]]
          
          ##############################################################################
          ######### FIT CLASSIFIERS WITH LDA/QDA/RF  ####################################
          ##############################################################################
          ###Note for Adriano:
          ###For the fit of the classifiers I believe a lot can be gained from trimming the shortest interactions from summary_AUTO, as these will introduce noise
          ####define to_keep variables to keep clearing memory between runs
          to_keep <- c(ls(),c("to_keep","trim_length_sec"))
          
          for (trim_length_sec in c(2,3)){##This level of the loop is AFTER extracting movement variables for training and test to avoid unnecessary repetition of this slow step.
            
            ###prepare output directories for Loop_ID
            subDir <- paste0("Loop_ID_",Loop_ID)
            if (file.exists(file.path(SAVEOUTPUT, "MachineLearning_outcomes",subDir))){
              unlink(file.path(SAVEOUTPUT, "MachineLearning_outcomes",subDir),recursive=T)
            }
            dir.create(file.path(SAVEOUTPUT, "MachineLearning_outcomes",subDir,"fits"),recursive = T)
            
            ###fit classifiers
            print(paste("Perform Classification for Loop_ID ",Loop_ID,"...",sep=""))
            classifiers <- fit_classifiers(
              trim_short_interactions ( ######trim_short_interactions: function that cuts interactions shorter than a certain duration
                training[["summary_AUTO"]] ######interaction table: summary_AUTO from training object 
                ,trim_length_sec ### maximum duration of interactions that will be trimmed (set to 0 to trim nothing)
                ,"duration_sec" ###name of the column that contains the duration of each interaction in seconds (could vary between object so specify it here)
                
              )
            )
            
            ###############################################################################
            ###### SAVING LOOP-RELEVANT OBJECTS        ####################################
            ###############################################################################
            #save the SIRUS rules object
            dput(classifiers[["SirusRules"]], file.path(SAVEOUTPUT, "MachineLearning_outcomes",subDir,"SirusRules.txt"))
            
            #save the BN_list object containing all information regarding the variables selected by RELIEF and the method to calculate them for new values
            dput(classifiers[["selected_variables_BNobject_list"]], file.path(SAVEOUTPUT, "MachineLearning_outcomes",subDir,"BN_object_list.dat"))
            
            ####################################################################################################################
            ###then for each classifier method, perform prediction on training and test data and get quality scores    #########
            ####################################################################################################################
            print(paste("Predict training and test data for Loop_ID ",Loop_ID,"...",sep=""))
            
            for (class_idx in 1:length(classifiers[["fits"]])){
              classifier <- list(classifiers[["fits"]][[class_idx]]); names(classifier) <- names(classifiers[["fits"]])[class_idx]
              
              ##check match with manual 
              for (beta in c(0.5)){
                quality_scores_pre_classifier       <- quality_scores(loop_interaction_detection_TPTNFPFN,beta)
                print(paste("Automatic interaction detectionin trunk loop ID ",Trunk_Loop_ID,"had a sensitivity of ",round(quality_scores_pre_classifier["sensitivity"],digits=2)," and a precision of ",round(quality_scores_pre_classifier["precision"],digits=2),".",sep=""))
                
                
                for (what in c("training","test")){
                  ###extract objects of interest
                  summary_AUTO   <- get(what)[["summary_AUTO"]]
                  summary_MANUAL <- get(what)[["summary_MANUAL"]]
                  list_IF_Frames <- get(what)[["list_IF_Frames"]]
                  
                  ###Perform prediction on training data
                  summary_AUTO["predicted_Hit"] <- predict_class  (summary_AUTO =    summary_AUTO
                                                                   ,BN_list      =  classifiers[["selected_variables_BNobject_list"]]
                                                                   ,classifier   = classifier
                  )
                  
                  
                  
                  
                  
                  if (sum(summary_AUTO$predicted_Hit,na.rm=T)>0){
                    aut_man_agreement <- auto_manual_agreement (summary_AUTO[which(summary_AUTO$predicted_Hit==1),]
                                                                , summary_MANUAL
                                                                , list_IF_Frames
                    )
                    assign (paste("true_false_positive_negatives",what,sep="_"),aut_man_agreement[["true_false_positive_negatives"]])
                    assign (paste("quality_scores",what,sep="_"),round(quality_scores(aut_man_agreement[["true_false_positive_negatives"]],beta),digits=3 ))
                  }else{
                    assign (paste("quality_scores",what,sep="_"),round(c(CSI=0,Fbeta=0,precision=0,sensitivity=0),digits=3 ))
                    assign (paste("true_false_positive_negatives",what,sep="_"),data.frame(true_negatives=NA,true_positives=0,false_negatives=NA,false_positives=0))
                  }
                  
                }
                #create a row per each Loop run 
                Grooming_LDA_eachRun  <- data.frame(Loop_ID=Loop_ID,
                                                    CAPSULE_FILE=CAPSULE_FILE,
                                                    DT_dist_THRESHOLD=DT_dist_THRESHOLD, 
                                                    MAX_INTERACTION_GAP=MAX_INTERACTION_GAP,
                                                    DISAGREEMENT_THRESH=DISAGREEMENT_THRESH,
                                                    DT_frame_THRESHOLD = DT_frame_THRESHOLD,
                                                    trim_length_sec=trim_length_sec,
                                                    beta=beta,
                                                    classifier=names(classifier),
                                                    CSI_pre_classifier = quality_scores_pre_classifier["CSI"],
                                                    Fbeta_pre_classifier  = quality_scores_pre_classifier["Fbeta"],
                                                    precision_pre_classifier = quality_scores_pre_classifier["precision"],
                                                    sensitivity_pre_classifier  = quality_scores_pre_classifier["sensitivity"],                                                  
                                                    CSI_training = quality_scores_training["CSI"],
                                                    Fbeta_training  = quality_scores_training["Fbeta"],
                                                    precision_training = quality_scores_training["precision"],
                                                    sensitivity_training  = quality_scores_training["sensitivity"],
                                                    CSI_test = quality_scores_test["CSI"],
                                                    Fbeta_test  = quality_scores_test["Fbeta"],
                                                    precision_test = quality_scores_test["precision"],
                                                    sensitivity_test  = quality_scores_test["sensitivity"],
                                                    stringsAsFactors = F,row.names = NULL)
                
                ###Add this line in general output table
                Grooming_LDA_output   <- rbind(Grooming_LDA_output,     Grooming_LDA_eachRun)
                
                ###############################################################################
                ###### SAVING CLASSIFIER-RELEVANT OBJECTS        ##############################
                ###############################################################################
                ###1. save new line in quality score object, then clear it
                if (file.exists(output_name)){
                  write.table(Grooming_LDA_eachRun,file=output_name,append=T,col.names=F,row.names=F,quote=T)
                }else{
                  write.table(Grooming_LDA_eachRun,file=output_name,append=F,col.names=T,row.names=F,quote=T)
                }
              }

              ###2. save fit
              saveRDS(classifier[[names(classifier)]], file.path(SAVEOUTPUT, "MachineLearning_outcomes",subDir,"fits",paste(names(classifier),".rds",sep="")))
              
              
            }#class_idx
            
            ###add 1 to Loop_ID counter
            Loop_ID <- Loop_ID +  1
            rm(   list   =  ls()[which(!ls()%in%to_keep)]    )
            gc()
          }#trim_length_sec
        }#DT_frame_THRESHOLD

    }#DT_dist_THRESHOLD
    Trunk_Loop_ID <- Trunk_Loop_ID +1
    trunk_loop_end_time <- Sys.time()
    print (paste("Trunk loop took ",trunk_loop_end_time-trunk_loop_start_time," seconds to complete"))
  }#MAX_INTERACTION_GAP
}#CAPSULE_FILE

