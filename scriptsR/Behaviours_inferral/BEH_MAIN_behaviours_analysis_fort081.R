##########################################################################################
############## BEH MAIN Behaviours Analysis ##############################################
##########################################################################################

#### THIS VERSION IS FORT 0.8.1 COMPATIBLE ####

#this should be the version of the script maintained for long term use.
#For previous versions of this script and the sourced one, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral

#NOTE: at the current time (24 March 2022), the computational part of the script works but the saving of the plots does not (needs to be fixed, possibly according to plots_structure.R)

#clean start
rm(list=(c("e")))
gc()
rm(list=ls())

###############################################################################
###### DEFINE MY FUNCTIONS ####################################################
###############################################################################

#Movement angle diff
Movement.angle.diff     <- function(x){
  if( length(x) > 0 & !is_null(x)){
  c( atan(x[-nrow(x), "x"] / x[-nrow(x), "y"]) - atan(x[-1, "x"] / x[-1, "y"]), NA)}
    else {vector()}}

#pi-wrap
pi_wrap <- function ( angle){ ####function that force any angle value to lie between -pi and + pi
  while(angle>pi){
    angle <- angle - 2*pi
  }
  while(angle<=-pi){
    angle <- angle + 2*pi
  }
  return(angle) 
}

####inclination_angle = function that converts an oriented angle difference between two vectors to an inclination angle (i.e. acute angle between the inclinations of the two lines on which the vectors lie)
#### values will range from 0 (parallel lines) to pi/2 (orthogonal lines)
inclination_angle <- function(oriented_angle){
  abs_angle <- abs(pi_wrap(oriented_angle))
  return (min (abs_angle, pi-abs_angle))
}

#change system of coords

#Orientation.diff 
Orientation.diff        <- function(x)  { c( x[-nrow(x), "angle"] - x[-1, "angle"], NA)}


## custom arrows function 
## quiver plot for vectors
arrows.az <- function(x, y, azimuth, rho, HeadWidth, ..., units=c("degrees", "radians"), Kol, Lwd)
{
  units <- match.arg(units)
  az.rad <- switch(units,
                   degrees=azimuth*pi/180,
                   radians=azimuth)
  arrows(x, y, x+cos(az.rad)*rho, y+sin(az.rad)*rho, length=HeadWidth,
         ..., col=Kol, lwd=Lwd)
}

###############################################################################
###### LOAD LIBRARIES #########################################################
###############################################################################
library(adehabitatHR) ####for home range and trajectory analyses
library(FortMyrmidon) ####R bindings https://formicidae-tracker.github.io/myrmidon/latest/index.html
library(igraph)       ####for network analysis
library(parsedate)
library (trajr)
library(plotrix) 
library(circular) #to work with circular data. objects not in circular class are coerced
library(tidyverse)
library(ggplot2)
library(reshape2) #to use melt and similar
library(bit64)
library(nanotime)
library(MALDIquant)
library(stringr)
library(data.table)
library(fields)
library(sp) #calculate convex hull area
#library(BAMBI) #angles wrapping

library(Matrix) #behavioural matrices subtraction

#LDA, PCA, etc
library(FactoMineR)
library(factoextra)
library(missMDA) #PCA with missing values
library(corrplot)
require(MASS)
require(ggplot2)
require(scales)
require(gridExtra)
library(reshape2)
library(ggbeeswarm)
library(GGally) #plot multicollinearity
library(bestNormalize)
library(car) # homogeinity of variance
library(MVN) # test multivariate normality ditribution
library(biotools) #BoxM test for homoscedasticity (maybe is better t use the Bartlett sphericity test)
library(CORElearn) # for RELIEF algorithm
#install.packages("remotes"); remotes::install_github("insilico/STIR")
library(stir) # calculation of statistical significance of features and adjustment for multiple testing of Relief-based scores
library(heplots) #Hypotesis testing plots to check Equal Covariance Ellipses in LDA vars assumptions check
library(devtools)
#RandForOld <- "https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz"
#install.packages(RandForOld, repos=NULL, type="source") 
library(randomForest)
#devtools::install_github("aflapan/sparseKOS")
library(sparseKOS) # sparseKOS, unlikely used for issues reported in the TASKs SHEET
#devtools::install_github("dongyuanwu/RSBID")
library(RSBID) # Resampling Strategies for Binary Imbalanced Datasets
library(caret) #to tune RandomForest Hyperparameters
library(sirus) # Stable and Interpretable RUle Set for RandomForests

###############################################################################
###### GLOSSARY ###############################################################
###############################################################################

# MAN         : manual interactions/trajectories deriving from the hand annotated data (from the file annotations)
# AUTO        : automatically extracted interactions/trajectories deriving from fmQueryComputeAntInteractions
# REPLICATE   : nest ID, either "R3SP" or "R9SP" (SP: small pathogen)
# PERIOD      : the treatment period, either "pre" or "post" pathogen exposure
# REP_PER     : each of the 4 blocks analised, crossing the REPLICATE and PERIOD

###############################################################################
###### SOURCES ON/OFF #########################################################
###############################################################################
#sources (calling outer scripts)
run_collisions                      <- TRUE #this controls the fmQueryCollideFrames, BEH_collisions_fort081.R and some extra calculations and plotting on Caspules distribution in MAN interactions
run_AUTO_MAN_agreement              <- TRUE # computes agreement matrix
create_interaction_AUTO_REP_PER     <- FALSE 
LDA_TP_FP_AUTO                      <- TRUE # calls the LDA analysis script, called BEH_PCA_fort081.R

#plots
run_Parameters_plots                <- FALSE # summary_MANUAL plots
plot_all_BEH                        <- FALSE # inside run_Parameters_plots, if FALSE plots only for Grooming (to run it an update in BEH statement is needed)


###############################################################################
###### PARAMETERS #############################################################
###############################################################################
USER <- "Adriano"
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

#currently unused
desired_step_length_time    <- 0.125 ###in seconds, the desired step length for the analysis

#trajectories jumps/gaps thresholds to avoid getting skewed means summary values (see their use in params extraction scripts: BEH_Traj_from_Man_annotations_fort081.R and BEH_Parameters_Auto_fort081.R)
DT_frame_THRESHOLD          <- 16  #arbitrary value
DT_dist_THRESHOLD           <- 0.3 # NOT HIGHER THAN 0.5 as it will cut a very large portion of data (most movements are on a very small scale) #tag length is 62 px approx (measured on full size pics in R9SP)

#FMmatchers specific for GROOMING used in script testing, some will be called in the  "OUTER PARAMETERS LOOP"
MAX_INTERACTION_GAP         <- 10 ## in SECONDS, the maximum gap in tracking before cutting the trajectory or interactions in two different object
AntDistanceSmallerThan      <- 300 #for higher accuracy, recalculate it from: max(interaction_MANUAL$straightline_dist_px,na.rm = T)
AntDistanceGreaterThan      <- 70 #for higher accuracy, recalculate it from: min(interaction_MANUAL$straightline_dist_px,na.rm = T)
ANT_LENGHT_PX               <- 153 #useful for matcher::meanAntDisplacement mean and median value are similar
MAX_DISPLACEMENT            <- fmSecond(0.5) ##check every X seconds if ant has displaced more than ANT_LENGHT_PX

#Assign Hit or Miss after matrix difference according to threshold
DISAGREEMENT_THRESH <- 0.4

###############################################################################
###### LOAD ANNOTATION AND MYRMIDON FILES #####################################
###############################################################################

####### navigate to folder containing myrmidon file
if (USER=="Adriano") {WORKDIR <- "/home/cf19810/Documents/Ants_behaviour_analysis"}
if (USER=="Tom")     {WORKDIR <- "/media/tom/MSATA/Dropbox/Ants_behaviour_analysis"}

DATADIR <- paste(WORKDIR,"Data",sep="/")
SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")

## perform analyses on annotation data
#source(paste(SCRIPTDIR,"BEH_Annotations_analysis.R",sep="/"))
# OR
#load only the training dataset, produced as output from the BEH_Annotations_analysis.R script (check the last few lines of it)
annotations <- read.csv(paste(DATADIR,"/annotations_TRAINING_DATASET.csv",sep = ""), sep = ",")
#transform zulu time in GMT
annotations$T_start_UNIX <- as.POSIXct(annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
annotations$T_stop_UNIX  <- as.POSIXct(annotations$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
#assign time in sec to avoid issues on time management and matching
annotations$T_start_sec <- round(as.numeric(annotations$T_start_UNIX),N_DECIMALS)
annotations$T_stop_sec <- round(as.numeric(annotations$T_stop_UNIX),N_DECIMALS)

#start fresh
Grooming_LDA_output             <- data.frame()


###############################################################################
###### OUTER PARAMETERS LOOP ##################################################
###############################################################################

CAPSULE_FILE <- NA #placeholder
Loop_ID <- 0
# #Varying capule shapes
# for (CAPSULE_FILE in vector) { #list of CAPUSLE FILES TO BE USED
#   ##THRESHOLD to exclude jitter in the individuals' movement (DISTANCE)
#   for (DT_dist_THRESHOLD in c(0.3,0.5)) { #NOT HIGHER THAN 0.5 # tag length is 62 px approx (measured on full size pics in R9SP)
#     #fmQueryComputeAntInteractions matcher for the max time interval in iteraction when the ant pair disengages the interaction. Specific for GROOMING
#     #Sequentially vary the interaction gap-filling to check what effect this has on the agreement between the MANUAL & AUTOMATIC interactions
#     for (MAX_INTERACTION_GAP in c(5,10)) { #IN SECONDS
#       # Assign Hit based on threshold
#       # maybe can be put somewhere better as it involves a later stage of the analysis 
#       for (DISAGREEMENT_THRESH in c(0.4,0.2)) {
#         
#       }}}} #these parenteses should go at the end of the script - here as placeholders-

Loop_ID <- Loop_ID +1

#set plots parameters (for plotting coords)
#needs to be fixed (should be included in the plotting structure outlined in plots_structure.R)
# plots_structure.R is a sketch on how plots can be distributed and activated.
pdf(file=paste(DATADIR,"Interactions_plots_8feb2022.pdf", sep = ""), width=6, height=4.5)
par(mfrow=c(2,3), mai=c(0.3,0.3,0.4,0.1), mgp=c(1.3,0.3,0), family="serif", tcl=-0.2)

### start fresh
interaction_MANUAL      <- NULL
interaction_MANUAL_COLL <- NULL
summary_MANUAL          <- NULL
interaction_AUTO        <- NULL
summary_AUTO            <- NULL

#check how long it takes to run the inner loop, consisting of the 4 REP_PER and the LDA
start.loop.time <- Sys.time()
for (REPLICATE in c("R3SP","R9SP")) 
  {
  ###############################################################################
  ###### OPEN EXPERIMENT INFORMATION ############################################
  ###############################################################################
  
  ## locate the ant info file for REPLICATE
  MyrmidonCapsuleFile <- list.files(path=DATADIR, pattern=REPLICATE, full.names=T)
  MyrmidonCapsuleFile <- MyrmidonCapsuleFile[grepl("BODY-HEAD_17March22_auto_oriented_meanlength.myrmidon",MyrmidonCapsuleFile)]

  e <- fmExperimentOpen(MyrmidonCapsuleFile)
  fmQueryGetDataInformations(e) # $details$tdd.path
  #attempts at extracting the file name to use it during plotting/file saves.
  # TO BE WORKED ON
  experiment_name <- unlist(strsplit(MyrmidonCapsuleFile,split="/"))[length(unlist(strsplit(MyrmidonCapsuleFile,split="/")))]
  experiment_name_end <-  str_sub(experiment_name,-30,-1)
  CapsuleDef <-  substr(experiment_name_end, 1, 7)
  ###tag statistics
  tag_stats <- fmQueryComputeTagStatistics(e)
  
  ################################################################################
  ########### START PERIOD LOOP ##################################################
  ################################################################################
  
  for (PERIOD in c("pre","post"))
    {
    print(paste("Replicate",REPLICATE, PERIOD))
    ## Prepare empty within-replicate/period data object
    interacts_MAN_REP_PER       <- NULL  
    summary_MAN_REP_PER         <- NULL
    interaction_AUTO_REP_PER    <- NULL
    summary_AUTO_REP_PER        <- NULL
    ## set experiment time window 
    time_start <- fmTimeCreate(min(annotations$T_start_UNIX[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD])) ###experiment start time
    #time_stop <- fmTimeCreate(min(annotations$T_start_UNIX[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD])+3) ###experiment start time
    #time_stop  <- fmTimeCreate(min(annotations$T_start_[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD]) + (34*60)  ) ###experiment stop time ####arbitrary time in the correct format + (N mins * N seconds)
    time_stop  <- fmTimeCreate(max(annotations$T_stop_UNIX[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD]) ) ###experiment stop time ####arbitrary time in the correct format + (N mins * N seconds)
  
    ###############################################################################
    ###### IDENTIFY FRAMES ########################################################
    ###############################################################################
    
    #Because of various issues raised, including the ones reported here https://github.com/formicidae-tracker/myrmidon/issues/240 ,
    #the analysis is performed not using UNIX_time but using frames. Every new queried file will then show a part including FRAME assignment.
    #Frames are then used for any trajectory cutting later on in the sub-scripts
    IdentifyFrames      <- fmQueryIdentifyFrames(e,start=time_start, end=time_stop)
    IF_frames           <- IdentifyFrames$frames
    # Assign a frame to each time since start and use it as baseline for all matching and computation
    IF_frames$frame_num <- as.numeric(seq.int(nrow(IF_frames)))
    
    # assign a time in sec to match annotation time (UNIX hard to use for class mismatch)
    IF_frames$time_sec <- round(as.numeric(IF_frames$time),N_DECIMALS)
 
    #assign frame numbering to annotations 
    annotations$frame_start <- match.closest(x = annotations$T_start_sec, table = IF_frames$time_sec, tolerance = 0.05)
    annotations$frame_stop  <- match.closest(x = annotations$T_stop_sec,  table = IF_frames$time_sec, tolerance = 0.05)
    
    # Creating a new zeroed-time since the start of the exp by  summing the cumulated differences between each pair of consecutive frames (to account for the time mismatch)
    IF_frames$cum_diff <- c(0, with(IF_frames, diff(time)))
    IF_frames$cum_diff <- cumsum(IF_frames$cum_diff)

    ###############################################################################
    ###### READING COLLISIONS #####################################################
    ###############################################################################
    if (run_collisions)
      {
      collisions <- fmQueryCollideFrames(e, start=time_start, end=time_stop)   ###collisions are for each frame, the list of ants whose shapes intersect one another. Normally not used
      collisions$frames$frame_num <- seq.int(nrow(IF_frames)) 
      }
    
    ###############################################################################
    ###### READING TRAJECTORIES ###################################################
    ###############################################################################
    #COMPUTE TRAJECTORIES  
    positions <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = TRUE) #set true to obtain the zone of the ant
    positions$trajectories_summary$frame_num <- NA
    #assign starting frame number
    positions$trajectories_summary["frame_num"] <- lapply(positions$trajectories_summary["start"], function(x) IF_frames$frame_num[match(x, IF_frames$time)])
    
    ## immediately after computing your positions object:
    positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
    names(positions$trajectories)       <- positions$trajectories_summary$antID_str ###and use the content of that column to rename the objects within trajectory list

    trajectories_summary <- positions$trajectories_summar
    
    ## Add ant_x names and times to the positions to convert from "FRAME since the start of the experiment", to FRAMES
    for (A in positions$trajectories_summary$antID)
    {
      AntID                         <- paste("ant_",A,sep="") 
      First_Obs_Frame               <- positions$trajectories_summary$frame_num [which(positions$trajectories_summary$antID_str==AntID)]
      IF_frames$new_zero_diff  <- IF_frames$cum_diff - IF_frames[IF_frames$frame_num==First_Obs_Frame,"cum_diff"] #subtracting the $frames zeroed-time  corresponding to the $start time from the zeroed-time column itself (New-Zeroed-time)
      print(paste("Adding first obs FRAME", First_Obs_Frame, "to the time-zeroed trajectory of ant", AntID))
      #assign corresponding frame N when the New-Zeroed-time and $time correspond, closest.match 0.05 (well inside the 0.125 frame length in sec)
      positions$trajectories[[AntID]]$frame <- match.closest(x = positions$trajectories[[AntID]]$time, table = IF_frames$new_zero_diff, tolerance = 0.05)
      IF_frames$new_zero_diff <- NA
    }
    AntID <- NULL; First_Obs_Frame <- NULL
    # the older version used UNIX_time
    # ## Add ant_x names and times to the positions to convert from time since the start of the experiment, to UNIX time
    # for (A in positions$trajectories_summary$antID)
    #   {
    #   AntID <- paste("ant_",A,sep="") 
    #   First_Obs_Time <- as.POSIXct(positions$trajectories_summary$start [which(positions$trajectories_summary$antID_str==AntID)], tz="GMT",origin = "1970-01-01 00:00:00") ## find the first time after the user defined time_start_ISO that this ant was seen
    #   print(paste("Adding first obs time", First_Obs_Time, "to the time-zeroed trajectory of ant", AntID))
    #   positions$trajectories[[AntID]] $UNIX_time <- as.POSIXct(positions$trajectories[[AntID]]$time, tz="GMT",origin = "1970-01-01 00:00:00")  + First_Obs_Time ##convert back to UNIX time  
    # }

    ######################################################################################
    ###### EXTRACT ANT TRAJECTORIES FROM MANUALLY-ANNOTATED DATA AND CALC PARAMETERS #####
    ######################################################################################
    source(paste(SCRIPTDIR,"BEH_Traj_from_Man_annotations_fort081.R",sep="/"))
    
    ## stack
    interaction_MANUAL <- rbind(interaction_MANUAL, interacts_MAN_REP_PER)
    summary_MANUAL     <- rbind(summary_MANUAL,       summary_MAN_REP_PER)
    
    #########################################################################################
    ###### READING AUTOMATIC INTERACTIONS ###################################################
    #########################################################################################
    
    #Get info on the capules shapes and use the relevant ones in the ant interaction query
    capsules  <- e$antShapeTypeNames
    names(capsules) <- as.character( 1:length(capsules))
    head_id <- as.numeric(names(capsules)[[which(capsules=="head")]])
    body_id <- as.numeric(names(capsules)[[which(capsules=="body")]])
    
    # body_id <- capsules[which(capsules$name=="body"),"typeID"]
    matcherCapType <- fmMatcherInteractionType(body_id,head_id)
    matcherCapTypeAntDists <- fmMatcherAnd(list(fmMatcherInteractionType(body_id,head_id),
                                 fmMatcherAntDistanceSmallerThan(AntDistanceSmallerThan),
                                 fmMatcherAntDistanceGreaterThan(AntDistanceGreaterThan),
                                 fmMatcherAntDisplacement(ANT_LENGHT_PX, MAX_DISPLACEMENT)))     #check every x seconds if ant has displaced more than ANT_LENGHT_PX
    # AN EXTRA 2 THAT CAN BE USED ARE #fmMatcherAntAngleSmallerThan(), fmMatcherAntAngleGreaterThan() 
    
    #note: adding AntDistances seems to reduce the false positives rate but has no effect on false negatives

      interacts_AUTO_REP_PER <- fmQueryComputeAntInteractions(e,
                                                              start=time_start, 
                                                              end=time_stop,
                                                              maximumGap =fmSecond(MAX_INTERACTION_GAP), ## WHEN A PAIR DISENGAGE, HOW LONG IS THE INTERVAL? 
                                                              reportFullTrajectories = T,
                                                              matcher = matcherCapTypeAntDists)
      

      #assign starting frame number to $trajectories
      interacts_AUTO_REP_PER$trajectories_summary["frame_num"]  <- lapply(interacts_AUTO_REP_PER$trajectories_summary["start"], function(x) IF_frames$frame_num[match(x, IF_frames$time)])
      ## immediately after computing your interacts_AUTO_REP_PER object:
      ## assign Ant_ID_str and Ant_ID_str+rowname_ID
      interacts_AUTO_REP_PER$trajectories_summary$antID_str     <- paste("ant_",interacts_AUTO_REP_PER$trajectories_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
      interacts_AUTO_REP_PER$trajectories_summary$ant.row.index <- paste0(interacts_AUTO_REP_PER$trajectories_summary$antID_str,"row",rownames(interacts_AUTO_REP_PER$trajectories_summary))
      names(interacts_AUTO_REP_PER$trajectories)                <- paste0(interacts_AUTO_REP_PER$trajectories_summary$antID_str,"row",rownames(interacts_AUTO_REP_PER$trajectories_summary)) ###and use the content of that column to rename the objects within trajectory list
      #names(interacts_AUTO_REP_PER$trajectories)       <- interacts_AUTO_REP_PER$trajectories_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
      ## Add ant_x names and times to the interacts_AUTO_REP_PER to convert from FRAME since the start of the experiment, to FRAMES
      ##creates a ID string for each ant in $interactions
      interacts_AUTO_REP_PER$interactions$ant1ID_str            <- paste("ant_",interacts_AUTO_REP_PER$interactions$ant1,sep="")
      interacts_AUTO_REP_PER$interactions$ant2ID_str            <- paste("ant_",interacts_AUTO_REP_PER$interactions$ant2,sep="")
      
      for (A_Row in interacts_AUTO_REP_PER$trajectories_summary$ant.row.index)
      {
        #ANTID_STR                         <- paste("ant_",interacts_AUTO_REP_PER$trajectories_summary$antID[i],sep="") 
        First_Obs_Frame               <- interacts_AUTO_REP_PER$trajectories_summary$frame_num [which(interacts_AUTO_REP_PER$trajectories_summary$ant.row.index==A_Row)]
        #min(interacts_AUTO_REP_PER$trajectories_summary$frame_num [which(interacts_AUTO_REP_PER$trajectories_summary$antID_str==ANTID_STR)])
        IF_frames$new_zero_diff  <- IF_frames$cum_diff - IF_frames[IF_frames$frame_num==First_Obs_Frame,"cum_diff"] #subtracting the $frames zeroed-time  corresponding to the $start time from the zeroed-time column itself (New-Zeroed-time)
        print(paste("Adding first obs FRAME", First_Obs_Frame, "to the time-zeroed trajectory of ant", A_Row))
        #assign corresponding frame N when the New-Zeroed-time and $time correspond, closest.match 0.05 (well inside the 0.125 frame length in sec)
        interacts_AUTO_REP_PER$trajectories[[A_Row]]$frame <- match.closest(x = interacts_AUTO_REP_PER$trajectories[[A_Row]]$time, table = IF_frames$new_zero_diff, tolerance = 0.05)
        IF_frames$new_zero_diff <- NA
      } 
      AntID <- NULL; First_Obs_Frame <- NULL
     
      # Assign frame number to $interactions
      interacts_AUTO_REP_PER$interactions["int_start_frame"]   <- lapply(interacts_AUTO_REP_PER$interactions["start"] , function(x) IF_frames$frame_num[match(x, IF_frames$time)])
      interacts_AUTO_REP_PER$interactions["int_end_frame"]   <- lapply(interacts_AUTO_REP_PER$interactions["end"] , function(x) IF_frames$frame_num[match(x, IF_frames$time)])
      # Assign interaction pair
      interacts_AUTO_REP_PER$interactions $pair <- paste(interacts_AUTO_REP_PER$interactions$ant1, interacts_AUTO_REP_PER$interactions$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...
      #calc duration (not including extremes,e.g. end frame not included)
      interacts_AUTO_REP_PER$interactions$Duration <- interacts_AUTO_REP_PER$interactions$int_end_frame - interacts_AUTO_REP_PER$interactions$int_start_frame
      #hist(interacts_AUTO_REP_PER$interactions$Duration,breaks = 30)


      ######################################################################################
      ###### CALC PARAMETERS FOR AUTOMATIC INTERACTIONS ####################################
      ######################################################################################
      source(paste(SCRIPTDIR,"BEH_Parameters_Auto_fort081.R",sep="/"))
      
      # #interesting plot for visualisation, maybe can be placed elsewhere...
      # #pdf(file=paste(DATADIR,"Interactions","AUTO_MAN_REP_PER","Gap",MAX_INTERACTION_GAP,"matcher_CapTypeAntDistsDispl.pdf", sep = "_"), width=10, height=60)
      # #plot the interactions by pair as timeline
      # #generate 1 plot per every iteration 
      # summary_MAN_REP_PER_sub <- summary_MAN_REP_PER[,c("REPLICATE","PERIOD","pair","int_start_frame","int_end_frame")] ; summary_MAN_REP_PER_sub$flag <- "manual"
      # interacts_AUTO_REP_PER_sub <- interacts_AUTO_REP_PER$interactions[,c("pair","int_start_frame","int_end_frame")] ; interacts_AUTO_REP_PER_sub$flag <- "auto"
      # AUTO_MAN_REP_PER <- dplyr::bind_rows(summary_MAN_REP_PER_sub,interacts_AUTO_REP_PER_sub) #use bind_rows to keep rep info https://stackoverflow.com/questions/42887217/difference-between-rbind-and-bind-rows-in-r
      # 
      # ggplot(AUTO_MAN_REP_PER) +
      #   geom_linerange(aes(y = pair, xmin = int_start_frame, xmax = int_end_frame,colour = flag),size=3,alpha = 0.5) +
      #   theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
      #   labs(title = paste("Grooming by pair" ,unique(REPLICATE),unique(PERIOD),"- MAN vs AUTO"),
      #        subtitle = paste("Nrows auto detected:" ,NROW(interacts_AUTO_REP_PER$interactions),"Nrows manual annotated:" ,NROW(summary_MAN_REP_PER),"\nMaxIntGap",MAX_INTERACTION_GAP, "s","Capsule file =",CapsuleDef )) +
      #   scale_color_manual(values = c("manual" = "red",
      #                                  "auto"="black"))

    ###############################################################################
    ###### AUTO_MAN AGREEMENT MATRIX ##############################################
    ###############################################################################
    if (run_AUTO_MAN_agreement){source(paste(SCRIPTDIR,"BEH_Auto_Man_agreement_matrix_fort081.R",sep="/"))}
    
    #stack output of BEH_Parameters_Auto_fort081.R + Hit information if run_AUTO_MAN_agreement is TRUE
    interaction_AUTO <- rbind(interaction_AUTO,interaction_AUTO_REP_PER)
    summary_AUTO <- rbind(summary_AUTO, summary_AUTO_REP_PER)
    
    
    ###############################################################################
    ###### COLLISIONS #############################################################
    ###############################################################################
    if (run_collisions){source(paste(SCRIPTDIR,"BEH_collisions_fort081.R",sep="/"))
      
      #stack output of collisions
      # THIS OUTPUT SHOULD BE PROBABLY SAVED INTO interaction_MANUAL DIRECTLY (AS IT IS JUST interaction_MANUAL + CAPSULE INFO)
      # BUT WHEN DOING SO CHECK THAT THIS DOES NOT HURT THE FOLLOWING interaction_AUTO FILE
      interaction_MANUAL_COLL <- rbind(interaction_MANUAL_COLL, interacts_MAN_REP_PER)
      
      }# run_collisions
    
    #SAVE AND KEEP THE VARIOUS PRODUCED DATAFRAMES (outside of the REPLICATE and PERIOD)
    # the following two objects are used in the CSI score determination
    assign(paste0("IF_frames","_", REPLICATE,PERIOD), IF_frames) 
    assign(paste0("int_mat_manual","_", REPLICATE,PERIOD), int_mat_manual)
    
  }##PERIOD
}##REPLICATE
dev.off()

##############################################################################
######### LDA ANALISYS  ######################################################
##############################################################################

if (LDA_TP_FP_AUTO){source(paste(SCRIPTDIR,"BEH_PCA_fort081.R",sep="/"))}

#######################################################################################
#### PLOTS FOR TOTAL VARS (are also in PCA_fort081.R but may be moved here for good) ##
#######################################################################################

if (run_Parameters_plots){source(paste(SCRIPTDIR,"BEH_Parameters_plots_fort081.R",sep="/"))}

###############################################################################
###### AUTO-MAN DISAGREEMENT PLOT #############################################
###############################################################################
if (run_AUTO_MAN_agreement)
  {
  ## visualise the cutoffs determined in the agreement_matrix script
  ## Make y-axis logarithmic in histogram and fix issues for 0 occurences with log
  par(mfrow=c(3,1))
  #all
  summary_AUTO$disagreement <- as.numeric(as.character(summary_AUTO$disagreement))
  hist.data = hist(summary_AUTO$disagreement, breaks=seq(-1,0,0.05), plot=F)
  hist.data$counts <- replace(log10(hist.data$counts), log10(hist.data$counts)==-Inf, 0)
  highestCount <- max(hist.data$counts)
  plot(hist.data, main="All disagreements (log)",
       sub = paste("N.AUTO TOT",length(summary_AUTO$Hit)),
       ylab='log10(Frequency)', xlab=" AUTO-MAN disagreement rate ", ylim=c(0,highestCount)); abline(v=-DISAGREEMENT_THRESH, lty=2)
  
  #True Positives
  hist.data1 = hist(summary_AUTO$disagreement[which(summary_AUTO$Hit==1)], breaks=seq(-1,0,0.05), plot=F)
  hist.data1$counts <- replace(log10(hist.data1$counts), log10(hist.data1$counts)==-Inf, 0)
  plot(hist.data1, main="True Positives",
       sub = paste("N.AUTO",length(summary_AUTO$disagreement[which(summary_AUTO$Hit==1)]),"/ N.MANUAL",length(summary_MAN_REP_PER)),
       ylab='log10(Frequency)',xlab=" AUTO-MAN agreement rate ", col=2, ylim=c(0,highestCount), xlim=c(-1,0)); abline(v=-DISAGREEMENT_THRESH, lty=2) #breaks=seq(-1,1,0.05)
  
  #False Positives
  hist.data1 = hist(summary_AUTO$disagreement[which(summary_AUTO$Hit==0)], breaks=seq(-1,0,0.05), plot=F)
  hist.data1$counts <- replace(log10(hist.data1$counts), log10(hist.data1$counts)==-Inf, 0)
  plot(hist.data1, main="False Positives",
       sub = paste("N.AUTO",length(summary_AUTO$disagreement[which(summary_AUTO$Hit==0)])),
       ylab='log10(Frequency)', xlab=" AUTO-MAN disagreement rate ",col=4, ylim=c(0,highestCount),xlim=c(-1,0)); abline(v=-DISAGREEMENT_THRESH, lty=2) #breaks=seq(-1,1,0.05)
  
  mtext(paste("VARS: MaxIntGap",MAX_INTERACTION_GAP, "s",", Capsule file =",CapsuleDef ), side = 3, line = -1.5, outer = TRUE)
  
}


##############################################################################
######### CAPSULES CALCULATIONS ##############################################
##############################################################################
#Understand which capsules are used during the manual interaction

############# CHANGE ALL TO interaction_MANUAL once interaction_MANUAL_COLL is removed
############# CARE MUST BE TAKEN WHEN PRODUCING THE FOLLOWING PLOTS AS THE CAPSULE N. SHOULD BE SUBSTITUTED WITH THE CAPSULE NAME!

#ACCESS THE CAPSULE INFO
interaction_MANUAL_COLL$REP_PER_R_B <- paste(interaction_MANUAL_COLL$REPLICATE,interaction_MANUAL_COLL$PERIOD,interaction_MANUAL_COLL$ROW,interaction_MANUAL_COLL$BEH,sep="_")

#Check there is ANY overlap 
if (length(interaction_MANUAL_COLL$types[!is.na(interaction_MANUAL_COLL$types)]) > 0) {
  
  split_types <- plyr::ldply(strsplit(interaction_MANUAL_COLL$types,","), rbind)
  split_types$REP_PER_R_B <- interaction_MANUAL_COLL$REP_PER_R_B
  
  #select the unique TYPES combinations and make them into a dataframe to later add the missing cases to the counts
  #uniq.split_types <- unique(interaction_MANUAL_COLL$types); uniq.split_types <- uniq.split_types[!is.na(uniq.split_types)]
  uniq.split_types_FREQ <- as.data.frame(table(unlist(split_types[,-which(names(split_types)=="REP_PER_R_B")]))); names(uniq.split_types_FREQ) <- c("types","Freq")
  uniq.split_types <- as.data.frame(uniq.split_types_FREQ[,which(names(uniq.split_types_FREQ)=="types")]); names(uniq.split_types) <- "types"
  
  #REP_PER_R_B <- as.data.frame(unique(split_types$REP_PER_R_B)); names(REP_PER_R_B) <- "ID"
  int_types_counts_MAN <- NULL
  for (ids in unique(split_types$REP_PER_R_B)) {
    int_types <- split_types[which(split_types$REP_PER_R_B==ids),]
    count <- table(unlist(int_types[,-which(names(int_types)=="REP_PER_R_B")]))
    result <- data.frame(types = (names(count)),count = as.integer(count))
    
    ## add the missing cases
    int_types_counts <- plyr::join (x = result , y=uniq.split_types, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            
    int_types_counts[is.na(int_types_counts)] <- 0
    #calc pecerntage
    #int_types_counts$perc <- round(int_types_counts$count/sum(int_types_counts$count),2)
    int_types_counts$REP_PER_R_B <- ids
    
    int_types_counts_MAN   <- rbind(int_types_counts_MAN,     int_types_counts)
    #the loop is looping over the same REP_PER_R_B and should be written to avoid it. A dirty fix is:
  }
  int_types_counts_MAN <- unique(int_types_counts_MAN)
  
  #int_types_count_TOT <- as.data.frame(tapply(int_types_counts_MAN$count, INDEX=list(int_types_counts_MAN$types),FUN=sum))
  
  #get total counts and percentages
  int_types_counts_MAN.DT <- data.table(int_types_counts_MAN)
  #for plotting purposes, let's make 2-3 and 3-2 the same
  # CAREFUL HERE, USE THE CAPSULES NAMES NOT THEIR NUMBERS AS THEY MAY CHANGE!
  int_types_counts_MAN.DT$types_inv <- int_types_counts_MAN.DT$types
  int_types_counts_MAN.DT$types_inv <- as.character( int_types_counts_MAN.DT$types_inv)
  int_types_counts_MAN.DT$types_inv[ int_types_counts_MAN.DT$types_inv == "1-2"] <- "2-1"
  
  # plot per int
  ggplot( int_types_counts_MAN.DT, aes(fill=types_inv, y=count, x=REP_PER_R_B)) + 
    geom_bar(position="fill", stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("(CHECK IF INPUT FILE IS CORRECT) perc. capsule type involved in each interaction",paste("body_id is type",body_id,", head_id is type",head_id,sep=" "))
  
  #Getting % of capsules from collisions to understand which capsules are interacting more
  int_types_count_TOT <- int_types_counts_MAN.DT[, list(fsum=sum(count)), by=types]
  int_types_count_TOT$perc <- round(int_types_count_TOT$fsum/sum(int_types_count_TOT$fsum),2)
  int_types_count_TOT
  
  #Frame by Frame Number of Manual annotations that do not include 2-3 or 3-2 capsules
  types_FREQ <- as.data.frame(table(interaction_MANUAL_COLL$types))
  types_FREQ$body_head <- NULL
  types_FREQ$body_head <- ifelse(grepl('2-1|1-2', types_FREQ$Var1), 'yes', "no")
  types_FREQ$body_head
  aggregate(types_FREQ$Freq, by=list(body_head=types_FREQ$body_head), FUN=sum)
  
  #Interactionn by interaction Number of Manual annotations that do not include 2-3 or 3-2 capsules
  interaction_MANUAL_COLL$body_head <- NULL
  interaction_MANUAL_COLL$body_head <- ifelse(grepl('2-1|1-2', interaction_MANUAL_COLL$types), 'yes', "no")
  interaction_MANUAL_COLL$body_head <- as.factor(interaction_MANUAL_COLL$body_head)
  interaction_MANUAL_COLL$Freq <- 1
  body_head_AGGREG <- aggregate(interaction_MANUAL_COLL$Freq, by=list(interaction_MANUAL_COLL$body_head,interaction_MANUAL_COLL$REP_PER_R_B), FUN=sum);names(body_head_AGGREG) <- c("body_head","REP_PER_R_B","Freq")
  
  #Plot Perc. of Manual annotations containing head-body capules
  ggplot( body_head_AGGREG, aes(fill=body_head, y=Freq, x=REP_PER_R_B)) + 
    geom_bar(position="fill", stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Perc. of Manual annotations that include 2-1 or 1-2 capsules")
  
  #get from long to wide format (this can be appended to SUMMARY file)
  int_types_per_INT <- reshape(int_types_counts_MAN.DT, idvar = "REP_PER_R_B", timevar = "types", direction = "wide")
}


###################################################
#### SUMMARY MANUAL ACT REC ID BY PARAMETERS ######
###################################################
#create a copy of the summary_MANUAL dataframe
summary_MANUAL_delta <- summary_MANUAL

### To decide how to assign the ACT-REC ID to the AUTOMATICALLY EXTRACTED INTERACTIONS, we can look at how
### parameters vary for the MANUAL INTERACTIONS

#calculate delta parameters ACT-REC to check when ACT is bigger then REC
#note: the output var will be named ..ACT_delta but it is the diff between the the two individuals
summary_MANUAL_delta[, paste0(grep('_ACT$', colnames(summary_MANUAL), value = TRUE), '_delta')] <- 
summary_MANUAL[, grep('_ACT$', colnames(summary_MANUAL))] - summary_MANUAL[, grep('_REC$', colnames(summary_MANUAL))]

#plot delta parameters 
indx <- grepl('ACT_delta', colnames(summary_MANUAL_delta))
par(mfrow=c(2,3))
for (deltavar in names(summary_MANUAL_delta[indx])) {
plot(summary_MANUAL_delta[,deltavar] ~ rep(1:length(summary_MANUAL_delta[,deltavar])),
     main=deltavar,xlab="delta value", ylab="int num.",
     col=ifelse(summary_MANUAL_delta[,deltavar]>0,"black","red")) + 
  abline(h=0)
}

# #calculate delta vars manually CAN BE REMOVED
# #currently not working, the plot calls the wrong object...
# summary_MANUAL$delta_speed_pxpersec <- summary_MANUAL$mean_speed_pxpersec_ACT - summary_MANUAL$mean_speed_pxpersec_REC
# 
# #reshape data for plotting. Split by REC and ACT
# summary_data_ACT <-summary_MANUAL %>% dplyr::select(contains(c("BEH", "ACT"), ignore.case = TRUE))
# summary_data_REC <- summary_MANUAL %>% dplyr::select(contains(c("BEH", "REC"), ignore.case = TRUE))
# summary_data_ACT$IntRow <- rep(1:nrow(summary_data_ACT))
# summary_data_REC$IntRow <- rep(1:nrow(summary_data_REC))
# #Rename columns to make them match and bind+ melt columns
# summa_data_ACT  <- summary_data_ACT %>% rename_with(~str_remove(., c("_ACT|Act_"))); summa_data_ACT$Role <- "ACT"
# summa_data_REC  <- summary_data_REC %>% rename_with(~str_remove(., c("_REC|Rec_"))); summa_data_REC$Role <- "REC"
# summa_data_bind <- rbind(summa_data_ACT,summa_data_REC)
# 
# # plot(summa_data_bind$mean_speed_pxpersec, as.factor(summa_data_bind$Role),type="o")
# ggplot(summa_data_bind, aes(x = Role, y = mean_speed_pxpersec,colour=delta_speed_pxpersec > 0)) +
#   geom_line(aes(group = IntRow)) +
#   geom_point()


#temporary plot, understand steplength for DT_dist_THRESHOLD
# par(mfrow=c(1,2))
# hist(interaction_MANUAL$ACT.distance,breaks = 60,main="ACT stepwise distance",sub="blue line = 2") + abline(v=2,col='blue', lwd=2) 
# hist(interaction_MANUAL$REC.distance,breaks = 60, main="REC stepwise distance",sub="blue line = 2") + abline(v=2,col='blue', lwd=2) 
# plot dt frame!!!!
#
# par(mfrow=c(1,2))
# hist(summary_AUTO_REP_PER$moved_distance_px_ACT,breaks = 600,main="AUTO ACT stepwise distance",sub="blue line = 2") + abline(v=2,col='blue', lwd=2)
# hist(summary_AUTO_REP_PER$moved_distance_px_REC,breaks = 600, main="AUTO REC stepwise distance",sub="blue line = 2") + abline(v=2,col='blue', lwd=2)
# 
# hist(summary_AUTO_REP_PER$moved_distance_px_ACT,breaks = 1000,main="AUTO ACT stepwise distance - truncated",sub="blue line = 2",xlim = c(0,50)) + abline(v=2,col='blue', lwd=2)
# hist(summary_AUTO_REP_PER$moved_distance_px_REC,breaks = 1000, main="AUTO REC stepwise distance - truncated",sub="blue line = 2",xlim = c(0,50)) + abline(v=2,col='blue', lwd=2)
# #plot dt frame!!!!

###############################################################################
###### FINAL DATASET BUILDING/SAVING  #########################################
###############################################################################

# produce Hit/Miss info
AUTO_Hit              <- nrow(summary_AUTO[which(summary_AUTO$Hit==1),])
AUTO_Miss             <- nrow(summary_AUTO[which(summary_AUTO$Hit==0),])

# N of interactions per each MANUAL REP_PER  
MAN_int               <- aggregate(BEH ~ REPLICATE + PERIOD, FUN=length, summary_MANUAL)
MAN_int_Count         <- data.frame( REP_PER =paste0("MAN_int_",MAN_int$REPLICATE, "-", MAN_int$PERIOD), Freq = MAN_int$BEH)

#create a row per each Loop run 
Grooming_LDA_eachRun  <- data.frame(Loop_ID,CAPSULE_FILE,DT_dist_THRESHOLD, MAX_INTERACTION_GAP,DISAGREEMENT_THRESH, #looping Variables
                                   t(column_to_rownames(MAN_int_Count,"REP_PER")), TOT_MAN_int = nrow(summary_MANUAL), #TOT and REP_PER info on MANUAL interactions
                                   TOT_AUTO_int =nrow(summary_AUTO),AUTO_Hit, AUTO_Miss, #TOT and Hit/Miss info on AUTO interactions
                                   t(column_to_rownames(CSI_scores,"REP_PER")),  perc_CSI, #TOT and REP_PER CSI score
                                   stringsAsFactors = F,row.names = NULL)
#report the CSI score here
CSI

#stack 
Grooming_LDA_output   <- rbind(Grooming_LDA_output,     Grooming_LDA_eachRun)

#clear cache before running next loop. TO BE TESTED 
# rm(list=(c("e")))
# gc() # clear cache

#check that a full run through REP and PER was performed.
unique(interaction_MANUAL$PERIOD)
unique(interaction_MANUAL$REPLICATE)

#evaluate time taken for the loop to work
end.loop.time <- Sys.time()
time.taken.loop <- end.loop.time - start.loop.time

## Save outputs with dput using:
#dput(positions, file = "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_positions.txt")
#positions_dget <- dget("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_positions.txt") # load file created with dput 

cat(paste0("**LOOP COMPLETED**" ,
           "\n\nNotes: \n -Activate outer loop with all varying vars
                       \n -Save all the plots per loop in a dedicated folder named as the Loop_ID
                       \n -Save in Loop_ID the output for all the main components for further tests (Decision Trees,Logistic Regression,Random Forests,Support Vector Machines,Neural Networks)"
))
time.taken.loop

# }}}} ### OUTER PARAMETERS LOOP end


