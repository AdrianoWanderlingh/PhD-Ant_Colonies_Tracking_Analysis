##########################################################################################
############## THIS VERSION IS FORT 0.8.1 COMPATIBLE #####################################
##########################################################################################
#this should be the version of the script maintained for long term use.
#"https://formicidae-tracker.github.io/myrmidon/docs/latest/api/index.html"

rm(list=(c("e")))
gc()
rm(list=ls())

############## DEFINE MY FUNCTIONS
Movement.angle.diff     <- function(x)  {
  if( length(x) > 0 & !is_null(x)){
  c( atan(x[-nrow(x), "x"] / x[-nrow(x), "y"]) - atan(x[-1, "x"] / x[-1, "y"]), NA)}
else {vector()}}

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

#change system of coords
Orientation.diff        <- function(x)  { c( x[-nrow(x), "angle"] - x[-1, "angle"], NA)}

######load necessary libraries
library(adehabitatHR) ####for home range and trajectory analyses
library(FortMyrmidon) ####R bindings
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

#SOURCES ON/OFF
run_collisions          <- FALSE
run_AUTO_MAN_agreement  <- FALSE
run_Parameters_plots    <- FALSE


## PARAMETERS
USER            <- "Adriano"
Xmin <- 2000
Xmax <- 7900
Ymin <- 0000
Ymax <- 5500

N_DECIMALS                  <- 3 ## number of decimals to preserve when rounding to match between interactions & collisions
FUZZY_MATCH                 <- TRUE  ## fuzzy matching between data frames

max_gap                     <- fmHour(24*365)   ## important parameter to set! Set a maximumGap high enough that no cutting will happen. For example, set it at an entire year: fmHour(24*365)
desired_step_length_time    <- 0.125 ###in seconds, the desired step length for the analysis

#matchers specific for GROOMING
MAX_INTERACTION_GAP         <- 10
ANT_LENGHT_PX               <- 153 #useful for matcher::meanAntDisplacement mean and median value are similar
AntDistanceSmallerThan      <- 264 #for higher accuracy, recalculate it from: max(interaction_MANUAL$straightline_dist_px,na.rm = T)
AntDistanceGreaterThan      <- 63 #for higher accuracy, recalculate it from: min(interaction_MANUAL$straightline_dist_px,na.rm = T)
MAX_DISPLACEMENT            <- fmSecond(0.5)
# AntDisplacement <- #max trajectory step lenght per each ant during interaction, use higher value for matcher


####### navigate to folder containing myrmidon file
if (USER=="Adriano") {WORKDIR <- "/home/cf19810/Documents/Ants_behaviour_analysis"}
if (USER=="Tom")     {WORKDIR <- "/media/tom/MSATA/Dropbox/Ants_behaviour_analysis"}

DATADIR <- paste(WORKDIR,"Data",sep="/")
SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")

## perform analyses on annotation data
#source(paste(SCRIPTDIR,"BEH_Annotations_analysis.R",sep="/"))
# OR
#load only the training dataset
annotations <- read.csv(paste(DATADIR,"/annotations_TRAINING_DATASET.csv",sep = ""), sep = ",")
#transform zulu time in GMT
annotations$T_start_UNIX <- as.POSIXct(annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
annotations$T_stop_UNIX  <- as.POSIXct(annotations$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
#assign time in sec 
annotations$T_start_sec <- round(as.numeric(annotations$T_start_UNIX),N_DECIMALS)
annotations$T_stop_sec <- round(as.numeric(annotations$T_stop_UNIX),N_DECIMALS)

### start fresh
interaction_MANUAL    <- NULL
summary_MANUAL        <- NULL
Sensitivity           <- data.frame()

for (REPLICATE in c("R3SP","R9SP")) 
  {
  ###############################################################################
  ###### OPEN EXPERIMENT INFORMATION ############################################
  ###############################################################################
  
  ## locate the ant info file for REPLICATE
  MyrmidonCapsuleFile <- list.files(path=DATADIR, pattern=REPLICATE, full.names=T)
  MyrmidonCapsuleFile <- MyrmidonCapsuleFile[grepl("Capsule_Zones_defined.myrmidon",MyrmidonCapsuleFile)]
  e <- fmExperimentOpen(MyrmidonCapsuleFile)
  fmQueryGetDataInformations(e) # $details$tdd.path
  experiment_name <- unlist(strsplit(MyrmidonCapsuleFile,split="/"))[length(unlist(strsplit(MyrmidonCapsuleFile,split="/")))]
  experiment_name_end <-  str_sub(experiment_name,-30,-1)
  CapsuleDef <-  substr(experiment_name_end, 1, 7)
  ###tag statistics
  tag_stats <- fmQueryComputeTagStatistics(e)
  
  
  for (PERIOD in c("pre","post"))
    {
    print(paste("Replicate",REPLICATE, PERIOD))
    ## Prepare empty within-replicate/period data object
    interacts_MAN_REP_PER <- NULL  
    summary_MAN_REP_PER  <- NULL
    ## set experiment time window 
    time_start <- fmTimeCreate(min(annotations$T_start_UNIX[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD])) ###experiment start time
    #time_stop  <- fmTimeCreate(min(annotations$T_start_[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD]) + (34*60)  ) ###experiment stop time ####arbitrary time in the correct format + (N mins * N seconds)
    time_stop  <- fmTimeCreate(max(annotations$T_stop_UNIX[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD]) ) ###experiment stop time ####arbitrary time in the correct format + (N mins * N seconds)
    
    ###############################################################################
    ###### IDENTIFY FRAMES ########################################################
    ###############################################################################
    IdentifyFrames      <- fmQueryIdentifyFrames(e,start=time_start, end=time_stop)
    IF_frames           <- IdentifyFrames$frames
    # Assign a frame to each time since start and use it as baseline for all matching and computation
    IF_frames$frame_num <- seq.int(nrow(IF_frames))
    
    # assign a time in sec to match annotation time (UNIX hard to use for class mismatch)
    IF_frames$time_sec <- round(as.numeric(IF_frames$time),N_DECIMALS)
 
    #assign frame numbering to annotations 
    annotations$frame_start <- match.closest(x = annotations$T_start_sec, table = IF_frames$time_sec, tolerance = 0.05)
    annotations$frame_stop  <- match.closest(x = annotations$T_stop_sec, table = IF_frames$time_sec, tolerance = 0.05)
    
    # Creating a new zeroed-time since the start of the exp by  summing the cumulated differences between each pair of consecutive frames (to account for the time mismatch)
    IF_frames$cum_diff <- c(0, with(IF_frames, diff(time)))
    IF_frames$cum_diff <- cumsum(IF_frames$cum_diff)

    ###############################################################################
    ###### READING COLLISIONS #####################################################
    ###############################################################################
    collisions <- fmQueryCollideFrames(e, start=time_start, end=time_stop)   ###collisions are for each frame, the list of ants whose shapes intersect one another. Normally not used
    collisions$frames$frame_num <- seq.int(nrow(IF_frames))
    
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

    trajectories_summary <- positions$trajectories_summary
    
    
    ## Add ant_x names and times to the positions to convert from FRAME since the start of the experiment, to FRAMES
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
    
    # ## Add ant_x names and times to the positions to convert from time since the start of the experiment, to UNIX time
    # for (A in positions$trajectories_summary$antID)
    #   {
    #   AntID <- paste("ant_",A,sep="") 
    #   First_Obs_Time <- as.POSIXct(positions$trajectories_summary$start [which(positions$trajectories_summary$antID_str==AntID)], tz="GMT",origin = "1970-01-01 00:00:00") ## find the first time after the user defined time_start_ISO that this ant was seen
    #   print(paste("Adding first obs time", First_Obs_Time, "to the time-zeroed trajectory of ant", AntID))
    #   positions$trajectories[[AntID]] $UNIX_time <- as.POSIXct(positions$trajectories[[AntID]]$time, tz="GMT",origin = "1970-01-01 00:00:00")  + First_Obs_Time ##convert back to UNIX time  
    # }

    ###############################################################################
    ###### EXTRACT ANT TRAJECTORIES FROM MANUALLY-ANNOTATED DATA ##################
    ###############################################################################
    source(paste(SCRIPTDIR,"BEH_Traj_from_Man_annotations_fort081.R",sep="/"))
    
    
    #create new variable by pasting ant numbers "low,high" for summary_MAN_REP_PER
    summary_MAN_REP_PER$ant1 <- as.numeric(gsub("ant_","", summary_MAN_REP_PER$Act_Name))
    summary_MAN_REP_PER$ant2 <- as.numeric(gsub("ant_","", summary_MAN_REP_PER$Rec_Name))
    ## the ant 1 & ant 2 labels are not strictly ascending, so sort them so ant 1 alwasy < ant 2
    summary_MAN_REP_PER$pair <- apply(summary_MAN_REP_PER[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
    
    ## stack
    interaction_MANUAL <- rbind(interaction_MANUAL, interacts_MAN_REP_PER)
    summary_MANUAL     <- rbind(summary_MANUAL,       summary_MAN_REP_PER)
    
    
    #########################################################################################
    ###### READING AUTOMATIC INTERACTIONS ###################################################
    #########################################################################################
    
    #fmMatcherAntAngleSmallerThan(), fmMatcherAntAngleGreaterThan() 
    capsules  <- e$antShapeTypeNames
    names(capsules) <- as.character( 1:length(capsules))
    antennae_id <- as.numeric(names(capsules)[[which(capsules=="antennae")]])
    body_id <- as.numeric(names(capsules)[[which(capsules=="body")]])
    
    # body_id <- capsules[which(capsules$name=="body"),"typeID"]
    matcherCapType <- fmMatcherInteractionType(body_id,antennae_id)
    matcherCapTypeAntDists <- fmMatcherAnd(list(fmMatcherInteractionType(body_id,antennae_id),
                                 fmMatcherAntDistanceSmallerThan(AntDistanceSmallerThan),
                                 fmMatcherAntDistanceGreaterThan(AntDistanceGreaterThan),
                                 fmMatcherAntDisplacement(ANT_LENGHT_PX, MAX_DISPLACEMENT)))     #check every 5 seconds if ant has displaced more than ANT_LENGHT_PX
    #different gap size doesn't seem to have an impact but maybe the length does!
    
    #NOTES: THESE MAY ALL BE WRONG!!!!
    # - adding AntDistances removes false positives but has no effect on false negatives
    # - with MaxGap=5 there is a reduction in f.neg. and f. positives compared to MaxGap=10
    # - with MaxGap=3 there is a reduction in f.neg. and f. positives compared to MaxGap=5
    # - with MaxGap=1 there is a reduction in f.neg. and f. positives compared to MaxGap=3
    
    ## sequentially vary the interaction gap-filling to check what effect this has on the agreement between the MANUAL & AUTOMATIC interactions
    
    # for (Buffer in seq(0,30,5))
    #   {
    #   for (MAX_INTERACTION_GAP in c(seq(1,9,2), seq(10,60,10)))
    #     {
       
        interacts_AUTO_REP_PER <- fmQueryComputeAntInteractions(e,
                                                                start=time_start, 
                                                                end=time_stop,
                                                                maximumGap =fmSecond(MAX_INTERACTION_GAP), ## WHEN A PAIR DISENGAGE, HOW LONG IS THE INTERVAL? 
                                                                reportFullTrajectories = T,
                                                                matcher = matcherCapTypeAntDists)
        
        ## Assign frame number
        interacts_AUTO_REP_PER$interactions["int_start_frame"]   <- lapply(interacts_AUTO_REP_PER$interactions["start"] , function(x) IF_frames$frame_num[match(x, IF_frames$time)])
        interacts_AUTO_REP_PER$interactions["int_end_frame"]   <- lapply(interacts_AUTO_REP_PER$interactions["end"] , function(x) IF_frames$frame_num[match(x, IF_frames$time)])
        
        interacts_AUTO_REP_PER$interactions $pair <- paste(interacts_AUTO_REP_PER$interactions$ant1, interacts_AUTO_REP_PER$interactions$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...
        
        interacts_AUTO_REP_PER$interactions$Duration <- interacts_AUTO_REP_PER$interactions$int_end_frame - interacts_AUTO_REP_PER$interactions$int_start_frame
        hist(interacts_AUTO_REP_PER$interactions$Duration,breaks = 30)
        nrow(interacts_AUTO_REP_PER$interactions)
        
        # ##  convert POSIX format to raw secs since 1.1.1970; brute force for match.closest with collisions
        # interacts_AUTO_REP_PER$interactions$int_start_secs  <- as.numeric(interacts_AUTO_REP_PER$interactions$start)  ## yes, it looks like the milisecs are gone, but they are there; traj_ACT$UNIX_secs
        # interacts_AUTO_REP_PER$interactions$int_end_secs     <- as.numeric(interacts_AUTO_REP_PER$interactions$end) 
        # ## now round these decimal seconds to 3 d.p.
        # interacts_AUTO_REP_PER$interactions$int_start_secs   <- round(interacts_AUTO_REP_PER$interactions$int_start_secs,   N_DECIMALS) ## eliminates the 'noise' below the 3rd d.p., and leaves each frame existing just once, see: table(traj_ACT$UNIX_secs)
        # interacts_AUTO_REP_PER$interactions$int_end_secs     <- round(interacts_AUTO_REP_PER$interactions$int_end_secs, N_DECIMALS)
      
        # ## CHECK WHAT MAX_INTERACTION_GAP is doing - mean interval between successive interactions for a pair {A,B}
        # PairInterInteractInterv <- data.frame(Pair=unique(interacts_AUTO_REP_PER$interactions$pair))  ## 1 element per pair
        # PairInterInteractInterv$Mean_Interval <- NA
        # PairInterInteractInterv$Min_Interval <- NA
        # for (PR in PairInterInteractInterv$Pair)
        #   {
        #   ## extract successive interactions for PR
        #   interacts_AUTO_REP_PER_PAIR <- interacts_AUTO_REP_PER$interactions [ which(interacts_AUTO_REP_PER$interactions$pair==PR),] 
        #   
        #   if (nrow(interacts_AUTO_REP_PER_PAIR)>2) ## need at least two inteactions to have an inter-interaction interval
        #     {
        #     ## measure inter-interaction interval
        #     Inter_Interact_Interval <- interacts_AUTO_REP_PER_PAIR$int_end_frame   [2:nrow(interacts_AUTO_REP_PER_PAIR)] - 
        #                                interacts_AUTO_REP_PER_PAIR$int_start_frame [1:(nrow(interacts_AUTO_REP_PER_PAIR)-1)]
        #     ## get mean interval for this pair
        #     Mean_Interval <- mean(Inter_Interact_Interval, na.rm=T)
        #     Min_Interval <- min(Inter_Interact_Interval, na.rm=T)  ## SHOUL DNEVER BE LOWER THAN MAX_INTERACTION_GAP
        #     PairInterInteractInterv$Mean_Interval [which(PairInterInteractInterv==PR)] <-  Mean_Interval
        #     PairInterInteractInterv$Min_Interval  [which(PairInterInteractInterv==PR)] <-  Min_Interval
        #     }
        #   }
        # GrandMeanInterval <- mean(PairInterInteractInterv$Min_Interval, na.rm=T)
        # GrandMinInterval  <- min(PairInterInteractInterv$Min_Interval, na.rm=T)
        # 

        # summaryTrajectory_all       <- interacts_AUTO_REP_PER$summaryTrajectory  ## same object as positions$trajectory_summary described above - however here it will have many more rows given that the trajectories will be broken whenever there is a 10 second non-detection gap for an ant.
        # ## therefore in this case my trick with ant_ID_str won't work - it only works if there's exactly one trajectory as ants in the data
        # ## so I would NOT use the output of this function to analyse trajectories
        # trajectories_all           <- interacts_AUTO_REP_PER$trajectories       ## same object as positions$trajectories described above - but containing more objects for the same reason as stated above
        # interacts_AUTO_REP_PER            <- interacts_AUTO_REP_PER$interactions       ## data frame containing all interactions
        # 
        
        ########################################################################################################################
        ################# cross reference manually-laeblled start & end times with auto interactions ###########################
        ########################################################################################################################
#         
#          ## TRUE POSITIVE RATE
#          interacts_AUTO_REP_PER_OVERLAP_ALL <- NULL
#          Buffer <- 0
#          for (I in 1:nrow(summary_MAN_REP_PER))  ## loop across each (time-aggregated) manually-labelled behaviour
#             {
#             ## Find start and end times of each *aggregated* **MANUALLY-LABELLED** interaction
#             Man_Start <- summary_MAN_REP_PER$int_start_frame [I] - Buffer
#             Man_Stop  <- summary_MAN_REP_PER$int_end_frame [I]  + Buffer
#             
#             ## extract manual participants for the Ith interaction in the auto interactions
#             Man_Pair <- summary_MAN_REP_PER$pair [I]
#           
#             
#             ## look for the Ith MANUAL interaction in the AUTOMATIC list
#             interacts_AUTO_REP_PER_OVERLAP <- unique(interacts_AUTO_REP_PER$interactions [ interacts_AUTO_REP_PER$interactions$int_start_frame >= Man_Start &
#                                                                                            interacts_AUTO_REP_PER$interactions$int_end_frame   <= Man_Stop & 
#                                                                                            interacts_AUTO_REP_PER$interactions $pair      == Man_Pair , ])
#             
#        
#             
#             ## stack the auto interactions when there's a match with the manual summary interaction
#             if (nrow(interacts_AUTO_REP_PER_OVERLAP)>0)
#               {
#               interacts_AUTO_REP_PER_OVERLAP_ALL <- rbind(interacts_AUTO_REP_PER_OVERLAP_ALL, 
#                                                           interacts_AUTO_REP_PER_OVERLAP)
#               
#               interacts_AUTO_REP_PER_OVERLAP_ALL <- unique(interacts_AUTO_REP_PER_OVERLAP_ALL)
#               }
#             }#I
#        
#         ## what is the overlap in seconds?
#         interacts_AUTO_REP_PER_OVERLAP_ALL$Duration <- interacts_AUTO_REP_PER_OVERLAP_ALL$int_end_frame - interacts_AUTO_REP_PER_OVERLAP_ALL$int_start_frame
#         
#         ## the observed overlap time
#         Overlap <- sum(interacts_AUTO_REP_PER_OVERLAP_ALL$Duration)
#         
#         ## what is the expected interaction time ?
#         Expectation <- sum(summary_MAN_REP_PER$interaction_length_secs)
#         ## stack
#         Sensitivity <- rbind(Sensitivity, data.frame(REPLICATE, PERIOD, Buffer, MAX_INTERACTION_GAP, GrandMeanInterval, GrandMinInterval, Overlap, Expectation, Hit_Rate=100 * (Overlap/Expectation)) )
# 
#         
# 
# #####################################################################################
#         
#         #oder dataframe to make it work with Fuzzy matcher which needs sorted dataframes
#         sum_MAN_REP_PER_ord <- summary_MAN_REP_PER[with(summary_MAN_REP_PER, order(int_start_frame)),]
#         sum_MAN_REP_PER_ord1 <- summary_MAN_REP_PER[with(summary_MAN_REP_PER, order(int_end_frame)),]
# 
#         TOLERANCE <- 20 #number of seconds of tolerance window, is like the previous buffer
#         FUZZY_INT_MATCH <- TRUE
#         A_in_M_OVERLAP_ALL <- NULL
#         if (FUZZY_INT_MATCH==TRUE)
#         {
#           for (I in 1:nrow(summary_MAN_REP_PER))  ## loop across each (time-aggregated) manually-labelled behaviour
#           {
#             ## Find start and end times of each *aggregated* **MANUALLY-LABELLED** interaction
#             Man_Start <- summary_MAN_REP_PER$int_start_frame [I]
#             Man_Stop  <- summary_MAN_REP_PER$int_end_frame [I]
# 
#             ## extract manual participants for the Ith interaction in the auto interactions
#             Man_Pair <- summary_MAN_REP_PER$pair [I]
# 
#           ### TRUE POSITIVES
#           # match.start and end
#           #start
#           A_in_M_rows_start <- unique(interacts_AUTO_REP_PER$interactions[ match.closest(x =     interacts_AUTO_REP_PER$interactions$int_start_frame,
#                                                                                          table = sum_MAN_REP_PER_ord$int_start_frame , tolerance = TOLERANCE)
#                                                                            & interacts_AUTO_REP_PER$interactions $pair == Man_Pair  , ])
#           #end
#           A_in_M_rows_end   <- unique(interacts_AUTO_REP_PER$interactions[ match.closest(x =     interacts_AUTO_REP_PER$interactions$int_end_frame,
#                                                                                          table = sum_MAN_REP_PER_ord1$int_end_frame , tolerance = TOLERANCE)
#                                                                            & interacts_AUTO_REP_PER$interactions $pair == Man_Pair  , ])
#           #issue raise with times, drop them
#           A_in_M_rows_end <-subset(A_in_M_rows_end, select = -c(start,end))
#           A_in_M_rows_start <-subset(A_in_M_rows_start, select = -c(start,end))
# 
#           A_in_M_OVERLAP <- generics::intersect(A_in_M_rows_start, A_in_M_rows_end)
#           #remove empty rows (all NA's)
#           A_in_M_OVERLAP <- A_in_M_OVERLAP[rowSums(is.na(A_in_M_OVERLAP)) != ncol(A_in_M_OVERLAP),]
#           
#           #### A_in_M ARE TRUE POSITIVES  -> ROWS IN BOTH AUTO E MAN
#         
#         ## stack the auto interactions when there's a match with the manual summary interaction
#         if (nrow(A_in_M_OVERLAP)>0)
#         {
#           A_in_M_OVERLAP_ALL <- rbind(A_in_M_OVERLAP_ALL, 
#                                                       A_in_M_OVERLAP)
#           A_in_M_OVERLAP_ALL <- unique(A_in_M_OVERLAP_ALL)
#         }
#       }#I
#     }#FUZZY_INT_MATCH    
# 
#  #FALSE NEGATIVE RATE
# 
#  ####ERRORRRRRR AUTO_IN_MANUAL SHOULD BE LESS OR EQUAL TO SUMMARY_MANUAL!!!!! DUPS AND NAS ARE REMOVED, WHAT'S HAPPENING???         
#  #cure momentarily by trimming auto$duration < manual$duration
#  nrow(A_in_M_OVERLAP_ALL) 
#  nrow(summary_MAN_REP_PER)
#  
#  summary_MAN_REP_PER$Duration <- summary_MAN_REP_PER$int_end_frame - summary_MAN_REP_PER$int_start_frame
#  A_in_M_OVERLAP_ALL$Duration <- A_in_M_OVERLAP_ALL$int_end_frame - A_in_M_OVERLAP_ALL$int_start_frame
#  
#  min(summary_MAN_REP_PER$Duration)
#  min(A_in_M_OVERLAP_ALL$Duration)
#  
#  min( summary_MAN_REP_PER$Duration[summary_MAN_REP_PER$Duration!=min(summary_MAN_REP_PER$Duration)] )
#  
#  MIN_TRIM <- min( summary_MAN_REP_PER$Duration[summary_MAN_REP_PER$Duration!=min(summary_MAN_REP_PER$Duration)] )
#  A_in_M_OVERLAP_ALL_trim <- A_in_M_OVERLAP_ALL[which(A_in_M_OVERLAP_ALL$Duration > MIN_TRIM),]
# 
#  
#  hist(summary_MAN_REP_PER$Duration,breaks = 30)
#  hist(A_in_M_OVERLAP_ALL$Duration,breaks = 30)
#  hist(A_in_M_OVERLAP_ALL_trim$Duration,breaks = 30)
#         
######################################################################################        

        pdf(file=paste(DATADIR,"Interactions","AUTO_MAN_REP_PER","Gap",MAX_INTERACTION_GAP,"matcher_CapTypeAntDistsDispl.pdf", sep = "_"), width=10, height=60)
        #plot the interactions by pair as timeline
        #generate 1 plot per every iteration 
        summary_MAN_REP_PER_sub <- summary_MAN_REP_PER[,c("REPLICATE","PERIOD","pair","int_start_frame","int_end_frame")] ; summary_MAN_REP_PER_sub$flag <- "manual"
        interacts_AUTO_REP_PER_sub <- interacts_AUTO_REP_PER$interactions[,c("pair","int_start_frame","int_end_frame")] ; interacts_AUTO_REP_PER_sub$flag <- "auto"
        AUTO_MAN_REP_PER <- dplyr::bind_rows(summary_MAN_REP_PER_sub,interacts_AUTO_REP_PER_sub) #use bind_rows to keep rep info https://stackoverflow.com/questions/42887217/difference-between-rbind-and-bind-rows-in-r
 
        
        ggplot(AUTO_MAN_REP_PER) +
          geom_linerange(aes(y = pair, xmin = int_start_frame, xmax = int_end_frame,colour = flag),size=3,alpha = 0.5) +
          theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
          labs(title = paste("Grooming by pair" ,unique(REPLICATE),unique(PERIOD),"- MAN vs AUTO"),
               subtitle = paste("Nrows auto detected:" ,NROW(interacts_AUTO_REP_PER$interactions),"Nrows manual annotated:" ,NROW(summary_MAN_REP_PER),"\nMaxIntGap",MAX_INTERACTION_GAP, "s","Capsule file =",CapsuleDef )) +
          scale_color_manual(values = c("manual" = "red",
                                         "auto"="black"))
       #dev.off()
        
      #}#MAX_INTERACTION_GAP
    #}#Buffer
    
    
    ## Select ONLY those AUTO interactions that are INSIDE the manual interactions
    # plot( Sensitivity[Sensitivity$Buffer==0 , c("MAX_INTERACTION_GAP","Buffer","GrandMinInterval","Overlap","Hit_Rate")])    
    
    ###############################################################################
    ###### AUTO_MAN AGREEMENT MATRIX ##############################################
    ###############################################################################
    if (run_AUTO_MAN_agreement){source(paste(SCRIPTDIR,"BEH_Auto_Man_agreement_matrix_fort081.R",sep="/"))}
  
    
    ###############################################################################
    ###### COLLISIONS #############################################################
    ###############################################################################
    if (run_collisions){source(paste(SCRIPTDIR,"BEH_collisions_fort081.R",sep="/"))}
      
  
  

       
    ## generate summary data
    # for (variable in names(summary_MAN_REP_PER)[!names(summary_MAN_REP_PER)%in%c("BEH","Act_Name","Rec_Name","PERIOD")])
    #   {
    #   summary_MAN_REP_PER[,"variable"] <- summary_MAN_REP_PER[,variable] #boxplot(variable~BEH,ylab=variable,data=summary_MAN_REP_PER)
    #   }##summary_MAN_REP_PER
    
    
  }##PERIOD
  #clear cache before opening following exp. TO BE TESTED 
  # rm(list=(c("e")))
  # gc() # clear cache
}##REPLICATE
dev.off()



###################### TO DOS ##############################
############################################################
# - position of heads is different between various behaviours, if you can measure distance between capsules you can use that as a variable
# - ask Mathias if he knows how to extract the geometry of capsules
# - angular momentum
# - acceleration of the angle
# - frame by frame rate of change of angle
# - list of measures for enrico
# - DELTA ANGLES formulas FIX
# Add collisions capsules
#other possible thing to do: Mixture model to separate between behaviours (but it is nott what I need to do: https://towardsdatascience.com/mixture-modelling-from-scratch-in-r-5ab7bfc83eef
# https://stats.stackexchange.com/questions/111145/how-to-fit-mixture-model-for-clustering)

#TO DO:
#Conversion of data in mm? Makes sense to do that on extracted trajectories files, using the tag size and box as reference...
# reuse this and tags size in pixel/mm as point of reference - traj[,which(grepl("x",names(traj))|grepl("y",names(traj)))] <-   traj[,which(grepl("x",names(traj))|grepl("y",names(traj)))] * petridish_diameter / range 


#check
unique(interaction_MANUAL$PERIOD)
unique(interaction_MANUAL$REPLICATE)
min(interaction_MANUAL$UNIX_time)
max(interaction_MANUAL$UNIX_time)

###############################################################################
###### PARAMETERS PLOTS #######################################################
###############################################################################
if (run_Parameters_plots){source(paste(SCRIPTDIR,"BEH_Parameters_plots_fort081.R",sep="/"))}





  rm(list=(c("e")))
  gc()
  
  
  
  
  
  
  
################# SCRAPS ##########################
  
  # #later on, deparse the capsule information to add it as an individual-linked row value in the interaction data.
  # #the structure for a 3-1, 1-1, could look like this:
  # capsules_example <- read.table(textConnection('
  # time  ACT_caps1  ACT_caps2  ACT_caps2  REC_caps1  REC_caps2 REC_caps3
  # 21211 1 0 1 2 0 0
  # 21212 1 0 1 2 0 0
  # '), header=TRUE)
  # #not very convinced that it would work, maybe is better to avoid full deparsing into multiple columns
  # 
  # some workaround has to be done to ensure that the ant (ACT-REC) is assigned the right capsule from the pair
  #unlist(strsplit(interacts_AUTO_REP_PER$types[4],","))[2]
  
  
  
  #----------------------------------------------------------------
  
  #FIND HOW CAPSULES ARE SPATIALLY ARRANGED COMPARED TO TAG POSITION TO MAKE POSSIBLE TO CALCULATE distance between interacting ants's capsules AND use that as a variable
  
  #THE FOLLOWING BIT OF CODE COMES FROM THE AUTOMATIC ANGLE DETERMINATION SCRIPT. IT HELPS TO UNDERSTAND HOW TO ACCESS
  # CAPSULE PARTS. USE BITS OF IT TO GET CAPSULE GEOMETRIES, THEIR POSITION IN SPACE AND TO CALCULATE DISTANCE AMONG 
  # VARIOUS CAPSULES DURING INTERACTIONS (IE. HEAD-ABDOMEN DISTANCE DURING GROOMING)
  # 
  # oriented_metadata <- NULL
  # capsule_list <- list()
  # #for (myrmidon_file in data_list){
  # experiment_name <- unlist(strsplit(MyrmidonCapsuleFile,split="/"))[length(unlist(strsplit(MyrmidonCapsuleFile,split="/")))]
  # oriented_data <- fmExperimentOpen(MyrmidonCapsuleFile) #this step is already performed at the beginning
  # oriented_ants <- oriented_data$ants
  # capsule_names <- oriented_data$antShapeTypeNames
  # for (ant in oriented_ants){
  #   ###extract ant length and capsules
  #   #ant_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_px)
  #   capsules      <- ant$capsules
  #   for (caps in 1:length(capsules)){
  #     capsule_name  <- capsule_names[[capsules[[caps]]$type]]
  #     capsule_coord <- capsules[[caps]]$capsule
  #     capsule_info <- data.frame(experiment = experiment_name,
  #                                antID      = ant$ID,
  #                                c1x = capsule_coord$c1[1],
  #                                c1y = capsule_coord$c1[2],
  #                                c2x = capsule_coord$c2[1],
  #                                c2y = capsule_coord$c2[2],
  #                                r1  = capsule_coord$r1[1],
  #                                r2   = capsule_coord$r2[1]
  #     )
  #     
  #     if (!capsule_name %in%names(capsule_list)){ ###if this is the first time we encounter this capsule, add it to capsule list...
  #       capsule_list <- c(capsule_list,list(capsule_info)) 
  #       if(length(names(capsule_list))==0){
  #         names(capsule_list) <- capsule_name
  #       }else{
  #         names(capsule_list)[length(capsule_list)] <- capsule_name
  #       }
  #     }else{###otherwise, add a line to the existing dataframe within capsule_list
  #       capsule_list[[capsule_name]] <- rbind(capsule_list[[capsule_name]] , capsule_info)
  #     }
  #   }
  # }
  # 
  # 
  #----------------------------------------------------------------
  
  
# 
# x <- 1
# y <- 1
# 
# atan(x/y) * 180/pi
# 
# #rad to deg
# rad * 180/pi

  

#dput(positions, file = "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_positions.txt")
#positions_dget <- dget("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_positions.txt") # load file created with dput 
