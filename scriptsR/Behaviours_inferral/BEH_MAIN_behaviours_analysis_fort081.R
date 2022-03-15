##########################################################################################
############## THIS VERSION IS FORT 0.8.1 COMPATIBLE #####################################
##########################################################################################
#this should be the version of the script maintained for long term use.
#"https://formicidae-tracker.github.io/myrmidon/docs/latest/api/index.html"

rm(list=(c("e")))
gc()
rm(list=ls())

############## DEFINE MY FUNCTIONS
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
library(fields)
library(sp) #calculate convex hull area
library(bestNormalize)
library(BAMBI) #angles wrapping

#SOURCES ON/OFF
run_collisions                      <- TRUE
run_AUTO_MAN_agreement              <- TRUE #computes agreement matrix
create_interaction_AUTO_REP_PER     <- FALSE 
#plots
run_Parameters_plots                <- FALSE
plot_all_BEH                        <- FALSE #inside run_Parameters_plots, if FALSE plots only for Grooming (to run it an update in BEH statement is needed)



## PARAMETERS
USER <- "Adriano"
Xmin <- 2000
Xmax <- 7900
Ymin <- 0000
Ymax <- 5500

N_DECIMALS                  <- 3 ## number of decimals to preserve when rounding to match between interactions & collisions
FUZZY_MATCH                 <- TRUE  ## fuzzy matching between data frames in collisions detection

max_gap                     <- fmHour(24*365)   ## important parameter to set! Set a maximumGap high enough that no cutting will happen. For example, set it at an entire year: fmHour(24*365)
desired_step_length_time    <- 0.125 ###in seconds, the desired step length for the analysis

#trajectories jumps/gaps thresholds to avoid getting skewed means (see their use in params extraction scripts)
DT_frame_THRESHOLD <- 16
DT_dist_THRESHOLD  <- 0.5 #tag length is 62 px approx (measured on full size pics in R9SP)

#matchers specific for GROOMING
MAX_INTERACTION_GAP         <- 10 ## in SECONDS 
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
interaction_MANUAL_COLL <- NULL
summary_MANUAL        <- NULL
interaction_AUTO      <- NULL
summary_AUTO          <- NULL
Sensitivity           <- data.frame()

#set plots parameters (for plotting coords)
pdf(file=paste(DATADIR,"Interactions_plots_8feb2022.pdf", sep = ""), width=6, height=4.5)
par(mfrow=c(2,3), mai=c(0.3,0.3,0.4,0.1), mgp=c(1.3,0.3,0), family="serif", tcl=-0.2)


for (REPLICATE in c("R3SP","R9SP")) 
  {
  ###############################################################################
  ###### OPEN EXPERIMENT INFORMATION ############################################
  ###############################################################################
  
  ## locate the ant info file for REPLICATE
  MyrmidonCapsuleFile <- list.files(path=DATADIR, pattern=REPLICATE, full.names=T)
  MyrmidonCapsuleFile <- MyrmidonCapsuleFile[grepl("Capsule_Zones_defined_BODY-HEAD.myrmidon",MyrmidonCapsuleFile)]

  e <- fmExperimentOpen(MyrmidonCapsuleFile)
  fmQueryGetDataInformations(e) # $details$tdd.path
  experiment_name <- unlist(strsplit(MyrmidonCapsuleFile,split="/"))[length(unlist(strsplit(MyrmidonCapsuleFile,split="/")))]
  experiment_name_end <-  str_sub(experiment_name,-30,-1)
  CapsuleDef <-  substr(experiment_name_end, 1, 7)
  ###tag statistics
  tag_stats <- fmQueryComputeTagStatistics(e)
  
  
  ################################
  #CHANGE BASE HEAD CAPSULE FROM LONG TO LARGE (see notebook notes 23Feb)
  #ASSIGN VARIOUS CAPSULES SHAPES
  
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
    AntID <- NULL; First_Obs_Frame <- NULL
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
    
    #fmMatcherAntAngleSmallerThan(), fmMatcherAntAngleGreaterThan() 
    capsules  <- e$antShapeTypeNames
    names(capsules) <- as.character( 1:length(capsules))
    head_id <- as.numeric(names(capsules)[[which(capsules=="head")]])
    body_id <- as.numeric(names(capsules)[[which(capsules=="body")]])
    
    # body_id <- capsules[which(capsules$name=="body"),"typeID"]
    matcherCapType <- fmMatcherInteractionType(body_id,head_id)
    matcherCapTypeAntDists <- fmMatcherAnd(list(fmMatcherInteractionType(body_id,head_id),
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
      

        # #pdf(file=paste(DATADIR,"Interactions","AUTO_MAN_REP_PER","Gap",MAX_INTERACTION_GAP,"matcher_CapTypeAntDistsDispl.pdf", sep = "_"), width=10, height=60)
        # #plot the interactions by pair as timeline
        # #generate 1 plot per every iteration 
        # summary_MAN_REP_PER_sub <- summary_MAN_REP_PER[,c("REPLICATE","PERIOD","pair","int_start_frame","int_end_frame")] ; summary_MAN_REP_PER_sub$flag <- "manual"
        # interacts_AUTO_REP_PER_sub <- interacts_AUTO_REP_PER$interactions[,c("pair","int_start_frame","int_end_frame")] ; interacts_AUTO_REP_PER_sub$flag <- "auto"
        # AUTO_MAN_REP_PER <- dplyr::bind_rows(summary_MAN_REP_PER_sub,interacts_AUTO_REP_PER_sub) #use bind_rows to keep rep info https://stackoverflow.com/questions/42887217/difference-between-rbind-and-bind-rows-in-r
        # 
        # 
        # ggplot(AUTO_MAN_REP_PER) +
        #   geom_linerange(aes(y = pair, xmin = int_start_frame, xmax = int_end_frame,colour = flag),size=3,alpha = 0.5) +
        #   theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
        #   labs(title = paste("Grooming by pair" ,unique(REPLICATE),unique(PERIOD),"- MAN vs AUTO"),
        #        subtitle = paste("Nrows auto detected:" ,NROW(interacts_AUTO_REP_PER$interactions),"Nrows manual annotated:" ,NROW(summary_MAN_REP_PER),"\nMaxIntGap",MAX_INTERACTION_GAP, "s","Capsule file =",CapsuleDef )) +
        #   scale_color_manual(values = c("manual" = "red",
        #                                  "auto"="black"))
        
      #}#MAX_INTERACTION_GAP
    #}#Buffer
    
    
    ## Select ONLY those AUTO interactions that are INSIDE the manual interactions
    # plot( Sensitivity[Sensitivity$Buffer==0 , c("MAX_INTERACTION_GAP","Buffer","GrandMinInterval","Overlap","Hit_Rate")])    
    
    ###############################################################################
    ###### AUTO_MAN AGREEMENT MATRIX ##############################################
    ###############################################################################
    if (run_AUTO_MAN_agreement){source(paste(SCRIPTDIR,"BEH_Auto_Man_agreement_matrix_fort081.R",sep="/"))}
    
    #stack ComputeAntInteractions with extra info calculated in agreement matrix
    
    #stack output of BEH_Parameters_Auto_fort081.R + Hit information if run_AUTO_MAN_agreement is TRUE
    interaction_AUTO <- rbind(interaction_AUTO,interaction_AUTO_REP_PER)
    summary_AUTO <- rbind(summary_AUTO, summary_AUTO_REP_PER)
    
    
    ###############################################################################
    ###### COLLISIONS #############################################################
    ###############################################################################
    if (run_collisions){source(paste(SCRIPTDIR,"BEH_collisions_fort081.R",sep="/"))}
      
    #stack output of collisions
    # THIS OUTPUT SHOULD BE PROBABLY SAVED INTO interaction_MANUAL DIRECTLY (AS IT IS JUST interaction_MANUAL + CAPSULE INFO)
    # BUT WHEN DOING SO CHECK THAT THIS DOES NOT HURT THE FOLLOWING interaction_AUTO FILE
    interaction_MANUAL_COLL <- rbind(interaction_MANUAL_COLL, interacts_MAN_REP_PER)
    
    
    ## PCA to check for natural differences in the behaviour of actor versus receiver during manually-defined grooming interactions
    interaction_MANUAL_observables <- interaction_MANUAL             %>% dplyr::select(contains(c("ACT","REC"), ignore.case = FALSE))
    interaction_MANUAL_observables <- interaction_MANUAL_observables %>% dplyr::select(!contains(c(".x",".y","ACT.angle","REC.angle"), ignore.case = FALSE))
    
    
    #transform to long format
    interaction_MANUAL_ACT <- interaction_MANUAL_observables[,grep("ACT",colnames(interaction_MANUAL_observables))]; colnames(interaction_MANUAL_ACT) <- gsub("ACT.","",colnames(interaction_MANUAL_ACT));  colnames(interaction_MANUAL_ACT) <- gsub("_ACT","",colnames(interaction_MANUAL_ACT))
    interaction_MANUAL_REC <- interaction_MANUAL_observables[,grep("REC",colnames(interaction_MANUAL_observables))]; colnames(interaction_MANUAL_REC) <- gsub("REC.","",colnames(interaction_MANUAL_REC));  colnames(interaction_MANUAL_REC) <- gsub("_REC","",colnames(interaction_MANUAL_REC))
    ## add actor/receiver labels to each
    interaction_MANUAL_ACT$ActRec_label <- "A"
    interaction_MANUAL_REC$ActRec_label <- "R"
    
    ## stack actor & receiver
    interaction_MANUAL_ACTREC <- rbind(interaction_MANUAL_REC, interaction_MANUAL_ACT)
    ##
    interaction_MANUAL_ACTREC_noNA <- na.omit(interaction_MANUAL_ACTREC)
    
    ## signs -> absolutes
    interaction_MANUAL_ACTREC_noNA[, sapply(interaction_MANUAL_ACTREC_noNA[1,], is.numeric)] <- abs(interaction_MANUAL_ACTREC_noNA[, sapply(interaction_MANUAL_ACTREC_noNA[1,], is.numeric)])
    
    ## scale the inputs
    par(mfrow=c(3,4))
    ObsNames <- colnames(interaction_MANUAL_ACTREC_noNA) [-match("ActRec_label",colnames(interaction_MANUAL_ACTREC_noNA))]
    for (OBS in 1:length(ObsNames))
      {
      hist(interaction_MANUAL_ACTREC_noNA[,OBS] , col=1, main=paste(ObsNames[OBS],"pre-transform"))
      interaction_MANUAL_ACTREC_noNA     [,OBS] <- scale((interaction_MANUAL_ACTREC_noNA[,OBS])^0.1)
      hist(interaction_MANUAL_ACTREC_noNA[,OBS] , col=2, main=paste(ObsNames[OBS],"post-transform"))
      }
    
    ## PCA is inappropriate as we know who is who ...
    PCA <- prcomp (x = interaction_MANUAL_ACTREC_noNA[,-match("ActRec_label",colnames(interaction_MANUAL_ACTREC_noNA))], scale=TRUE, center=TRUE)
    ##  add point colour labels (same dimensions)
    Eigenvalues <- as.data.frame(PCA$x)
    Eigenvalues$Colour <- as.numeric(as.factor(interaction_MANUAL_ACTREC_noNA$ActRec_label))
    ## THERE IS A DIFFERENCE!!
    plot(PCA$x[,1:2], pch=1, col= Eigenvalues$Colour, bg= Eigenvalues$Colour)



    
    ## LDA
    LDA <- lda(ActRec_label ~ ., interaction_MANUAL_ACTREC_noNA)
    #get / compute LDA scores from LDA coefficients / loadings
    plda <- predict(object = LDA,
                    newdata = interaction_MANUAL_ACTREC_noNA)

    par(mai=c(0.4,0.4,0.1,0.1))
    ldahist(data = plda$x[,1], g=interaction_MANUAL_ACTREC_noNA$ActRec_label)

    ## TO DO: APPLY THE LDA TO THE TEST DATA SETS...!! (see 'predict' example in ?lda help file)    
    
    
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

##############################################################################
######### CAPSULES CALCULATIONS ##############################################
##############################################################################

############# CHANGE ALL TO interaction_MANUAL

#ACCESS THE CAPSULE INFO
interaction_MANUAL_COLL$REP_PER_R_B <- paste(interaction_MANUAL_COLL$REPLICATE,interaction_MANUAL_COLL$PERIOD,interaction_MANUAL_COLL$ROW,interaction_MANUAL_COLL$BEH,sep="_")
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
int_types_counts_MAN.DT$types_inv <- int_types_counts_MAN.DT$types
int_types_counts_MAN.DT$types_inv <- as.character( int_types_counts_MAN.DT$types_inv)
int_types_counts_MAN.DT$types_inv[ int_types_counts_MAN.DT$types_inv == "3-2"] <- "2-3"

# plot per int
ggplot( int_types_counts_MAN.DT, aes(fill=types_inv, y=count, x=REP_PER_R_B)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("perc. capsule type involved in each interaction",paste("body_id is type",body_id,", head_id is type",head_id,sep=" "))

#Getting % of capsules from collisions to understand which capsules are interacting more
int_types_count_TOT <- int_types_counts_MAN.DT[, list(fsum=sum(count)), by=types]
int_types_count_TOT$perc <- round(int_types_count_TOT$fsum/sum(int_types_count_TOT$fsum),2)
int_types_count_TOT

#Frame by Frame Number of Manual annotations that do not include 2-3 or 3-2 capsules
types_FREQ <- as.data.frame(table(interaction_MANUAL_COLL$types))
types_FREQ$body_head <- NULL
types_FREQ$body_head <- ifelse(grepl('2-3|3-2', types_FREQ$Var1), 'yes', "no")
types_FREQ$body_head
aggregate(types_FREQ$Freq, by=list(body_head=types_FREQ$body_head), FUN=sum)

#Interactionn by interaction Number of Manual annotations that do not include 2-3 or 3-2 capsules
interaction_MANUAL_COLL$body_head <- NULL
interaction_MANUAL_COLL$body_head <- ifelse(grepl('2-3|3-2', interaction_MANUAL_COLL$types), 'yes', "no")
interaction_MANUAL_COLL$body_head <- as.factor(interaction_MANUAL_COLL$body_head)
interaction_MANUAL_COLL$Freq <- 1
body_head_AGGREG <- aggregate(interaction_MANUAL_COLL$Freq, by=list(interaction_MANUAL_COLL$body_head,interaction_MANUAL_COLL$REP_PER_R_B), FUN=sum);names(body_head_AGGREG) <- c("body_head","REP_PER_R_B","Freq")


ggplot( body_head_AGGREG, aes(fill=body_head, y=Freq, x=REP_PER_R_B)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Perc. of Manual annotations that include 2-3 or 3-2 capsules")



#get from long to wide format (this can be appended to SUMMARY file)
int_types_per_INT <- reshape(int_types_counts_MAN.DT, idvar = "REP_PER_R_B", timevar = "types", direction = "wide")


##############################################################################
#### PLOTS FOR TOTAL VARS (are also in PCA but can be moved here for good) ###
##############################################################################

#################################MANUAL VARS###################################
if (run_Parameters_plots){source(paste(SCRIPTDIR,"BEH_Parameters_plots_fort081.R",sep="/"))}




#summary_MANUAL_vars <- summary_MANUAL[, -match(c("REPLICATE", "PERIOD","BEH","ROW","Act_Name","Rec_Name","int_start_frame","int_end_frame","prop_time_undetected_REC","prop_time_undetected_ACT","ant1","ant2","pair"), names(summary_MANUAL))] 
#
#transform to long format BY DIVIDING BY ACTOR AND RECEIVER!! (DONE ALREADY ELSEWHERE)
# summary_MANUAL_vars_long <- melt(summary_MANUAL_vars,id.vars=c("Hit")) #explanation on the warning message https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message
# 
# 
# ###plot divided by variable and Hit for Grooming
# par(mfrow=c(2,3), family="serif" , mar = c(0.1, 0.1, 2.2, 0))
# 
# summ_vars_plot <- ggplot(summary_MANUAL_vars_long, aes(value, fill = Hit)) +
#   facet_wrap(variable ~ .,scales="free") +
#   theme_bw() +
#   theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) #,legend.position="bottom",legend.justification='right'
# summ_vars_plot + geom_density(alpha = 0.2) +
#   labs(title = "Density plot for movement variables by Hit rate")#,
# #subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))


###############################################################################
###### AUTO-MAN DISAGREEMENT PLOT #############################################
###############################################################################
if (run_AUTO_MAN_agreement)
  {
  ## visualise the cutoffs
  ## Make y-axis logarithmic in histogram and fix issues for 0 occurences with log
  par(mfrow=c(3,1))
  #all
  hist.data = hist(summary_AUTO$disagreement, breaks=seq(-1,0,0.05), plot=F)
  hist.data$counts <- replace(log10(hist.data$counts), log10(hist.data$counts)==-Inf, 0)
  highestCount <- max(hist.data$counts)
  plot(hist.data, main="All disagreements (log)",
       sub = paste("N.AUTO TOT",length(summary_AUTO$Hit)),
       ylab='log10(Frequency)', xlab=" AUTO-MAN disagreement rate ", ylim=c(0,highestCount)); abline(v=-THRESH, lty=2)
  
  #True Positives
  hist.data1 = hist(summary_AUTO$disagreement[which(summary_AUTO$Hit==1)], breaks=seq(-1,0,0.05), plot=F)
  hist.data1$counts <- replace(log10(hist.data1$counts), log10(hist.data1$counts)==-Inf, 0)
  plot(hist.data1, main="True Positives",
       sub = paste("N.AUTO",length(summary_AUTO$disagreement[which(summary_AUTO$Hit==1)]),"/ N.MANUAL",length(summary_MAN_REP_PER)),
       ylab='log10(Frequency)',xlab=" AUTO-MAN agreement rate ", col=2, ylim=c(0,highestCount), xlim=c(-1,0)); abline(v=-THRESH, lty=2) #breaks=seq(-1,1,0.05)
  
  #False Positives
  hist.data1 = hist(summary_AUTO$disagreement[which(summary_AUTO$Hit==0)], breaks=seq(-1,0,0.05), plot=F)
  hist.data1$counts <- replace(log10(hist.data1$counts), log10(hist.data1$counts)==-Inf, 0)
  plot(hist.data1, main="False Positives",
       sub = paste("N.AUTO",length(summary_AUTO$disagreement[which(summary_AUTO$Hit==0)])),
       ylab='log10(Frequency)', xlab=" AUTO-MAN disagreement rate ",col=4, ylim=c(0,highestCount),xlim=c(-1,0)); abline(v=-THRESH, lty=2) #breaks=seq(-1,1,0.05)
  
  mtext(paste("VARS: MaxIntGap",MAX_INTERACTION_GAP, "s",", Capsule file =",CapsuleDef ), side = 3, line = -1.5, outer = TRUE)
  
}

#Save the uber-large output of all cut trajectories as computing takes minutes
#dput(interaction_AUTO_REP_PER, file = "/home/cf19810/Documents/Ants_behaviour_analysis/Data/interaction_AUTO_REP_PER_16feb22.txt")
#dput(summary_AUTO_REP_PER, file = "/home/cf19810/Documents/Ants_behaviour_analysis/Data/summary_AUTO_REP_PER_16feb22.txt")



# TO COMPARE PLOTS FOR THE FALSE POSITIVES (but also for the F negs) YOU WILL NEED TO USE vars from THE AUTO FILE
# AFTER THE MATRIX CALCULATION HAS APPENDED AGREEMENT COLUMNS, AND TO GENERATE PCAs (2 - FALSE POS, FALSE NEG, AGREEing - per every particular parameter of ComputeANTInteraction ) with subsets of the AGREEMENT COL = 0 or 1.
# ATTENTION:  FALSE NEGATIVES CAN BE GENERATED BY CREATING INVERSE DISAGREEMENT FOR LOOP (SEE ROW 77-84 PF AGREEMENT MATRIX SCRIPT)


#ADD ANOTHER CAPSULE DEF FOR TESTING
# DECIDE A RANGE OF PARAMS FOR COMPUTEANTINTERACTS


###################################################
#### SUMMARY MANUAL  ACT REC ID BY PARAMETERS #####
###################################################
summary_MANUAL_delta<- summary_MANUAL

#calculate delta parameters ACT-REC to check when ACT is bigger then REC
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

#calculate delta vars manually CAN BE REMOVED
summary_MANUAL$delta_speed_pxpersec <- summary_MANUAL$mean_speed_pxpersec_ACT - summary_MANUAL$mean_speed_pxpersec_REC

#reshape data for plotting. Split by REC and ACT
summary_data_ACT <-summary_MANUAL %>% dplyr::select(contains(c("BEH", "ACT"), ignore.case = TRUE))
summary_data_REC <- summary_MANUAL %>% dplyr::select(contains(c("BEH", "REC"), ignore.case = TRUE))
summary_data_ACT$IntRow <- rep(1:nrow(summary_data_ACT))
summary_data_REC$IntRow <- rep(1:nrow(summary_data_REC))
#Rename columns to make them match and bind+ melt columns
summa_data_ACT  <- summary_data_ACT %>% rename_with(~str_remove(., c("_ACT|Act_"))); summa_data_ACT$Role <- "ACT"
summa_data_REC  <- summary_data_REC %>% rename_with(~str_remove(., c("_REC|Rec_"))); summa_data_REC$Role <- "REC"
summa_data_bind <- rbind(summa_data_ACT,summa_data_REC)

# plot(summa_data_bind$mean_speed_pxpersec, as.factor(summa_data_bind$Role),type="o")
ggplot(summa_data_bind, aes(x = Role, y = mean_speed_pxpersec,colour=delta_speed_pxpersec > 0)) +
  geom_line(aes(group = IntRow)) +
  geom_point()


#temporary plot, understand steplength for DT_dist_THRESHOLD
# par(mfrow=c(1,2))
# hist(interaction_MANUAL$ACT.distance,breaks = 60,main="ACT stepwise distance",sub="blue line = 2") + abline(v=2,col='blue', lwd=2) 
# hist(interaction_MANUAL$REC.distance,breaks = 60, main="REC stepwise distance",sub="blue line = 2") + abline(v=2,col='blue', lwd=2) 
# plot dt frame!!!!



###################### TO DOS ##############################
############################################################
# - acceleration of the angle
# - frame by frame rate of change of angle
# - DELTA ANGLES formulas FIX


#check
unique(interaction_MANUAL$PERIOD)
unique(interaction_MANUAL$REPLICATE)


  # rm(list=(c("e")))
  # gc()
  
  
  

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
