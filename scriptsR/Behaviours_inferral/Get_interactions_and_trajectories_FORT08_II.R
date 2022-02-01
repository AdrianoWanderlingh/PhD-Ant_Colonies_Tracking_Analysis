##########################################################################################
############## THIS VERSION IS FORT 0.8.1 COMPATIBLE #####################################
##########################################################################################
#this should be the version of the script maintained for long term use.
#"https://formicidae-tracker.github.io/myrmidon/docs/latest/api/index.html"

rm(list=(c("e")))
gc()
rm(list=ls())

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

#TO DO:
#Conversion of data in mm. Makes sense to do that on extracted trajectories files, using the tag size and box as reference...
# reuse this and tags size in pixel/mm as point of reference - traj[,which(grepl("x",names(traj))|grepl("y",names(traj)))] <-   traj[,which(grepl("x",names(traj))|grepl("y",names(traj)))] * petridish_diameter / range 

############## DEFINE MY FUNCTIONS
Movement.angle.diff     <- function(x)  {
  if( length(x) > 0 & !is_null(x)){
  c( atan(x[-nrow(x), "x"] / x[-nrow(x), "y"]) - atan(x[-1, "x"] / x[-1, "y"]), NA)}
else {vector()}}

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

## PARAMETERS
USER            <- "Adriano"
Xmin <- 2000
Xmax <- 7900
Ymin <- 0000
Ymax <- 5500

N_DECIMALS <- 3 ## number of decimals to preserve when rounding to match between interactions & collisions
FUZZY_MATCH <- TRUE  ## fuzzy matching between data frames
MAX_INTERACTION_GAP <- 1

max_gap         <- fmHour(24*365)   ## important parameter to set! Set a maximumGap high enough that no cutting will happen. For example, set it at an entire year: fmHour(24*365)
desired_step_length_time    <- 0.125 ###in seconds, the desired step length for the analysis
ANT_LENGHT_PX <- 153 #useful for matcher::meanAntDisplacement mean and median value are similar

####### navigate to folder containing myrmidon file
if (USER=="Adriano") {WORKDIR <- "/home/cf19810/Documents/Ants_behaviour_analysis"}
if (USER=="Tom")     {WORKDIR <- "/media/tom/MSATA/Dropbox/Ants_behaviour_analysis"}

DATADIR <- paste(WORKDIR,"Data",sep="/")
SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")

## perform analyses on annotation data
#source(paste(SCRIPTDIR,"Annotations_analysis.R",sep="/"))

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
  experiment_name_end <-  str_sub(experiment_name,-16,-1)
  CapsuleDef <-  substr(experiment_name_end, 1, 7)
  ###tag statistics
  tag_stats <- fmQueryComputeTagStatistics(e)
  
  pdf(file=paste(DATADIR,"Interactions_Collsions_plots_25Jan2021.pdf", sep = ""), width=6, height=4.5)
  par(mfrow=c(2,3), mai=c(0.3,0.3,0.4,0.1), mgp=c(1.3,0.3,0), family="serif", tcl=-0.2)
  
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
    ###### extract ant trajectories from hand-annotated data ######################
    ###############################################################################
    ####First let's extract ant's trajectories
    #set plots parameters (for plotting coords)
  
    for (BEH in c("G"))#,"T","FR","CR"))
      {
      ## subset all hand-labelled bahavs for this behaviour type in this colony
      annot_BEH <- annotations[which(annotations$treatment_rep==REPLICATE  & annotations$period==PERIOD & annotations$Behaviour==BEH ),]
      ## loop through each event in annot_BEH
      for (ROW in 1:nrow(annot_BEH))
        { 
  
        ## extract actor, receiver IDs & start & end times from the hand-annotated data
        ACT <- annot_BEH$Actor[ROW]
        REC <- annot_BEH$Receiver[ROW]
     
        ENC_FRAME_start <- annot_BEH$frame_start[ROW]
        ENC_FRAME_stop  <- annot_BEH$frame_stop[ROW]
        
        Act_Name <- paste("ant",ACT,sep="_")
        Rec_Name <- paste("ant",REC,sep="_")
        print(paste("Behaviour:",BEH,"number",ROW,"Actor:",Act_Name,"Receiver:",Rec_Name))
        
        ## extract the trajectory for ACT
        traj_ACT <-  positions$trajectories[[Act_Name]]
        traj_REC <-  positions$trajectories[[Rec_Name]]
        
        # ##  convert POSIX format to raw secs since 1.1.1970; brute force for match.closest with collisions
        # traj_ACT$UNIX_secs <- as.numeric(traj_ACT$UNIX_time)  ## yes, it looks like the milisecs are gone, but they are there; traj_ACT$UNIX_secs
        # traj_REC$UNIX_secs <- as.numeric(traj_REC$UNIX_time) 
        # ## now round these decimal seconds to 3 d.p.
        # traj_ACT$UNIX_secs <- round(traj_ACT$UNIX_secs,N_DECIMALS) ## eliminates the 'noise' below the 3rd d.p., and leaves each frame existing just once, see: table(traj_ACT$UNIX_secs)
        # traj_REC$UNIX_secs <- round(traj_REC$UNIX_secs,N_DECIMALS)
        # 
        ## remove 'time' column as it is confusing - it's not a common time
        traj_ACT$time <- NULL
        traj_REC$time <- NULL
          
        # ## Plot trajectories of both actor & receiver, show on the same panel
        # plot  (y ~ x, traj_ACT, pch=".", col=rgb(0,0,1,0.3,1), main=Title, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
        # points(y ~ x, traj_REC, pch=".", col=rgb(1,0,0,0.3,1))
         
        ## subset the trajectories of both actor & receiver using the start & end times
        traj_ACT <- traj_ACT [ which(traj_ACT$frame >= ENC_FRAME_start & traj_ACT$frame <= ENC_FRAME_stop),]
        traj_REC <- traj_REC [ which(traj_REC$frame >= ENC_FRAME_start & traj_REC$frame <= ENC_FRAME_stop),]
        
        ## For some reason, the 4th decimal place of the UNIX times (thousandths of a sec) for two individuals seen in a given frame are not identical; cbind(format(traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS4"), format(traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS4"))
        ## so eliminate the 4th-nth decimal place ; retains accuract to a thousandth-of-a-second
        # traj_ACT$UNIX_time <- format(traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
        # traj_REC$UNIX_time <- format(traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
        # print(table(substr(x =  format( traj_ACT$UNIX_time, "%OS6"), start = 7, stop = 7)))
        # print(table(substr(x =  format( traj_REC$UNIX_time, "%OS6"), start = 7, stop = 7)))

      #   ## .. and convert back to POSIX format, durr
      #   traj_ACT$UNIX_time <- as.POSIXct(traj_ACT$UNIX_time , tz = "GMT", origin = "1970-01-01 00:00:00")
      #   traj_REC$UNIX_time <- as.POSIXct(traj_REC$UNIX_time , tz = "GMT", origin = "1970-01-01 00:00:00")
      #   ## use posixlt - more accurate than posixct, but still not perfect
      #   # traj_ACT$UNIX_time <- as.POSIXlt(traj_ACT$UNIX_time , tz = "GMT", origin = "1970-01-01 00:00:00")
      #   # traj_REC$UNIX_time <- as.POSIXlt(traj_REC$UNIX_time , tz = "GMT", origin = "1970-01-01 00:00:00")
      #   
      #   ## https://stackoverflow.com/questions/7726034/how-r-formats-posixct-with-fractional-seconds
      #   myformat.POSIXct <- function(x, digits=0) 
      #   {
      #     x2 <- round(unclass(x), digits)
      #     attributes(x2) <- attributes(x)
      #     x <- as.POSIXlt(x2)
      #     x$sec <- round(x$sec, digits)
      #     format.POSIXlt(x, paste("%Y-%m-%d %H:%M:%OS",digits,sep=""))
      #     }
      #   
      #   traj_ACT$UNIX_time <- as.POSIXct( myformat.POSIXct(traj_ACT$UNIX_time,3), tz = "GMT", origin = "1970-01-01 00:00:00")
      #   traj_REC$UNIX_time <- as.POSIXct( myformat.POSIXct(traj_REC$UNIX_time,3), tz = "GMT", origin = "1970-01-01 00:00:00")
      #   
      #   
      #   ## ERROR ERROR ! For some reason, as.POSIXct converts the 3.d.p. trimmed times to have e.g. '999' from the 4th-nth decimals ??!?!?! 
      #   
      #   format( traj_ACT$UNIX_time, "%OS3")
      #   format( traj_REC$UNIX_time, "%OS3")
      #   
      # data.frame( format( traj_ACT$UNIX_time, "%OS3"), 
      #             format( traj_ACT$UNIX_time, "%OS6"))
      #   
      #   print(table(substr(x =  format( traj_ACT$UNIX_time, "%OS6"), start = 7, stop = 7)))
      #   print(table(substr(x =  format( traj_REC$UNIX_time, "%OS6"), start = 7, stop = 7)))
      #   
        
        
        #check start time correspondance
        #print(paste("Behaviour:",BEH,"ACT:",Act_Name,"REC:",Rec_Name, "annot_start", ENC_FRAME_start, "traj_start", min(traj_ACT$frame,na.rm = TRUE)))
        # # ## Plot trajectories of both actor & receiver, show on the same panel
        Title <- paste(REPLICATE, ", ", PERIOD, ", ", BEH, ROW,", ", "Act:",ACT, ", ", "Rec:",REC, "\nframes ", ENC_FRAME_start, "-", ENC_FRAME_stop, sep="")
        plot   (y ~ x, traj_ACT, type="l", lwd=4, col="blue4",asp=1, main=Title,cex.main=0.9 ,xlim=c(min(traj_ACT$x,traj_REC$x),max(traj_ACT$x,traj_REC$x)),ylim=c(min(traj_ACT$y,traj_REC$y),max(traj_ACT$y,traj_REC$y))) #, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
        points (y ~ x, traj_REC, type="l", lwd=4,  col="red4",asp=1)
  
        #plot angles of both actor & receiver
        #plot(angle ~ time, traj_ACT, col=rgb(0,0,1,0.3,1), main=paste("angle_rad",BEH,", Act:",ACT, "Rec:",REC, ENC_TIME_start, "-", ENC_TIME_stop)) #,xlim = c(0,500),ylim = c(0,2*pi))
        #points(angle ~ time, traj_REC, col=rgb(1,0,0,0.3,1))
        
        ##################
        ## INDIVIDUAL TRAJECTORY MEASURES
        # uncorrelated vars used for Constance Tests:  median_acceleration_mmpersec2, rmse_mm, distance_mm, 
        # mean_acceleration_mmpersec2, straightness_index , periodicity_sec (last two not calculated here as they need rediscretisation)
        
        ### angular st.dev 
        StDev_angle_ACT                    <-  angular.deviation(traj_ACT$angle, na.rm = TRUE)
        StDev_angle_REC                    <-  angular.deviation(traj_REC$angle, na.rm = TRUE)
        
        # mean direction of a vector of circular data
        
        # ADRIANO to double check  wheTHER THE circular average is based on the wrong coordinate systEM - E.G. 0-360 CLOCKWISE !!!!!!!
        
        angle_mean_ACT                    <-  mean.circular(traj_ACT$angle, na.rm = TRUE) 
        angle_mean_REC                    <-  mean.circular(traj_REC$angle, na.rm = TRUE)
        # 
        # ##Delta_angles: differential between body angle and movement angle
        # # movement_angle: as the body_angle (orientation) remains the same, the movement_angle can change without if the movement is perpendicular to the body_angle
        # # change in movement angle calculated as atan(x2/y2) - atan(x1/y1)
        # traj_ACT  <- data.frame(time=double(), x=double(), y=double(), angle=double())#, UNIX_time=)
        # 
        # 
        # ## #ADRIANO TO FILL THESE IN
        # # if (nrow(traj_ACT)>1)
        # #   {
        #   
        #   # variation of movement angle frame by frame
        #   traj_ACT$Movement_angle_difference <- Movement.angle.diff(x = traj_ACT) ## To check that this function works, run this: traj_ACT$Movement_angle_difference_CHECK <- NA; for (i in 1: (nrow(traj_ACT)-1)) {traj_ACT$Movement_angle_difference_CHECK[i]   <- atan(traj_ACT[i, "x"] / traj_ACT[i, "y"]) -  atan(traj_ACT[i+1, "x"] / traj_ACT[i+1, "y"])}
        #   
        # # }else{print(paste("NO DATA FOR", ACT))}
        # 
        #   
        # 
        # # if (nrow(traj_REC)>1)
        # # {
        #   # variation of movement angle frame by frame
        #   
        #   traj_REC$Movement_angle_difference <- Movement.angle.diff(x = traj_REC)## anticlockwise turning ants have POSTIVE values of Movement.angle.diff & vice-versa - see sanity check above
        #   
        # # }else{print(paste("NO DATA FOR", REC))}
        # 
        # 
        #   
        # ## SANITY CHECK: create toy data to check the MEANING OF THE +/-'ve signs generated by Movement.angle.diff:
        # # traj_ACT <- traj_ACT[1:4,]
        # # ## anticlockwise turning ants have POSTIVE values of Movement.angle.diff
        # # traj_ACT$x <- c(0,1,1,0)
        # # traj_ACT$y <- c(0,0,1,1)
        # # ## clockwise turning ants have NEGATIVE values of Movement.angle.diff
        # # traj_ACT$x <- c(0,0,1,1)
        # # traj_ACT$y <- c(0,1,1,0)
        # 
        # #.
        # #... CONTINUE FROM HERE: 
        # 
        # ## and take the absolute value of the 'movement angle' 
        # traj_ACT$Movement_angle_difference_ABS <- abs(traj_ACT$Movement_angle_difference)
        # traj_REC$Movement_angle_difference_ABS <- abs(traj_REC$Movement_angle_difference)
        # 
        # # variation of Orientation angle frame by frame
        # traj_ACT$Orientation_diff <- Orientation.diff(x = traj_ACT)
        # traj_REC$Orientation_diff <- Orientation.diff(x = traj_REC)
        # 
        # 
        # 
        #traj_BOTH$angle_diff <- abs((traj_BOTH$REC.angle - pi) - (traj_BOTH$ACT.angle -pi)) %% (2*pi)
        # Orientation.diff_RE        <- function(x)  {   c(  abs( (x[-nrow(x), "angle"]-pi) - 
        #                                                         (x[-1,      "angle"])-pi) , NA)}
        # 
        # traj_ACT$Orientation_diff_RE <- Orientation.diff_RE(x = traj_ACT)
        # 
        # traj_ACT$Orientation_diff_RE_Normalised <- traj_ACT$Orientation_diff_RE %% (2*pi)
        # ## DOUBLE TRIPLE CHECK THIS!!
        # 
        
        ## SANITY CHECK : 
        # traj_ACT <- traj_ACT[1:4,]
        # traj_ACT$angle <- NA
        # traj_ACT$angle <- c(pi/2, pi - 0.001, -pi + 0.001, -pi/2)# ant moves from N to WNW, to WSW, to S
        # traj_ACT$Orientation_diff <- Orientation.diff(x = traj_ACT)
        # traj_ACT$Orientation_diff %% 2*pi
        # 
        ## sanity check 2: use circular package to handle everything
        # traj_ACT$angle_mirror_negative <- NA
        # traj_ACT$angle_mirror_negative <- pi - traj_ACT$angle[which(traj_ACT$angle<0)] 
        # traj_ACT$angle <- circular(traj_ACT$angle_mirror_negative, type="angles", units="radians", modulo = "pi", rotation="clock")
        
        # # delta_angles
        # traj_ACT$delta_angles <- traj_ACT$Movement_angle_difference - traj_ACT$Orientation_diff
        # traj_REC$delta_angles <- traj_REC$Movement_angle_difference - traj_REC$Orientation_diff
        # # mean delta_angles
        # mean_delta_angles_ACT <- mean(traj_ACT$delta_angles,na.rm=TRUE)
        # mean_delta_angles_REC <- mean(traj_REC$delta_angles,na.rm=TRUE)
        # 
        ##define trajectory
        #CHECK THAT THIS IS RIGHT OR DISCARD
        trajectory_ACT <- TrajFromCoords(data.frame(x=traj_ACT$x,y=traj_ACT$y, frames=traj_ACT$frame), spatialUnits = "px",timeUnits="frames")
        trajectory_REC <- TrajFromCoords(data.frame(x=traj_REC$x,y=traj_REC$y, frames=traj_REC$frame), spatialUnits = "px",timeUnits="frames")

        #trajectory                      <- TrajResampleTime (trajectory_ori,stepTime =desired_step_length_time )
        ### total distance moved
        moved_distance_px_ACT                  <- TrajLength(trajectory_ACT)
        moved_distance_px_REC                  <- TrajLength(trajectory_REC)
        ### Calculate trajectory derivatives
        deriv_traj_ACT                <- TrajDerivatives(trajectory_ACT)
        deriv_traj_REC                <- TrajDerivatives(trajectory_REC)
        #mean_speed_mmpersec_ACT        <- mean(deriv_traj_ACT$speed,na.rm=T)
        #  median_speed_mmpersec           <- median(deriv_traj$speed,na.rm=T)
        #  Mean and median acceleration
        mean_accel_pxpersec2_ACT     <- mean(deriv_traj_ACT$acceleration,na.rm=T) #units: pxpersec2
        mean_accel_pxpersec2_REC     <- mean(deriv_traj_REC$acceleration,na.rm=T) #units: pxpersec2
        #median_accel_ACT  <- median(deriv_traj_ACT$acceleration,na.rm=T) #units: pxpersec2
        #median_accel_REC  <- median(deriv_traj_REC$acceleration,na.rm=T) #units: pxpersec2
        ### Mean and median turn angle, in radians
        #  mean_turnangle_radians          <- mean(abs( TrajAngles(trajectory)),na.rm=T)
        #  median_turnangle_radians        <- median(abs( TrajAngles(trajectory)),na.rm=T)
        ### Straightness
        #  straightness_index              <- Mod(TrajMeanVectorOfTurningAngles(deriv_traj_ACT)) #requires redisretisation of the traj https://www.rdocumentation.org/packages/trajr/versions/1.4.0/topics/TrajMeanVectorOfTurningAngles
        #  sinuosity                       <- TrajSinuosity(trajectory)
        #  sinuosity_corrected             <- TrajSinuosity2(trajectory)
        ### Expected Displacement
        #  expected_displacement_mm        <- TrajEmax(trajectory,eMaxB =T)
        ### Autocorrelations
        # all_Acs                         <- TrajDirectionAutocorrelations(deriv_traj_ACT) #requires redisretisation of the traj
        # periodicity_sec                 <- TrajDAFindFirstMinimum(all_Acs)["deltaS"]*desired_step_length_time
        ### root mean square displacement
        rmsd_px_ACT                            <-  sqrt(sum( (traj_ACT$x-mean(traj_ACT$x))^2 + (traj_ACT$y-mean(traj_ACT$y))^2 )/length(na.omit(traj_ACT$x)))
        rmsd_px_REC                            <-  sqrt(sum( (traj_REC$x-mean(traj_REC$x))^2 + (traj_REC$y-mean(traj_REC$y))^2 )/length(na.omit(traj_REC$x)))
        
        
        #rename trajectories columns NOT TO BE MERGED for Act and Rec
        names(traj_ACT)[names(traj_ACT) == 'x'] <- 'ACT.x'; names(traj_ACT)[names(traj_ACT) == 'y'] <- 'ACT.y'; names(traj_ACT)[names(traj_ACT) == 'angle'] <- 'ACT.angle'
        names(traj_REC)[names(traj_REC) == 'x'] <- 'REC.x'; names(traj_REC)[names(traj_REC) == 'y'] <- 'REC.y'; names(traj_REC)[names(traj_REC) == 'angle'] <- 'REC.angle'
     
        #merge trajectories matching by time
        # traj_BOTH <- merge(traj_ACT, traj_REC,all=T, by=c("UNIX_time"))  ## CAREFUL NOT TO INCLUDE 'time' AS IF THE FIRST OBSERVATION TIME OF ACT & REC IS NOT IDENTICAL, THEN MERGE WILL TREAT THE SAME UNIX TIME AS DIFFERENT ...
        traj_BOTH <- merge(traj_ACT, traj_REC,all=T, by=c("frame"))  ## 21 Jan 2022: Changed to sue raw unix seconds to allow simpler matching below (posix formats cause problems with the matching...)
        
        
        ## merge drops the column attributes - add back the POSIX time format 
        # traj_BOTH$UNIX_time <- as.POSIXct(traj_BOTH$UNIX_time, tz="GMT",origin="1970-01-01 00:00:00")
        
        ## measure the length *in seconds* of the interaction between ACT & REC
        # interaction_length_secs <- as.numeric(difftime ( max(traj_BOTH$UNIX_time, na.rm=T), min(traj_BOTH$UNIX_time, na.rm=T), units="secs"))
        int_start_frame <- min(traj_BOTH$frame, na.rm=T)
        int_end_frame <- max(traj_BOTH$frame, na.rm=T)
        interaction_length_secs <-  (int_end_frame - int_start_frame)/8  ## 21 Jan 2022
        
        prop_time_undetected_ACT <- (sum(is.na(traj_BOTH$ACT.x)) / 8) / interaction_length_secs  ## the prop of the interaction in which ACT was seen 
        prop_time_undetected_REC <- (sum(is.na(traj_BOTH$REC.x)) / 8)  / interaction_length_secs ## UNITS are secs / secs --> proportion; MUST BE STRICTLY < 1 !!
       
        # FRAME BY FRAME PARAMETERS
        
        # distance walked
        # traj_BOTH$ACT.distance        <- c(NA, with(traj_BOTH, (sqrt(diff(ACT.x)^2 + diff(ACT.y)^2))) # euclidean distance
        # traj_BOTH$REC.distance        <- c(NA, with(traj_BOTH, sqrt(diff(REC.x)^2 + diff(REC.y)^2))) # euclidean distance
        
        ## ADRIAO TO FIX THIS IN OWN TIME!
        
        # traj_BOTH$UNIX_interval <- c(NA, as.numeric(diff(traj_BOTH$UNIX_time), units="secs"))
        # #speed
        # traj_BOTH$ACT.speed_PxPerSec        <- c( with(traj_BOTH, (sqrt(diff(ACT.x)^2 + diff(ACT.y)^2)) / UNIX_interval)) # it works as it corresponds to mean(deriv_traj_ACT$speed,na.rm=T)
        # traj_BOTH$REC.speed_PxPerSec        <- c(NA, with(traj_BOTH, (sqrt(diff(REC.x)^2 + diff(REC.y)^2)) / UNIX_interval)) # it works as it corresponds to mean(deriv_traj_ACT$speed,na.rm=T)
        # 
        # #  acceleration
        # traj_BOTH$ACT.accel_PxPerSec2 <- c(NA, with(traj_BOTH, diff(ACT.speed_PxPerSec)/UNIX_interval))
        # traj_BOTH$REC.accel_PxPerSec2 <- c(NA, with(traj_BOTH, diff(REC.speed_PxPerSec)/UNIX_interval))
        # 
        # # jerk (diff in accelerations)
        # traj_BOTH$ACT.jerk_PxPerSec2 <- c(NA, with(traj_BOTH, diff(ACT.accel_PxPerSec2)/UNIX_interval))
        # traj_BOTH$REC.jerk_PxPerSec2 <- c(NA, with(traj_BOTH, diff(REC.accel_PxPerSec2)/UNIX_interval))
        
        ##################
        ## INTERACTING PAIR TRAJECTORY MEASURES
  
        #straight line - euclidean distance
        traj_BOTH$straightline_dist_px <-  sqrt((traj_BOTH$ACT.x-traj_BOTH$REC.x)^2+(traj_BOTH$ACT.y-traj_BOTH$REC.y)^2)
        strghtline_dist_px <- mean(traj_BOTH$straightline_dist_px, na.rm=TRUE)
        #angular difference
        traj_BOTH$angle_diff <- abs((traj_BOTH$REC.angle - pi) - (traj_BOTH$ACT.angle -pi)) %% (2*pi)
        
        ###############
      
        #calculate acceleration row by row
        # traj_BOTH$accel_pxpersec2_ACT     <- deriv_traj_ACT$acceleration #units: pxpersec2
        # traj_BOTH$accel_pxpersec2_REC     <- deriv_traj_REC$acceleration #units: pxpersec2
        
        summary_MAN_ROW <- data.frame(REPLICATE, PERIOD, BEH, ROW, Act_Name, Rec_Name, 
                   # StDev_angle_ACT=StDev_angle_ACT, StDev_angle_REC=StDev_angle_REC,
                   # angle_mean_ACT=angle_mean_ACT, angle_mean_REC=angle_mean_REC,
                   # #mean_delta_angles_ACT=mean_delta_angles_ACT, mean_delta_angles_REC=mean_delta_angles_REC,
                   # moved_distance_px_ACT=moved_distance_px_ACT, moved_distance_px_REC=moved_distance_px_REC,
                   # mean_accel_pxpersec2_ACT=mean_accel_pxpersec2_ACT, mean_accel_pxpersec2_REC=mean_accel_pxpersec2_REC,
                   # #median_accel_ACT=median_accel_ACT,median_accel_REC=median_accel_REC,
                   # rmsd_px_ACT=rmsd_px_ACT,rmsd_px_REC=rmsd_px_REC,
                   int_start_frame , int_end_frame, interaction_length_secs,
                   prop_time_undetected_ACT, prop_time_undetected_REC,
                   # strghtline_dist_px=strghtline_dist_px,
                   #when adding a new variable, it must be included in the reshape rule for data plotting
                   stringsAsFactors = F)
         
        #NO UNIX TIME BUT FRAMES. interaction AREA= NEST, FORAGING delete x y coords
        interacts_MAN_ROW <- data.frame(REPLICATE=REPLICATE,PERIOD=PERIOD,ROW=ROW,BEH=BEH,Act_Name=Act_Name,Rec_Name=Rec_Name,
                                           traj_BOTH,
                                           stringsAsFactors = F)
        
        ## stack -by the end of the loop, this should just contain all events for replicate X, period Y
        interacts_MAN_REP_PER <- rbind(interacts_MAN_REP_PER, interacts_MAN_ROW)
        summary_MAN_REP_PER   <- rbind(summary_MAN_REP_PER,     summary_MAN_ROW)
        
        }##ROW
      }##BEH
    
    ## stack
    interaction_MANUAL <- rbind(interaction_MANUAL, interacts_MAN_REP_PER)
    summary_MANUAL     <- rbind(summary_MANUAL,       summary_MAN_REP_PER)
    
    ##}
    #}
    
    #store these elsewhere as they are generated by the last iteration of the outer loop.
    #ONCE ANNOTATION FILE IS DONE FOR GOOD, PUT THEM AMONG THE VARIABLES ON TOP AS NUMERIC VALUE KEEPING THE FORMULA COMMENTED 
    AntDistanceSmallerThan <- max(interaction_MANUAL$straightline_dist_px,na.rm = T)
    AntDistanceGreaterThan <- min(interaction_MANUAL$straightline_dist_px,na.rm = T)
   # AntDisplacement <- #max trajectory step lenght per each ant during interaction, use higher value for matcher
    
    #create new variable by pasting ant numbers "low,high" for summary_MAN_REP_PER
    summary_MAN_REP_PER$ant1 <- as.numeric(gsub("ant_","", summary_MAN_REP_PER$Act_Name))
    summary_MAN_REP_PER$ant2 <- as.numeric(gsub("ant_","", summary_MAN_REP_PER$Rec_Name))
    ## the ant 1 & ant 2 labels are not strictly ascending, so sort them so ant 1 alwasy < ant 2
    summary_MAN_REP_PER$pair <- apply(summary_MAN_REP_PER[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
    
    
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
                                 fmMatcherAntDistanceGreaterThan(AntDistanceGreaterThan)
                                  ))     #check every 5 seconds if ant has displaced more than ANT_LENGHT_PX 
    #fmMatcherAntDisplacement(ANT_LENGHT_PX, 52)
    
    #NOTES: THESE MAY ALL BE WRONG!!!!
    # - adding AntDistances removes false positives but has no effect on false negatives
    # - with MaxGap=5 there is a reduction in f.neg. and f. positives compared to MaxGap=10
    # - with MaxGap=3 there is a reduction in f.neg. and f. positives compared to MaxGap=5
    # - with MaxGap=1 there is a reduction in f.neg. and f. positives compared to MaxGap=3
    
    ## sequentially vary the interaction gap-filling to check what effect this has on the agreement between the MANUAL & AUTOMATIC interactions
    for (Buffer in seq(0,30,5))
      {
      for (MAX_INTERACTION_GAP in c(seq(1,9,2), seq(10,60,10)))
        {
       
        interacts_AUTO_REP_PER <- fmQueryComputeAntInteractions(e,
                                                                start=time_start, 
                                                                end=time_stop,
                                                                maximumGap =fmSecond(MAX_INTERACTION_GAP), ## WHEN A PAIR DISENGAGE, HOW LONG IS THE INTERVAL? 
                                                                reportFullTrajectories = T,
                                                                matcher = matcherCapTypeAntDists)
        
        # ##  convert POSIX format to raw secs since 1.1.1970; brute force for match.closest with collisions
        # interacts_AUTO_REP_PER$interactions$int_start_secs  <- as.numeric(interacts_AUTO_REP_PER$interactions$start)  ## yes, it looks like the milisecs are gone, but they are there; traj_ACT$UNIX_secs
        # interacts_AUTO_REP_PER$interactions$int_end_secs     <- as.numeric(interacts_AUTO_REP_PER$interactions$end) 
        # ## now round these decimal seconds to 3 d.p.
        # interacts_AUTO_REP_PER$interactions$int_start_secs   <- round(interacts_AUTO_REP_PER$interactions$int_start_secs,   N_DECIMALS) ## eliminates the 'noise' below the 3rd d.p., and leaves each frame existing just once, see: table(traj_ACT$UNIX_secs)
        # interacts_AUTO_REP_PER$interactions$int_end_secs     <- round(interacts_AUTO_REP_PER$interactions$int_end_secs, N_DECIMALS)
        
        ## Assign frame number
        interacts_AUTO_REP_PER$interactions["int_start_frame"]   <- lapply(interacts_AUTO_REP_PER$interactions["start"] , function(x) IF_frames$frame_num[match(x, IF_frames$time)])
        interacts_AUTO_REP_PER$interactions["int_end_frame"]   <- lapply(interacts_AUTO_REP_PER$interactions["end"] , function(x) IF_frames$frame_num[match(x, IF_frames$time)])
        
        interacts_AUTO_REP_PER$interactions $pair <- paste(interacts_AUTO_REP_PER$interactions$ant1, interacts_AUTO_REP_PER$interactions$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...
        
        interacts_AUTO_REP_PER$interactions$Duration <- interacts_AUTO_REP_PER$interactions$int_end_frame - interacts_AUTO_REP_PER$interactions$int_start_frame
        hist(interacts_AUTO_REP_PER$interactions$Duration,breaks = 30)
        nrow(interacts_AUTO_REP_PER$interactions)
        
        ## CHECK WHAT MAX_INTERACTION_GAP is doing - mean interval between successive interactions for a pair {A,B}
        PairInterInteractInterv <- data.frame(Pair=unique(interacts_AUTO_REP_PER$interactions$pair))  ## 1 element per pair
        PairInterInteractInterv$Mean_Interval <- NA
        PairInterInteractInterv$Min_Interval <- NA
        for (PR in PairInterInteractInterv$Pair)
          {
          ## extract successive interactions for PR
          interacts_AUTO_REP_PER_PAIR <- interacts_AUTO_REP_PER$interactions [ which(interacts_AUTO_REP_PER$interactions$pair==PR),] 
          
          if (nrow(interacts_AUTO_REP_PER_PAIR)>2) ## need at least two inteactions to have an inter-interaction interval
            {
            ## measure inter-interaction interval
            Inter_Interact_Interval <- interacts_AUTO_REP_PER_PAIR$int_end_frame   [2:nrow(interacts_AUTO_REP_PER_PAIR)] - 
                                       interacts_AUTO_REP_PER_PAIR$int_start_frame [1:(nrow(interacts_AUTO_REP_PER_PAIR)-1)]
            ## get mean interval for this pair
            Mean_Interval <- mean(Inter_Interact_Interval, na.rm=T)
            Min_Interval <- min(Inter_Interact_Interval, na.rm=T)  ## SHOUL DNEVER BE LOWER THAN MAX_INTERACTION_GAP
            PairInterInteractInterv$Mean_Interval [which(PairInterInteractInterv==PR)] <-  Mean_Interval
            PairInterInteractInterv$Min_Interval  [which(PairInterInteractInterv==PR)] <-  Min_Interval
            }
          }
        GrandMeanInterval <- mean(PairInterInteractInterv$Min_Interval, na.rm=T)
        GrandMinInterval  <- min(PairInterInteractInterv$Min_Interval, na.rm=T)
        

        # summaryTrajectory_all       <- interacts_AUTO_REP_PER$summaryTrajectory  ## same object as positions$trajectory_summary described above - however here it will have many more rows given that the trajectories will be broken whenever there is a 10 second non-detection gap for an ant.
        # ## therefore in this case my trick with ant_ID_str won't work - it only works if there's exactly one trajectory as ants in the data
        # ## so I would NOT use the output of this function to analyse trajectories
        # trajectories_all           <- interacts_AUTO_REP_PER$trajectories       ## same object as positions$trajectories described above - but containing more objects for the same reason as stated above
        # interacts_AUTO_REP_PER            <- interacts_AUTO_REP_PER$interactions       ## data frame containing all interactions
        # 
        
        ########################################################################################################################
        ################# cross reference manually-laeblled start & end times with auto interactions ###########################
        ########################################################################################################################
        
         ## TRUE POSITIVE RATE
         interacts_AUTO_REP_PER_OVERLAP_ALL <- NULL
         Buffer <- 0
         for (I in 1:nrow(summary_MAN_REP_PER))  ## loop across each (time-aggregated) manually-labelled behaviour
            {
            ## Find start and end times of each *aggregated* **MANUALLY-LABELLED** interaction
            Man_Start <- summary_MAN_REP_PER$int_start_frame [I] - Buffer
            Man_Stop  <- summary_MAN_REP_PER$int_end_frame [I]  + Buffer
            
            ## extract manual participants for the Ith interaction in the auto interactions
            Man_Pair <- summary_MAN_REP_PER$pair [I]
          
            
            ## look for the Ith MANUAL interaction in the AUTOMATIC list
            interacts_AUTO_REP_PER_OVERLAP <- unique(interacts_AUTO_REP_PER$interactions [ interacts_AUTO_REP_PER$interactions$int_start_frame >= Man_Start &
                                                                                           interacts_AUTO_REP_PER$interactions$int_end_frame   <= Man_Stop & 
                                                                                           interacts_AUTO_REP_PER$interactions $pair      == Man_Pair , ])
            
       
            
            ## stack the auto interactions when there's a match with the manual summary interaction
            if (nrow(interacts_AUTO_REP_PER_OVERLAP)>0)
              {
              interacts_AUTO_REP_PER_OVERLAP_ALL <- rbind(interacts_AUTO_REP_PER_OVERLAP_ALL, 
                                                          interacts_AUTO_REP_PER_OVERLAP)
              
              interacts_AUTO_REP_PER_OVERLAP_ALL <- unique(interacts_AUTO_REP_PER_OVERLAP_ALL)
              }
            }#I
       
        ## what is the overlap in seconds?
        interacts_AUTO_REP_PER_OVERLAP_ALL$Duration <- interacts_AUTO_REP_PER_OVERLAP_ALL$int_end_frame - interacts_AUTO_REP_PER_OVERLAP_ALL$int_start_frame
        
        ## the observed overlap time
        Overlap <- sum(interacts_AUTO_REP_PER_OVERLAP_ALL$Duration)
        
        ## what is the expected interaction time ?
        Expectation <- sum(summary_MAN_REP_PER$interaction_length_secs)
        ## stack
        Sensitivity <- rbind(Sensitivity, data.frame(REPLICATE, PERIOD, Buffer, MAX_INTERACTION_GAP, GrandMeanInterval, GrandMinInterval, Overlap, Expectation, Hit_Rate=100 * (Overlap/Expectation)) )

        

#####################################################################################
        
        #oder dataframe to make it work with Fuzzy matcher which needs sorted dataframes
        sum_MAN_REP_PER_ord <- summary_MAN_REP_PER[with(summary_MAN_REP_PER, order(int_start_frame)),]
        sum_MAN_REP_PER_ord1 <- summary_MAN_REP_PER[with(summary_MAN_REP_PER, order(int_end_frame)),]

        TOLERANCE <- 20 #number of seconds of tolerance window, is like the previous buffer
        FUZZY_INT_MATCH <- TRUE
        A_in_M_OVERLAP_ALL <- NULL
        if (FUZZY_INT_MATCH==TRUE)
        {
          for (I in 1:nrow(summary_MAN_REP_PER))  ## loop across each (time-aggregated) manually-labelled behaviour
          {
            ## Find start and end times of each *aggregated* **MANUALLY-LABELLED** interaction
            Man_Start <- summary_MAN_REP_PER$int_start_frame [I]
            Man_Stop  <- summary_MAN_REP_PER$int_end_frame [I]

            ## extract manual participants for the Ith interaction in the auto interactions
            Man_Pair <- summary_MAN_REP_PER$pair [I]

          ### TRUE POSITIVES
          # match.start and end
          #start
          A_in_M_rows_start <- unique(interacts_AUTO_REP_PER$interactions[ match.closest(x =     interacts_AUTO_REP_PER$interactions$int_start_frame,
                                                                                         table = sum_MAN_REP_PER_ord$int_start_frame , tolerance = TOLERANCE)
                                                                           & interacts_AUTO_REP_PER$interactions $pair == Man_Pair  , ])
          #end
          A_in_M_rows_end   <- unique(interacts_AUTO_REP_PER$interactions[ match.closest(x =     interacts_AUTO_REP_PER$interactions$int_end_frame,
                                                                                         table = sum_MAN_REP_PER_ord1$int_end_frame , tolerance = TOLERANCE)
                                                                           & interacts_AUTO_REP_PER$interactions $pair == Man_Pair  , ])
          #issue raise with times, drop them
          A_in_M_rows_end <-subset(A_in_M_rows_end, select = -c(start,end))
          A_in_M_rows_start <-subset(A_in_M_rows_start, select = -c(start,end))

          A_in_M_OVERLAP <- generics::intersect(A_in_M_rows_start, A_in_M_rows_end)
          #remove empty rows (all NA's)
          A_in_M_OVERLAP <- A_in_M_OVERLAP[rowSums(is.na(A_in_M_OVERLAP)) != ncol(A_in_M_OVERLAP),]
          
          #### A_in_M ARE TRUE POSITIVES  -> ROWS IN BOTH AUTO E MAN
        
        ## stack the auto interactions when there's a match with the manual summary interaction
        if (nrow(A_in_M_OVERLAP)>0)
        {
          A_in_M_OVERLAP_ALL <- rbind(A_in_M_OVERLAP_ALL, 
                                                      A_in_M_OVERLAP)
          A_in_M_OVERLAP_ALL <- unique(A_in_M_OVERLAP_ALL)
        }
      }#I
    }#FUZZY_INT_MATCH    

 #FALSE NEGATIVE RATE

 ####ERRORRRRRR AUTO_IN_MANUAL SHOULD BE LESS OR EQUAL TO SUMMARY_MANUAL!!!!! DUPS AND NAS ARE REMOVED, WHAT'S HAPPENING???         
 #cure momentarily by trimming auto$duration < manual$duration
 nrow(A_in_M_OVERLAP_ALL) 
 nrow(summary_MAN_REP_PER)
 
 summary_MAN_REP_PER$Duration <- summary_MAN_REP_PER$int_end_frame - summary_MAN_REP_PER$int_start_frame
 A_in_M_OVERLAP_ALL$Duration <- A_in_M_OVERLAP_ALL$int_end_frame - A_in_M_OVERLAP_ALL$int_start_frame
 
 min(summary_MAN_REP_PER$Duration)
 min(A_in_M_OVERLAP_ALL$Duration)
 
 min( summary_MAN_REP_PER$Duration[summary_MAN_REP_PER$Duration!=min(summary_MAN_REP_PER$Duration)] )
 
 MIN_TRIM <- min( summary_MAN_REP_PER$Duration[summary_MAN_REP_PER$Duration!=min(summary_MAN_REP_PER$Duration)] )
 A_in_M_OVERLAP_ALL_trim <- A_in_M_OVERLAP_ALL[which(A_in_M_OVERLAP_ALL$Duration > MIN_TRIM),]

 
 hist(summary_MAN_REP_PER$Duration,breaks = 30)
 hist(A_in_M_OVERLAP_ALL$Duration,breaks = 30)
 hist(A_in_M_OVERLAP_ALL_trim$Duration,breaks = 30)
        
######################################################################################        

        #plot the interactions by pair as timeline
        #generate 1 plot per every iteration 
        summary_MAN_REP_PER_sub <- summary_MAN_REP_PER[,c("REPLICATE","PERIOD","pair","int_start_frame","int_end_frame")] ; summary_MAN_REP_PER_sub$flag <- "manual"
        interacts_AUTO_REP_PER_sub <- A_in_M_OVERLAP_ALL_trim[,c("pair","int_start_frame","int_end_frame")] ; interacts_AUTO_REP_PER_sub$flag <- "auto"
        AUTO_MAN_REP_PER <- dplyr::bind_rows(summary_MAN_REP_PER_sub,interacts_AUTO_REP_PER_sub) #use bind_rows to keep rep info https://stackoverflow.com/questions/42887217/difference-between-rbind-and-bind-rows-in-r
        # AUTO_MAN_REP_PER$int_end_frame <- AUTO_MAN_REP_PER$int_end_frame-min(AUTO_MAN_REP_PER$int_start_frame)
        # AUTO_MAN_REP_PER$int_start_frame <- AUTO_MAN_REP_PER$int_start_frame-min(AUTO_MAN_REP_PER$int_start_frame)
        
        # ggplot(summary_MAN_REP_PER) +
        #   geom_linerange(aes(y = pair, xmin = int_start_frame, xmax = int_end_frame)) +#,
        #   labs(title = "Grooming by pair",
        #        subtitle = paste(unique(REPLICATE),unique(PERIOD))) +
        #   theme_minimal()
 
        #pdf(file=paste(DATADIR,"Interactions","AUTO_MAN_REP_PER","Gap",MAX_INTERACTION_GAP,"matcherCapTypeAntDists.pdf", sep = "_"), width=10, height=60)
        #ADD CAPSULES SIZE USED BY MODIFING THE OUTPUT OF CapsuleDef
        #CapsuleDef
        ggplot(AUTO_MAN_REP_PER) +
          geom_linerange(aes(y = pair, xmin = int_start_frame, xmax = int_end_frame,colour = flag),size=3,alpha = 0.5) +
          theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
          labs(title = paste("Grooming by pair" ,unique(REPLICATE),unique(PERIOD),"- MAN vs AUTO"),
               subtitle = paste("Nrows auto detected:" ,NROW(summary_MAN_REP_PER),"Nrows manual annotated:" ,NROW(A_in_M_OVERLAP_ALL_trim),"\nMaxIntGap",MAX_INTERACTION_GAP, "s" )) +
          scale_color_manual(values = c("manual" = "red",
                                         "auto"="black"))
       # dev.off()
        
      }
    }
    
    ## Select ONLY those AUTO interactions that are INSIDE the manual interactions
    plot( Sensitivity[Sensitivity$Buffer==0 , c("MAX_INTERACTION_GAP","Buffer","GrandMinInterval","Overlap","Hit_Rate")])    
    
    
    
    ########################################################################################################################
    ########################################################################################################################
    # COLLISIONS
    
    ##creates a ID string for each ant: ant_1, ant_2,...
    collisions$collisions$ant1_str <- paste("ant_",collisions$collisions$ant1,sep="") 
    collisions$collisions$ant2_str <- paste("ant_",collisions$collisions$ant2,sep="")
    
    ###collisions contains 3 objects: frames, positions and collisions
    collisions_frames    <- collisions$frames      ###data frames giving detail (time, space, width,height) about each frame
    collisions_positions <- collisions$positions   ###list containing as many objects as there are frames in collisions_frames. For each one, lists the positions of all detected ants on that frame
    ###(those two are the exactly the same as the output from fmQueryIdentifyFrames)
    collisions_coll <- collisions$collisions      #### dataframe containing all detected collisions
    
    ## trim the UNIX time to 3 decimals as we did for the hand-labelled interactions above _ DOESN@T WORK :-()
    # collisions_frames$UNIX_time <- format(collisions_frames$time, "%Y-%m-%d %H:%M:%OS3"); collisions_frames$time <- NULL
    ## .. and convert back to POSIX format, dur
    # collisions_frames$UNIX_time <- as.POSIXct(collisions_frames$UNIX_time , tz = "GMT", origin = "1970-01-01 00:00:00")
    
    ## for reliable matching with the events in interactions, above (don't trust POSIX times!)
    #collisions_frames$UNIX_secs <- round(as.numeric(collisions_frames$time),N_DECIMALS)

    ###let's view the first few lines of collisions
    head(collisions_frames) 
    head(collisions_coll) ###columns giving the ID of the two colliding ants, the zone where they are, the types of collisions (all types), and the frames_row_index referring to which frame that collision was detected in (matches the list indices in collisions_positions)
    head(collisions_positions)
    #use rownames of collision_frames to filter frames of collisions according to collisions-coll (which acts as a Look Up Table)
    #rename frame_num for the merging
    collisions_frames$frames_row_index <- collisions_frames$frame_num
    
    #merge collisions_frame and coll_coll based on index
    collisions_merge <- merge(collisions_coll, collisions_frames[,-match(c("height","width"),colnames(collisions_frames))], by="frames_row_index")
  
    
    #check that to the same frame corresponds the same time
    # format(collisions_merge$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # format(collisions_frames$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # format(    traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # format(    traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # 
    # format(    traj_BOTH$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # format(interacts_MAN_ROW$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # 
    # table(substr(x =  format( traj_BOTH$UNIX_time, "%OS6"), start = 7, stop = 7))
    
    
    ##if(collisions_coll$frames_row_index[125] == collisions_merge$frames_row_index[125]) {print("TRUE")}
    
    # add new column to interaction data when conditions are met (eg. time + ant ID)
    #RENAME COLUMN TRAJBOTHREC IN time 
    #merge with equal time and any order of the couple Act_Name-Rec_Name=ant1_str-ant2_str (any order). HOW?
    nrow(collisions_merge)
    str(collisions_merge)
    nrow(interacts_MAN_REP_PER)
    str(interacts_MAN_REP_PER)
    
    #create new variable by pasting ant numbers "low,high" for collisions_merge
    # collisions_merge$pair <- apply(collisions_merge[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
    collisions_merge $pair <- paste(collisions_merge$ant1, collisions_merge$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...
    
    #create new variable by pasting ant numbers "low,high" for interacts_MAN_REP_PER
    interacts_MAN_REP_PER$ant1 <- gsub("ant_","", interacts_MAN_REP_PER$Act_Name)
    interacts_MAN_REP_PER$ant2 <- gsub("ant_","", interacts_MAN_REP_PER$Rec_Name)
    ## the ant 1 & ant 2 labels are not strictly ascending, so sort them so ant 1 alwasy < ant 2
    interacts_MAN_REP_PER$pair <- apply(interacts_MAN_REP_PER[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })

    # check that the time formats are the same
    attributes(interacts_MAN_REP_PER$UNIX_time[1])
    attributes(collisions_merge$UNIX_time[1])
    
    ## check that the pair-time combinations in interacts_MAN_REP_PER are in collisions_merge
    table(  paste(interacts_MAN_REP_PER$pair) %in% 
              paste(collisions_merge$pair) )
    table(  paste(interacts_MAN_REP_PER$UNIX_secs) %in% 
              paste(collisions_merge$UNIX_secs) )
    table(  paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$UNIX_secs) %in% 
              paste(collisions_merge$pair,collisions_merge$UNIX_secs) )
    ## get a list of the missing pair-time combinations:
    paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$UNIX_secs) [!paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$UNIX_secs) %in% paste(collisions_merge$pair,collisions_merge$UNIX_secs) ]
    
    #join dataframes to have collisions on interacts_MAN_REP_PER  
    #FUNCTION: WHEN pair=pair AND time=UNIX_time -> ASSIGN collisions_merge$types TO interacts_MAN_REP_PER
    if (FUZZY_MATCH==FALSE)
      {
      interacts_MAN_REP_PER <- plyr::join(x = interacts_MAN_REP_PER, 
                                             y = collisions_merge[,c("UNIX_secs","pair","types")], 
                                             by=c("UNIX_secs","pair"), type="left")
      }
    ## alternative: fuzzy match
    if (FUZZY_MATCH==TRUE)
      {
      I_in_C_rows <- match.closest(x = interacts_MAN_REP_PER$UNIX_secs,
                                   table = collisions_merge$UNIX_secs, tolerance = 0.125)
      ## copy tpes across
      interacts_MAN_REP_PER$types  <- collisions_merge$types [I_in_C_rows]
      }
    
    ## PROBLEM: *MANY* rows in interactions aren't present in collisions; plot the spatial interactions to see whether this makes sense or not:
    par(mfrow=c(3,4), mai=c(0.3,0.3,0.4,0.1))
    for(BH in unique(interacts_MAN_REP_PER$BEH))
    {
      for (RW in unique(interacts_MAN_REP_PER$ROW))
      {
        IntPair <- interacts_MAN_REP_PER[which(interacts_MAN_REP_PER$BEH==BH & interacts_MAN_REP_PER$ROW==RW),]
        ## check if this interaction is present in collisions
        paste(IntPair$pair[1],IntPair)
        ## plot it
        #to fix as some of the vars are defined only in the previous function #TitleInt <- paste(REPLICATE, ", ", PERIOD, ", ", BH, RW,", ", "Act:",ACT, ", ", "Rec:",REC, "\n", ENC_TIME_start, "-", ENC_TIME_stop, sep="")
        plot(NA, xlim=range(c(IntPair$ACT.x,IntPair$REC.x),na.rm=T),  ylim=range(c(IntPair$ACT.y,IntPair$REC.y),na.rm=T), type="n", main=paste(RW,BH), xlab="x", ylab="y")
        points (ACT.y ~ ACT.x, IntPair, type="p", col="blue4"); lines (ACT.y ~ ACT.x, IntPair, col="blue4")
        points (REC.y ~ REC.x, IntPair, type="p", col="red4"); lines (REC.y ~ REC.x, IntPair, col="red4")
        ## show the headings of each ACT
        arrows.az (x = IntPair$ACT.x, 
                   y = IntPair$ACT.y, 
                   azimuth = IntPair$ACT.angle, 
                   rho = 10,
                   HeadWidth=0.1,
                   units="radians", Kol="blue2", Lwd=1)
        ## show the headings of each REC
        arrows.az (x = IntPair$REC.x, 
                   y = IntPair$REC.y, 
                   azimuth = IntPair$REC.angle, 
                   rho = 10,
                   HeadWidth=0.1,
                   units="radians", Kol="red2", Lwd=1)
        
        
        ## add lines connecting ACT & REC when there is a capsule overlap 
        IntPairCapOverlap  <- IntPair[which(!is.na(IntPair$types)),]
        if (nrow(IntPairCapOverlap)>0)
           {
           segments(x0 = IntPairCapOverlap$ACT.x, y0 = IntPairCapOverlap$ACT.y,
                   x1 = IntPairCapOverlap$REC.x, y1 = IntPairCapOverlap$REC.y)
           }
        
      }
    }
    
    ## if the missing rows occur when the distance between Act & Rec is large, check:
    interacts_MAN_REP_PER$types_present_absent <- "absent"
    interacts_MAN_REP_PER$types_present_absent[which(!is.na(interacts_MAN_REP_PER$types))] <- "present"
    boxplot(straightline_dist_px ~ types_present_absent, interacts_MAN_REP_PER)
    ### ... so rows in interactions that don't match collisions is not  a function of distance ...
    
    

    #cut collisions for the specific G 1 case

    # ROW2 <- 1 #first behaviour

    # interactio_data_SELECTED <- interaction_MANUAL[which(interaction_MANUAL$ROW==ROW2  & interaction_MANUAL$BEH=="G" & interaction_MANUAL$PERIOD=="post"),]
    # 
    # COLL_TIME_start <- min(interactio_data_SELECTED$UNIX_secs)
    # COLL_TIME_stop  <- max(interactio_data_SELECTED$UNIX_secs)
    # PAIR <- "1_5"
    # 
    # ## subset the collisions using the start & end times
    # collisions_merge_SELECTED <- collisions_merge[ which(collisions_merge$UNIX_secs >= COLL_TIME_start & collisions_merge$UNIX_secs <= COLL_TIME_stop & collisions_merge$pair==PAIR),]
    # min(collisions_merge_SELECTED$UNIX_secs)
    # max(collisions_merge_SELECTED$UNIX_secs)
    # 
    # format( collisions_merge_SELECTED$UNIX_secs, nsmall=5)
    # format( interactio_data_SELECTED$UNIX_secs, nsmall=5)
    # 
    # 
    # tail(format(collisions_merge_SELECTED$UNIX_secs, "%Y-%m-%d %H:%M:%OS6"),10)
    # tail(format(interactio_data_SELECTED$UNIX_secs, "%Y-%m-%d %H:%M:%OS6"),10)

    ########################################################################################################################
    ########################################################################################################################
    
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
dev.off() # save to pds!




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



#check
unique(interaction_MANUAL$PERIOD)
unique(interaction_MANUAL$REPLICATE)
min(interaction_MANUAL$UNIX_time)
max(interaction_MANUAL$UNIX_time)



###############################################################################
###### PARAMETERS PLOTS #######################################################
###############################################################################

#SUMMARY_DATA
#descriptive analysis
str(summary_MAN_REP_PER)

#reshape data for plotting. Split by REC and ACT
summary_data_ACT <-summary_MAN_REP_PER %>% dplyr::select(contains(c("BEH", "ACT","interaction_length","strghtline"), ignore.case = TRUE))
summary_data_REC <- summary_MAN_REP_PER %>% dplyr::select(contains(c("BEH", "REC","interaction_length","strghtline"), ignore.case = TRUE))
#Rename columns to make them match and bind+ melt columns
summ_data_ACT  <- summary_data_ACT %>% rename_with(~str_remove(., c("_ACT|Act_"))); summ_data_ACT$Role <- "ACT"
summ_data_REC  <- summary_data_REC %>% rename_with(~str_remove(., c("_REC|Rec_"))); summ_data_REC$Role <- "REC"
summ_data_bind <- rbind(summ_data_ACT,summ_data_REC)
summ_data_long <- melt(summ_data_bind,id.vars=c("BEH","Role","Name")) #explanation on the warning message https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message

#subset data by BEH for plotting
summary_data_G <- summ_data_long[which(summ_data_long$BEH == "G"),]
summary_data_T <- summ_data_long[which(summ_data_long$BEH == "T"),]
summary_data_FR <- summ_data_long[which(summ_data_long$BEH == "FR"),]
summary_data_CR <- summ_data_long[which(summ_data_long$BEH == "CR"),]
#keep Trophallaxis, Cross Rest and Front REst togheter and fill by BEH when plotting
summary_data_T_FR_CR <- summ_data_long[which(!summ_data_long$BEH == "G"),]


#set up plot
pdf(file=paste(DATADIR,"Parameters_plots_post_30minWindow_CHECK_NAME_IS_CORRECT.pdf", sep = ""), width=6, height=6)
par(mfrow=c(2,3), family="serif" , mar = c(0.1, 0.1, 2.2, 0))

###plot divided by variable and Role for Grooming
vars_plot_G <- ggplot(summary_data_G, aes(value, fill = Role)) +
  facet_wrap(variable ~ .,scales="free") +
  theme_bw() +
  theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) #,legend.position="bottom",legend.justification='right'
vars_plot_G + geom_histogram(colour='black',alpha = 0.2,position="identity")+
  labs(title = "Histogram plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming",
       subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))
vars_plot_G + geom_density(alpha = 0.2) +
  labs(title = "Density plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming",
       subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))


###plot divided by variable and Role for Trophallaxis
vars_plot_T <- ggplot(summary_data_T, aes(value)) +
  facet_wrap(variable ~ .,scales="free") +
  theme_bw() +
  theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) 
vars_plot_T + geom_histogram(colour='black',alpha = 0.2,position="identity") +
  labs(title = "Histogram for movement variables calculated from coordinates \n for interacting ants during Trophallaxis",
       subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))
vars_plot_T + geom_density(alpha = 0.2) +
  labs(title = "Density plot for movement variables calculated from coordinates \n for interacting ants during Trophallaxis",
       subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))

###plot divided by variable and Role for Trophallaxis compared to FR and CR
vars_plot_T_FR_CR <- ggplot(summary_data_T_FR_CR, aes(value, color = forcats::fct_inorder(BEH))) +
  facet_wrap(variable ~ .,scales="free") +
  theme_bw() +
  theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm'))
#vars_plot_T_FR_CR + geom_histogram(colour='black',alpha = 0.2,position="identity") +
#  labs(title = "Histogram for movement variables calculated from coordinates \n for interacting ants during Trophallaxis and for non-interacting ants during Front Rest and Cross Rest")
vars_plot_T_FR_CR + geom_density(alpha = 0.2,aes(linetype=forcats::fct_inorder(BEH))) +
  labs(title = "Density plot for movement variables calculated from coordinates \n for interacting ants during Trophallaxis and for non-interacting ants during Front Rest and Cross Rest",
       subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))


#INTERACTION_DATA
#-------------------------------------------------------------
#use same plots as before but for interactions
#THIS CASE WILL BE JUST CUTTING FOR 1 INTERACTION
interaction_data_19 <- interaction_MANUAL[which(interaction_MANUAL$ROW == 19),]


#reshape data for plotting. Split by REC and ACT
interaction_data_19$frame <- seq.int(nrow(interaction_data_19)) 
interaction_data_ACT <-interaction_data_19 %>% dplyr::select(contains(c("ROW","BEH", "Act_Name","traj_BOTH.ACT","angle_diff","straightline","frame"), ignore.case = TRUE))
interaction_data_REC <- interaction_data_19 %>% dplyr::select(contains(c("ROW","BEH", "Rec_Name","traj_BOTH.REC","angle_diff","straightline","frame"), ignore.case = TRUE))
#Rename columns to make them match and bind+ melt columns
int_data_ACT  <- interaction_data_ACT %>% rename_with(~str_remove(., c("_ACT|Act_|traj_BOTH.ACT.|traj_BOTH."))); int_data_ACT$Role <- "ACT"
int_data_REC  <- interaction_data_REC %>% rename_with(~str_remove(., c("_REC|Rec_|traj_BOTH.REC.|traj_BOTH."))); int_data_REC$Role <- "REC"
int_data_bind <- rbind(int_data_ACT,int_data_REC)
#int_data_bind$frame <- seq.int(nrow(int_data_bind)) 
int_data_long <- melt(int_data_bind,id.vars=c("BEH","Role","Name","frame")) #explanation on the warning message https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message

#subset data by BEH for plotting and #remove ROW
interaction_data_G <- int_data_long[which(int_data_long$BEH == "G" & !int_data_long$variable == "ROW" ),]
interaction_data_G <- interaction_data_G[which(!is.na(interaction_data_G$value)),]
#cut the receiver as it will be the same for common parameters and add row number as time 
interaction_data_G <- interaction_data_G[which(!interaction_data_G$Role == "REC"),]
interaction_data_G <- interaction_data_G[which(!interaction_data_G$variable == "x"),]
interaction_data_G <- interaction_data_G[which(!interaction_data_G$variable == "y"),]
interaction_data_G <- interaction_data_G[which(!interaction_data_G$variable == "angle"),]
interaction_data_G <- interaction_data_G[which(!is.na(interaction_data_G$value)),]


#set up plot

pdf(file=paste(DATADIR,"TEST_TEST_2.pdf", sep = ""), width=1.7, height=3.1)
par(mfrow=c(2,6), family="serif" , mar = c(0.1, 0.1, 2.2, 0))
#par(mfrow=c(2,3), family="serif" , mar = c(0.1, 0.1, 2.2, 0))


ggplot(data=interaction_data_G,
       aes(x=frame, y=value, colour=variable)) + facet_wrap(variable ~ .,scales="free",ncol = 1) + theme(legend.position = "none") +
  geom_line(size=0.8)

dev.off()

###plot divided by variable and Role for Grooming
# vars_plot_G <- ggplot(interaction_data_G, aes(frame, fill = BEH)) +
#   facet_wrap(variable ~ .,scales="free") +
#   theme_bw() +
#   theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) #,legend.position="bottom",legend.justification='right'
# # vars_plot_G + geom_histogram(colour='black',alpha = 0.2,position="identity",bins = 500)+
# #   labs(title = "Histogram plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming"#,
# #        #subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO)
# #        )
# vars_plot_G + geom_density(alpha = 0.2) +
#   labs(title = "Density plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming",
#        #subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO)
#        )


#------------------------------------------------------------

interaction_data_circ <- circular(interaction_MANUAL$traj_BOTH.angle_diff)

##angular_differences plot per interaction
interaction_data_G <- subset(interaction_MANUAL[!is.na(interaction_MANUAL$traj_BOTH.angle_diff),], BEH == "G")
for (row in unique(interaction_data_G$ROW)) {
  single_interaction <-subset(interaction_data_G,ROW == row)
  p <- plot.circular(single_interaction$traj_BOTH.angle_diff, pch = 16, cex = 0.8, stack=TRUE, bins=100, tcl.text	 = 0.5, tol= 0, shrink = 2.5, sep= 0.04, #xlim=c(-1,1), ylim=c(-2.5,2.5), 
                     main = paste("Behaviour: G, \n", "interaction N: ",row, sep=""), sub=NULL) #REPLICATE, ", ", PERIOD, ", \n", "behaviour:", BEH,
  arrows.circular(mean.circular(single_interaction$traj_BOTH.angle_diff))
  arrows.circular(0, col = "red")
  }

interaction_data_T <- subset(interaction_MANUAL[!is.na(interaction_MANUAL$traj_BOTH.angle_diff),], BEH == "T")
for (row in unique(interaction_data_T$ROW)) {
  single_interaction <-subset(interaction_data_T,ROW == row)
  p <- plot.circular(single_interaction$traj_BOTH.angle_diff, pch = 16, cex = 0.8, stack=TRUE, bins=100, tcl.text	 = 0.5, tol= 0, shrink = 2.5, sep= 0.04, #xlim=c(-1,1), ylim=c(-2.5,2.5), 
                     main = paste("Behaviour: T, \n", "interaction N: ",row, sep=""), sub=NULL)
  arrows.circular(mean.circular(single_interaction$traj_BOTH.angle_diff))
  arrows.circular(0, col = "red")
  }


#means
## calculate circular mean angles for each interaction Row
#interaction_data_mean <- aggregate(traj_BOTH.angle_diff ~ ROW + BEH + REPLICATE + PERIOD, FUN=mean, na.action=na.omit, interaction_MANUAL)
interaction_data_circ_mean <- aggregate(circular(traj_BOTH.angle_diff) ~ ROW + BEH + REPLICATE + PERIOD, FUN=mean, na.action=na.omit, interaction_MANUAL)


for (beh in unique(interaction_data_circ_mean$BEH)) {
  single_interaction <-subset(interaction_data_circ_mean,BEH == beh)
  plot.circular(single_interaction$`circular(traj_BOTH.angle_diff)`, pch = 16, cex = 0.8,stack=TRUE, bins=100, tcl.text	 = 0.5, tol= 0, shrink = 2.5, sep= 0.08, xlim=c(-0.4,0.4), ylim=c(-0.4,0.4), 
                     main = paste("Mean interaction angle for ", beh, "\n", "Tot N interactions: ",NROW(subset(interaction_data_circ_mean,BEH == beh)), sep=""), sub=NULL)
  #densityline <- density(single_interaction$`circular(traj_BOTH.angle_diff)`, bw=30)
  #lines(densityline, col=2)
  arrows.circular(mean(single_interaction$`circular(traj_BOTH.angle_diff)`))
  arrows.circular(0, col = "red")
}


#close the pdf
dev.off()


#polar coordinates plot alternative
# p <- ggplot(interaction_data_G, aes(x=as.factor(ROW), y=traj_BOTH.angle_diff)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
# geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
# # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
# ylim(-10,300) +
# # Custom the theme: no axis title and no cartesian grid
# theme_minimal() +
# # theme(
# #   axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank()# , plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
# # ) +
# labs(title = paste(REPLICATE, ", ", PERIOD, ", ", BEH,", ", sep="")) +
# coord_polar(start = 0) # This makes the coordinate polar instead of cartesian.
# print(p)






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
  
  oriented_metadata <- NULL
  capsule_list <- list()
  #for (myrmidon_file in data_list){
  experiment_name <- unlist(strsplit(MyrmidonCapsuleFile,split="/"))[length(unlist(strsplit(MyrmidonCapsuleFile,split="/")))]
  oriented_data <- fmExperimentOpen(MyrmidonCapsuleFile) #this step is already performed at the beginning
  oriented_ants <- oriented_data$ants
  capsule_names <- oriented_data$antShapeTypeNames
  for (ant in oriented_ants){
    ###extract ant length and capsules
    #ant_length_px <- mean(fmQueryComputeMeasurementFor(oriented_data,antID=ant$ID)$length_px)
    capsules      <- ant$capsules
    for (caps in 1:length(capsules)){
      capsule_name  <- capsule_names[[capsules[[caps]]$type]]
      capsule_coord <- capsules[[caps]]$capsule
      capsule_info <- data.frame(experiment = experiment_name,
                                 antID      = ant$ID,
                                 c1x = capsule_coord$c1[1],
                                 c1y = capsule_coord$c1[2],
                                 c2x = capsule_coord$c2[1],
                                 c2y = capsule_coord$c2[2],
                                 r1  = capsule_coord$r1[1],
                                 r2   = capsule_coord$r2[1]
      )
      
      if (!capsule_name %in%names(capsule_list)){ ###if this is the first time we encounter this capsule, add it to capsule list...
        capsule_list <- c(capsule_list,list(capsule_info)) 
        if(length(names(capsule_list))==0){
          names(capsule_list) <- capsule_name
        }else{
          names(capsule_list)[length(capsule_list)] <- capsule_name
        }
      }else{###otherwise, add a line to the existing dataframe within capsule_list
        capsule_list[[capsule_name]] <- rbind(capsule_list[[capsule_name]] , capsule_info)
      }
    }
  }
  
  
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
