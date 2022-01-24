##########################################################################################
############## THIS VERSION IS FORT 0.8.1 COMPATIBLE #####################################
##########################################################################################
#this should be the version of the script maintained for long term use.

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



#"https://formicidae-tracker.github.io/myrmidon/docs/latest/api/index.html"

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

# ## TEMPORARY - Hard code plotrix function
# std.error <- function (x, na.rm) 
# {
#   vn <- function(x) return(sum(!is.na(x)))
#   dimx <- dim(x)
#   if (is.null(dimx)) {
#     stderr <- sd(x, na.rm = TRUE)
#     vnx <- vn(x)
#   }
#   else {
#     if (is.data.frame(x)) {
#       vnx <- unlist(sapply(x, vn))
#       stderr <- unlist(sapply(x, sd, na.rm = TRUE))
#     }
#     else {
#       vnx <- unlist(apply(x, 2, vn))
#       stderr <- unlist(apply(x, 2, sd, na.rm = TRUE))
#     }
#   }
#   return(stderr/sqrt(vnx))
# }

## PARAMETERS
USER            <- "Adriano"
Xmin <- 2000
Xmax <- 7900
Ymin <- 0000
Ymax <- 5500

N_DECIMALS <- 3 ## number of decimals to preserve when rounding to match between interactions & collisions
FUZZY_MATCH <- TRUE  ## fuzzy matching between data frames


max_gap         <- fmHour(24*365)   ## important parameter to set! Set a maximumGap high enough that no cutting will happen. For example, set it at an entire year: fmHour(24*365)
desired_step_length_time    <- 0.125 ###in seconds, the desired step length for the analysis

####### navigate to folder containing myrmidon file
if (USER=="Adriano") {WORKDIR <- "/home/cf19810/Dropbox/Ants_behaviour_analysis"}
if (USER=="Tom")     {WORKDIR <- "/media/tom/MSATA/Dropbox/Ants_behaviour_analysis"}

DATADIR <- paste(WORKDIR,"Data",sep="/")
SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")

## perform analyses on annotation data
source(paste(SCRIPTDIR,"Annotations_analysis.R",sep="/"))


#######################
interaction_data_ALL <- NULL
summary_data         <- NULL

for (REPLICATE in c("R3SP"))#,"R9SP")) 
  {
  ###############################################################################
  ###### OPEN EXPERIMENT INFORMATION ############################################
  ###############################################################################
  
  ## locate the ant info file for REPLICATE
  MyrmidonCapsuleFile <- list.files(path=DATADIR, pattern=REPLICATE, full.names=T)
  MyrmidonCapsuleFile <- MyrmidonCapsuleFile[grepl("Capsule_Zones_defined.myrmidon",MyrmidonCapsuleFile)]
  e <- fmExperimentOpen(MyrmidonCapsuleFile)
  fmQueryGetDataInformations(e)
  ###tag statistics
  tag_stats <- fmQueryComputeTagStatistics(e)
  

  for (PERIOD in c("pre","post"))
    {
    print(paste("Replicate",REPLICATE, PERIOD))
    interaction_data_REP_PER <- NULL   ############Prepare empty within-replicate/period data object

    
    ###############################################################################
    ###### READING TRAJECTORIES ###################################################
    ###############################################################################
      
    time_start <- fmTimeCreate(min(annotations$T_start_[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD])) ###experiment start time
    time_stop  <- fmTimeCreate(min(annotations$T_start_[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD]) + 40*60  ) ###experiment start time ####arbitrary time in the correct format + (N mins * N seconds)
  
    ###############################################################################
    ###### READING COLLISIONS #####################################################
    ###############################################################################
    
    collisions <- fmQueryCollideFrames(e, start=time_start, end=time_stop)   ###collisions are for each frame, the list of ants whose shapes intersect one another. Normally not used

    
    #   
    # StartTime      <- min(annotations$T_start_[annotations$treatment_rep==REPLICATE & annotations$period==PERIOD])
    # time_start_ISO <- parse_iso_8601(StartTime)
    # time_start     <- fmTimeCPtrFromAnySEXP(time_start_ISO)
    # #time_stop      <- max(annotations$T_start[annotations$treatment_rep==REPLICATE & annotations$period=="post"])
    # time_stop      <- fmTimeCPtrFromAnySEXP(time_start_ISO + (40*60) ) ####arbitrary time in the correct format + (N mins * N seconds)
    # 
    #keep end time as date and not environment for plotting purposes
    #time_stop_ISO <- time_start_ISO + (40*60)
  
    #COMPUTE TRAJECTORIES  
    positions <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = TRUE) #set true to obtain the zone of the ant
    #print(head(postions))
    
    ######### no subsetting - take the entire tracking period
    # positions <- fmQueryComputeAntTrajectories(e,start = e$getDataInformations()$details$tdd.start[1],end = e$getDataInformations()$details$tdd.end[3],maximumGap = max_gap,computeZones = FALSE) #set true to obtain the zone of the ant
    # ?fmQueryComputeAntTrajectories

    ###WARNING!!!!!
    ###in this particular example the values in antID are consecutive (1 to 22), and so positions$trajectories[[1]] will correspond to the trajectory of ant 1
    ### but imagine you had to delete an ant in fort-studio and the antID list "jumped" from 10 to 12 (so it would be 1,...10 then 12...23)
    ### then the 22nd object of trajectory list would not be the trajectory of ant 22 but of ant 23. That is not fool proof, but very risky!
    ##so what I suggest you do immediately after computing your positions object:
    positions$trajectories_summary$antID_str <- paste("ant_",positions$trajectories_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
    names(positions$trajectories)       <- positions$trajectories_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
    ##and now we can extract each ant trajectory using the character ID stored in antID_str - it's much safer!
    
    # #unlisting the trajectories 
    # library(data.table)
    # trajectories_unlist <- as.data.frame(rbindlist(positions$trajectories,use.names=TRUE,idcol=TRUE))
    # write.csv(trajectories_unlist, file=paste(e$getDataInformations()[["details"]][["tdd.URI"]],'trajectories_30min.csv'),row.names = FALSE) #specify name of the block better!
    
    trajectories_summary <- positions$trajectories_summary
    
    #generate reproducible example
    #dput(positions, file = "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_positions.txt")
    #positions_dget <- dget("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_positions.txt") # load file created with dput 
  
    # because time gaps are not fixed (not exactly 0.125) every x frames (~800/1000) there is a time mismatch at the third decimal place
    #[787]  98.250995  98.375997  98.500998  98.625999  *98.751000*  98.876002
    
    ## Add ant_x names and times to the positions to convert from time since the start of the experiment, to UNIX time
    for (A in positions$trajectories_summary$antID)
      {
      AntID <- paste("ant_",A,sep="") 
      First_Obs_Time <- as.POSIXct(positions$trajectories_summary$start [which(positions$trajectories_summary$antID_str==AntID)], tz="GMT",origin = "1970-01-01 00:00:00") ## find the first time after the user defined time_start_ISO that this ant was seen
      print(paste("Adding first obs time", First_Obs_Time, "to the time-zeroed trajectory of ant", AntID))
      positions$trajectories[[AntID]] $UNIX_time <- as.POSIXct(positions$trajectories[[AntID]]$time, tz="GMT",origin = "1970-01-01 00:00:00")  + First_Obs_Time ##convert back to UNIX time  
      }
   
    #check that milliseconds are preserved after adding the ant's time_start_ISO #positions$trajectories[["ant_27"]][[2,"UNIX_time"]]-positions$trajectories[["ant_27"]][[1,"UNIX_time"]]
    
    ###############################################################################
    ###### extract ant trajectories from hand-annotated data ######################
    ###############################################################################
    ####First let's extract ant's trajectories
    #set plots parameters (for plotting coords)
    #pdf(file=paste(DATADIR,"Interactions_Coordinates_plots.pdf", sep = ""), width=6, height=4.5)
    par(mfrow=c(2,3), mai=c(0.3,0.3,0.4,0.1), mgp=c(1.3,0.3,0), family="serif", tcl=-0.2)
   
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
     
        ENC_TIME_start <- annot_BEH$T_start_UNIX[ROW]
        ENC_TIME_stop  <- annot_BEH$T_stop_UNIX[ROW]
        
        Act_Name <- paste("ant",ACT,sep="_")
        Rec_Name <- paste("ant",REC,sep="_")
        print(paste("Behaviour:",BEH,"number",ROW,"Actor:",Act_Name,"Receiver:",Rec_Name))
        
        ## extract the trajectory for ACT
        traj_ACT <-  positions$trajectories[[Act_Name]]
        traj_REC <-  positions$trajectories[[Rec_Name]]
        
        ##  convert POSIX format to raw secs since 1.1.1970; brute force for match.closest with collisions
        traj_ACT$UNIX_secs <- as.numeric(traj_ACT$UNIX_time)  ## yes, it looks like the milisecs are gone, but they are there; traj_ACT$UNIX_secs
        traj_REC$UNIX_secs <- as.numeric(traj_REC$UNIX_time) 
        ## now round these decimal seconds to 3 d.p.
        traj_ACT$UNIX_secs <- round(traj_ACT$UNIX_secs,N_DECIMALS) ## eliminates the 'noise' below the 3rd d.p., and leaves each frame existing just once, see: table(traj_ACT$UNIX_secs)
        traj_REC$UNIX_secs <- round(traj_REC$UNIX_secs,N_DECIMALS)
        
        ## remove 'time' column as it is confusing - it's not a common time
        traj_ACT$time <- NULL
        traj_REC$time <- NULL
          
        # ## Plot trajectories of both actor & receiver, show on the same panel
        # plot  (y ~ x, traj_ACT, pch=".", col=rgb(0,0,1,0.3,1), main=Title, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
        # points(y ~ x, traj_REC, pch=".", col=rgb(1,0,0,0.3,1))
         
        ## subset the trajectories of both actor & receiver using the start & end times
        traj_ACT <- traj_ACT [ which(traj_ACT$UNIX_time >= ENC_TIME_start & traj_ACT$UNIX_time <= ENC_TIME_stop),]
        traj_REC <- traj_REC [ which(traj_REC$UNIX_time >= ENC_TIME_start & traj_REC$UNIX_time <= ENC_TIME_stop),]
        
        ## For some reason, the 4th decimal place of the UNIX times (thousandths of a sec) for two individuals seen in a given frame are not identical; cbind(format(traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS4"), format(traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS4"))
        ## so eliminate the 4th-nth decimal place ; retains accuract to a thousandth-of-a-second
        # traj_ACT$UNIX_time <- format(traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
        # traj_REC$UNIX_time <- format(traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
        # print(table(substr(x =  format( traj_ACT$UNIX_time, "%OS6"), start = 7, stop = 7)))
        # print(table(substr(x =  format( traj_REC$UNIX_time, "%OS6"), start = 7, stop = 7)))
        
      #   
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
        # print(paste("Behaviour:",BEH,"ACT:",Act_Name,"REC:",Rec_Name, "annot_start", ENC_TIME_start, "traj_start", min(traj_ACT$UNIX_time,na.rm = TRUE)))
        # 
        # # ## Plot trajectories of both actor & receiver, show on the same panel
        Title <- paste(REPLICATE, ", ", PERIOD, ", ", BEH, ROW,", ", "Act:",ACT, ", ", "Rec:",REC, "\n", ENC_TIME_start, "-", ENC_TIME_stop, sep="")
        plot   (y ~ x, traj_ACT, type="l", lwd=6, col="blue4",asp=1, main=Title,cex.main=0.9 ,xlim=c(min(traj_ACT$x,traj_REC$x),max(traj_ACT$x,traj_REC$x)),ylim=c(min(traj_ACT$y,traj_REC$y),max(traj_ACT$y,traj_REC$y))) #, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
        points (y ~ x, traj_REC, type="l", lwd=6,  col="red4",asp=1)
  
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
        trajectory_ACT <- TrajFromCoords(data.frame(x=traj_ACT$x,y=traj_ACT$y, UNIX_time=as.numeric(traj_ACT$UNIX_time)), spatialUnits = "px",timeUnits="s")
        trajectory_REC <- TrajFromCoords(data.frame(x=traj_REC$x,y=traj_REC$y, UNIX_time=as.numeric(traj_REC$UNIX_time)), spatialUnits = "px",timeUnits="s")

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
        traj_BOTH <- merge(traj_ACT, traj_REC,all=T, by=c("UNIX_secs"))  ## 21 Jan 2022: Changed to sue raw unix seconds to allow simpler matching below (posix formats cause problems with the matching...)
        
        
        ## merge drops the column attributes - add back the POSIX time format 
        # traj_BOTH$UNIX_time <- as.POSIXct(traj_BOTH$UNIX_time, tz="GMT",origin="1970-01-01 00:00:00")
        
        ## measure the length *in seconds* of the interaction between ACT & REC
        # interaction_length_secs <- as.numeric(difftime ( max(traj_BOTH$UNIX_time, na.rm=T), min(traj_BOTH$UNIX_time, na.rm=T), units="secs"))
        interaction_length_secs <-  max(traj_BOTH$UNIX_secs, na.rm=T) - min(traj_BOTH$UNIX_secs, na.rm=T)  ## 21 Jan 2022
        
        ## CAREFUL - this is wrong, as traj_BOTH does NOT CONTAIN blank rows for frames whenneither Act or REC were seen - need to use interaction_length_secs
        # prop_time_undetected_ACT <- (sum(is.na(traj_BOTH$ACT.x))) / length(traj_BOTH$ACT.x)
        # prop_time_undetected_REC <- (sum(is.na(traj_BOTH$REC.x))) / length(traj_BOTH$REC.x)
        prop_time_undetected_ACT <- (sum(is.na(traj_BOTH$ACT.x)) / 8) / interaction_length_secs  ## the prop of the interaction in which ACT was seen 
        prop_time_undetected_REC <- (sum(is.na(traj_BOTH$REC.x)) / 8)  / interaction_length_secs ## UNITS are secs / secs --> proportion; MUST BE STRICTLY < 1 !!
       
        # FRAME BY FRAME PARAMETERS
        
        #  distance walked
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
        # 
  
        # #SUMMARY_DATA is a summary, INTERACTION_DATA reports interactions frame by frame
        # summary_data <- rbind(summary_data,
        #                       data.frame(ROW=ROW, BEH=BEH, Act_Name=Act_Name, Rec_Name=Rec_Name, PERIOD=PERIOD,
        #                                  StDev_angle_ACT=StDev_angle_ACT, StDev_angle_REC=StDev_angle_REC,
        #                                  angle_mean_ACT=angle_mean_ACT, angle_mean_REC=angle_mean_REC,
        #                                  #mean_delta_angles_ACT=mean_delta_angles_ACT, mean_delta_angles_REC=mean_delta_angles_REC,
        #                                  moved_distance_px_ACT=moved_distance_px_ACT, moved_distance_px_REC=moved_distance_px_REC,
        #                                  mean_accel_pxpersec2_ACT=mean_accel_pxpersec2_ACT, mean_accel_pxpersec2_REC=mean_accel_pxpersec2_REC,
        #                                  #median_accel_ACT=median_accel_ACT,median_accel_REC=median_accel_REC,
        #                                  rmsd_px_ACT=rmsd_px_ACT,rmsd_px_REC=rmsd_px_REC,
        #                                  interaction_length_frames=interaction_length_frames,
        #                                  prop_time_undetected_ACT=prop_time_undetected_ACT, prop_time_undetected_REC=prop_time_undetected_REC,
        #                                  strghtline_dist_px=strghtline_dist_px,
        #                                  #when adding a new variable, it must be included in the reshape rule for data plotting
        #                                  stringsAsFactors = F))
        # 
        #NO UNIX TIME BUT TIME SINCE START OF INTERACTION
        #interaction AREA= NEST, FORAGING delete x y coords
        interaction_data_ROW <- data.frame(REPLICATE=REPLICATE,PERIOD=PERIOD,ROW=ROW,BEH=BEH,Act_Name=Act_Name,Rec_Name=Rec_Name,
                                           traj_BOTH,
                                           stringsAsFactors = F)
        ## stack -by the end of the loop, this should just contain all events for replicate X, period Y
        interaction_data_REP_PER <- rbind(interaction_data_REP_PER, interaction_data_ROW)
        
        
        
       
        
        }##ROW
      }##BEH
    
    ## stack
    interaction_data_ALL <- rbind(interaction_data_ALL, interaction_data_REP_PER)
    
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
    collisions_frames$UNIX_secs <- round(as.numeric(collisions_frames$time),N_DECIMALS)
    
    ###let's view the first few lines of collisions
    head(collisions_frames) 
    head(collisions_coll) ###columns giving the ID of the two colliding ants, the zone where they are, the types of collisions (all types), and the frames_row_index referring to which frame that collision was detected in (matches the list indices in collisions_positions)
    head(collisions_positions)
    ?fmQueryCollideFrames
    #use rownames of collision_frames to filter frames of collisions according to collisions-coll (which acts as a Look Up Table)
    #add frames_row_index
    collisions_frames$frames_row_index <- rownames(collisions_frames)
    
    #merge collisions_frame and coll_coll based on index
    collisions_merge <- merge(collisions_coll, collisions_frames[,-match(c("height","width"),colnames(collisions_frames))], by="frames_row_index")
  
    
    #check that to the same frame corresponds the same time
    # format(collisions_merge$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # format(collisions_frames$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # format(    traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # format(    traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # 
    # format(    traj_BOTH$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # format(interaction_data_ROW$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
    # 
    # table(substr(x =  format( traj_BOTH$UNIX_time, "%OS6"), start = 7, stop = 7))
    
    
    ##if(collisions_coll$frames_row_index[125] == collisions_merge$frames_row_index[125]) {print("TRUE")}
    
    # add new column to interaction data when conditions are met (eg. time + ant ID)
    #RENAME COLUMN TRAJBOTHREC IN time 
    #merge with equal time and any order of the couple Act_Name-Rec_Name=ant1_str-ant2_str (any order). HOW?
    nrow(collisions_merge)
    str(collisions_merge)
    nrow(interaction_data_REP_PER)
    str(interaction_data_REP_PER)
    
    #create new variable by pasting ant numbers "low,high" for collisions_merge
    # collisions_merge$pair <- apply(collisions_merge[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
    collisions_merge $pair <- paste(collisions_merge$ant1, collisions_merge$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...
    
    #create new variable by pasting ant numbers "low,high" for interaction_data_REP_PER
    interaction_data_REP_PER$ant1 <- gsub("ant_","", interaction_data_REP_PER$Act_Name)
    interaction_data_REP_PER$ant2 <- gsub("ant_","", interaction_data_REP_PER$Rec_Name)
    ## the ant 1 & ant 2 labels are not strictly ascending, so sort them so ant 1 alwasy < ant 2
    interaction_data_REP_PER$pair <- apply(interaction_data_REP_PER[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })

    # check that the time formats are the same
    attributes(interaction_data_REP_PER$UNIX_time[1])
    attributes(collisions_merge$UNIX_time[1])
    
    ## check that the pair-time combinations in interaction_data_REP_PER are in collisions_merge
    table(  paste(interaction_data_REP_PER$pair) %in% 
              paste(collisions_merge$pair) )
    table(  paste(interaction_data_REP_PER$UNIX_secs) %in% 
              paste(collisions_merge$UNIX_secs) )
    table(  paste(interaction_data_REP_PER$pair,interaction_data_REP_PER$UNIX_secs) %in% 
              paste(collisions_merge$pair,collisions_merge$UNIX_secs) )
    ## get a list of the missing pair-time combinations:
    paste(interaction_data_REP_PER$pair,interaction_data_REP_PER$UNIX_secs) [!paste(interaction_data_REP_PER$pair,interaction_data_REP_PER$UNIX_secs) %in% paste(collisions_merge$pair,collisions_merge$UNIX_secs) ]
    
    #join dataframes to have collisions on interaction_data_REP_PER  
    #FUNCTION: WHEN pair=pair AND time=UNIX_time -> ASSIGN collisions_merge$types TO interaction_data_REP_PER
    if (FUZZY_MATCH==FALSE)
      {
      interaction_data_REP_PER <- plyr::join(x = interaction_data_REP_PER, 
                                             y = collisions_merge[,c("UNIX_secs","pair","types")], 
                                             by=c("UNIX_secs","pair"), type="left")
      }
    ## alternative: fuzzy match
    if (FUZZY_MATCH==TRUE)
      {
      I_in_C_rows <- match.closest(x = interaction_data_REP_PER$UNIX_secs,
                                   table = collisions_merge$UNIX_secs, tolerance = 0.125)
      ## copy tpes across
      interaction_data_REP_PER$types  <- collisions_merge$types [I_in_C_rows]
      }
    
    ## PROBLEM: *MANY* rows in interactions aren't present in collisions; plot the spatial interactions to see whether this makes sense or not:
    par(mfrow=c(3,4), mai=c(0.3,0.3,0.4,0.1))
    for(BH in unique(interaction_data_REP_PER$BEH))
    {
      for (RW in unique(interaction_data_REP_PER$ROW))
      {
        IntPair <- interaction_data_REP_PER[which(interaction_data_REP_PER$BEH==BH & interaction_data_REP_PER$ROW==RW),]
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
    interaction_data_REP_PER$types_present_absent <- "absent"
    interaction_data_REP_PER$types_present_absent[which(!is.na(interaction_data_REP_PER$types))] <- "present"
    boxplot(straightline_dist_px ~ types_present_absent, interaction_data_REP_PER)
    ### ... so rows in interactions that don't match collisions is not  a function of distance ...
    
    

    #cut collisions for the specific G 1 case

    ROW2 <- 1 #first behaviour

    interactio_data_SELECTED <- interaction_data_ALL[which(interaction_data_ALL$ROW==ROW2  & interaction_data_ALL$BEH=="G" & interaction_data_ALL$PERIOD=="post"),]

    COLL_TIME_start <- min(interactio_data_SELECTED$UNIX_secs)
    COLL_TIME_stop  <- max(interactio_data_SELECTED$UNIX_secs)
    PAIR <- "1_5"

    ## subset the collisions using the start & end times
    collisions_merge_SELECTED <- collisions_merge[ which(collisions_merge$UNIX_secs >= COLL_TIME_start & collisions_merge$UNIX_secs <= COLL_TIME_stop & collisions_merge$pair==PAIR),]
    min(collisions_merge_SELECTED$UNIX_secs)
    max(collisions_merge_SELECTED$UNIX_secs)

    format( collisions_merge_SELECTED$UNIX_secs, nsmall=5)
    format( interactio_data_SELECTED$UNIX_secs, nsmall=5)
    

    tail(format(collisions_merge_SELECTED$UNIX_secs, "%Y-%m-%d %H:%M:%OS6"),10)
    tail(format(interactio_data_SELECTED$UNIX_secs, "%Y-%m-%d %H:%M:%OS6"),10)

    ########################################################################################################################
    ########################################################################################################################
    
    ## generate summary data
    for (variable in names(summary_data)[!names(summary_data)%in%c("BEH","Act_Name","Rec_Name","PERIOD")])
      {
      summary_data[,"variable"] <- summary_data[,variable] #boxplot(variable~BEH,ylab=variable,data=summary_data)
      }##summary_data
    
    }##PERIOD
    #clear cache before opening following exp. TO BE TESTED 
    # rm(list=(c("e")))
    # gc() # clear cache
  }##REPLICATE
# dev.off() # save to pds!


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
unique(interaction_data_ALL$PERIOD)
min(interaction_data_ALL$UNIX_time)
max(interaction_data_ALL$UNIX_time)

plot(interaction_data_ALL$ACT.y ~ interaction_data_ALL$ACT.x, asp=1, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
###############################################################################
###### PARAMETERS PLOTS #######################################################
###############################################################################

#SUMMARY_DATA
#descriptive analysis
str(summary_data)

#reshape data for plotting. Split by REC and ACT
summary_data_ACT <-summary_data %>% dplyr::select(contains(c("BEH", "ACT","interaction_length","strghtline"), ignore.case = TRUE))
summary_data_REC <- summary_data %>% dplyr::select(contains(c("BEH", "REC","interaction_length","strghtline"), ignore.case = TRUE))
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
       subtitle = paste( "Periods:",unique(interaction_data_ALL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))
vars_plot_G + geom_density(alpha = 0.2) +
  labs(title = "Density plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming",
       subtitle = paste( "Periods:",unique(interaction_data_ALL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))


###plot divided by variable and Role for Trophallaxis
vars_plot_T <- ggplot(summary_data_T, aes(value)) +
  facet_wrap(variable ~ .,scales="free") +
  theme_bw() +
  theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) 
vars_plot_T + geom_histogram(colour='black',alpha = 0.2,position="identity") +
  labs(title = "Histogram for movement variables calculated from coordinates \n for interacting ants during Trophallaxis",
       subtitle = paste( "Periods:",unique(interaction_data_ALL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))
vars_plot_T + geom_density(alpha = 0.2) +
  labs(title = "Density plot for movement variables calculated from coordinates \n for interacting ants during Trophallaxis",
       subtitle = paste( "Periods:",unique(interaction_data_ALL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))

###plot divided by variable and Role for Trophallaxis compared to FR and CR
vars_plot_T_FR_CR <- ggplot(summary_data_T_FR_CR, aes(value, color = forcats::fct_inorder(BEH))) +
  facet_wrap(variable ~ .,scales="free") +
  theme_bw() +
  theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm'))
#vars_plot_T_FR_CR + geom_histogram(colour='black',alpha = 0.2,position="identity") +
#  labs(title = "Histogram for movement variables calculated from coordinates \n for interacting ants during Trophallaxis and for non-interacting ants during Front Rest and Cross Rest")
vars_plot_T_FR_CR + geom_density(alpha = 0.2,aes(linetype=forcats::fct_inorder(BEH))) +
  labs(title = "Density plot for movement variables calculated from coordinates \n for interacting ants during Trophallaxis and for non-interacting ants during Front Rest and Cross Rest",
       subtitle = paste( "Periods:",unique(interaction_data_ALL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))


#INTERACTION_DATA
#-------------------------------------------------------------
#use same plots as before but for interactions
#THIS CASE WILL BE JUST CUTTING FOR 1 INTERACTION
interaction_data_19 <- interaction_data_ALL[which(interaction_data_ALL$ROW == 19),]


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
# #        #subtitle = paste( "Periods:",unique(interaction_data_ALL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO)
# #        )
# vars_plot_G + geom_density(alpha = 0.2) +
#   labs(title = "Density plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming",
#        #subtitle = paste( "Periods:",unique(interaction_data_ALL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO)
#        )


#------------------------------------------------------------

interaction_data_circ <- circular(interaction_data_ALL$traj_BOTH.angle_diff)

##angular_differences plot per interaction
interaction_data_G <- subset(interaction_data_ALL[!is.na(interaction_data_ALL$traj_BOTH.angle_diff),], BEH == "G")
for (row in unique(interaction_data_G$ROW)) {
  single_interaction <-subset(interaction_data_G,ROW == row)
  p <- plot.circular(single_interaction$traj_BOTH.angle_diff, pch = 16, cex = 0.8, stack=TRUE, bins=100, tcl.text	 = 0.5, tol= 0, shrink = 2.5, sep= 0.04, #xlim=c(-1,1), ylim=c(-2.5,2.5), 
                     main = paste("Behaviour: G, \n", "interaction N: ",row, sep=""), sub=NULL) #REPLICATE, ", ", PERIOD, ", \n", "behaviour:", BEH,
  arrows.circular(mean.circular(single_interaction$traj_BOTH.angle_diff))
  arrows.circular(0, col = "red")
  }

interaction_data_T <- subset(interaction_data_ALL[!is.na(interaction_data_ALL$traj_BOTH.angle_diff),], BEH == "T")
for (row in unique(interaction_data_T$ROW)) {
  single_interaction <-subset(interaction_data_T,ROW == row)
  p <- plot.circular(single_interaction$traj_BOTH.angle_diff, pch = 16, cex = 0.8, stack=TRUE, bins=100, tcl.text	 = 0.5, tol= 0, shrink = 2.5, sep= 0.04, #xlim=c(-1,1), ylim=c(-2.5,2.5), 
                     main = paste("Behaviour: T, \n", "interaction N: ",row, sep=""), sub=NULL)
  arrows.circular(mean.circular(single_interaction$traj_BOTH.angle_diff))
  arrows.circular(0, col = "red")
  }


#means
## calculate circular mean angles for each interaction Row
#interaction_data_mean <- aggregate(traj_BOTH.angle_diff ~ ROW + BEH + REPLICATE + PERIOD, FUN=mean, na.action=na.omit, interaction_data_ALL)
interaction_data_circ_mean <- aggregate(circular(traj_BOTH.angle_diff) ~ ROW + BEH + REPLICATE + PERIOD, FUN=mean, na.action=na.omit, interaction_data_ALL)


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
#unlist(strsplit(interactions_all$types[4],","))[2]



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




################# SCRAPS ##########################

# 
# ######################################################
# 
# #define capsules in fortstudio, then look at the data
# capsules  <- e$antShapeTypeNames()
# body_id <- capsules[which(capsules$name=="body"),"typeID"]
# 
# ###############################################################################
# ###### 9. READING INTERACTIONS ################################################
# ###############################################################################
# ###This provides more synthetic and  complete information than collisions, as if the same pair of ants interacts over successive frames, it will be reported as a single interactions with a start and end time
# ### The function to extract interactions is fmQueryComputeAntInteractions
# ###Let's have a look at the arguments needed:
# ?fmQueryComputeAntInteractions
# ######### experiment, start, end
# ######################## As in trajectories
# 
# ######### maximumGap
# ######################## EXTREMELY IMPORTANT PARAMETER - THE CHOICE WILL BE DIFFERENT FROM THAT MADE FOR THE TRAJECTORY QUERY ABOVE!
# ######################## The basic idea is the following:
# ######################## Within a single biological, real interaction, there are likely to be gaps
# ######################## i.e. frame in which one or both of te ants are not detected, and so no collision is detected between these two ants for that frame
# ######################## This is expected because the tracking system is not perfect - or the ants may tilt slightly during the interaction, making their tag undetectable
# ######################## Yet it would be a mistake to say that whenever one or both ants disappear, the interaction stops and a new starts
# ######################## So we need to define an "allowance", i.e. a maximum period of continuous interruption of the interaction for which we will still be happy to say that it is the same interaction that continues
# ######################## If we were to set that at one year, like we did for the trajectories above, then each pair of ant would only ever be as having a single interaction
# ######################## In the past we have used something akin to 10 seconds. If the interruption is longer than this, we consider that a new interaction starts - and a new line will be written 
# 
# ######### reportTrajectories
# ######################## Argument taking either T or F value, determining whether trajectories will be output in parallel to the interactions
# ######################## Useful if we want to know the x-y coordinates of each ant in each interaction
# ######################## At the moment the function only works if this is set as TRUE (in my case I get a segmentation fault otherwise)
# 
# ###so:
interactions_all       <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportFullTrajectories = T)
# 
# ###let's View the object: 
str(interactions_all)
# ###it's a list containing 3 objects: summaryTrajectory, trajectories and interactions
# ###let's extract them for simplicity:
summaryTrajectory_all       <- interactions_all$trajectories_summary  ## same object as positions$trajectories_summary described above - however here it will have many more rows given that the trajectories will be broken whenever there is a 10 second non-detection gap for an ant.
# ## therefore in this case my trick with ant_ID_str won't work - it only works if there's exactly one trajectory as ants in the data
# ## so I would NOT use the output of this function to analyse trajectories
trajectories_all           <- interactions_all$trajectories      ## same object as positions$trajectories described above - but containing more objects for the same reason as stated above
interactions_all            <- interactions_all$interactions       ## data frame containing all interactions
# 
# str(interactions_all)
# ###let's have a look at the interactions object
# head(interactions_all)
# ####dataframe comtaining the following info:
# ######ant1: ID of the first ant in the interaction
# ######ant2: ID of the second ant in the interaction
# ######start: time at which interaction started
# ######end: time at which interaction ended
# ######space: space in which the interaction took place
# ######types: all types of capsule intersection observed in this interaction. Commas separate different types of intersections, and each intersection type lists which capsule of ant1 interesected with which capsule in at2, separated by a "-"
# ###############thus 1-5, 5-1,5-5: means that during this interaction, capsule 1 in ant 1 intersected with capsule 5 in ant 2; capsule 5 in ant 1 intersected with capsule 1 in ant 2, and capsule 5 in ant1 intersected with capsule 5 in ant2
# ######ant1.trajectory.row: gives the index to use within summaryTrajectory and trajectories to extract the appropriate trajectory segment for ant1 during this interaction
# ######ant1.trajectory.start: within trajectory segment trajectoriessummaryTrajectory[ant1.trajectory.row], ant1.trajectory.start gives the row corresponding to the start of this interaction
# ######ant1.trajectory.end: within trajectory segment trajectories[ant1.trajectory.row], ant1.trajectory.end gives the row corresponding to the end of this interaction. BUT THIS SEEMS WRONG IN CURRENT VERSION. NEED TO REPORT BUG
# ######ant2.trajectory.row, ant2.trajectory.start, ant2.trajectory.end: same for ant2
# 
# ####In the example above we have not specified a matcher, i.e. all intersections between all types of capsules have been considered
# ####But in many cases we will be interested in a specific type of interactions. 
# #### For example: we may be interested only in interactions involving the intersection between body shapes
# ###First we need to figure out the index of the shape corresponding to body shape
# capsules  <- e$antShapeTypeNames()
# body_id <- capsules[which(capsules$name=="body"),"typeID"]
# ###Then we will specify a matcher in which we are interested in interactions involving capsule 1 for both anta: fmMatcherInteractionType(body_id,body_id)
# ##for more info, type:
# ?fmMatcherInteractionType
# 
# ####################################################
# #################################CONTINUE FROM HERE!
# 
# nrow(interactions_all)    
# 
# #with this function we do access the trajectories information specifying the interactions parameters
# ###so for example, if we wanted to know the x-y coordinates of ant1 (7) and ant2 (9) on the first frame of the first interaction listed int his table:
# coord_ant1_start <- trajectories_all[[   interactions_all[1,"ant1.trajectory.row"]    ]][interactions_all[1,"ant1.trajectory.start"],c("x","y")]
# coord_ant2_start <- trajectories_all[[   interactions_all[1,"ant2.trajectory.row"]    ]][interactions_all[1,"ant2.trajectory.start"],c("x","y")]
# 
# ###if we wanted to know the mean x-y coordinates of ant1 (7) and ant2 (9) in first interaction listed int his table:
# coord_ant1_mean <- colMeans(trajectories[[   interactions[1,"ant1.trajectory.row"]    ]][interactions[1,"ant1.trajectory.start"]:interactions[1,"ant1.trajectory.end"],c("x","y")])
# coord_ant2_mean <- colMeans(trajectories[[   interactions[1,"ant2.trajectory.row"]    ]][interactions[1,"ant2.trajectory.start"]:interactions[1,"ant2.trajectory.end"],c("x","y")])
# 
# ###and if you wanted to compute this for all interactions, you could either loop over all interactions or write a function and use apply
# ####but that would be a bit more complicated to get to work than the examples I gave above
# ####here I don't use match but I need to be sure that trajectories remained non-shuffled
# mean_coord <- function (x,which_ant,trajectories_all){
#   coord_mean <- colMeans(trajectories_all[[as.numeric(x[paste(which_ant,".trajectory.row",sep="")] )]]  [as.numeric(x[paste(which_ant,".trajectory.start",sep="")]):as.numeric(x[paste(which_ant,".trajectory.end",sep="")]),c("x","y")] )
#   return(data.frame(x=coord_mean["x"],y=coord_mean["y"]))
# }
# mean_coordinates <- unlist(apply(interactions_all, 1, FUN=mean_coord,which_ant="ant1",trajectories_all=trajectories_all))
# interactions_all$mean_x_ant1 <- mean_coordinates[grepl("x",names(mean_coordinates))]
# interactions_all$mean_y_ant1 <- mean_coordinates[grepl("y",names(mean_coordinates))]
# 
# mean_coordinates <- unlist(apply(interactions_all, 1, FUN=mean_coord,which_ant="ant2",trajectories_all=trajectories_all))
# interactions_all$mean_x_ant2 <- mean_coordinates[grepl("x",names(mean_coordinates))]
# interactions_all$mean_y_ant2 <- mean_coordinates[grepl("y",names(mean_coordinates))]
# 
# 
# interactions_all
# 
# plot(interactions_all$mean_x_ant2,interactions_all$mean_y_ant2 )
# plot(interactions_all$mean_x_ant1,interactions_all$mean_y_ant1 )
# 
# #trying to calculate the StDev results in errors but should not bee too hard! 
# #errors: 
# # number 1 - keeping the same structure as for colMeans
# # > StDev_coord <- function (x,which_ant,trajectories_all){
# #   +   coord_StDev <- sd(trajectories_all[[as.numeric(x[paste(which_ant,".trajectory.row",sep="")] )]]  [as.numeric(x[paste(which_ant,".trajectory.start",sep="")]):as.numeric(x[paste(which_ant,".trajectory.end",sep="")]),c("x","y")] )
# #   +   return(data.frame(x=coord_StDev["x"],y=coord_StDev["y"]))
# #   + }
# # > StDev_coordinates <- unlist(apply(interactions_all, 1, FUN=StDev_coord,which_ant="ant1",trajectories_all=trajectories_all))
# # Error in is.data.frame(x) : 
# #   'list' object cannot be coerced to type 'double'
# 
# # number 2 - removing as.numeric from StDev_coord function
# #it results into all NAs
# 
# 
# StDev_coord <- function (x,which_ant,trajectories_all){
#   coord_StDev <- sd(trajectories_all[[as.numeric(x[paste(which_ant,".trajectory.row",sep="")] )]]  [as.numeric(x[paste(which_ant,".trajectory.start",sep="")]):as.numeric(x[paste(which_ant,".trajectory.end",sep="")]),c("x","y")])
#   return(data.frame(x=coord_StDev["x"],y=coord_StDev["y"]))
# }
# StDev_coordinates <- unlist(apply(interactions_all, 1, FUN=StDev_coord,which_ant="ant1",trajectories_all=trajectories_all))
# 
# interactions_all$mean_x_ant1 <- StDev_coordinates[grepl("x",names(StDev_coordinates))]
# interactions_all$mean_y_ant1 <- StDev_coordinates[grepl("y",names(StDev_coordinates))]
# 
# StDev_coordinates <- unlist(apply(interactions_all, 1, FUN=StDev_coord,which_ant="ant2",trajectories_all=trajectories_all))
# interactions_all$mean_x_ant2 <- StDev_coordinates[grepl("x",names(StDev_coordinates))]
# interactions_all$mean_y_ant2 <- StDev_coordinates[grepl("y",names(StDev_coordinates))]
# 
# 
# 
# 
# 
# 
# #extract first interaction with list of coordinates
# 
# 
# ###so for example, if we wanted to know the x-y coordinates of ant1 (7) and ant2 (9) on the first frame of the first interaction listed int his table:
# coord_ant1_start <- trajectories[[   interactions[1,"ant1.trajectory.row"]    ]][interactions[1,"ant1.trajectory.start"],c("x","y")]
# coord_ant2_start <- trajectories[[   interactions[1,"ant2.trajectory.row"]    ]][interactions[1,"ant2.trajectory.start"],c("x","y")]
# 
# ###if we wanted to know the x-y coordinates of ant1 (7) and ant2 (9) on the last frame of the first interaction listed int his table:
# coord_ant1_end <- trajectories[[   interactions[1,"ant1.trajectory.row"]    ]][interactions[1,"ant1.trajectory.end"],c("x","y")]
# coord_ant2_end <- trajectories[[   interactions[1,"ant2.trajectory.row"]    ]][interactions[1,"ant2.trajectory.end"],c("x","y")]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ########################################
# ##################################################################################################################################
# ###### 6. EXAMPLE USE OF TRAJECTORIES: using R libraries to calculate turn angles and home ranges#################################
# ##################################################################################################################################
# 
# #......................
# 
# ###Second - let's convert the relative time (in seconds since start) into absolute times
# ###For this we need the starting time of that ant's trajectory
# ###Which we find in object positions$trajectories_summary : positions$trajectories_summary[which(positions$trajectories_summary$antID_str=="ant1"),"start"]
# trajectory$time_abs <- positions$trajectories_summary[which(positions$trajectories_summary$antID_str=="ant_1"),"start"] + trajectory$time
# 
# 
# trajectories_unlist$time_abs <- positions$trajectories_summary[which(positions$trajectories_summary$antID_str),"start"] + trajectory$time
# positions$trajectories <- unlist(lapply(positions$trajectories[c(match(positions$trajectories_summary$antID_str,names(positions$trajectories)))],FUN=trajectory_abs_time)) ###the match is once again to ensure we will extract the right information for the right ant
# 
# ##################################################################################################################################
# ###### 7. EXAMPLE USE OF TRAJECTORIES: how to apply a function to all trajectories simultaneously#################################
# ##################################################################################################################################
# ###let's assume we want to know the duration of each trajectory and fill in the information into the positions$trajectories_summary file 
# funTest <- function(trajectory) {
#   if(positions$trajectories_summary$antID_str==trajectories_unlist$.id) {
#     return (positions$trajectories_summary$start + trajectories_unlist$time)
#   } 
# }
# 
# trajectories_unlist$time_abs <- lapply(positions$trajectories[c(match(positions$trajectories_summary$antID_str,names(positions$trajectories)))],FUN=funTest)
# 
# ###First we would define our own function:
# trajectory_absol_time                 <- function(trajectory){ return (positions$trajectories_summary$start + trajectory$time)}
# positions$trajectories_summary$start
# trajectory$time
# trajectory_duration                   <- function(trajectory){ return (max(trajectory$time,na.rm=T)-min(trajectory$time,na.rm=T))}
# 
# 
# ###Second let's apply that function to all trajectories and fill in the results
# positions$trajectories <- unlist(lapply(positions$trajectories[c(match(positions$trajectories_summary$antID_str,names(positions$trajectories)))],FUN=trajectory_abs_time)) ###the match is once again to ensure we will extract the right information for the right ant
# 
# ###Another example: number of coordinates per trajectory
# positions$trajectories_summary$nb_coordinates <- unlist(lapply(positions$trajectories[c(match(positions$trajectories_summary$antID_str,names(positions$trajectories)))],FUN=nrow)) ###the match is once again to ensure we will extract the right information for the right ant
# 
# ###etc. Generally this will save you time compared to a loop
# 


##########################################################################################################################  
##########considerations about interactions from fmQueryComputeAntInteractions reported here for convenience:#############
##########################################################################################################################  
# ####.....
# ####But in many cases we will be interested in a specific type of interactions.
# #### For example: we may be interested only in interactions involving the intersection between body shapes
# ###First we need to figure out the index of the shape corresponding to body shape
# capsules  <- e$antShapeTypeNames()
# body_id <- capsules[which(capsules$name=="Body"),"typeID"]
# ###Then we will specify a matcher in which we are interested in interactions involving capsule 1 for both anta: fmMatcherInteractionType(body_id,body_id)
# ##Let's re-run interactions with that matcher:
# interactions_body       <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportTrajectories = T,matcher = fmMatcherInteractionType(body_id,body_id))
# summaryTrajectory_body      <- interactions_body$summaryTrajectory
# trajectories_body            <- interactions_body$trajectories
# interactions_body            <- interactions_body$interactions
# ###we could also be interested in antennations, i.e. interactions where the antenna of one ant touches the body of the others
# antenna_id <- capsules[which(capsules$name=="Antenna"),"typeID"]
# interactions_antennations            <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportTrajectories = T,matcher = fmMatcherInteractionType(body_id,antenna_id))
# 
# ###we can also combine several matchers
# ###for example, if we want interactions that involve either body/body or antenna/body intersections:
# ? fmMatcherOr
# interactions_body_or_antennations           <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportTrajectories = T,matcher = fmMatcherOr(list(fmMatcherInteractionType(body_id,body_id),fmMatcherInteractionType(body_id,antenna_id))))
# summaryTrajectory_body_or_antennations      <- interactions_body_or_antennations$summaryTrajectory
# trajectories_body_or_antennations           <- interactions_body_or_antennations$trajectories
# interactions_body_or_antennations           <- interactions_body_or_antennations$interactions
# nrow(interactions_body_or_antennations) ###49886 antennations or body/body contacts - indicating that most bioy/body interactions were included within the antennations only
# 
# 
# ###for the rest, let's focus on the body-body interactions - for simplicity, I will rename them without suffix
# summaryTrajectory    <- summaryTrajectory_body
# trajectories           <- trajectories_body
# interactions           <- interactions_body
# 
# 
# ###if we wanted to know the mean x-y coordinates of all interactions, you could either loop over all interactions or write a function and use apply
# ####here I don't use match but I need to be sure that trajectories remained non-shuffled
# mean_coord <- function (x,which_ant,trajectories){
#   coord_mean <- colMeans(trajectories[[as.numeric(x[paste(which_ant,".trajectory.row",sep="")] )]]  [as.numeric(x[paste(which_ant,".trajectory.start",sep="")]):as.numeric(x[paste(which_ant,".trajectory.end",sep="")]),c("x","y")] )
#   return(data.frame(x=coord_mean["x"],y=coord_mean["y"]))
# }
# mean_coordinates <- unlist(apply(interactions, 1, FUN=mean_coord,which_ant="ant1",trajectories=trajectories))
# interactions$mean_x_ant1 <- mean_coordinates[grepl("x",names(mean_coordinates))]
# interactions$mean_y_ant1 <- mean_coordinates[grepl("y",names(mean_coordinates))]
# 
# 
# mean_coordinates <- unlist(apply(interactions, 1, FUN=mean_coord,which_ant="ant2",trajectories=trajectories))
# interactions$mean_x_ant2 <- mean_coordinates[grepl("x",names(mean_coordinates))]
# interactions$mean_y_ant2 <- mean_coordinates[grepl("y",names(mean_coordinates))]


rm(list=(c("e")))
gc()




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
