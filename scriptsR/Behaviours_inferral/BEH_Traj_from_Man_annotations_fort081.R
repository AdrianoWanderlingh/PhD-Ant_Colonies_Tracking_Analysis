###############################################################################
###### extract ant trajectories from hand-annotated data ######################
###############################################################################
print(paste("EXTRACT ANT TRAJ FROM HAND ANNOTATED DATA ",REPLICATE, PERIOD))
####First let's extract ant's trajectories

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
    #Title <- paste(REPLICATE, ", ", PERIOD, ", ", BEH, ROW,", ", "Act:",ACT, ", ", "Rec:",REC, "\nframes ", ENC_FRAME_start, "-", ENC_FRAME_stop, sep="")
    #plot   (y ~ x, traj_ACT, type="l", lwd=4, col="blue4",asp=1, main=Title,cex.main=0.9 ,xlim=c(min(traj_ACT$x,traj_REC$x),max(traj_ACT$x,traj_REC$x)),ylim=c(min(traj_ACT$y,traj_REC$y),max(traj_ACT$y,traj_REC$y))) #, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
    #points (y ~ x, traj_REC, type="l", lwd=4,  col="red4",asp=1)
    # 
    # ## show the headings of each ACT
    # arrows.az (x = traj_ACT$x, 
    #            y = traj_ACT$y, 
    #            azimuth = traj_ACT$angle, 
    #            rho = 10,
    #            HeadWidth=0.1,
    #            units="radians", Kol="blue2", Lwd=1)
    # ## show the headings of each REC
    # arrows.az (x = traj_REC$x, 
    #            y = traj_REC$y, 
    #            azimuth = traj_REC$angle, 
    #            rho = 10,
    #            HeadWidth=0.1,
    #            units="radians", Kol="red2", Lwd=1)
    # 

    
    ##################
    ## INDIVIDUAL TRAJECTORY MEASURES
    # uncorrelated vars used for Constance Tests:  median_acceleration_mmpersec2, rmse_mm, distance_mm, 
    # mean_acceleration_mmpersec2, straightness_index , periodicity_sec (last two not calculated here as they need rediscretisation)
    
    ### angular st.dev 
    StDev_angle_ACT                    <-  angular.deviation(traj_ACT$angle, na.rm = TRUE)
    StDev_angle_REC                    <-  angular.deviation(traj_REC$angle, na.rm = TRUE)
    
    # mean direction of a vector of circular data
    
    # ADRIANO to double check  wheTHER THE circular average is based on the wrong coordinate systEM - E.G. 0-360 CLOCKWISE !!!!!!!
    
    mean_angle_ACT                    <-  mean.circular(traj_ACT$angle, na.rm = TRUE) 
    mean_angle_REC                    <-  mean.circular(traj_REC$angle, na.rm = TRUE)
    
    #---------------------------------------------------------------------------------

    ##Delta_angles: differential between orientation_angle and movement_angle
    # movement_angle: as the orientation_angle remains the same, the movement_angle can change if the movement is perpendicular to the orientation_angle
   
  #movement_angle diff & orientation_angle diff PER ANT PER FRAME
    # change in movement angle calculated as atan(x2/y2) - atan(x1/y1)
    # # SANITY CHECK: create toy data to check the MEANING OF THE +/-'ve signs generated by Movement.angle.diff:
    # traj_ACT <- traj_ACT[1:4,]
    # ## anticlockwise turning ants have POSTIVE values of Movement.angle.diff
    # traj_ACT$x <- c(0,1,1,0); traj_ACT$y <- c(0,0,1,1)
    # ## clockwise turning ants have NEGATIVE values of Movement.angle.diff
    # traj_ACT$x <- c(0,0,1,1); traj_ACT$y <- c(0,1,1,0); plot(y ~ x, traj_ACT )
    
  if (nrow(traj_ACT)>1)
    {
      # variation of movement angle frame by frame
      traj_ACT$Movement_angle_difference_ACT <- Movement.angle.diff(x = traj_ACT) ## To check that this function works, run this: traj_ACT$Movement_angle_difference_CHECK <- NA; for (i in 1: (nrow(traj_ACT)-1)) {traj_ACT$Movement_angle_difference_CHECK[i]   <- atan(traj_ACT[i, "x"] / traj_ACT[i, "y"]) -  atan(traj_ACT[i+1, "x"] / traj_ACT[i+1, "y"])}
      ## and take the absolute value of the 'movement angle'
      #traj_ACT$Movement_angle_difference_ABS <- abs(traj_ACT$Movement_angle_difference)
      # variation of Orientation angle frame by frame
      traj_ACT$Orientation_diff <- Orientation.diff(x = traj_ACT)
      
     }else{print(paste("NO DATA FOR", ACT))}


    if (nrow(traj_REC)>1)
    {
      # variation of movement angle frame by frame
      traj_REC$Movement_angle_difference_REC <- Movement.angle.diff(x = traj_REC)## anticlockwise turning ants have POSTIVE values of Movement.angle.diff & vice-versa - see sanity check above
      #traj_REC$Movement_angle_difference_ABS <- abs(traj_REC$Movement_angle_difference)
      traj_REC$Orientation_diff <- Orientation.diff(x = traj_REC)
      
    }else{print(paste("NO DATA FOR", REC))}
    
    # #test alternative function with abs values, not sure it really works....
    # Orientation.diff_RE        <- function(x)  {   c(  abs( (x[-nrow(x), "angle"]-pi) - (x[-1,      "angle"])-pi) , NA)}
    # traj_ACT$Orientation_diff_RE <- Orientation.diff_RE(x = traj_ACT)
    # traj_ACT$Orientation_diff_RE_Normalised <- traj_ACT$Orientation_diff_RE %% (2*pi)
    # ## DOUBLE TRIPLE CHECK THIS!!
    
    
    # #movement_angle diff & orientation_angle diff PER INTERACTING COUPLE PER FRAME
    # #Calculate Orientation_angle difference
    # Orient_angle_diff <- abs((traj_REC$angle - pi) - (traj_ACT$angle -pi)) %% (2*pi)
    # #traj_BOTH$Orient_angle_diff <- abs((traj_BOTH$REC.angle - pi) - (traj_BOTH$ACT.angle -pi)) %% (2*pi)
    # deg(angle_diff)
    
    # # SANITY CHECK :
    # traj_ACT <- traj_ACT[1:4,]
    # traj_ACT$angle <- NA
    # traj_ACT$angle <- c(pi/2, pi - 0.001, -pi + 0.001, -pi/2)# ant moves from N to WNW, to WSW, to S
    # traj_ACT$Orientation_diff <- Orientation.diff(x = traj_ACT)
    # #traj_ACT$Orientation_diff %% (2*pi)
    # abs(Orientation.diff(x = traj_ACT)-pi) %% (2*pi)
    # #next test:
    # - test what the difference in angles should be in deg and rads 
    # - have a modified version of orientation_diff where -pi is subtracted from each coord as done in traj_BOTH$Orient_angle_diff (why is that done? test if it works beforehand!)

    
    # # sanity check 2: use circular package to handle everything
    # traj_ACT$angle_mirror_negative <- NA
    # traj_ACT$angle_mirror_negative <- pi - traj_ACT$angle[which(traj_ACT$angle<0)]
    # traj_ACT$angle2 <- circular(traj_ACT$angle_mirror_negative, type="angles", units="radians", modulo = "pi", rotation="clock")
    # # delta_angles
    traj_ACT$delta_angles <- traj_ACT$Movement_angle_difference - traj_ACT$Orientation_diff
    traj_REC$delta_angles <- traj_REC$Movement_angle_difference - traj_REC$Orientation_diff
    # mean delta_angles
    mean_delta_angles_ACT <- mean.circular(traj_ACT$delta_angles,na.rm=TRUE)
    mean_delta_angles_REC <- mean.circular(traj_REC$delta_angles,na.rm=TRUE)
    #---------------------------------------------------------------------------------
    
    # ##define trajectory
    # #CHECK THAT THIS IS RIGHT OR DISCARD
    # trajectory_ACT <- TrajFromCoords(data.frame(x=traj_ACT$x,y=traj_ACT$y, frames=traj_ACT$frame), spatialUnits = "px",timeUnits="frames")
    # trajectory_REC <- TrajFromCoords(data.frame(x=traj_REC$x,y=traj_REC$y, frames=traj_REC$frame), spatialUnits = "px",timeUnits="frames")
    # 
    # #trajectory                      <- TrajResampleTime (trajectory_ori,stepTime =desired_step_length_time )
    # ### total distance moved
    # moved_distance_px_ACT                  <- TrajLength(trajectory_ACT)
    # moved_distance_px_REC                  <- TrajLength(trajectory_REC)
    # ### Calculate trajectory derivatives
    # deriv_traj_ACT                <- TrajDerivatives(trajectory_ACT)
    # deriv_traj_REC                <- TrajDerivatives(trajectory_REC)
    # #mean_speed_mmpersec_ACT        <- mean(deriv_traj_ACT$speed,na.rm=T)
    # #  median_speed_mmpersec           <- median(deriv_traj$speed,na.rm=T)
    # #  Mean and median acceleration
    # mean_accel_pxpersec2_ACT     <- mean(deriv_traj_ACT$acceleration,na.rm=T) #units: pxpersec2
    # mean_accel_pxpersec2_REC     <- mean(deriv_traj_REC$acceleration,na.rm=T) #units: pxpersec2
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
    interaction_length_secs <-  ((int_end_frame - int_start_frame)+1)/8  ##includes end frame with +1 . 16 Feb 2022
    
    prop_time_undetected_ACT <- (sum(is.na(traj_BOTH$ACT.x)) / 8) / interaction_length_secs  ## the prop of the interaction in which ACT was seen 
    prop_time_undetected_REC <- (sum(is.na(traj_BOTH$REC.x)) / 8)  / interaction_length_secs ## UNITS are secs / secs --> proportion; MUST BE STRICTLY < 1 !!
    
    # FRAME BY FRAME PARAMETERS
    
    # distance walked
    traj_BOTH$ACT.distance        <- c(NA, with(traj_BOTH, (sqrt(diff(ACT.x)^2 + diff(ACT.y)^2)))) # euclidean distance
    traj_BOTH$REC.distance        <- c(NA, with(traj_BOTH, sqrt(diff(REC.x)^2 + diff(REC.y)^2))) # euclidean distance
    moved_distance_px_ACT         <- sum(traj_BOTH$ACT.distance,na.rm = T)
    moved_distance_px_REC         <- sum(traj_BOTH$REC.distance,na.rm = T)
    

    #traj_BOTH$UNIX_interval <- c(NA, as.numeric(diff(traj_BOTH$UNIX_time), units="secs"))
    traj_BOTH$time_interval <- c(NA, diff(traj_BOTH$frame)/8)
    #speed
    traj_BOTH$ACT.speed_PxPerSec  <- c( with(traj_BOTH, c(NA,(sqrt(diff(ACT.x)^2 + diff(ACT.y)^2))) / time_interval))
    traj_BOTH$REC.speed_PxPerSec  <- c( with(traj_BOTH, c(NA,(sqrt(diff(REC.x)^2 + diff(REC.y)^2))) / time_interval))
    mean_speed_pxpersec_ACT       <- mean(traj_BOTH$ACT.speed_PxPerSec, na.rm=T) 
    mean_speed_pxpersec_REC       <- mean(traj_BOTH$REC.speed_PxPerSec, na.rm=T) 
    #  acceleration
    traj_BOTH$ACT.accel_PxPerSec2 <- c( with(traj_BOTH, c(NA,(diff(ACT.speed_PxPerSec))) / time_interval))
    traj_BOTH$REC.accel_PxPerSec2 <- c( with(traj_BOTH, c(NA,(diff(REC.speed_PxPerSec))) / time_interval))
    mean_accel_pxpersec2_ACT      <- mean(traj_BOTH$ACT.accel_PxPerSec2, na.rm=T)
    mean_accel_pxpersec2_REC      <- mean(traj_BOTH$REC.accel_PxPerSec2, na.rm=T) 

    # jerk (diff in accelerations)
    traj_BOTH$ACT.jerk_PxPerSec3 <- c( with(traj_BOTH, c(NA,(diff(ACT.accel_PxPerSec2))) / time_interval))
    traj_BOTH$REC.jerk_PxPerSec3 <- c( with(traj_BOTH, c(NA,(diff(REC.accel_PxPerSec2))) / time_interval))
    mean_jerk_PxPerSec3_ACT      <- mean( traj_BOTH$ACT.jerk_PxPerSec3, na.rm=T)
    mean_jerk_PxPerSec3_REC      <- mean( traj_BOTH$REC.jerk_PxPerSec3, na.rm=T)
    ##################
    ## INTERACTING PAIR TRAJECTORY MEASURES
    
    #straight line - euclidean distance
    traj_BOTH$straightline_dist_px <-  sqrt((traj_BOTH$ACT.x-traj_BOTH$REC.x)^2+(traj_BOTH$ACT.y-traj_BOTH$REC.y)^2)
    mean_strghtline_dist_px <- mean(traj_BOTH$straightline_dist_px, na.rm=TRUE)
    #angular difference
    traj_BOTH$Orient_angle_diff <- abs((traj_BOTH$REC.angle - pi) - (traj_BOTH$ACT.angle -pi)) %% (2*pi)
    mean_orient_angle_diff <-  mean(traj_BOTH$Orient_angle_diff, na.rm=TRUE)
    
    traj_BOTH$Movement_angle_diff <- abs((traj_BOTH$Movement_angle_difference_REC - pi) - (traj_BOTH$Movement_angle_difference_ACT -pi)) %% (2*pi)
    mean_movement_angle_diff <-  mean(traj_BOTH$Movement_angle_diff, na.rm=TRUE)
    
    
    ###############
    
    summary_MAN_ROW <- data.frame(REPLICATE, PERIOD, BEH, ROW, Act_Name, Rec_Name, 
                                  StDev_angle_ACT, StDev_angle_REC,
                                  mean_angle_ACT, mean_angle_REC,
                                  mean_delta_angles_ACT, mean_delta_angles_REC,
                                  moved_distance_px_ACT, moved_distance_px_REC,
                                  mean_speed_pxpersec_ACT, mean_speed_pxpersec_REC,
                                  mean_accel_pxpersec2_ACT, mean_accel_pxpersec2_REC,
                                  mean_jerk_PxPerSec3_ACT, mean_jerk_PxPerSec3_REC,
                                  rmsd_px_ACT,rmsd_px_REC,
                                  int_start_frame , int_end_frame, interaction_length_secs,
                                  prop_time_undetected_ACT, prop_time_undetected_REC,
                                  mean_strghtline_dist_px, 
                                  mean_orient_angle_diff, mean_movement_angle_diff,
                                  #when adding a new variable, it must be included in the reshape rule for data plotting
                                  stringsAsFactors = F)
    
    
    
    #create new variable by pasting ant numbers "low,high" for summary_MAN_ROW
    summary_MAN_ROW$ant1 <- as.numeric(gsub("ant_","", summary_MAN_ROW$Act_Name))
    summary_MAN_ROW$ant2 <- as.numeric(gsub("ant_","", summary_MAN_ROW$Rec_Name))
    ## the ant 1 & ant 2 labels are not strictly ascending, so sort them so ant 1 alwasy < ant 2
    summary_MAN_ROW$pair <- apply(summary_MAN_ROW[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
    
    
    ## lets be cautious:
    StDev_angle_ACT<-NULL 
    StDev_angle_REC<-NULL 
    mean_angle_ACT<-NULL 
    mean_angle_REC<-NULL 
    mean_delta_angles_ACT<-NULL 
    mean_delta_angles_REC<-NULL 
    moved_distance_px_ACT<-NULL 
    moved_distance_px_REC<-NULL 
    mean_speed_pxpersec_ACT<-NULL 
    mean_speed_pxpersec_REC<-NULL 
    mean_accel_pxpersec2_ACT<-NULL 
    mean_accel_pxpersec2_REC<-NULL 
    mean_jerk_PxPerSec3_ACT<-NULL 
    mean_jerk_PxPerSec3_REC<-NULL 
    rmsd_px_ACT<-NULL 
    rmsd_px_REC<-NULL 
    int_start_frame <-NULL 
    int_end_frame<-NULL 
    interaction_length_secs<-NULL 
    prop_time_undetected_ACT<-NULL 
    prop_time_undetected_REC<-NULL 
    mean_strghtline_dist_px<-NULL 
    mean_orient_angle_diff<-NULL 
    mean_movement_angle_diff <- NULL
    
    #NO UNIX TIME BUT FRAMES. interaction AREA= NEST, FORAGING delete x y coords
    interacts_MAN_ROW <- data.frame(REPLICATE=REPLICATE,PERIOD=PERIOD,ROW=ROW,BEH=BEH,Act_Name=Act_Name,Rec_Name=Rec_Name,
                                    traj_BOTH,
                                    stringsAsFactors = F)
    
    ## stack -by the end of the loop, this should just contain all events for replicate X, period Y
    interacts_MAN_REP_PER <- rbind(interacts_MAN_REP_PER, interacts_MAN_ROW)
    summary_MAN_REP_PER   <- rbind(summary_MAN_REP_PER,     summary_MAN_ROW)
    
  }##ROW
}##BEH
