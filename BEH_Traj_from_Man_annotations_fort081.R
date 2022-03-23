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
    # ## Plot trajectories of both actor & receiver, show on the same panel
    # Title <- paste(REPLICATE, ", ", PERIOD, ", ", BEH, ROW,", ", "Act:",ACT, ", ", "Rec:",REC, "\nframes ", ENC_FRAME_start, "-", ENC_FRAME_stop, sep="")
    # plot   (y ~ x, traj_ACT, type="l", lwd=4, col="blue4",asp=1, main=Title,cex.main=0.9 ,xlim=c(min(traj_ACT$x,traj_REC$x),max(traj_ACT$x,traj_REC$x)),ylim=c(min(traj_ACT$y,traj_REC$y),max(traj_ACT$y,traj_REC$y))) #, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
    # points (y ~ x, traj_REC, type="l", lwd=4,  col="red4",asp=1)
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

    ################################################################    
    ######### JUMPS & JITTER #######################################
    ################################################################
    
    ##  jumps detection by frame
    traj_ACT$dt_FRAME            <- c(1, diff(traj_ACT$frame)) ## time difference in frames
    traj_REC$dt_FRAME            <- c(1, diff(traj_REC$frame))
    #traj_ACT$dt_TIME             <- c(diff(traj_ACT$Frame)/8,0.125) ## time difference in seconds
    
    #jumps detection by distance is reported below in traj_BOTH$ACT.distance and traj_BOTH$REC.distance   
    
    ##################
    ## INDIVIDUAL TRAJECTORY MEASURES
    
    #turnangles
    trjACT <- TrajFromCoords(data.frame(traj_ACT$x,traj_ACT$y,traj_ACT$frame))
    trjREC <- TrajFromCoords(data.frame(traj_REC$x,traj_REC$y,traj_REC$frame))
    traj_ACT$TurnAngle <- c(NA,NA,TrajAngles(trjACT))
    traj_REC$TurnAngle <- c(NA,NA,TrajAngles(trjREC))

    ##Mov_Orient_delta_angle: differential between orientation_angle and movement_angle
    # movement_angle: as the orientation_angle remains the same, the movement_angle can change if the movement is perpendicular to the orientation_angle
   
    # # SANITY CHECK 
    # # movement_angle diff & orientation_angle diff PER ANT PER FRAME
    # # change in movement angle calculated as atan(x2/y2) - atan(x1/y1)
    # # SANITY CHECK: create toy data to check the MEANING OF THE +/-'ve signs generated by Movement.angle.diff:
    # traj_ACT <- traj_ACT[1:4,]
    # ## anticlockwise turning ants have POSTIVE values of Movement.angle.diff
    # traj_ACT$x <- c(0,1,1,0); traj_ACT$y <- c(0,0,1,1)
    # ## clockwise turning ants have NEGATIVE values of Movement.angle.diff
    # traj_ACT$x <- c(0,0,1,1); traj_ACT$y <- c(0,1,1,0); plot(y ~ x, traj_ACT )
    
    if (nrow(traj_ACT)>1)
    {
      # variation of movement angle frame by frame
      traj_ACT$Movement_angle_difference <- Movement.angle.diff(x = traj_ACT) ## To check that this function works, run this: traj_ACT$Movement_angle_difference_CHECK <- NA; for (i in 1: (nrow(traj_ACT)-1)) {traj_ACT$Movement_angle_difference_CHECK[i]   <- atan(traj_ACT[i, "x"] / traj_ACT[i, "y"]) -  atan(traj_ACT[i+1, "x"] / traj_ACT[i+1, "y"])}
      # variation of Orientation angle frame by frame
      traj_ACT$Orientation_diff <- Orientation.diff(x = traj_ACT)
    }else{print(paste("NO DATA FOR", ACT))}
    
    if (nrow(traj_REC)>1)
    {
      # variation of movement angle frame by frame
      traj_REC$Movement_angle_difference <- Movement.angle.diff(x = traj_REC)## anticlockwise turning ants have POSTIVE values of Movement.angle.diff & vice-versa - see sanity check above
      traj_REC$Orientation_diff <- Orientation.diff(x = traj_REC)
    }else{print(paste("NO DATA FOR", REC))}
    
    # # SANITY CHECK :
    # traj_ACT <- traj_ACT[1:4,]
    # traj_ACT$angle <- NA
    # traj_ACT$angle <- c(pi/2, pi - 0.001, -pi + 0.001, -pi/2)# ant moves from N to WNW, to WSW, to S
    # traj_ACT$Orientation_diff <- Orientation.diff(x = traj_ACT)
    # #traj_ACT$Orientation_diff %% (2*pi)

    # delta_angle = difference between movement angle and orientation angle of each ant
    # between 0 - 90 deg
    Mov_Orient_delta_angle_ACT_raw <-  traj_ACT$Movement_angle_difference - traj_ACT$Orientation_diff
    traj_ACT$Mov_Orient_delta_angle <- c(unlist(lapply(Mov_Orient_delta_angle_ACT_raw[!is.na(Mov_Orient_delta_angle_ACT_raw)],FUN=inclination_angle)),NA)
    Mov_Orient_delta_angle_REC_raw <-  traj_REC$Movement_angle_difference - traj_REC$Orientation_diff
    traj_REC$Mov_Orient_delta_angle <- c(unlist(lapply(Mov_Orient_delta_angle_REC_raw[!is.na(Mov_Orient_delta_angle_REC_raw)],FUN=inclination_angle)),NA)
    
    ### mean square displacement
    mean_sqrt_err_px_ACT              <-  sum( (traj_ACT$x-mean(traj_ACT$x))^2 + (traj_ACT$y-mean(traj_ACT$y))^2 )/length(na.omit(traj_ACT$x))
    mean_sqrt_err_px_REC              <-  sum( (traj_REC$x-mean(traj_REC$x))^2 + (traj_REC$y-mean(traj_REC$y))^2 )/length(na.omit(traj_REC$x))
    
    # Convex Hull
    box.coords <- as.matrix(traj_ACT[, c("x", "y")]); box.hpts <- chull(x = traj_ACT$x, y = traj_ACT$y) # calculate convex hull for x and y columns
    box.hpts <- c(box.hpts, box.hpts[1]) #add first coord at the bottom to close area
    box.chull.coords <- box.coords[box.hpts,]
    chull.poly <- Polygon(box.chull.coords, hole=F); chull_area_ACT <- chull.poly@area  #calculate area
    
    box.coords <- as.matrix(traj_REC[, c("x", "y")]); box.hpts <- chull(x = traj_REC$x, y = traj_REC$y) # calculate convex hull for x and y columns
    box.hpts <- c(box.hpts, box.hpts[1]) #add first coord at the bottom to close area
    box.chull.coords <- box.coords[box.hpts,]
    chull.poly <- Polygon(box.chull.coords, hole=F); chull_area_REC <- chull.poly@area  #calculate area
    
    #rename trajectories columns NOT TO BE MERGED for Act and Rec
    colnames(traj_ACT) <- paste("ACT", colnames(traj_ACT), sep = ".")
    colnames(traj_REC) <- paste("REC", colnames(traj_REC), sep = ".")
    #except frame
    names(traj_ACT)[names(traj_ACT) == 'ACT.frame'] <- 'frame'
    names(traj_REC)[names(traj_REC) == 'REC.frame'] <- 'frame'
    
    #merge trajectories matching by time
    # traj_BOTH <- merge(traj_ACT, traj_REC,all=T, by=c("UNIX_time"))  ## CAREFUL NOT TO INCLUDE 'time' AS IF THE FIRST OBSERVATION TIME OF ACT & REC IS NOT IDENTICAL, THEN MERGE WILL TREAT THE SAME UNIX TIME AS DIFFERENT ...
    traj_BOTH <- merge(traj_ACT, traj_REC,all=T, by=c("frame"))  ## 21 Jan 2022: Changed to sue raw unix seconds to allow simpler matching below (posix formats cause problems with the matching...)
    
    
    ## merge drops the column attributes - add back the POSIX time format 
    # traj_BOTH$UNIX_time <- as.POSIXct(traj_BOTH$UNIX_time, tz="GMT",origin="1970-01-01 00:00:00")
    
    ## measure the length *in seconds* of the interaction between ACT & REC
    # int_length_secs <- as.numeric(difftime ( max(traj_BOTH$UNIX_time, na.rm=T), min(traj_BOTH$UNIX_time, na.rm=T), units="secs"))
    int_start_frame <- min(traj_BOTH$frame, na.rm=T)
    int_end_frame <- max(traj_BOTH$frame, na.rm=T)
    int_length_secs <-  ((int_end_frame - int_start_frame)+1)/8  ##includes end frame with +1 . 16 Feb 2022
    
    prop_time_undetected_ACT <- (sum(is.na(traj_BOTH$ACT.x)) / 8) / int_length_secs  ## the prop of the interaction in which ACT was seen 
    prop_time_undetected_REC <- (sum(is.na(traj_BOTH$REC.x)) / 8)  / int_length_secs ## UNITS are secs / secs --> proportion; MUST BE STRICTLY < 1 !!
    
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

    #  acceleration
    traj_BOTH$ACT.accel_PxPerSec2 <- c( with(traj_BOTH, c(NA,(diff(ACT.speed_PxPerSec))) / time_interval))
    traj_BOTH$REC.accel_PxPerSec2 <- c( with(traj_BOTH, c(NA,(diff(REC.speed_PxPerSec))) / time_interval))

    # jerk (diff in accelerations)
    traj_BOTH$ACT.jerk_PxPerSec3  <- c( with(traj_BOTH, c(NA,(diff(ACT.accel_PxPerSec2))) / time_interval))
    traj_BOTH$REC.jerk_PxPerSec3  <- c( with(traj_BOTH, c(NA,(diff(REC.accel_PxPerSec2))) / time_interval))
    ##################
    
    ## INTERACTING PAIR TRAJECTORY MEASURES
    
    #straight line - euclidean distance
    traj_BOTH$straightline_dist_px<-  sqrt((traj_BOTH$ACT.x-traj_BOTH$REC.x)^2+(traj_BOTH$ACT.y-traj_BOTH$REC.y)^2)
    #Orientation difference of the pair (0-180deg)
    # traj_BOTH$Orient_angle_diff   <- abs((traj_BOTH$REC.angle - pi) - (traj_BOTH$ACT.angle -pi)) %% (2*pi)
    #traj_BOTH$pair_orient_diff <- abs(zero_to_2pi(traj_BOTH$ACT.angle  - traj_BOTH$REC.angle))
    traj_BOTH$pair_orient_diff <- abs((traj_BOTH$ACT.angle  - traj_BOTH$REC.angle) %% pi)
    

    # angular velocity: ω = (α₂ - α₁) / t = Δα / t
    traj_BOTH$ang_Velocity        <- (traj_BOTH$ACT.angle - traj_BOTH$REC.angle) / traj_BOTH$time_interval
    
    ##################################
    ## mean values for summary #######
    ##################################
    
    ## Apply THRESHOLDs to exclude gaps in which the individuals were not detected for very long (FRAME)
    # remove traj points when Trajectory$dt_FRAME > DT_frame_THRESHOLD
    # dt_FRAME corresponds to frame between T and T+1. This should be applied to all movement characteristics
    traj_BOTH[which(traj_BOTH$ACT.dt_FRAME>DT_frame_THRESHOLD), c("ACT.speed_PxPerSec","ACT.accel_PxPerSec2","ACT.jerk_PxPerSec3","ACT.TurnAngle","ang_Velocity","pair_orient_diff")] <- NA
    traj_BOTH[which(traj_BOTH$REC.dt_FRAME>DT_frame_THRESHOLD), c("REC.speed_PxPerSec","REC.accel_PxPerSec2","REC.jerk_PxPerSec3","REC.TurnAngle","ang_Velocity","pair_orient_diff")] <- NA
    ## Apply THRESHOLDs to exclude jitter in the individuals' movement (DISTANCE)
    # remove traj points when traj_BOTH$ACT.distance$distance < DT_dist_THRESHOLD
    traj_BOTH[which(traj_BOTH$ACT.distance<DT_dist_THRESHOLD), c("ACT.speed_PxPerSec","ACT.accel_PxPerSec2","ACT.jerk_PxPerSec3","ACT.TurnAngle","ang_Velocity","pair_orient_diff")] <- NA
    traj_BOTH[which(traj_BOTH$REC.distance<DT_dist_THRESHOLD), c("REC.speed_PxPerSec","REC.accel_PxPerSec2","REC.jerk_PxPerSec3","REC.TurnAngle","ang_Velocity","pair_orient_diff")] <- NA

    ### angular st.dev 
    # THIS SHOULD BE SOMEHOW INCORPORATE THE JUMPS CUTTING PROCEDURE (but I can't really just cut the coords x and y)
    StDev_angle_ACT                    <-  angular.deviation(traj_BOTH$ACT.angle, na.rm = TRUE)
    StDev_angle_REC                    <-  angular.deviation(traj_BOTH$REC.angle, na.rm = TRUE)
    
    #turnangle stDev
    stDev_turnAngle_ACT                 <-  angular.deviation(traj_BOTH$ACT.TurnAngle, na.rm = TRUE)
    stDev_turnAngle_REC                 <-  angular.deviation(traj_BOTH$REC.TurnAngle, na.rm = TRUE)
    
    #means are then calculated on the reduced dataset
    #individual
    mean_abs_turnAngle_ACT                  <- mean.circular(abs(traj_BOTH$ACT.TurnAngle),na.rm=T)
    mean_abs_turnAngle_REC                  <- mean.circular(abs(traj_BOTH$REC.TurnAngle),na.rm=T)
    mean_Mov_Orient_delta_angle_ACT         <- mean.circular(traj_BOTH$ACT.Mov_Orient_delta_angle,na.rm=T)
    mean_Mov_Orient_delta_angle_REC         <- mean.circular(traj_BOTH$REC.Mov_Orient_delta_angle,na.rm=T)
    
    mean_speed_pxpersec_ACT                 <- mean(traj_BOTH$ACT.speed_PxPerSec, na.rm=T) 
    mean_speed_pxpersec_REC                 <- mean(traj_BOTH$REC.speed_PxPerSec, na.rm=T)
    mean_accel_pxpersec2_ACT                <- mean(traj_BOTH$ACT.accel_PxPerSec2, na.rm=T)
    mean_accel_pxpersec2_REC                <- mean(traj_BOTH$REC.accel_PxPerSec2, na.rm=T) 
    mean_jerk_PxPerSec3_ACT                 <- mean(traj_BOTH$ACT.jerk_PxPerSec3, na.rm=T)
    mean_jerk_PxPerSec3_REC                 <- mean(traj_BOTH$REC.jerk_PxPerSec3, na.rm=T)
  
    #pair
    mean_strghtline_dist_px                 <- mean(traj_BOTH$straightline_dist_px, na.rm=T)
    mean_ang_Velocity                       <- mean.circular(abs(traj_BOTH$ang_Velocity), na.rm=T) #abs value to avoid a mean around 0 for back and forth movement
    mean_pair_orient_diff                   <- mean.circular(traj_BOTH$pair_orient_diff, na.rm=T) #between 0-180 deg
    #mean_movement_angle_diff             <- mean.circular(traj_BOTH$Movement_angle_diff, na.rm=T)
    
    
    ###############
    
    summary_MAN_ROW <- data.frame(REPLICATE, PERIOD, BEH, ROW, Act_Name, Rec_Name, 
                                  StDev_angle_ACT, StDev_angle_REC,
                                  mean_abs_turnAngle_ACT,mean_abs_turnAngle_REC,
                                  stDev_turnAngle_ACT,stDev_turnAngle_REC,
                                  mean_Mov_Orient_delta_angle_ACT, mean_Mov_Orient_delta_angle_REC,
                                  moved_distance_px_ACT, moved_distance_px_REC,
                                  mean_speed_pxpersec_ACT, mean_speed_pxpersec_REC,
                                  mean_accel_pxpersec2_ACT, mean_accel_pxpersec2_REC,
                                  mean_jerk_PxPerSec3_ACT, mean_jerk_PxPerSec3_REC,
                                  mean_sqrt_err_px_ACT,mean_sqrt_err_px_REC,
                                  chull_area_ACT,chull_area_REC,
                                  int_start_frame , int_end_frame, 
                                  int_length_secs,
                                  prop_time_undetected_ACT, prop_time_undetected_REC,
                                  mean_strghtline_dist_px, 
                                  mean_pair_orient_diff,
                                  mean_ang_Velocity,
                                  #when adding a new variable, it must be included in the reshape rule for data plotting
                                  stringsAsFactors = F)
    
    summary_MAN_ROW$ROW <- as.factor(summary_MAN_ROW$ROW )
    
    #create new variable by pasting ant numbers "low,high" for summary_MAN_ROW
    summary_MAN_ROW$ant1 <- as.numeric(gsub("ant_","", summary_MAN_ROW$Act_Name))
    summary_MAN_ROW$ant2 <- as.numeric(gsub("ant_","", summary_MAN_ROW$Rec_Name))
    ## the ant 1 & ant 2 labels are not strictly ascending, so sort them so ant 1 alwasy < ant 2
    summary_MAN_ROW$pair <- apply(summary_MAN_ROW[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
    
    
    ## lets be cautious:
    StDev_angle_ACT<-NULL 
    StDev_angle_REC<-NULL 
    mean_abs_angle_ACT<-NULL 
    mean_abs_angle_REC<-NULL 
    mean_abs_turnAngle_ACT<-NULL
    mean_abs_turnAngle_REC<-NULL
    mean_Mov_Orient_delta_angle_ACT<-NULL 
    mean_Mov_Orient_delta_angle_REC<-NULL 
    moved_distance_px_ACT<-NULL 
    moved_distance_px_REC<-NULL 
    mean_speed_pxpersec_ACT<-NULL 
    mean_speed_pxpersec_REC<-NULL 
    mean_accel_pxpersec2_ACT<-NULL 
    mean_accel_pxpersec2_REC<-NULL 
    mean_jerk_PxPerSec3_ACT<-NULL 
    mean_jerk_PxPerSec3_REC<-NULL 
    mean_sqrt_err_px_ACT<-NULL 
    mean_sqrt_err_px_REC<-NULL 
    int_start_frame <-NULL 
    int_end_frame<-NULL 
    int_length_secs<-NULL 
    prop_time_undetected_ACT<-NULL 
    prop_time_undetected_REC<-NULL 
    mean_strghtline_dist_px<-NULL 
    mean_pair_orient_diff<-NULL 
    mean_movement_angle_diff <- NULL
    stDev_turnAngle_ACT<- NULL
    stDev_turnAngle_REC<- NULL
    chull_area_ACT<- NULL
    chull_area_REC<- NULL
    mean_ang_Velocity<- NULL
    
    #NO UNIX TIME BUT FRAMES. interaction AREA= NEST, FORAGING delete x y coords
    interacts_MAN_ROW <- data.frame(REPLICATE=REPLICATE,PERIOD=PERIOD,ROW=ROW,BEH=BEH,Act_Name=Act_Name,Rec_Name=Rec_Name,
                                    traj_BOTH,
                                    stringsAsFactors = F)
    interacts_MAN_ROW$ROW <- as.factor(interacts_MAN_ROW$ROW )
    
    ## stack -by the end of the loop, this should just contain all events for replicate X, period Y
    interacts_MAN_REP_PER <- rbind(interacts_MAN_REP_PER, interacts_MAN_ROW)
    summary_MAN_REP_PER   <- rbind(summary_MAN_REP_PER,     summary_MAN_ROW)
    
  }##ROW
}##BEH



# plot.circular(summary_AUvcTO_REP_PER$mean_ang_Velocity, pch = 16, cex = 0.8, stack=TRUE, bins=100)
# 
# plot.circular(summary_AUvcTO_REP_PER$mean_pair_orient_diff, pch = 16, cex = 0.8, stack=TRUE, bins=100)

# x <- 1; y <- 1
# atan(x/y) * 180/pi
# #rad to deg
# rad * 180/pi