
######  GETTING VARS OFF INTERACTIONS ######

#IT HAS TO BE DONE INSIDE THE REP-PER LOOP BEFORE THE interaction_AUTO STACKING

#fmQueryComputeInteractions $interactions
#interacts_AUTO_REP_PER$interactions

#fmQueryComputeInteractions $trajectories
int_traj_AUTO_REP_PER <- interacts_AUTO_REP_PER$trajectories


# CONSIDER THAT THE REC/ACT IS NOT KNOWN FOR AUTO_INTERACTS SO THE VARIABLES CONSIDERED SHOULD ONLY BE THE ONES NOT LINKED TO 1 SPECIFIC INDIVIDUAL. 
# USE FOR LOOP  TO CREATE, INTERACTION BY INTERACTION, THE LIST OF PARAMS X Y ANGLE WITH EXP INFOS ATTACHED

#let's create an empy data.frame that will be the interacts_AUTO_REP_PER but with trajectory info

# FINAL RESULT SHOULD LOOK LIKE interacts_MAN_REP_PER
for (INT in 1:nrow(interacts_AUTO_REP_PER$interactions)) {
    TRAJ_ANT_list <- list()
    for (which_ant in c("ant1","ant2")) {
      
      ## extract actor, receiver IDs & start & end times from the hand-annotated data
      #INTER <- interacts_AUTO_REP_PER$interactions[INT,]
      ANT <- interacts_AUTO_REP_PER$interactions[,which_ant][INT]
      
      #extract frame lenght between start and end (int_start_frame ? int_end_frame)
      INT_start_frame_ANT <- interacts_AUTO_REP_PER$interactions[,"int_start_frame"][INT]
      INT_end_frame_ANT <- interacts_AUTO_REP_PER$interactions[,"int_end_frame"][INT]
      
      ## extract the trajectory for ANT
      #trajectory.row
      INT_TRAJ_ROW_ANT <- interacts_AUTO_REP_PER$interactions[,paste0(which_ant,".trajectory.row")][INT]
      #trajectory.start
      INT_FRAME_start  <- as.numeric(interacts_AUTO_REP_PER$interactions[,paste0(which_ant,".trajectory.start")][INT])
      #trajectory.end
      INT_FRAME_end  <- as.numeric(interacts_AUTO_REP_PER$interactions[,paste0(which_ant,".trajectory.end")][INT])
      
      print(paste("Interaction number",INT,which_ant,ANT,"INT_TRAJ_ROW_ANT",INT_TRAJ_ROW_ANT,"INT_FRAME_start",INT_FRAME_start,"INT_FRAME_end",INT_FRAME_end))
     
      # which_ant.trajectory.row the corresponding index in $trajectory.
      TRAJ_ANT <- interacts_AUTO_REP_PER$trajectories[[INT_TRAJ_ROW_ANT]]
      TRAJ_ANT_INT <- TRAJ_ANT [ which(as.numeric(rownames(TRAJ_ANT)) >= INT_FRAME_start & as.numeric(rownames(TRAJ_ANT)) <= INT_FRAME_end),]
      #add frame info
      #NOTE: ANT1 AND ANT2 TRAJ ARE NOT THE SAME LENGHT (MAYBE FOR MISSED FRAMES?). cHECK FOR GAPS IN TRAJS AND EXPAND GRID FOR MISSIN VALS (TIME >0.125)
      
      TRAJ_ANT_INT$frame <- as.numeric(seq(INT_start_frame_ANT:INT_end_frame_ANT))
      
      ##################
      ## INDIVIDUAL TRAJECTORY MEASURES
      
      ### angular st.dev 
      StDev_angle                  <-  angular.deviation(TRAJ_ANT_INT$angle, na.rm = TRUE)
      #circular average
      mean_angle                   <-  mean.circular(TRAJ_ANT_INT$angle, na.rm = TRUE) 
      
      ## Delta_angles: differential between orientation_angle and movement_angle
      # movement_angle: as the orientation_angle remains the same, the movement_angle can change if the movement is perpendicular to the orientation_angle
      
      #movement_angle diff & orientation_angle diff PER ANT PER FRAME
      if (nrow(TRAJ_ANT_INT)>1)
      {
        # variation of movement angle frame by frame
        TRAJ_ANT_INT$Movement_angle_difference <- Movement.angle.diff(x = TRAJ_ANT_INT) ## To check that this function works, run this: TRAJ_ANT_INT$Movement_angle_difference_CHECK <- NA; for (i in 1: (nrow(TRAJ_ANT_INT)-1)) {TRAJ_ANT_INT$Movement_angle_difference_CHECK[i]   <- atan(TRAJ_ANT_INT[i, "x"] / TRAJ_ANT_INT[i, "y"]) -  atan(TRAJ_ANT_INT[i+1, "x"] / TRAJ_ANT_INT[i+1, "y"])}
        ## and take the absolute value of the 'movement angle'
        #TRAJ_ANT_INT$Movement_angle_difference_ABS <- abs(TRAJ_ANT_INT$Movement_angle_difference)
        # variation of Orientation angle frame by frame
        TRAJ_ANT_INT$Orientation_diff <- Orientation.diff(x = TRAJ_ANT_INT)
        
      }else{print(paste("NO DATA FOR", ANT))}
      
      # # SANITY CHECK :
      # TRAJ_ANT_INT <- TRAJ_ANT_INT[1:4,]
      # TRAJ_ANT_INT$angle <- NA
      # TRAJ_ANT_INT$angle <- c(pi/2, pi - 0.001, -pi + 0.001, -pi/2)# ant moves from N to WNW, to WSW, to S
      # TRAJ_ANT_INT$Orientation_diff <- Orientation.diff(x = TRAJ_ANT_INT)
      # #TRAJ_ANT_INT$Orientation_diff %% (2*pi)
      # abs(Orientation.diff(x = TRAJ_ANT_INT)-pi) %% (2*pi)
      # #next test:
      # - test what the difference in angles should be in deg and rads 
      # - have a modified version of orientation_diff where -pi is subtracted from each coord as done in traj_BOTH$Orient_angle_diff (why is that done? test if it works beforehand!)
      
      ## delta_angles
      TRAJ_ANT_INT$delta_angles <- TRAJ_ANT_INT$Movement_angle_difference - TRAJ_ANT_INT$Orientation_diff
      ## mean delta_angles
      mean_delta_angles <- mean.circular(TRAJ_ANT_INT$delta_angles,na.rm=TRUE)
      ## root mean square displacement
      rmsd_px                            <-  sqrt(sum( (TRAJ_ANT_INT$x-mean(TRAJ_ANT_INT$x))^2 + (TRAJ_ANT_INT$y-mean(TRAJ_ANT_INT$y))^2 )/length(na.omit(TRAJ_ANT_INT$x)))

      # ## measure the length *in seconds* of the interaction between ACT & REC
      # # interaction_length_secs <- as.numeric(difftime ( max(traj_BOTH$UNIX_time, na.rm=T), min(traj_BOTH$UNIX_time, na.rm=T), units="secs"))
      # int_start_frame <- min(traj_BOTH$frame, na.rm=T)
      # int_end_frame <- max(traj_BOTH$frame, na.rm=T)
      # interaction_length_secs <-  (int_end_frame - int_start_frame)/8  ## 21 Jan 2022
      # 
      # prop_time_undetected_ANT <- (sum(is.na(TRAJ_ANT_INT$x)) / 8) / interaction_length_secs  ## the prop of the interaction in which ACT was seen 

      ##### FRAME BY FRAME PARAMETERS
      
      # distance walked
      TRAJ_ANT_INT$distance        <- c(NA, with(TRAJ_ANT_INT, (sqrt(diff(x)^2 + diff(y)^2)))) # euclidean distance
      moved_distance_px            <- sum(TRAJ_ANT_INT$distance,na.rm = T)
      
      TRAJ_ANT_INT$time_interval   <- c(NA, diff(TRAJ_ANT_INT$frame)/8)
      #speed
      TRAJ_ANT_INT$speed_PxPerSec  <- c( with(TRAJ_ANT_INT, c(NA,(sqrt(diff(x)^2 + diff(y)^2))) / time_interval))
      mean_speed_pxpersec         <- mean(TRAJ_ANT_INT$speed_PxPerSec, na.rm=T) 
      #  acceleration
      TRAJ_ANT_INT$accel_PxPerSec2 <- c( with(TRAJ_ANT_INT, c(NA,(diff(speed_PxPerSec))) / time_interval))
      mean_accel_pxpersec2         <- mean(TRAJ_ANT_INT$accel_PxPerSec2, na.rm=T)
      # jerk (diff in accelerations)
      TRAJ_ANT_INT$jerk_PxPerSec3  <- c( with(TRAJ_ANT_INT, c(NA,(diff(accel_PxPerSec2))) / time_interval))
      mean_jerk_PxPerSec3          <- mean( TRAJ_ANT_INT$jerk_PxPerSec3, na.rm=T)
  
        
      
      #names(TRAJ_ANT_INT)[names(TRAJ_ANT_INT) == 'x'] <- paste0(which_ant,".x"); names(TRAJ_ANT_INT)[names(TRAJ_ANT_INT) == 'y'] <- paste0(which_ant,".y"); names(TRAJ_ANT_INT)[names(TRAJ_ANT_INT) == 'angle'] <- paste0(which_ant,".angle")
      #assign each which_ant to list
      TRAJ_ANT_list[[which_ant]] <- TRAJ_ANT_INT

  }#which_ant
  
    # turn your list into a dataframe
    TRAJ_AUTO_BOTH <- data.frame(TRAJ_ANT_list)
    TRAJ_AUTO_BOTH[c("ant1.time","ant2.time","ant2.zone","ant2.frame")] <- list(NULL)
    names(TRAJ_AUTO_BOTH)[names(TRAJ_AUTO_BOTH) == 'ant1.zone'] <- "zone"; names(TRAJ_AUTO_BOTH)[names(TRAJ_AUTO_BOTH) == 'ant1.frame'] <- "frame"
    #add info from INT
    TRAJ_AUTO_BOTH <- merge(interacts_AUTO_REP_PER$interactions[INT,c("ant1","ant2","pair")],TRAJ_AUTO_BOTH)
    #col_order <- c("frame","ant1","ant2","pair","ant1.x","ant1.y","ant1.angle","ant2.x","ant2.y","ant2.angle","zone")
    #TRAJ_AUTO_BOTH <- TRAJ_AUTO_BOTH[, col_order]

    ##################
    ## INTERACTING PAIR TRAJECTORY MEASURES
    
    #straight line - euclidean distance
    TRAJ_AUTO_BOTH$straightline_dist_px <-  sqrt((TRAJ_AUTO_BOTH$ant1.x-TRAJ_AUTO_BOTH$ant2.x)^2+(TRAJ_AUTO_BOTH$ant1.y-TRAJ_AUTO_BOTH$ant1.y)^2)
    mean_strghtline_dist_px <- mean(TRAJ_AUTO_BOTH$straightline_dist_px, na.rm=TRUE)
    #angular difference
    TRAJ_AUTO_BOTH$Orient_angle_diff <- abs((TRAJ_AUTO_BOTH$ant2.angle - pi) - (TRAJ_AUTO_BOTH$ant1.angle -pi)) %% (2*pi)
    mean_orient_angle_diff <-  mean(TRAJ_AUTO_BOTH$Orient_angle_diff, na.rm=TRUE)
    
    TRAJ_AUTO_BOTH$Movement_angle_diff <- abs((TRAJ_AUTO_BOTH$ant2.Movement_angle_difference - pi) - (TRAJ_AUTO_BOTH$ant1.Movement_angle_difference -pi)) %% (2*pi)
    mean_movement_angle_diff <-  mean(TRAJ_AUTO_BOTH$Movement_angle_diff, na.rm=TRUE)
    
    ###############
    summary_AUTO_INT <- data.frame(REPLICATE, PERIOD, INT, unique(TRAJ_AUTO_BOTH$ant1), unique(TRAJ_AUTO_BOTH$ant2), unique(TRAJ_AUTO_BOTH$pair),
                                  StDev_angle,
                                  mean_angle,
                                  mean_delta_angles, 
                                  moved_distance_px, 
                                  mean_speed_pxpersec, 
                                  mean_accel_pxpersec2, 
                                  mean_jerk_PxPerSec3, 
                                  rmsd_px,
                                  INT_start_frame_ANT, INT_end_frame_ANT, interaction_length_secs,
                                  #prop_time_undetected_ANT, prop_time_undetected_REC,
                                  mean_strghtline_dist_px,
                                  mean_orient_angle_diff, mean_movement_angle_diff,
                                  #when adding a new variable, it must be included in the reshape rule for data plotting
                                  stringsAsFactors = F)

    names(summary_AUTO_INT)[names(summary_AUTO_INT) == 'unique.TRAJ_AUTO_BOTH.ant1.'] <- "ant1"
    names(summary_AUTO_INT)[names(summary_AUTO_INT) == 'unique.TRAJ_AUTO_BOTH.ant2.'] <- "ant2"
    names(summary_AUTO_INT)[names(summary_AUTO_INT) == 'unique.TRAJ_AUTO_BOTH.pair.'] <- "pair"
    
    #stack per TRAJ_AUTO_BOTH
    interacts_AUTO_INT <- data.frame(REPLICATE,PERIOD,INT,TRAJ_AUTO_BOTH,
                                     stringsAsFactors = F)
    
    ## Plot trajectories of both actor & receiver, show on the same panel
    Title <- "Automatic Interactions"
    plot   (ant1.y ~ ant1.x, TRAJ_AUTO_BOTH, type="l", lwd=4, col="blue4",asp=1, main=Title,cex.main=0.9 ,xlim=c(min(TRAJ_AUTO_BOTH$ant1.x,TRAJ_AUTO_BOTH$ant2.x),max(TRAJ_AUTO_BOTH$ant1.x,TRAJ_AUTO_BOTH$ant2.x)),ylim=c(min(TRAJ_AUTO_BOTH$ant1.y,TRAJ_AUTO_BOTH$ant2.y),max(TRAJ_AUTO_BOTH$ant1.y,TRAJ_AUTO_BOTH$ant2.y))) #, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
    points (ant2.y~ ant2.x, TRAJ_AUTO_BOTH, type="l", lwd=4,  col="red4",asp=1)
    

    ## stack -by the end of the loop, this should just contain all events for replicate X, period Y
    interacts_AUTO_REP_PER_FULL <- rbind(interacts_AUTO_REP_PER_FULL, interacts_AUTO_INT)
    summary_AUTO_REP_PER   <- rbind(summary_AUTO_REP_PER,     summary_AUTO_INT)
    
}










