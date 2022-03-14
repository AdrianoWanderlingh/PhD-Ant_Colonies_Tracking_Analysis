######  GETTING VARS OFF  AUTO INTERACTIONS ######
print(paste("GETTING VARS OFF AUTO INTERACTIONS ",REPLICATE, PERIOD))

# WHAT TO DO WHEN IT IS NOT POSSIBLE TO DISCRIMINATE BETWEEN ACTOR AND RECEIVER?
#USE ACT - REC definition based on speed


#fmQueryComputeInteractions $interactions
#interacts_AUTO_REP_PER$interactions

dim(interacts_AUTO_REP_PER$interactions)

#fmQueryComputeInteractions $trajectories

# CONSIDER THAT THE REC/ACT IS NOT KNOWN FOR AUTO_INTERACTS SO THE VARIABLES CONSIDERED SHOULD ONLY BE THE ONES NOT LINKED TO 1 SPECIFIC INDIVIDUAL. 
# USE FOR LOOP  TO CREATE, INTERACTION BY INTERACTION, THE LIST OF PARAMS X Y ANGLE WITH EXP INFOS ATTACHED

#let's create an empy data.frame that will be the interacts_AUTO_REP_PER but with trajectory info

# FINAL RESULT SHOULD LOOK LIKE interacts_MAN_REP_PER
for (INT in 1:nrow(interacts_AUTO_REP_PER$interactions)) {
    #TRAJ_ANT_list <- list()
    capsule_ANT         <- NULL
    ROLE                <- NULL
    StDev_angle         <- NULL
    mean_abs_angle      <- NULL
    mean_abs_turnAngle  <- NULL
    mean_sqrt_err_px    <- NULL
    chull_area          <- NULL
    moved_distance_px   <- NULL
    mean_speed_PxPerSec <- NULL
    mean_accel_PxPerSec2<- NULL
    mean_jerk_PxPerSec3 <- NULL
    mean_turnAngle      <- NULL
    for (which_ant in c("ant1","ant2")) {
      
      ## extract actor, receiver IDs & start & end times from the hand-annotated data
      #INTER <- interacts_AUTO_REP_PER$interactions[INT,]
      ANT <- interacts_AUTO_REP_PER$interactions[INT,which_ant]
      
      #extract frame length between start and end (int_start_frame ? int_end_frame)
      #INT_start_frame_ANT <- interacts_AUTO_REP_PER$interactions[,"int_start_frame"][INT]
      #INT_end_frame_ANT   <- interacts_AUTO_REP_PER$interactions[,"int_end_frame"]  [INT]
      
      # interacts_AUTO_REP_PER$trajectories_summary$ant.row.index <- paste0(interacts_AUTO_REP_PER$trajectories_summary$antID_str,"row",rownames(interacts_AUTO_REP_PER$trajectories_summary))
      # names(interacts_AUTO_REP_PER$trajectories)                <- paste0(interacts_AUTO_REP_PER$trajectories_summary$antID_str,"row",rownames(interacts_AUTO_REP_PER$trajectories_summary)) ###and use the content of that column to rename the objects within trajectory list
      # 
      
      ## extract the trajectory for ANT
      #trajectory.row
      INT_TRAJ_ROW_ANT <- interacts_AUTO_REP_PER$interactions[INT,paste0(which_ant,".trajectory.row")]
      ANT.ROW.INDEX <- paste0("ant_",ANT,"row",INT_TRAJ_ROW_ANT)
      #trajectory.start
      INT_FRAME_start  <- as.numeric(interacts_AUTO_REP_PER$interactions[INT,paste0(which_ant,".trajectory.start")])
      #trajectory.end
      INT_FRAME_end  <- as.numeric(interacts_AUTO_REP_PER$interactions[INT,paste0(which_ant,".trajectory.end")])
      
      print(paste("Interaction number",INT,which_ant,ANT,"ANT.ROW.INDEX",ANT.ROW.INDEX,"INT_FRAME_start",INT_FRAME_start,"INT_FRAME_end",INT_FRAME_end))
     
      # which_ant.trajectory.row the corresponding index in $trajectory.
      TRAJ_ANT     <- as.data.table(interacts_AUTO_REP_PER$trajectories[[ANT.ROW.INDEX]])
      
      TRAJ_ANT_INT <- TRAJ_ANT [ which(as.numeric(rownames(TRAJ_ANT)) >= INT_FRAME_start & as.numeric(rownames(TRAJ_ANT)) <= INT_FRAME_end),]
      
      TRAJ_ANT_INT<- cbind(ANT,TRAJ_ANT_INT)

      #interacting capsules
      capsules  <- e$antShapeTypeNames
      names(capsules) <- as.character( 1:length(capsules))
      head_id <- names(capsules)[[which(capsules=="head")]]
      body_id <- names(capsules)[[which(capsules=="body")]]
      INT_capsules  <- interacts_AUTO_REP_PER$interactions[INT,"types"]
      #get info only when there are only 2 interacting capsules
      if (!grepl(",",INT_capsules)) {
        capsule_ANT <- unlist(strsplit(INT_capsules,"-"))
        if (which_ant=="ant1") {
          #when body, assign REC_Role. when head, assign ACT_role
          if (capsule_ANT[1]==body_id) {
            ROLE[which_ant]  <- "REC"
          } else { ROLE[which_ant]  <- "ACT"}
        }
      
      
        if (which_ant=="ant2"){
          #when body, assign REC_Role. when head, assign ACT_role
          if (capsule_ANT[2]==body_id) {
            ROLE[which_ant]  <- "REC"
          } else { ROLE[which_ant]  <- "ACT"}
        }
      }#select only type with a single capsule interacting


      TRAJ_ANT_INT <- as.data.frame(TRAJ_ANT_INT)
      
      ################################################################
      ######### JUMPS ################################################
      ################################################################

      ##  jumps detection by frame
      TRAJ_ANT_INT$dt_FRAME            <- c(1, diff(TRAJ_ANT_INT$frame)) ## time difference in frames
      
      #jumps detection by distance is reported below in TRAJ_ANT_INT$distance

      # ##################
      # ## INDIVIDUAL TRAJECTORY MEASURES

      #turnangles
      trj <- TrajFromCoords(data.frame(TRAJ_ANT_INT$x,TRAJ_ANT_INT$y,TRAJ_ANT_INT$frame))
      TRAJ_ANT_INT$TurnAngle <- c(NA,NA,TrajAngles(trj))
      
      
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
      # - have a modified version of orientation_diff where -pi is subtracted from each coord as done in TRAJ_AUTO_BOTH$Orient_angle_diff (why is that done? test if it works beforehand!)
      
      ## delta_angles
      TRAJ_ANT_INT$delta_angles <- TRAJ_ANT_INT$Movement_angle_difference - TRAJ_ANT_INT$Orientation_diff
     
      ## root mean square displacement
      mean_sqrt_err_px[which_ant]                            <-  sum( (TRAJ_ANT_INT$x-mean(TRAJ_ANT_INT$x))^2 + (TRAJ_ANT_INT$y-mean(TRAJ_ANT_INT$y))^2 )/length(na.omit(TRAJ_ANT_INT$x))

      #Convex Hull
      box.coords <- as.matrix(TRAJ_ANT_INT[, c("x", "y")]); box.hpts <- chull(x = TRAJ_ANT_INT$x, y = TRAJ_ANT_INT$y) # calculate convex hull for x and y columns
      box.hpts <- c(box.hpts, box.hpts[1]) #add first coord at the bottom to close area
      box.chull.coords <- box.coords[box.hpts,]
      chull.poly <- Polygon(box.chull.coords, hole=F); chull_area[which_ant]   <- chull.poly@area  #calculate area
      
      ##### FRAME BY FRAME PARAMETERS

      # distance walked
      TRAJ_ANT_INT$distance           <- c(NA, with(TRAJ_ANT_INT, (sqrt(diff(x)^2 + diff(y)^2)))) # euclidean distance
      moved_distance_px[which_ant]    <- sum(TRAJ_ANT_INT$distance,na.rm = T)
      
      TRAJ_ANT_INT$time_interval      <- c(NA, diff(TRAJ_ANT_INT$frame)/8)
      #speed
      TRAJ_ANT_INT$speed_PxPerSec     <- c( with(TRAJ_ANT_INT, c(NA,(sqrt(diff(x)^2 + diff(y)^2))) / time_interval))
      #  acceleration
      TRAJ_ANT_INT$accel_PxPerSec2    <- c( with(TRAJ_ANT_INT, c(NA,(diff(speed_PxPerSec))) / time_interval))
      # jerk (diff in accelerations)
      TRAJ_ANT_INT$jerk_PxPerSec3     <- c( with(TRAJ_ANT_INT, c(NA,(diff(accel_PxPerSec2))) / time_interval))
  
      #names(TRAJ_ANT_INT)[names(TRAJ_ANT_INT) == 'x'] <- paste0(which_ant,".x"); names(TRAJ_ANT_INT)[names(TRAJ_ANT_INT) == 'y'] <- paste0(which_ant,".y"); names(TRAJ_ANT_INT)[names(TRAJ_ANT_INT) == 'angle'] <- paste0(which_ant,".angle")
      #assign each which_ant to list
      if (which_ant=="ant1") {
        TRAJ_ANT_ant1 <-TRAJ_ANT_INT
      }
      if (which_ant=="ant2"){
        TRAJ_ANT_ant2 <-TRAJ_ANT_INT
      #TRAJ_ANT_list[[which_ant]] <- TRAJ_ANT_INT
      }
  }#which_ant
  
    #-----------------------------------------------------------
    # #assign basing on ROLE
    # if (which_ant=="ant1") {
    #   TRAJ_ANT_ant1 <-TRAJ_ANT_INT
    # }
    # if (which_ant=="ant2"){
    #   TRAJ_ANT_ant2 <-TRAJ_ANT_INT
    # }

    
    #--------------------
    #assign ACT REC id, if 1 is faster, it is the ACT
    if ((mean(TRAJ_ANT_ant1$speed_PxPerSec, na.rm=T) - mean(TRAJ_ANT_ant2$speed_PxPerSec, na.rm=T))>0 ) {
      TRAJ_ANT_ACT <- TRAJ_ANT_ant1
      TRAJ_ANT_REC <- TRAJ_ANT_ant2
    }else{
      TRAJ_ANT_ACT <- TRAJ_ANT_ant2
      TRAJ_ANT_REC <- TRAJ_ANT_ant1
    }
    
    #--------------------
    
    
    #rename trajectories columns NOT TO BE MERGED for Act and Rec, all except frame 
    names(TRAJ_ANT_ACT) <- paste("ACT" ,names(TRAJ_ANT_ACT), sep = ".")
    names(TRAJ_ANT_REC) <- paste("REC" ,names(TRAJ_ANT_REC), sep = ".")
    names(TRAJ_ANT_ACT)[names(TRAJ_ANT_ACT) == 'ACT.frame'] <-"frame"
    names(TRAJ_ANT_REC)[names(TRAJ_ANT_REC) == 'REC.frame'] <-"frame"
  
    #merge trajectories matching by frame
    TRAJ_AUTO_BOTH <- merge(TRAJ_ANT_ACT, TRAJ_ANT_REC,all=T, by=c("frame"))
    TRAJ_ANT_ACT <- NULL
    TRAJ_ANT_REC <- NULL
    
    

    #TRAJ_AUTO_BOTH <- data.frame(TRAJ_ANT_list)
    TRAJ_AUTO_BOTH[c("ACT.time","REC.time","REC.zone","REC.frame")] <- list(NULL)
    names(TRAJ_AUTO_BOTH)[names(TRAJ_AUTO_BOTH) == 'ACT.zone'] <- "zone"; names(TRAJ_AUTO_BOTH)[names(TRAJ_AUTO_BOTH) == 'ACT.frame'] <- "frame"
    #add info from INT
    TRAJ_AUTO_BOTH <- merge(interacts_AUTO_REP_PER$interactions[INT,c("ant1","ant2","pair")],TRAJ_AUTO_BOTH)
    #col_order <- c("frame","ant1","ant2","pair","ant1.x","ant1.y","ant1.angle","ant2.x","ant2.y","ant2.angle","zone")
    #TRAJ_AUTO_BOTH <- TRAJ_AUTO_BOTH[, col_order]

    ##################
    ## INTERACTING PAIR TRAJECTORY MEASURES
    
    #straight line - euclidean distance
    TRAJ_AUTO_BOTH$straightline_dist_px <-  sqrt((TRAJ_AUTO_BOTH$ant1.x-TRAJ_AUTO_BOTH$ant2.x)^2+(TRAJ_AUTO_BOTH$ant1.y-TRAJ_AUTO_BOTH$ant1.y)^2)
    #angular difference
    TRAJ_AUTO_BOTH$Orient_angle_diff <- abs((TRAJ_AUTO_BOTH$ant2.angle - pi) - (TRAJ_AUTO_BOTH$ant1.angle -pi)) %% (2*pi)
    
    TRAJ_AUTO_BOTH$Movement_angle_diff <- abs((TRAJ_AUTO_BOTH$ant2.Movement_angle_difference - pi) - (TRAJ_AUTO_BOTH$ant1.Movement_angle_difference -pi)) %% (2*pi)
    
    # angular velocity: ω = (α₂ - α₁) / t = Δα / t
    #  FIX FIX  the time interval should be only 1 per interacting pair, correct last term of equation
    TRAJ_AUTO_BOTH$ang_Velocity        <- (TRAJ_AUTO_BOTH$ant2.angle - TRAJ_AUTO_BOTH$ant1.angle) / TRAJ_AUTO_BOTH$ant1.time_interval
    
    ## measure the length *in seconds* of the interaction
    int_length_secs <- ((max(TRAJ_AUTO_BOTH$frame) - min(TRAJ_AUTO_BOTH$frame))+1)/8  ## includes end frame with +1
    prop_time_undetected_ant1 <- (sum(is.na(TRAJ_AUTO_BOTH$ant1.x)) / 8) / int_length_secs  ## the prop of the interaction in which ACT was seen 
    prop_time_undetected_ant2 <- (sum(is.na(TRAJ_AUTO_BOTH$ant2.x)) / 8)  / int_length_secs ## UNITS are secs / secs --> proportion; MUST BE STRICTLY < 1 !!
    #mean_prop_time_undetected <- (prop_time_undetected_ant1+prop_time_undetected_ant2)/2

    ##################################
    ## mean values for summary #######
    ##################################
    
    ## Apply THRESHOLDs to exclude gaps in which the individuals were not detected for very long
    
    # # remove traj points when Trajectory$dt_FRAME > DT_frame_THRESHOLD
    # # dt_FRAME corresponds to frame between T and T+1. This should be applied to all movement characteristics
    # TRAJ_AUTO_BOTH[which(TRAJ_AUTO_BOTH$ant1.dt_FRAME>DT_frame_THRESHOLD), c("ACT.speed_PxPerSec","ACT.accel_PxPerSec2","ACT.jerk_PxPerSec3","ACT.TurnAngle","ang_Velocity")] <- NA
    # TRAJ_AUTO_BOTH[which(TRAJ_AUTO_BOTH$ant2.dt_FRAME>DT_frame_THRESHOLD), c("REC.speed_PxPerSec","REC.accel_PxPerSec2","REC.jerk_PxPerSec3","REC.TurnAngle","ang_Velocity")] <- NA
    # # remove traj points when TRAJ_AUTO_BOTH$ACT.distance$distance < DT_dist_THRESHOLD
    # TRAJ_AUTO_BOTH[which(TRAJ_AUTO_BOTH$ant1.distance<DT_dist_THRESHOLD), c("ACT.speed_PxPerSec","ACT.accel_PxPerSec2","ACT.jerk_PxPerSec3","ACT.TurnAngle","ang_Velocity")] <- NA
    # TRAJ_AUTO_BOTH[which(TRAJ_AUTO_BOTH$ant2.distance<DT_dist_THRESHOLD), c("REC.speed_PxPerSec","REC.accel_PxPerSec2","REC.jerk_PxPerSec3","REC.TurnAngle","ang_Velocity")] <- NA
    # 

    ### angular st.dev 
    # THIS SHOULD BE SOMEHOW INCORPORATE THE JUMPS CUTTING PROCEDURE (but I can't really just cut the coords x and y)
    StDev_angle_ant1    <-  angular.deviation(TRAJ_AUTO_BOTH$ant1.angle, na.rm = TRUE)
    StDev_angle_ant2    <-  angular.deviation(TRAJ_AUTO_BOTH$ant2.angle, na.rm = TRUE)
    
    #turnangle stDev
    stDev_turnAngle_ant1                 <-  angular.deviation(TRAJ_AUTO_BOTH$ant1.TurnAngle, na.rm = TRUE)
    stDev_turnAngle_ant2                 <-  angular.deviation(TRAJ_AUTO_BOTH$ant2.TurnAngle, na.rm = TRUE)
    
    #means are then calculated on the reduced dataset
    #individual
    #  mean_abs_turnAngle[which_ant]              <- mean.circular(abs(TRAJ_ANT_INT$ant1.TurnAngle))
    mean_abs_turnAngle_ant1   <- mean.circular(abs(TRAJ_AUTO_BOTH$ant1.TurnAngle),na.rm=T)
    mean_abs_turnAngle_ant2   <- mean.circular(abs(TRAJ_AUTO_BOTH$ant2.TurnAngle),na.rm=T)
    ##UPDATE ALSO THE FOLLOWING ONES....
    # mean_speed_PxPerSec[which_ant]  <- mean(TRAJ_ANT_INT$speed_PxPerSec, na.rm=T) 
    # mean_accel_PxPerSec2[which_ant] <- mean(TRAJ_ANT_INT$accel_PxPerSec2, na.rm=T)
    # mean_jerk_PxPerSec3[which_ant]  <- mean( TRAJ_ANT_INT$jerk_PxPerSec3, na.rm=T)
    # 
    
    #pair
    mean_strghtline_dist_px <- mean(TRAJ_AUTO_BOTH$straightline_dist_px, na.rm=TRUE)
    mean_ang_Velocity             <- mean(abs(TRAJ_AUTO_BOTH$ang_Velocity), na.rm=T)
    
    #TO FIXXXXX!!!!!!!!!!!!!!!!
    mean_orient_angle_diff <-  mean(TRAJ_AUTO_BOTH$Orient_angle_diff, na.rm=TRUE)
    mean_movement_angle_diff <-  mean(TRAJ_AUTO_BOTH$Movement_angle_diff, na.rm=TRUE)
  
    #circular average
    #mean_abs_angle[which_ant]   <-  mean.circular(abs(TRAJ_ANT_INT$angle), na.rm = TRUE)
    
    #Remove this var?? there is also a bit of confusion between Orientation_diff and Orient_angle_diff
    ## mean delta_angles
    mean_delta_angles <- mean.circular(TRAJ_ANT_INT$delta_angles,na.rm=TRUE)
    
    
    
    
    
    
    ###############
    summary_AUTO_INT <- data.frame(REPLICATE, PERIOD, INT, unique(TRAJ_AUTO_BOTH$ant1), unique(TRAJ_AUTO_BOTH$ant2), unique(TRAJ_AUTO_BOTH$pair),
                                  StDev_angle_ant1,StDev_angle_ant2, # StDev_angle[["ant1"]],StDev_angle[["ant2"]],
                                  mean_abs_turnAngle_ant1,mean_abs_turnAngle_ant2, #mean_abs_turnAngle[["ant1"]],mean_abs_turnAngle[["ant2"]],
                                  stDev_turnAngle_ant1,stDev_turnAngle_ant1,
                                ## mean_delta_angles_ACT, mean_delta_angles_REC,  
                                  mean_delta_angles, 
                                  moved_distance_px[["ant1"]], moved_distance_px[["ant2"]],
                                  # mean_speed_PxPerSec[["ant1"]],mean_speed_PxPerSec[["ant2"]],
                                  # mean_accel_PxPerSec2[["ant1"]], mean_accel_PxPerSec2[["ant2"]], 
                                  # mean_jerk_PxPerSec3[["ant1"]], mean_jerk_PxPerSec3[["ant2"]], 
                                  mean_sqrt_err_px[["ant1"]],mean_sqrt_err_px[["ant2"]],
                                  chull_area[["ant1"]],chull_area[["ant2"]],
                                  int_start_frame=min(TRAJ_AUTO_BOTH$frame), int_end_frame =max(TRAJ_AUTO_BOTH$frame),
                                  int_length_secs,  
                                  prop_time_undetected_ant1, prop_time_undetected_ant2, #mean_prop_time_undetected,
                                  mean_strghtline_dist_px,
                                  mean_orient_angle_diff, mean_movement_angle_diff,
                                  mean_ang_Velocity,
                                  #when adding a new variable, it must be included in the reshape rule for data plotting
                                  stringsAsFactors = F)

    names(summary_AUTO_INT)[names(summary_AUTO_INT) == 'unique.TRAJ_AUTO_BOTH.ant1.'] <- "ant1"
    names(summary_AUTO_INT)[names(summary_AUTO_INT) == 'unique.TRAJ_AUTO_BOTH.ant2.'] <- "ant2"
    names(summary_AUTO_INT)[names(summary_AUTO_INT) == 'unique.TRAJ_AUTO_BOTH.pair.'] <- "pair"
    colnames(summary_AUTO_INT) <- sub("...ant1...", "_ant1", colnames(summary_AUTO_INT))
    colnames(summary_AUTO_INT) <- sub("...ant2...", "_ant2", colnames(summary_AUTO_INT))
    
    
    #stack per TRAJ_AUTO_BOTH
    interaction_AUTO_INT <- data.frame(REPLICATE,PERIOD,INT,TRAJ_AUTO_BOTH,
                                     stringsAsFactors = F)
    
    ## Plot trajectories of both actor & receiver, show on the same panel
    #Title <- "Automatic Interactions"
    #plot   (ant1.y ~ ant1.x, TRAJ_AUTO_BOTH, type="l", lwd=4, col="blue4",asp=1, main=Title,cex.main=0.9 ,xlim=c(min(TRAJ_AUTO_BOTH$ant1.x,TRAJ_AUTO_BOTH$ant2.x),max(TRAJ_AUTO_BOTH$ant1.x,TRAJ_AUTO_BOTH$ant2.x)),ylim=c(min(TRAJ_AUTO_BOTH$ant1.y,TRAJ_AUTO_BOTH$ant2.y),max(TRAJ_AUTO_BOTH$ant1.y,TRAJ_AUTO_BOTH$ant2.y))) #, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
    #points (ant2.y~ ant2.x, TRAJ_AUTO_BOTH, type="l", lwd=4,  col="red4",asp=1)
    

    ## stack -by the end of the loop, this should just contain all events for replicate X, period Y
    if (create_interaction_AUTO_REP_PER) {
      interactions_AUTO_REP_PER <- rbind(interactions_AUTO_REP_PER, interaction_AUTO_INT)
    }
    summary_AUTO_REP_PER   <- rbind(summary_AUTO_REP_PER,     summary_AUTO_INT)
    
}










