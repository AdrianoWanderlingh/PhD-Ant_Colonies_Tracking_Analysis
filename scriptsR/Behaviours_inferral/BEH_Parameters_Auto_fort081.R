######  GETTING VARS OFF  AUTO INTERACTIONS ######
print(paste("GETTING VARS OFF AUTO INTERACTIONS ",REPLICATE, PERIOD))
#IT HAS TO BE DONE INSIDE THE REP-PER LOOP BEFORE THE interaction_AUTO STACKING


# WHAT TO DO WHEN IT IS NOT POSSIBLE TO DISCRIMINATE BETWEEN ACTOR AND RECEIVER?

#es.
#Y <- it is possible to distinguish ACT and REC (CONTAINS EITHER 2-3 OR 3-2)
#N <- it is NOT possible to distinguish ACT and REC as the int is bidirectional
# > table(interacts_AUTO_REP_PER$interactions$types)
# 2-2,2-3          13 Y
# 2-2,2-3,3-2      10 N
# 2-2,2-3,3-2,3-3   7 N
# 2-2,2-3,3-3       4 Y
# 2-2,3-2          20 Y
# 2-2,3-2,3-3       1 Y
# 2-3             226 Y
# 2-3,3-2          12 N
# 2-3,3-2,3-3      19 N
# 2-3,3-3          94 Y
# 3-2             272 Y
# 3-2,3-3          78 Y
# 
# 
# 










#fmQueryComputeInteractions $interactions
#interacts_AUTO_REP_PER$interactions

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
    rmsd_px             <- NULL
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
      #ISSUE HERE FOR ANT 10 -> NO FRAME INFO
      #DEPENDS ON THE WAY FRAME INFO IS ASSIGNED
      TRAJ_ANT     <- interacts_AUTO_REP_PER$trajectories[[ANT.ROW.INDEX]]
      
      TRAJ_ANT_INT <- TRAJ_ANT [ which(as.numeric(rownames(TRAJ_ANT)) >= INT_FRAME_start & as.numeric(rownames(TRAJ_ANT)) <= INT_FRAME_end),]


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


      
      ##################
      ## INDIVIDUAL TRAJECTORY MEASURES
      
      ### angular st.dev 
      StDev_angle[which_ant]      <-  angular.deviation(TRAJ_ANT_INT$angle, na.rm = TRUE)
      #circular average
      mean_abs_angle[which_ant]   <-  mean.circular(abs(TRAJ_ANT_INT$angle), na.rm = TRUE)
      #turn angles
      trj <- TrajFromCoords(data.frame(TRAJ_ANT_INT$x,TRAJ_ANT_INT$y,TRAJ_ANT_INT$frame))
      mean_abs_turnAngle[which_ant]              <- mean.circular(abs(TrajAngles(trj)))

      
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
      rmsd_px[which_ant]                            <-  sqrt(sum( (TRAJ_ANT_INT$x-mean(TRAJ_ANT_INT$x))^2 + (TRAJ_ANT_INT$y-mean(TRAJ_ANT_INT$y))^2 )/length(na.omit(TRAJ_ANT_INT$x)))

      #Convex Hull
      box.coords <- as.matrix(TRAJ_ANT_INT[, c("x", "y")]); box.hpts <- chull(x = TRAJ_ANT_INT$x, y = TRAJ_ANT_INT$y) # calculate convex hull for x and y columns
      box.hpts <- c(box.hpts, box.hpts[1]) #add first coord at the bottom to close area
      box.chull.coords <- box.coords[box.hpts,]
      chull.poly <- Polygon(box.chull.coords, hole=F); chull_area[which_ant]   <- chull.poly@area  #calculate area
      
      
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
      moved_distance_px[which_ant]            <- sum(TRAJ_ANT_INT$distance,na.rm = T)
      
      TRAJ_ANT_INT$time_interval   <- c(NA, diff(TRAJ_ANT_INT$frame)/8)
      #speed
      TRAJ_ANT_INT$speed_PxPerSec  <- c( with(TRAJ_ANT_INT, c(NA,(sqrt(diff(x)^2 + diff(y)^2))) / time_interval))
      mean_speed_PxPerSec[which_ant]         <- mean(TRAJ_ANT_INT$speed_PxPerSec, na.rm=T) 
      #  acceleration
      TRAJ_ANT_INT$accel_PxPerSec2 <- c( with(TRAJ_ANT_INT, c(NA,(diff(speed_PxPerSec))) / time_interval))
      mean_accel_PxPerSec2[which_ant]         <- mean(TRAJ_ANT_INT$accel_PxPerSec2, na.rm=T)
      # jerk (diff in accelerations)
      TRAJ_ANT_INT$jerk_PxPerSec3  <- c( with(TRAJ_ANT_INT, c(NA,(diff(accel_PxPerSec2))) / time_interval))
      mean_jerk_PxPerSec3[which_ant]          <- mean( TRAJ_ANT_INT$jerk_PxPerSec3, na.rm=T)
  
  
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

    #rename trajectories columns NOT TO BE MERGED for Act and Rec, all except frame 
    names(TRAJ_ANT_ant1) <- paste0("ant1." ,names(TRAJ_ANT_ant1))
    names(TRAJ_ANT_ant2) <- paste0("ant2." ,names(TRAJ_ANT_ant2))
    names(TRAJ_ANT_ant1)[names(TRAJ_ANT_ant1) == 'ant1.frame'] <-"frame"
    names(TRAJ_ANT_ant2)[names(TRAJ_ANT_ant2) == 'ant2.frame'] <-"frame"
  
    #merge trajectories matching by time
    TRAJ_AUTO_BOTH <- merge(TRAJ_ANT_ant1, TRAJ_ANT_ant2,all=T, by=c("frame"))
    TRAJ_ANT_ant1 <- NULL
    TRAJ_ANT_ant2 <- NULL

    #TRAJ_AUTO_BOTH <- data.frame(TRAJ_ANT_list)
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
    
    interaction_length_secs <- ((max(TRAJ_AUTO_BOTH$frame) - min(TRAJ_AUTO_BOTH$frame))+1)/8  ## includes end frame with +1
    prop_time_undetected_ant1 <- (sum(is.na(TRAJ_AUTO_BOTH$ant1.x)) / 8) / interaction_length_secs  ## the prop of the interaction in which ACT was seen 
    prop_time_undetected_ant2 <- (sum(is.na(TRAJ_AUTO_BOTH$ant2.x)) / 8)  / interaction_length_secs ## UNITS are secs / secs --> proportion; MUST BE STRICTLY < 1 !!
    mean_prop_time_undetected <- (prop_time_undetected_ant1+prop_time_undetected_ant2)/2
    
    ###############
    summary_AUTO_INT <- data.frame(REPLICATE, PERIOD, INT, unique(TRAJ_AUTO_BOTH$ant1), unique(TRAJ_AUTO_BOTH$ant2), unique(TRAJ_AUTO_BOTH$pair),
                                  StDev_angle[["ant1"]],StDev_angle[["ant2"]],
                                  mean_abs_angle[["ant1"]],mean_abs_angle[["ant2"]],
                                  mean_abs_turnAngle[["ant1"]],mean_abs_turnAngle[["ant2"]],
                                  mean_delta_angles, 
                                  moved_distance_px[["ant1"]], moved_distance_px[["ant2"]],
                                  mean_speed_PxPerSec[["ant1"]],mean_speed_PxPerSec[["ant2"]],
                                  mean_accel_PxPerSec2[["ant1"]], mean_accel_PxPerSec2[["ant2"]], 
                                  mean_jerk_PxPerSec3[["ant1"]], mean_jerk_PxPerSec3[["ant2"]], 
                                  rmsd_px[["ant1"]],rmsd_px[["ant2"]],
                                  chull_area[["ant1"]],chull_area[["ant2"]],
                                  int_start_frame=min(TRAJ_AUTO_BOTH$frame), int_end_frame =max(TRAJ_AUTO_BOTH$frame),
                                  interaction_length_secs,  
                                  prop_time_undetected_ant1, prop_time_undetected_ant2, mean_prop_time_undetected,
                                  mean_strghtline_dist_px,
                                  mean_orient_angle_diff, mean_movement_angle_diff,
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










