##########################################################################################
############## BEH COLLISIONS ############################################################
##########################################################################################

#script dependant on BEH_MAIN_behaviours_analysis_fort081.R

#For previous versions of this script, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral

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
#rename frame for the merging
collisions_frames$frame <- collisions_frames$frame_num        ; collisions_frames$frame_num <- NULL
collisions_coll  $frame <- collisions_coll  $frames_row_index ; collisions_coll  $frames_row_index <- NULL
  
#merge collisions_frame and coll_coll based on index
collisions_merge <- merge(collisions_coll, collisions_frames[,-match(c("height","width"),colnames(collisions_frames))], by="frame")

#check that to the same frame corresponds the same time
# format(collisions_merge$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
# format(collisions_frames$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
# format(    traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
# format(    traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
# format(    traj_BOTH$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
# format(interacts_MAN_ROW$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
# 
# table(substr(x =  format( traj_BOTH$UNIX_time, "%OS6"), start = 7, stop = 7))


##if(collisions_coll$frames_row_index[125] == collisions_merge$frames_row_index[125]) {print("TRUE")}

# # add new column to interaction data when conditions are met (eg. time + ant ID)
# #RENAME COLUMN TRAJBOTHREC IN time 
# #merge with equal time and any order of the couple Act_Name-Rec_Name=ant1_str-ant2_str (any order). HOW?
# nrow(collisions_merge)
# str(collisions_merge)
# nrow(interacts_MAN_REP_PER)
# str(interacts_MAN_REP_PER)

#create new variable by pasting ant numbers "low,high" for collisions_merge
# collisions_merge$pair <- apply(collisions_merge[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })
collisions_merge $pair <- paste(collisions_merge$ant1, collisions_merge$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...

#create new variable by pasting ant numbers "low,high" for interacts_MAN_REP_PER
interacts_MAN_REP_PER$ant1 <- as.numeric(gsub("ant_","", interacts_MAN_REP_PER$Act_Name)) ## CRUCIAL TO NSURE THESE ARE CONVERTED TO NUMERIC BEFORE THE NEXT OPERATION TO PRODUCE interacts_MAN_REP_PER$pair, WHICH RELIES ON A **SORT** function!!!
interacts_MAN_REP_PER$ant2 <- as.numeric(gsub("ant_","", interacts_MAN_REP_PER$Rec_Name))
## the ant 1 & ant 2 labels are not strictly ascending, so sort them so ant 1 alwasy < ant 2
interacts_MAN_REP_PER$pair <- apply(interacts_MAN_REP_PER[,c("ant1","ant2")],1,function(x){paste(sort(x),collapse = "_") })

# # check that the time formats are the same
# attributes(interacts_MAN_REP_PER$UNIX_time[1])
# attributes(collisions_merge$UNIX_time[1])

# check that the frame formats are the same
interacts_MAN_REP_PER$frame[1]
collisions_merge$frame     [1]


## check that the pair-time combinations in interacts_MAN_REP_PER are in collisions_merge
table(  paste(interacts_MAN_REP_PER$pair) %in% 
          paste(collisions_merge$pair) )
table(  paste(interacts_MAN_REP_PER$frame) %in% 
          paste(collisions_merge$frame) )
table(  paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$frame) %in% 
          paste(collisions_merge$pair,collisions_merge$frame) )
# table(  paste(interacts_MAN_REP_PER$UNIX_secs) %in% 
#           paste(collisions_merge$UNIX_secs) )
# table(  paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$UNIX_secs) %in% 
#           paste(collisions_merge$pair,collisions_merge$UNIX_secs) )

## get a list of the missing pair-frame combinations:
paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$frame) [!paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$frame) %in% paste(collisions_merge$pair,collisions_merge$frame) ]

## get a list of the missing pair-time combinations:
# paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$UNIX_secs) [!paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$UNIX_secs) %in% paste(collisions_merge$pair,collisions_merge$UNIX_secs) ]

## assign the collision types to interacts_MAN_REP_PER, using a loop to subset interacts_MAN_REP_PER$pair by each uique pair to ensure we are matching within pairs 
interacts_MAN_REP_PER$types <- NA
PLOT_OVERLAP <- TRUE
if (PLOT_OVERLAP) {par(mai=c(0.1,0.1,0.17,0), mfrow=c(3,4))}
for (PR in unique(interacts_MAN_REP_PER$pair))
  { print(paste("Adding collision capsule types from collisions_merge to interacts_MAN_REP_PER for pair",PAIR))
  
  ##  sort the entire data frame subset for pair
  interacts_MAN_REP_PER_PAIR <- interacts_MAN_REP_PER [which(interacts_MAN_REP_PER$pair==PR),]
  collisions_merge_PAIR      <- collisions_merge      [which(collisions_merge$pair==PR),]

  ## alternative: fuzzy match to one frame
  if (FUZZY_MATCH==TRUE)
    {
    I_in_C_rows1 <- match.closest(x =    interacts_MAN_REP_PER_PAIR$frame,
                                 table = collisions_merge_PAIR     $frame, tolerance = 1)
    ## and add back to interacts_MAN_REP_PER
    interacts_MAN_REP_PER$types [which(interacts_MAN_REP_PER$pair==PR)] <- collisions_merge_PAIR$types [I_in_C_rows1]
    ## show the auto (BLACK) - manual (RED) overlap
    ## OBservation: overall, the black-red corespondence is good, but the auto collisions are too sparse!
    if (PLOT_OVERLAP) 
      {
      plot  (x=collisions_merge_PAIR$frame,      y=rep(0,nrow(collisions_merge_PAIR)), pch="+", main=PR, ylab="", xlab="") ## auto = black
      points(x=interacts_MAN_REP_PER_PAIR$frame, y=rep(0.5,nrow(interacts_MAN_REP_PER_PAIR)), pch="+", col=2) ## man=red
      }
    }
  }##PR


## ADRIANO TO DO: WRAP EVERYTHIN IN AN OUTER LOOP TO PROGRESSIVELY INCREASE THE SIZE OF THE CAPSULES , 
## SO THAT THERE ARE FEWER MISSING (NA) MATCHES: 
#table(interacts_MAN_REP_PER$types, useNA="al")
  

  # rowIndex <- I_in_C_rows; type_match <- collisions_merge$types[I_in_C_rows]
  # type_DF <- as.data.frame(rowIndex,type_match); type_DF$type_match <- rownames(type_DF)
  # rownames(type_DF) <- rowIndex
  # ## copy types across
  # #paste(interacts_MAN_REP_PER$pair,interacts_MAN_REP_PER$frame) %in% paste(collisions_merge$pair,collisions_merge$frame)
  # 
  # merge(traj_ACT, traj_REC,all=T, by=c("UNIX_time"))
  # 
  # 
  # 
  # collisions_merge$types_match <- ifelse(rownames(collisions_merge) == I_in_C_rows,  collisions_merge$types, NA)
  # 
  # collisions_merge$types_match <- match(rownames(collisions_merge), I_in_C_rows)
  # 
  # collisions_merge$types_match <- NA
  # collisions_merge$types_match[I_in_C_rows] <- collisions_merge$types [I_in_C_rows]
  # 
  # collisions_merge$types_match <- collisions_merge$types [I_in_C_rows]
  #collisions_merge_join <- collisions_merge[,c("frame","pair","types")]

  # 
  # interacts_MAN_REP_PER <- plyr::join(x = interacts_MAN_REP_PER, 
  #                                     y = collisions_merge[,c("frame","pair","types_match")], 
  #                                     by=c("frame","pair"), type="left")
#}

# #FUNCTION: WHEN pair=pair AND time=UNIX_time -> ASSIGN collisions_merge$types TO interacts_MAN_REP_PER
# if (FUZZY_MATCH==FALSE)
# {
#   interacts_MAN_REP_PER <- plyr::join(x = interacts_MAN_REP_PER,
#                                       y = collisions_merge[,c("UNIX_secs","pair","types")],
#                                       by=c("UNIX_secs","pair"), type="left")
# }
# # alternative: fuzzy match
# if (FUZZY_MATCH==TRUE)
# {
#   I_in_C_rows <- match.closest(x = interacts_MAN_REP_PER$UNIX_secs,
#                                table = collisions_merge$UNIX_secs, tolerance = 0.125)
#   ## copy tpes across
#   interacts_MAN_REP_PER$types  <- collisions_merge$types [I_in_C_rows]
# }

## PROBLEM: *MANY* rows in interactions aren't present in collisions; plot the spatial interactions to see whether this makes sense or not:
# par(mfrow=c(3,4), mai=c(0.3,0.3,0.4,0.1))
# for(BH in unique(interacts_MAN_REP_PER$BEH))
# {
#   for (RW in unique(interacts_MAN_REP_PER$ROW))
#   {
#     IntPair <- interacts_MAN_REP_PER[which(interacts_MAN_REP_PER$BEH==BH & interacts_MAN_REP_PER$ROW==RW),]
#     ## check if this interaction is present in collisions
#     paste(IntPair$pair[1],IntPair)
#     ## plot it
#     #to fix as some of the vars are defined only in the previous function #TitleInt <- paste(REPLICATE, ", ", PERIOD, ", ", BH, RW,", ", "Act:",ACT, ", ", "Rec:",REC, "\n", ENC_TIME_start, "-", ENC_TIME_stop, sep="")
#     plot(NA, xlim=range(c(IntPair$ACT.x,IntPair$REC.x),na.rm=T),  ylim=range(c(IntPair$ACT.y,IntPair$REC.y),na.rm=T), type="n", main=paste(RW,BH), xlab="x", ylab="y")
#     points (ACT.y ~ ACT.x, IntPair, type="p", col="blue4"); lines (ACT.y ~ ACT.x, IntPair, col="blue4")
#     points (REC.y ~ REC.x, IntPair, type="p", col="red4"); lines (REC.y ~ REC.x, IntPair, col="red4")
#     ## show the headings of each ACT
#     arrows.az (x = IntPair$ACT.x, 
#                y = IntPair$ACT.y, 
#                azimuth = IntPair$ACT.angle, 
#                rho = 10,
#                HeadWidth=0.1,
#                units="radians", Kol="blue2", Lwd=1)
#     ## show the headings of each REC
#     arrows.az (x = IntPair$REC.x, 
#                y = IntPair$REC.y, 
#                azimuth = IntPair$REC.angle, 
#                rho = 10,
#                HeadWidth=0.1,
#                units="radians", Kol="red2", Lwd=1)
#     
#     ## add lines connecting ACT & REC when there is a capsule overlap 
#     IntPairCapOverlap  <- IntPair[which(!is.na(IntPair$types)),]
#     if (nrow(IntPairCapOverlap)>0)
#     {
#       segments(x0 = IntPairCapOverlap$ACT.x, y0 = IntPairCapOverlap$ACT.y,
#                x1 = IntPairCapOverlap$REC.x, y1 = IntPairCapOverlap$REC.y)
#     }
#     
#   }
# }
#dev.off()

## if the missing rows occur when the distance between Act & Rec is large, check:
# interacts_MAN_REP_PER$types_present_absent <- "absent"
# interacts_MAN_REP_PER$types_present_absent[which(!is.na(interacts_MAN_REP_PER$types))] <- "present"
# boxplot(straightline_dist_px ~ types_present_absent, interacts_MAN_REP_PER)
### ... so rows in interactions that don't match collisions is not  a function of distance ...






#get info only when there are only 2 interacting capsules
#(!grepl(",",INT_capsules)) 
 # capsule_ANT <- unlist(strsplit(INT_capsules,"-"))

  
  
  
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
