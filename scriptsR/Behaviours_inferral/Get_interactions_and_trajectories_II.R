rm(list=ls())

## TO DO ; change 'PP' column name - to E.G. PERIOD & make alterations throughout!
## TO DO ; start to calculate trajectory statistics! 

#"https://formicidae-tracker.github.io/myrmidon/docs/latest/api/index.html"

######load necessary libraries
library(adehabitatHR) ####for home range and trajectory analyses
library(FortMyrmidon) ####R bindings
library(igraph)       ####for network analysis
library(parsedate)
library (trajr)
library(plotrix) 


# ## TEMPORARY - Hard code plotrix function - Adriano to do; install plotrix!
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
Ymin <- 2000
Ymax <- 5500

max_gap         <- fmHour(24*365)   ## important parameter to set! Set a maximumGap high enough that no cutting will happen. For example, set it at an entire year: fmHour(24*365)
desired_step_length_time    <- 0.125 ###in seconds, the desired step length for the analysis

##################################################
###### 1. OPENING AN EXPERIMENT ##################
##################################################
####### navigate to folder containing myrmidon file
if (USER=="Adriano") {WORKDIR <- "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/Subset_for_Tom"}
if (USER=="Tom")     {WORKDIR <- "/media/tom/MSATA/Dropbox/Ants_behaviour_analysis"}

DATADIR <- paste(WORKDIR,"Data",sep="/")
SCRIPTDIR <- paste(WORKDIR,"ScriptsR",sep="/")

##
source(paste(SCRIPTDIR,"Annotations_analysis.R",sep="/"))

#######################

for (REPLICATE in c("R3SP"))#,"R9SP")) TEMPORARY
  {
  ###############################################################################
  ###### 3. QUERYING GENERAL EXPERIMENT/ANT INFORMATION #########################
  ###############################################################################
  
  ## locate the ant info file for REPLICATE
  MyrmidonCapsuleFile <- list.files(path=DATADIR, pattern=REPLICATE, full.names=T); MyrmidonCapsuleFile <- MyrmidonCapsuleFile[grepl("Capsule_defined.myrmidon",MyrmidonCapsuleFile)]
  e <- fmExperimentOpenReadOnly(MyrmidonCapsuleFile)
  e$getDataInformations()
  ###tag statistics
  tag_stats <- fmQueryComputeTagStatistics(e)
  

  for (PP in c("pre","post"))
    {
    print(paste("Replicate",REPLICATE, PP))


  ###############################################################################
  ###### 4. READING TRAJECTORIES ################################################
  ###############################################################################

  ## TO DO : CHANGE THIS START TIME DEPENDING ON THE REPLICATE - NOT SAME FOR R3 & R9 !!!
  
  # time_start_ISO <- parse_iso_8601("2021-03-16T12:13:21.670072093Z")  ### OLD: TO DO - CHECK THIS ADRIANO
  StartTime      <- min(annotations$T_start[annotations$treatment_rep==REPLICATE & annotations$period==PP])
  time_start_ISO <- parse_iso_8601(StartTime)
  time_start     <- fmTimeCPtrFromAnySEXP(time_start_ISO)
  time_stop      <- fmTimeCPtrFromAnySEXP(time_start_ISO + (32*60) ) ####arbitrary time in the correct format + 16min 
  
  ######### Subset 15 min block
  positions <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = FALSE) #set true to obtain the zone of the ant
  ######### no subsetting - take the entire tracking period
  # positions <- fmQueryComputeAntTrajectories(e,start = e$getDataInformations()$details$tdd.start[1],end = e$getDataInformations()$details$tdd.end[3],maximumGap = max_gap,computeZones = FALSE) #set true to obtain the zone of the ant
  # ?fmQueryComputeAntTrajectories
  # 
  # positions$trajectory_summary
  # positions$trajectories
  ###Let's have a look at the first element of the trajectories list within positions:
  # str(positions$trajectories)
  # str(positions$trajectory_summary)
  #plot(positions$trajectories[[20]]$x,positions$trajectories[[20]]$y)
  
  ###that is a dataframe with columns "time" (in seconds since start), x, y (coordinates in pixels), angle (orientation of the ant in radians), zone (which zone the ant is)
  
  ###WARNING!!!!!
  ###in this particular example the values in antID are consecutive (1 to 22), and so positions$trajectories[[1]] will correspond to the trajectory of ant 1
  ### but imagine you had to delete an ant in fort-studio and the antID list "jumped" from 10 to 12 (so it would be 1,...10 then 12...23)
  ### then the 22nd object of trajectory list would not be the trajectory of ant 22 but of ant 23. That is not fool proof, but very risky!
  ##so what I suggest you do immediately after computing your positions object:
  positions$trajectory_summary$antID_str <- paste("ant_",positions$trajectory_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
  names(positions$trajectories)       <- positions$trajectory_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
  ##and now we can extract each ant trajectory using the character ID stored in antID_str - it's much safer!
  
  # #unlisting the trajectories 
  # library(data.table)
  # trajectories_unlist <- as.data.frame(rbindlist(positions$trajectories,use.names=TRUE,idcol=TRUE))
  # write.csv(trajectories_unlist, file=paste(e$getDataInformations()[["details"]][["tdd.URI"]],'trajectories_30min.csv'),row.names = FALSE) #specify name of the block better!
  
  trajectory_summary <- positions$trajectory_summary
  
  #generate reproducible example
  #dput(positions, file = "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_positions.txt")
  #positions_dget <- dget("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_positions.txt") # load file created with dput 
  
  
  ## Add ant_x names and times to the positions to convert from time since the start of the experiment, to UNIX time
  for (A in positions$trajectory_summary$antID)
    {
    AntID <- paste("ant_",A,sep="") 
    First_Obs_Time <- positions$trajectory_summary$start [which(positions$trajectory_summary$antID_str==AntID)] ## find the first time after the user defined time_start_ISO that this ant was seen
    print(paste("Adding first obs time", First_Obs_Time, "to the time-zeroed trajectory of ant", AntID))
    positions$trajectories[[AntID]] $UNIX_time <- positions$trajectories[[AntID]] $time + First_Obs_Time ##convert back to UNIX time  
    }
  
  #check that milliseconds are preserved after adding the ant's time_start_ISO
  positions$trajectories[["ant_27"]][[2,"UNIX_time"]]-positions$trajectories[["ant_27"]][[1,"UNIX_time"]]
  
  #Conversion of data in mm, probably takes ages over the full dataset
  # reuse this and tags size in pixel/mm as point of reference - traj[,which(grepl("x",names(traj))|grepl("y",names(traj)))] <-   traj[,which(grepl("x",names(traj))|grepl("y",names(traj)))] * petridish_diameter / range 
  
  
  # subset of annotations to run tests
  # annotations <- annotations[which(annotations$period=="post" & annotations$treatment_rep == "R3SP"),]
  
  ############Prepare overall data object
  summary_data <- NULL
  
  ####First let's extract ant's trajectories
  par(mfrow=c(1,5), mai=c(0.3,0.3,0.2,0.1), mgp=c(1.3,0.3,0), family="serif", tcl=-0.2)
  
  for (BEH in c("G","T"))
    {
    ## subset all hand-labelled bahavs for this behaviour type in this colony
    annot_BEH <- annotations[which(annotations$treatment_rep==REPLICATE  & annotations$period==PP & annotations$Behaviour==BEH ),]
    ## remove doubled allo-grooming interactions
    if (BEH=="G") {print("duplicates removed")}#{annot_BEH <- annot_BEH[!duplicated(annot_BEH),]}  ## leave NOT to catch possible un-matched rows
    if (BEH=="T") {print("duplicates removed")}
    ## loop through each event in annot_BEH
    for (ROW in 1:nrow(annot_BEH))
      { 
      ## extract actor, receiver IDs & start & end times from the hand-annotated data
      ACT <- annot_BEH$Actor[ROW]
      REC <- annot_BEH$Receiver[ROW]
      ENC_TIME_start <- annot_BEH$T_start[ROW]
      ENC_TIME_stop  <- annot_BEH$T_stop[ROW]
      
      Act_Name <- paste("ant",ACT,sep="_")
      Rec_Name <- paste("ant",REC,sep="_")
      print(paste("Behaviour:",BEH,"number",ROW,"Actor:",Act_Name,"Receiver:",Rec_Name))
      
      ## extract the trajectory for ACT
      traj_ACT <-  positions$trajectories[[Act_Name]]
      traj_REC <-  positions$trajectories[[Rec_Name]]
      
      Title <- paste(REPLICATE, ", ", PP, ", ", BEH,", ", "Act:",ACT, ", ", "Rec:",REC, ", ", ENC_TIME_start, "-", ENC_TIME_stop, sep="")
      
      ## Plot trajectories of both actor & receiver, show on the same panel
      plot  (y ~ x, traj_ACT, pch=".", col=rgb(0,0,1,0.3,1), main=Title, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
      points(y ~ x, traj_REC, pch=".", col=rgb(1,0,0,0.3,1))
      
      ## subset the trajectories of both actor & receiver using the start & end times
      traj_ACT <- traj_ACT [ which(traj_ACT$UNIX_time >= ENC_TIME_start & traj_ACT$UNIX_time <= ENC_TIME_stop),]
      traj_REC <- traj_REC [ which(traj_REC$UNIX_time >= ENC_TIME_start & traj_REC$UNIX_time <= ENC_TIME_stop),]
      
      print(paste("For row",ROW,"traj_ACT has",nrow(traj_ACT),"traj_REC has", nrow(traj_REC)))
      
      #check start time correspondance
      # print(paste("Behaviour:",BEH,"ACT:",Act_Name,"REC:",Rec_Name, "annot_start", ENC_TIME_start, "traj_start", min(traj_ACT$UNIX_time,na.rm = TRUE)))
      # 
      # ## Plot trajectories of both actor & receiver, show on the same panel
      plot   (y ~ x, traj_ACT, type="l", lwd=6, col="blue4", main=Title, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
      points (y ~ x, traj_REC, type="l", lwd=6,  col="red4")
      
      #plot angles of both actor & receiver
      #plot(angle ~ time, traj_ACT, col=rgb(0,0,1,0.3,1), main=paste("angle_rad",BEH,", Act:",ACT, "Rec:",REC, ENC_TIME_start, "-", ENC_TIME_stop)) #,xlim = c(0,500),ylim = c(0,2*pi))
      #points(angle ~ time, traj_REC, col=rgb(1,0,0,0.3,1))
      
      
      ########################################################
      #CHECK: HOW IT DEALS WITH MISSING COORDINATES? should the trajectories be resampled/rediscretized?
      
      ##define trajectory
      trajectory_ACT                 <- TrajFromCoords(data.frame(x=traj_ACT$x,y=traj_ACT$y,time=as.numeric(traj_ACT$UNIX_time)), spatialUnits = "px",timeUnits="s")
      trajectory_REC                <- TrajFromCoords(data.frame(x=traj_REC$x,y=traj_REC$y,time=as.numeric(traj_REC$UNIX_time)), spatialUnits = "px",timeUnits="s")
      
       #trajectory                      <- TrajResampleTime (trajectory_ori,stepTime =desired_step_length_time )
  
      ###1. total distance moved 
      distance_ACT                   <- TrajLength(trajectory_ACT)
      distance_REC                   <- TrajLength(trajectory_REC)
  
  #   mean angle difference (data angle -> orientation of the ant in radians)
      
     ###2. Mean and median speed
     #  deriv_traj                      <- TrajDerivatives(trajectory_ACT)
     #  mean_speed_mmpersec_ACT             <- mean(deriv_traj$speed,na.rm=T)
     #  median_speed_mmpersec           <- median(deriv_traj$speed,na.rm=T)
     #  ###3. Mean and median acceleration
     #  mean_acceleration_mmpersec2     <- mean(deriv_traj$acceleration,na.rm=T)
     #  median_acceleration_mmpersec2   <- median(deriv_traj$acceleration,na.rm=T)
     #  ###4. Mean and median turn angle, in radians
     #  mean_turnangle_radians          <- mean(abs( TrajAngles(trajectory)),na.rm=T)
     #  median_turnangle_radians        <- median(abs( TrajAngles(trajectory)),na.rm=T)
     #  ###5. Straightness
     #  straightness_index              <- Mod(TrajMeanVectorOfTurningAngles(trajectory))
     #  sinuosity                       <- TrajSinuosity(trajectory)
     #  sinuosity_corrected             <- TrajSinuosity2(trajectory)
     #  ###6. Expected Displacement
     #  expected_displacement_mm        <- TrajEmax(trajectory,eMaxB =T)
     #  ###7. Autocorrelations
     #  all_Acs                         <- TrajDirectionAutocorrelations(trajectory)
     #  periodiciRUE) :ty_sec                 <- TrajDAFindFirstMinimum(all_Acs)["deltaS"]*desired_step_length_time
     #  ###8. root mean square displacement
     # # rmse_mm                            <-  sqrt(sum( (traj[,c(paste("x",ACT,sep=""))]-mean(traj[,c(paste("x",ACT,sep=""))]))^2 + (traj[,c(paste("y",ACT,sep=""))]-mean(traj[,c(paste("y",ACT,sep=""))]))^2 )/length(na.omit(traj[,c(paste("x",ACT,sep=""))])))
  
     # summary_data <- rbind(summary_data,data.frame(BEH=BEH,Act_Name=Act_Name,Rec_Name=Rec_Name,distance_mm=distance_mm,mean_speed_mmpersec=mean_speed_mmpersec,median_speed_mmpersec=median_speed_mmpersec,mean_acceleration_mmpersec2=mean_acceleration_mmpersec2,median_acceleration_mmpersec2=median_acceleration_mmpersec2,mean_turnangle_radians=mean_turnangle_radians,median_turnangle_radians=median_turnangle_radians,straightness_index=straightness_index,sinuosity=sinuosity,sinuosity_corrected=sinuosity_corrected,expected_displacement_mm=expected_displacement_mm,periodicity_sec=periodicity_sec,stringsAsFactors = F))
      summary_data <- rbind(summary_data,data.frame(BEH=BEH,Act_Name=Act_Name,Rec_Name=Rec_Name,distance_ACT=distance_ACT, stringsAsFactors = F))
      
      }##ACT
    }##BEH

  ########################################################
    # for (variable in names( summary_data)[!names(summary_data)%in%c("BEH","Act_Name","Rec_Name")]){
    #   summary_data[,"variable"] <- summary_data[,variable]
    #   boxplot(variable~BEH,ylab=variable,data=summary_data)
    # }##var_PLOT
    
    ########################################################
  }##PERIOD
}##REPLICATE

rm(list=(c("e")))
gc() # clear cache








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
# interactions_all       <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportTrajectories = T)
# 
# ###let's View the object: 
# str(interactions_all)
# ###it's a list containing 3 objects: summaryTrajectory, trajectories and interactions
# ###let's extract them for simplicity:
# summaryTrajectory_all       <- interactions_all$summaryTrajectory  ## same object as positions$trajectory_summary described above - however here it will have many more rows given that the trajectories will be broken whenever there is a 10 second non-detection gap for an ant.
# ## therefore in this case my trick with ant_ID_str won't work - it only works if there's exactly one trajectory as ants in the data
# ## so I would NOT use the output of this function to analyse trajectories
# trajectories_all           <- interactions_all$trajectories       ## same object as positions$trajectories described above - but containing more objects for the same reason as stated above
# interacions_all            <- interactions_all$interactions       ## data frame containing all interactions
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
# ###Which we find in object positions$trajectory_summary : positions$trajectory_summary[which(positions$trajectory_summary$antID_str=="ant1"),"start"]
# trajectory$time_abs <- positions$trajectory_summary[which(positions$trajectory_summary$antID_str=="ant_1"),"start"] + trajectory$time
# 
# 
# trajectories_unlist$time_abs <- positions$trajectory_summary[which(positions$trajectory_summary$antID_str),"start"] + trajectory$time
# positions$trajectories <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=trajectory_abs_time)) ###the match is once again to ensure we will extract the right information for the right ant
# 
# ##################################################################################################################################
# ###### 7. EXAMPLE USE OF TRAJECTORIES: how to apply a function to all trajectories simultaneously#################################
# ##################################################################################################################################
# ###let's assume we want to know the duration of each trajectory and fill in the information into the positions$trajectory_summary file 
# funTest <- function(trajectory) {
#   if(positions$trajectory_summary$antID_str==trajectories_unlist$.id) {
#     return (positions$trajectory_summary$start + trajectories_unlist$time)
#   } 
# }
# 
# trajectories_unlist$time_abs <- lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=funTest)
# 
# ###First we would define our own function:
# trajectory_absol_time                 <- function(trajectory){ return (positions$trajectory_summary$start + trajectory$time)}
# positions$trajectory_summary$start
# trajectory$time
# trajectory_duration                   <- function(trajectory){ return (max(trajectory$time,na.rm=T)-min(trajectory$time,na.rm=T))}
# 
# 
# ###Second let's apply that function to all trajectories and fill in the results
# positions$trajectories <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=trajectory_abs_time)) ###the match is once again to ensure we will extract the right information for the right ant
# 
# ###Another example: number of coordinates per trajectory
# positions$trajectory_summary$nb_coordinates <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=nrow)) ###the match is once again to ensure we will extract the right information for the right ant
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
