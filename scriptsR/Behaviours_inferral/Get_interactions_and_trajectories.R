rm(list=ls())

#"https://formicidae-tracker.github.io/myrmidon/docs/latest/api/index.html"

#install.packages("adehabitatHR")
#install.packages("igraph")

######load necessary libraries
library(adehabitatHR) ####for home range and trajectory analyses
library(FortMyrmidon) ####R bindings
library(igraph)       ####for network analysis
library(parsedate)

##################################################
###### 1. OPENING AN EXPERIMENT #####################
##################################################
####### navigate to folder containing myrmidon file
setwd("/home/cf19810/Documents/Postprocessing_training")

#######################
### Note 2: if you now want to switch to the read/write because you would like to fill in some information into the metadata file, you will need to delete the object containing the experiment from R and clear the cache using gc()
# rm(list=(c("e")))
# gc()
# e <- fmExperimentOpen("training.myrmidon")
# rm(list=(c("e")))
# gc()


###############################################################################
###### 3. QUERYING GENERAL EXPERIMENT/ANT INFORMATION #############################
###############################################################################
e <- fmExperimentOpenReadOnly("training.myrmidon")
e$getDataInformations()

objectss <- fmQueryComputeMeasurementFor(e,1,1) #blindly acess an element, in this cas: measures taken for the body lenghts

###tag statistics
tag_stats <- fmQueryComputeTagStatistics(e)

###Detail of all ant coordinates on each frame
# all_tracking <- fmQueryIdentifyFrames(e) ##computationally heavy

###Visualise all_tracking object:
# View(all_tracking) ###list containing a dataframe (frames) and a list of positions (positions), which can be extracted as follows
# Frames_summary <- all_tracking$frames
# Frames_detail <- all_tracking$positions

####CHECK IF USELESS!
# ###equivalence between tagID and antID at a certain time
# ###now
# e$identificationsAt(fmTimeNow(),FALSE)
# 
# ###at start of experiment
# e$identificationsAt(fmTimeCPtrFromAnySEXP(e$getDataInformations()$start),FALSE)
# ###at specific time
# e$identificationsAt(fmTimeParse("2020-02-11T11:53:26.790610Z"),FALSE) 
# ###at a specific frame
# e$identificationsAt(fmTimeCPtrFromAnySEXP(Frames_summary[2,"time"]),FALSE)


###############################################################################
###### 4. READING TRAJECTORIES ################################################
###############################################################################
###To obtain trajectories, and trajectories only, the function to use is fmQueryComputeAntTrajectories
###To get more detail about the function, type the following:
?fmQueryComputeAntTrajectories
###we see we have to specify a few arguments.
##Those of note are the following:
######### start/end: 
time_start <- parse_iso_8601("2021-03-22T16:35:44.149918831Z")
time_stop  <- fmTimeCPtrFromAnySEXP(time_start + (3600*24)/48) ####arbitrary time in the correct format + 30mins
time_start <- fmTimeCPtrFromAnySEXP(time_start)

######### maximumGap: 
######################## A very important parameter to set!
######################## he trick is to set a maximumGap high enough that no cutting will happen
######################## For example, set it at an entire year: fmHour(24*365)
max_gap <- fmHour(24*365)

######### matcher: 
######################## Allows to ask for trajectories to be output only for ants with a specific criterion
######################## See section 8 below for examples

######### computeZones
######################## Set to true to output which zone each coordinate is. Very useful to define space use!

######### So overall:
positions <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = FALSE) #set true to obtain the zone of the ant
?fmQueryComputeAntTrajectories
###Let's have a look at the first element of the trajectories list within positions:
View(positions$trajectories)
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






View(positions$trajectories)
library(data.table)
trajectories_unlist <- as.data.frame(rbindlist(positions$trajectories,use.names=TRUE,idcol=TRUE))
#output to send to Henry Reeves
#
#
write.csv(trajectories_unlist, file=paste(e$getDataInformations()[["details"]][["tdd.URI"]],'trajectories_30min.csv'),row.names = FALSE) #specify name of the block better!
#
#

trajectory_summary <- positions$trajectory_summary
########################################
##################################################################################################################################
###### 6. EXAMPLE USE OF TRAJECTORIES: using R libraries to calculate turn angles and home ranges#################################
##################################################################################################################################

###Third: create a trajectory object in R
####IT MAY BE USEFUL FOR ANALYES OF BEHAVIOUR?
#R_traj <- as.ltraj(trajectory[c("x","y")],date=as.POSIXct(trajectory$time_abs ),id="ant_1")[[1]]
###This object is a data.frame containing basic information, such as the distance moved at each step, rel.angle the relative turn angle at each step (with negative and positive values depending on whether the ant turns left and right ), abs.angle the direction of each move relative to the x-axis



######FOR 1 ANT
####First let's extract a single ant's trajectory
trajectory <-  positions$trajectories[["ant_1"]]

###Second - let's convert the relative time (in seconds since start) into absolute times
###For this we need the starting time of that ant's trajectory
###Which we find in object positions$trajectory_summary : positions$trajectory_summary[which(positions$trajectory_summary$antID_str=="ant1"),"start"]
trajectory$time_abs <- positions$trajectory_summary[which(positions$trajectory_summary$antID_str=="ant_1"),"start"] + trajectory$time


trajectories_unlist$time_abs <- positions$trajectory_summary[which(positions$trajectory_summary$antID_str),"start"] + trajectory$time
positions$trajectories <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=trajectory_abs_time)) ###the match is once again to ensure we will extract the right information for the right ant

##################################################################################################################################
###### 7. EXAMPLE USE OF TRAJECTORIES: how to apply a function to all trajectories simultaneously#################################
##################################################################################################################################
###let's assume we want to know the duration of each trajectory and fill in the information into the positions$trajectory_summary file 
funTest <- function(trajectory) {
  if(positions$trajectory_summary$antID_str==trajectories_unlist$.id) {
    return (positions$trajectory_summary$start + trajectories_unlist$time)
  } 
}

trajectories_unlist$time_abs <- lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=funTest)

###First we would define our own function:
trajectory_absol_time                 <- function(trajectory){ return (positions$trajectory_summary$start + trajectory$time)}
positions$trajectory_summary$start
trajectory$time
trajectory_duration                   <- function(trajectory){ return (max(trajectory$time,na.rm=T)-min(trajectory$time,na.rm=T))}


###Second let's apply that function to all trajectories and fill in the results
positions$trajectories <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=trajectory_abs_time)) ###the match is once again to ensure we will extract the right information for the right ant

###Another example: number of coordinates per trajectory
positions$trajectory_summary$nb_coordinates <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=nrow)) ###the match is once again to ensure we will extract the right information for the right ant

###etc. Generally this will save you time compared to a loop

###############################################################################
###### 9. READING INTERACTIONS ################################################
###############################################################################
###This provides more synthetic and  complete information than collisions, as if the same pair of ants interacts over successive frames, it will be reported as a single interactions with a start and end time
### The function to extract interactions is fmQueryComputeAntInteractions
###Let's have a look at the arguments needed:
?fmQueryComputeAntInteractions
######### experiment, start, end
######################## As in trajectories

######### maximumGap
######################## EXTREMELY IMPORTANT PARAMETER - THE CHOICE WILL BE DIFFERENT FROM THAT MADE FOR THE TRAJECTORY QUERY ABOVE!
######################## The basic idea is the following:
######################## Within a single biological, real interaction, there are likely to be gaps
######################## i.e. frame in which one or both of te ants are not detected, and so no collision is detected between these two ants for that frame
######################## This is expected because the tracking system is not perfect - or the ants may tilt slightly during the interaction, making their tag undetectable
######################## Yet it would be a mistake to say that whenever one or both ants disappear, the interaction stops and a new starts
######################## So we need to define an "allowance", i.e. a maximum period of continuous interruption of the interaction for which we will still be happy to say that it is the same interaction that continues
######################## If we were to set that at one year, like we did for the trajectories above, then each pair of ant would only ever be as having a single interaction
######################## In the past we have used something akin to 10 seconds. If the interruption is longer than this, we consider that a new interaction starts - and a new line will be written 

######### reportTrajectories
######################## Argument taking either T or F value, determining whether trajectories will be output in parallel to the interactions
######################## Useful if we want to know the x-y coordinates of each ant in each interaction
######################## At the moment the function only works if this is set as TRUE (in my case I get a segmentation fault otherwise)

###so:
interactions_all       <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportTrajectories = T)

###let's View the object: 
View(interactions_all)
###it's a list containing 3 objects: summaryTrajectory, trajectories and interactions
###let's extract them for simplicity:
summaryTrajectory_all       <- interactions_all$summaryTrajectory  ## same object as positions$trajectory_summary described above - however here it will have many more rows given that the trajectories will be broken whenever there is a 10 second non-detection gap for an ant.
## therefore in this case my trick with ant_ID_str won't work - it only works if there's exactly one trajectory as ants in the data
## so I would NOT use the output of this function to analyse trajectories
trajectories_all           <- interactions_all$trajectories       ## same object as positions$trajectories described above - but containing more objects for the same reason as stated above
interactions_all            <- interactions_all$interactions       ## data frame containing all interactions

###let's have a look at the interactions object
head(interactions_all)
####dataframe comtaining the following info:
######ant1: ID of the first ant in the interaction
######ant2: ID of the second ant in the interaction
######start: time at which interaction started
######end: time at which interaction ended
######space: space in which the interaction took place
######types: all types of capsule intersection observed in this interaction. Commas separate different types of intersections, and each intersection type lists which capsule of ant1 interesected with which capsule in at2, separated by a "-"
###############thus 1-5, 5-1,5-5: means that during this interaction, capsule 1 in ant 1 intersected with capsule 5 in ant 2; capsule 5 in ant 1 intersected with capsule 1 in ant 2, and capsule 5 in ant1 intersected with capsule 5 in ant2
######ant1.trajectory.row: gives the index to use within summaryTrajectory and trajectories to extract the appropriate trajectory segment for ant1 during this interaction
######ant1.trajectory.start: within trajectory segment trajectoriessummaryTrajectory[ant1.trajectory.row], ant1.trajectory.start gives the row corresponding to the start of this interaction
######ant1.trajectory.end: within trajectory segment trajectories[ant1.trajectory.row], ant1.trajectory.end gives the row corresponding to the end of this interaction. BUT THIS SEEMS WRONG IN CURRENT VERSION. NEED TO REPORT BUG
######ant2.trajectory.row, ant2.trajectory.start, ant2.trajectory.end: same for ant2

####In the example above we have not specified a matcher, i.e. all intersections between all types of capsules have been considered
####But in many cases we will be interested in a specific type of interactions. 
#### For example: we may be interested only in interactions involving the intersection between body shapes
###First we need to figure out the index of the shape corresponding to body shape
capsules  <- e$antShapeTypeNames()
body_id <- capsules[which(capsules$name=="Body"),"typeID"]
###Then we will specify a matcher in which we are interested in interactions involving capsule 1 for both anta: fmMatcherInteractionType(body_id,body_id)
##for more info, type:
?fmMatcherInteractionType

##Let's re-run interactions with that matcher:
interactions_body       <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportTrajectories = T,matcher = fmMatcherInteractionType(body_id,body_id))
summaryTrajectory_body      <- interactions_body$summaryTrajectory  
trajectories_body            <- interactions_body$trajectories       
interactions_body            <- interactions_body$interactions     

###we could also be interested in antennations, i.e. interactions where the antenna of one ant touches the body of the others
antenna_id <- capsules[which(capsules$name=="Antenna"),"typeID"]
interactions_antennations            <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportTrajectories = T,matcher = fmMatcherInteractionType(body_id,antenna_id))
summaryTrajectory_antennations       <- interactions_antennations$summaryTrajectory  
trajectories_antennations            <- interactions_antennations$trajectories       
interactions_antennations            <- interactions_antennations$interactions     

###note that the number of interactions produced for these different matcher types is different:
nrow(interactions_all)          ###60173 interactions of any type
nrow(interactions_body)         ###29717 interactions body-body
nrow(interactions_antennations) ###49124 antennations

###we can also combine several matchers
###for example, if we want interactions that involve either body/body or antenna/body intersections:
? fmMatcherOr 
interactions_body_or_antennations           <- fmQueryComputeAntInteractions(e,start=time_start, end=time_stop,maximumGap =fmSecond(10),reportTrajectories = T,matcher = fmMatcherOr(list(fmMatcherInteractionType(body_id,body_id),fmMatcherInteractionType(body_id,antenna_id))))
summaryTrajectory_body_or_antennations      <- interactions_body_or_antennations$summaryTrajectory  
trajectories_body_or_antennations           <- interactions_body_or_antennations$trajectories       
interactions_body_or_antennations           <- interactions_body_or_antennations$interactions     
nrow(interactions_body_or_antennations) ###49886 antennations or body/body contacts - indicating that most bioy/body interactions were included within the antennations only


###for the rest, let's focus on the body-body interactions - for simplicity, I will rename them without suffix
summaryTrajectory    <- summaryTrajectory_body
trajectories           <- trajectories_body
interactions           <- interactions_body

###so for example, if we wanted to know the x-y coordinates of ant1 (7) and ant2 (9) on the first frame of the first interaction listed int his table:
coord_ant1_start <- trajectories[[   interactions[1,"ant1.trajectory.row"]    ]][interactions[1,"ant1.trajectory.start"],c("x","y")]
coord_ant2_start <- trajectories[[   interactions[1,"ant2.trajectory.row"]    ]][interactions[1,"ant2.trajectory.start"],c("x","y")]

###if we wanted to know the x-y coordinates of ant1 (7) and ant2 (9) on the last frame of the first interaction listed int his table:
coord_ant1_end <- trajectories[[   interactions[1,"ant1.trajectory.row"]    ]][interactions[1,"ant1.trajectory.end"],c("x","y")]
coord_ant2_end <- trajectories[[   interactions[1,"ant2.trajectory.row"]    ]][interactions[1,"ant2.trajectory.end"],c("x","y")]

###if we wanted to know the mean x-y coordinates of ant1 (7) and ant2 (9) in first interaction listed int his table:
coord_ant1_mean <- colMeans(trajectories[[   interactions[1,"ant1.trajectory.row"]    ]][interactions[1,"ant1.trajectory.start"]:interactions[1,"ant1.trajectory.end"],c("x","y")])
coord_ant2_mean <- colMeans(trajectories[[   interactions[1,"ant2.trajectory.row"]    ]][interactions[1,"ant2.trajectory.start"]:interactions[1,"ant2.trajectory.end"],c("x","y")])

###and if you wanted to compute this for all interactions, you could either loop over all interactions or write a function and use apply
####but that would be a bit more complicated to get to work than the examples I gave above
####here I don't use match but I need to be sure that trajectories remained non-shuffled
mean_coord <- function (x,which_ant,trajectories){
  coord_mean <- colMeans(trajectories[[as.numeric(x[paste(which_ant,".trajectory.row",sep="")] )]]  [as.numeric(x[paste(which_ant,".trajectory.start",sep="")]):as.numeric(x[paste(which_ant,".trajectory.end",sep="")]),c("x","y")] )
  return(data.frame(x=coord_mean["x"],y=coord_mean["y"]))
}
mean_coordinates <- unlist(apply(interactions, 1, FUN=mean_coord,which_ant="ant1",trajectories=trajectories))
interactions$mean_x_ant1 <- mean_coordinates[grepl("x",names(mean_coordinates))]
interactions$mean_y_ant1 <- mean_coordinates[grepl("y",names(mean_coordinates))]


mean_coordinates <- unlist(apply(interactions, 1, FUN=mean_coord,which_ant="ant2",trajectories=trajectories))
interactions$mean_x_ant2 <- mean_coordinates[grepl("x",names(mean_coordinates))]
interactions$mean_y_ant2 <- mean_coordinates[grepl("y",names(mean_coordinates))]

###############################################################################
###### 13. READING COLLISIONS ##################################################
###############################################################################
###collisions are FOR EACH FRAME the list of ants whose shapes intersect one another. Normally not used
collisions <- fmQueryCollideFrames(e, start=time_start, end=time_stop)


###visualise collision object
View(collisions)
###collisions contains 3 objects: frames, positions and collisions
###let's extract them#
collisions_frames    <- collisions$frames      ###data frames giving detail (time, space, width,height) about each frame
collisions_positions <- collisions$positions   ###list containing as many objects as there are frames in collisions_frames. For each one, lists the positions of all detected ants on that frame
###(those two are the exactly the same as the output from fmQueryIdentifyFrames)
collisions_collisions <- collisions$collisions            #### dataframe containing all detected collisions

collisions_positions[[1200]] #all the collisions taking place in a specific frame 

###let's view the first few lines of collisions
head(collisions_collisions) ###columns giving the ID of the two colliding ants, the zone where they are, the types of collisions (all types), and the frames_row_index referring to which frame that collision was detected in (matches the list indices in collisions_positions)
head(collisions_collisions)




# rm(list=(c("e")))
# gc() 