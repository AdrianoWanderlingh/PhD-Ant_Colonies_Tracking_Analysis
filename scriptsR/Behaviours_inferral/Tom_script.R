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
setwd("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/")

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
e <- fmExperimentOpenReadOnly("R3SP_13-03-21_Capsule_defined.myrmidon")
e$getDataInformations()

objectss <- fmQueryComputeMeasurementFor(e,1,1) #blindly acess an element, in this cas: measures taken for the body lenghts

###tag statistics
tag_stats <- fmQueryComputeTagStatistics(e)



###############################################################################
###### 4. READING TRAJECTORIES ################################################
###############################################################################
###To obtain trajectories, and trajectories only, the function to use is fmQueryComputeAntTrajectories
###To get more detail about the function, type the following:
?fmQueryComputeAntTrajectories
###we see we have to specify a few arguments.
##Those of note are the following:
######### start/end: 
time_start_ISO <- parse_iso_8601("2021-03-16T12:13:21.670072093Z")
time_stop  <- fmTimeCPtrFromAnySEXP(time_start_ISO + ((3600*24)/90)) ####arbitrary time in the correct format + 16min
time_start <- fmTimeCPtrFromAnySEXP(time_start_ISO)

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

positions <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = FALSE) #set true to obtain the zone of the ant

###WARNING!!!!!
###in this particular example the values in antID are consecutive (1 to 22), and so positions$trajectories[[1]] will correspond to the trajectory of ant 1
### but imagine you had to delete an ant in fort-studio and the antID list "jumped" from 10 to 12 (so it would be 1,...10 then 12...23)
### then the 22nd object of trajectory list would not be the trajectory of ant 22 but of ant 23. That is not fool proof, but very risky!
##so what I suggest you do immediately after computing your positions object:
positions$trajectory_summary$antID_str <- paste("ant_",positions$trajectory_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
names(positions$trajectories)       <- positions$trajectory_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
##and now we can extract each ant trajectory using the character ID stored in antID_str - it's much safer!


######################## TEMPORARY - ADRIANO TO VALIDATE AGREEMENT BETWEEN 'CORRECTED' UNIX TIMES & TIMES ON THE VIDEOS
## Add times to the positions to convert from time since the start of the experiment, to UNIX time


for (A in positions$trajectory_summary$antID)
  {
  AntID <- paste("ant_",A,sep="") 
  First_Obs_Time <- positions$trajectory_summary$start [which(positions$trajectory_summary$antID_str==AntID)] ## find the first time after the user defined time_start_ISO that this ant was seen
  print(paste("Adding first obs time", First_Obs_Time, "to the time-zeroed trajectory of ant", AntID))
  positions$trajectories[[AntID]] $UNIX_time <- positions$trajectories[[AntID]] $time + First_Obs_Time ##convert back to UNIX time  
  }

########################  LOAD MANUAL-ANNOTATIONS

annotations <- read.csv("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/reproducible_example_Adriano/R3SP_Post1_updated.csv", sep = ",")
annotations$Actor <- as.character(annotations$Actor)

#convert Zulu time to UTC
annotations$T_start <- parse_iso_8601(annotations$T_start)
annotations$T_stop <- parse_iso_8601(annotations$T_stop)

table(annotations$Behaviour)## with(annotations, tapply(frequency, list(Behaviour), sum))
Counts_by_Actor <- aggregate(treatment_rep ~ Behaviour + Actor, FUN=length, annotations); colnames(Counts_by_Actor) [match("treatment_rep",colnames(Counts_by_Actor))] <- "Count"


######FOR 1 ANT
####First let's extract a single ant's trajectory
par(mfrow=c(2,3), mai=c(0.3,0.3,0.1,0.1), mgp=c(1.3,0.3,0), family="serif", tcl=-0.2)

for (BEH in c("G","T"))
  {
  annot_BEH <- annotations[which(annotations$Behaviour==BEH),]
  ## remove doubled allo-grooming interactions
  if (BEH=="G") {annot_BEH <- annot_BEH[!duplicated(annot_BEH),]}  ## leave NOT to catch possible un-matched rows
  if (BEH=="T") {print("STILL TO DO")}
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
    
    ## Plot trajectories of both actor & receiver, show on the same panel
    plot  (y ~ x, traj_ACT, pch=".", col=rgb(0,0,1,0.2,1), main=paste("Beh",BEH,", Act:",ACT, "Rec:",REC, ENC_TIME_start, "-", ENC_TIME_stop))#, xlim=c(2500,8000),ylim=c(500,5500))
    points(y ~ x, traj_REC, pch=".", col=rgb(1,0,0,0.2,1))
    
    ## subset the trajectories of both actor & receiver using the start & end times
    traj_ACT <- traj_ACT [ which(traj_ACT$UNIX_time >= ENC_TIME_start & traj_ACT$UNIX_time <= ENC_TIME_stop),]
    traj_REC <- traj_REC [ which(traj_REC$UNIX_time >= ENC_TIME_start & traj_REC$UNIX_time <= ENC_TIME_stop),]
    
    ## Plot trajectories of both actor & receiver, show on the same panel
    points  (y ~ x, traj_ACT, type="l", lwd=2, col="blue4")
    points(y ~ x, traj_REC, type="l", lwd=2,  col="red4")
    
    }##ACT
  }##BEH



rm(list=(c("e")))
gc()
