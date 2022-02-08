rm(list=ls())

#"https://formicidae-tracker.github.io/myrmidon/docs/latest/api/index.html"

######load necessary libraries
library(adehabitatHR) ####for home range and trajectory analyses
library(FortMyrmidon) ####R bindings
library(igraph)       ####for network analysis


##################################################
###### 1. OPENING AN EXPERIMENT #####################
##################################################
####### navigate to folder containing myrmidon file
setwd("/media/bzniks/SAMSUNG/Postprocessing_training")

######There are two ways to open a myrmidon file
######Either as read only (then the myrmidon file is locked and cannot be edited)
######Of as read/write (then the myrmidon file can be edited)
###Read-only version:
e <- fmExperimentOpenReadOnly("training.myrmidon")
#########################
### Note 1: a myrmidon file can only be open in a single program instance, by design
### So if you have fort-studio open in the background, this step will fail - and you will need to go and stop fort-studio first
### Similarly, if R is open and you have the experiment open in there, you won't be able to load the file in fort-studio and you will need to flush the memory first (see Note 2)
#########################
### Note 2: if you now want to switch to the read/write because you would like to fill in some information into the metadata file, you will need to delete the object containing the experiment from R and clear the cache using gc()
rm(list=(c("e")))
gc()
e <- fmExperimentOpen("training.myrmidon")
rm(list=(c("e")))
gc()

#####################################################
###### 2. CREATING ANTS #############################
#####################################################
##When you create ants directly in fort-studio, antIDs will be created in the order in which you orient the ants
##This means that if two users were to orient the same experiment, but did it in a slightly different orders, they would end up with different antID lists where different antID would correspond to different tagID
##It also means that if for any reasons you wanted to replicate your results / reorient your ants from scratch, and wanted to compare results to your previous analysis, you may end up with a different antID list and this would make the comparison difficult
##One way to avoid that is to create the ants from R using a fixed code
###To do so, you first need to create a myrmidon file and fill in the first page (defning the space, author, tag size and comment)
###Then, BEFORE DOING ANYTHING ELSE, run the following code in RStudio to create the ants in a set order
###At this step you may already want to apply a filter so that tags that had very few detections compard to others wll be discarded as false detections / low-quality ants

###For this step, you will need to open the experiment as read/write as you will want to overwrite your "blank" myrmidon file with a file that contains antIDs
setwd("/media/bzniks/SAMSUNG/Postprocessing_training")
tracking_data <- fmExperimentOpen("training_blank.myrmidon")

###next you will want to extract the tag statistics to know how many times each tag was detected using the following function:
ts <- fmQueryComputeTagStatistics(tracking_data)
###to understand the content of the object ts, you can use the following functions
class(ts) ### tells you what sort of object ts is - in this case, a data.frame 
str(ts)   ### gives an overview of the different components of ts. Quite useful if ts was a list rather than a data.frame
head(ts)  ### displays the first few lines of ts
####We see here that:
####the row names of this object correspond to the tagID, in hexadecimal format 
####the first column contains the same tagID, but in decimal values. Ignore this! They are different form the antID you will create.
###then for each tag there is a variety of information, including "count" which corresponds to the number of detections

###creating the ants
for ( i in 1:nrow(ts)) {  #####loop over each tag
  if ( ts[i,"count"] >= 0.001*max(ts[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
    a <- tracking_data$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
    tracking_data$addIdentification(tracking_data$cAnts()$summary[nrow(tracking_data$cAnts()$summary),"antID"],ts[i,"tagDecimalValue"],fmTimeInf(),fmTimeInf()) ####this lines create an ant, valid from -infinte to + infinite
  }
}

###overwriting the myrmidon file
###For some reason, if we just give the relative path name of the myrmidon file (as in the next line), R will throw an error. try it:
tracking_data$save("training_blank.myrmidon")  #### "Error in tracking_data$save("training_blank.myrmidon") : Could not acquire exclusive lock on 'training_blank.myrmidon':  another program has write or read access on it"
###But this does not happen if you specify the full path instead:
tracking_data$save("/media/bzniks/SAMSUNG/Postprocessing_training/training_blank.myrmidon")  
###Alternatively, you could of course save a new myrmidon file with a new name instead of overwriting!

###Advantages of creating ants in R
######1. Repeatability (as explained above)
######2. Useful when you want to identify and display ants with specific properties halfway through the experiment, and you don't have the time to orient them manually (e.g. Adriano's experiment requiring to identify nurses)
###Disadvantages: 
####Does not allow to have several tags pointing to the same ant (retag). 
####BUT this can be altered manually in fort-studio afterwards by deleting the second ant and then adding the second tagID as identification to the first ant

rm(list=c("tracking_data"))
gc()


###############################################################################
###### 3. QUERYING GENERAL EXPERIMENT/ANT INFORMATION #############################
###############################################################################
e <- fmExperimentOpenReadOnly("training.myrmidon")
e$getDataInformations()

###the following query allows you to access the measurements made in fort-studio for a particular ant:
fmQueryComputeMeasurementFor(e,antID=1,mTypeID=1)
###where antID is the ID of the ant ( not the tag) and mTypeId the ID of the measurement type
###By default fort-studio always compute the length, which is the distance between the two points defined when orienting the ant, normally from tip of gaster to mid-mandible points; and that would be type 1
###So that would correspond to the body length of the ant.
###But you could specify more types of measurements within fort-studio (e.g. headwidth, eye-to-eye distance) - then you'd need to know what each mTypeID corresponds and query the correct one
###Note: if you oriented the ant multiple times, you will get one row for each pose where an orientation was defined, together with the time corresponded to that pose

###tag statistics
tag_stats <- fmQueryComputeTagStatistics(e)

###Detail of all ant coordinates on each frame
all_tracking <- fmQueryIdentifyFrames(e)
###Visualise all_tracking object:
View(all_tracking) ###list containing a dataframe (frames) and a list of positions (positions), which can be extracted as follows
Frames_summary <- all_tracking$frames
Frames_detail <- all_tracking$positions
####single frame:
fmQueryIdentifyFrames(e,start=fmTimeCPtrFromAnySEXP(Frames_summary[2,"time"]),end=fmTimeCPtrFromAnySEXP(Frames_summary[3,"time"]))
fmQueryIdentifyFrames(e,start=Frames_summary[2,"time"],end=Frames_summary[3,"time"])


###equivalence between tagID and antID at a certain time
###now
e$identificationsAt(fmTimeNow(),FALSE)
###at start of experiment
e$identificationsAt(fmTimeCPtrFromAnySEXP(e$getDataInformations()$start),FALSE)
###at specific time
e$identificationsAt(fmTimeParse("2020-02-11T11:53:26.790610Z"),FALSE) 
###at a specific frame
e$identificationsAt(fmTimeCPtrFromAnySEXP(Frames_summary[2,"time"]),FALSE)


###############################################################################
###### 4. READING TRAJECTORIES ################################################
###############################################################################
###To obtain trajectories, and trajectories only, the function to use is fmQueryComputeAntTrajectories
###To get more detail about the function, type the following:
?fmQueryComputeAntTrajectories
###we see we have to specify a few arguments.
##Those of note are the following:
######### start/end: 
######################## if not specified / set as NULL, the trajectories will be computed over the entire duration of the experiment
######################## alternatively you may want to specify specific times ranges
######################## This is a little tricky - the handling of times is very obscure and not intuitive at all. it is one of the things that will be improved in the future
######################## Here are some examples that work:
time_start <- fmTimeParse("2021-03-22T16:35:44.149918831Z") ####time copied from fort-studio Visualization
time_stop  <- fmTimeParse("2021-03-22T18:35:16.203359796Z")  ####time copied from fort-studio Visualization

time_start <- fmTimeCPtrFromAnySEXP(e$getDataInformations()$end - 48*3600)####last time in tracking minus 48 hours 
time_stop  <- fmTimeCPtrFromAnySEXP(e$getDataInformations()$end)          ####last time in experiment

arbitrary_GMT_time   <- as.POSIXct("2021-03-25 16:00:00 GMT")     ####arbitrary time 
time_start           <- fmTimeCPtrFromAnySEXP(arbitrary_GMT_time)
time_stop            <- fmTimeCPtrFromAnySEXP(arbitrary_GMT_time + 24*3600) ####arbitrary time in the correct format + 24 hours

time_start <- fmTimeCPtrFromAnySEXP(e$getDataInformations()$start )                   ####first time in tracking 
time_stop  <- fmTimeCPtrFromAnySEXP(e$getDataInformations()$start + 24*3600)          ####first time in tracking plus 24 hours

######### maximumGap: 
######################## A very important parameter to set!
######################## Basically the program will "cut" each ant trajectories into segment of maximum duration maximumGap
######################## That would be an absolute nightmare to analyse
######################## So the trick is to set a maximumGap high enough that no cutting will happen
######################## For example, set it at an entire year: fmHour(24*365)
max_gap <- fmHour(24*365)

######### matcher: 
######################## Allows to ask for trajectories to be output only for ants with a specific criterion
######################## See section 8 below for examples

######### computeZones
######################## Set to true to output which zone each coordinate is. Very useful to define space use!

######### So overall:
positions <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = TRUE)

###Now let's have a look at the structure of this object
View(positions)
###positions is a list of two objects: trajectories_summaries (a dataframe with 22 rows, one row per ant) and trajectories (a list of length 22, one element per ant)

###Let's have a look at the trajectories_summaries data.frame within positions:
View(positions$trajectory_summary)
###it's a very simple data.frame with 3 columns: antID, start and space
###antID is the ID of the ant (NOT the tag ID! Forget all about tagIDs from now on)
###space is the index of the space the ant is in (in our example we have a single space, esterhase)
###start is the time at which that particular ant's trajectory starts (remember, not all ants are detected in the first frame so this may vary)

###Let's have a look at the first element of the trajectories list within positions:
View(positions$trajectories[[1]])
###that is a dataframe with columns "time" (in seconds since start), x, y (coordinates in pixels), angle (orientation of the ant in radians), zone (which zone the ant is)

###WARNING!!!!!
###in this particular example the values in antID are consecutive (1 to 22), and so positions$trajectories[[1]] will correspond to the trajectory of ant 1
### but imagine you had to delete an ant in fort-studio and the antID list "jumped" from 10 to 12 (so it would be 1,...10 then 12...23)
### then the 22nd object of trajectory list would not be the trajectory of ant 22 but of ant 23. That is not fool proof, but very risky!
##so what I suggest you do immediately after computing your positions object:
positions$trajectory_summary$antID_str <- paste("ant_",positions$trajectory_summary$antID,sep="") ##creates a ID string for each ant: ant1, ant2,...
names(positions$trajectories)       <- positions$trajectory_summary$antID_str ###and use the content of that column to rename the objects within trajectory list
##and now we can extract each ant trajectory using the character ID storesd in antID_str - it's much safer!


##################################################################################################################################
###### 5. EXAMPLE USE OF TRAJECTORIES: distinguishing between nurses and foragers ################################################
##################################################################################################################################
###Here we will use the zones to find out which ants always stay inside the nest (nurses), and which ants occasionally leave the nest (foragers)

###First, let's have a look at the zones we created. It's again a bit of a tricky query - just copy-paste it in the future
zones <- e$cSpaces()$objects[[1]]$cZones()$summary #function to show the Zones present in the Space
print(zones)

###Here we have four zones: Nest, Arena, Water, Honeywater; and the dataframe zone gives their correspondence in terms of zone ID (1,2,3,4)
###This may vary from experiment to experiment so ypu will always need to refer to that dataframe to extract the correct correspondence
###Here what we will do is look through all the trajectories and measured the percentage of the time an ant spent in a zone that is not the nest zone
##So first let's extract the nest_zone ID
nest_zone <- zones[which(zones$name=="Nest"),"ID"]

###Then let's prepare a column within the trajectories_summaries to hold the proportion of time spent in nest
positions$trajectory_summary$prop_time_in_nest <- NA

###Then let's loop over the trajectories_summaries dataframe:
for (i in 1:nrow(positions$trajectory_summary)){
  ####to be fool proof, and be sure you extract the trajectory corresponding the correct ant, make sure you make use of the antID_str column!
  traj <- positions$trajectories [[   positions$trajectory_summary[i,"antID_str"]    ]]
  
  ###then let's fill in the prop_time_in_nest for that ant in the trajectory_summary object:
  positions$trajectory_summary[i, "prop_time_in_nest"] <- length (which(traj$zone==nest_zone)) / nrow(traj)
}

###Now we can use this information to define who is a nurse and who is a forager.
### You may want to define your own criterion in the future, but for simplicity here we will say that nurses are those which stayed 100% of time in the nest
### Let's define a new column in positions$trajectory_summary that will contain that information (e.g. column group)
##By default we will give it the value "nurse"
positions$trajectory_summary$group <- "nurse"

###and now let's use the content of column prop_time_in_nest to identify foragers:
positions$trajectory_summary[which(positions$trajectory_summary$prop_time_in_nest<0.95),"group"] <- "forager"

##################################################################################################################################
###### 5b. Exporting list of foragers to be read by artemis and siaplayed in real time ##########################################
##################################################################################################################################

###Now let's assume we want to output the list of TagID (hexadecimal) corresponding to foragers, in the correct format so we can highlight them in real time in the tracking
#First, let's obtain the table that matches tagID to antID
IDs <- e$identificationsAt(fmTimeNow(),FALSE)

#Second, as the hexadecimal tagID are contained in the rownames rather than in a devoted column, let's create a column that will contain the tag hexadecimal IDs for each ant
IDs$tag_hex_ID=rownames(IDs)

# Third copy that information into positions$trajectory_summary - using function match in case one the two data frames have different row orders
positions$trajectory_summary$tag_hex_ID <- IDs[ match(positions$trajectory_summary$antID,IDs$antID)     , "tag_hex_ID"]

# Fourth export forager list in a format that can be read directly by artemis: hexadecimal tagIDs corresponding to foragers in data frame. 
forager_list <- data.frame(tag_hex_id=paste("- ", as.character(positions$trajectory_summary[which(positions$trajectory_summary$group=="forager"),"tag_hex_ID"] ),sep=""))
write.table(forager_list, file = paste(e$getDataInformations()[["details"]][["tdd.URI"]],"_forager_list.txt",sep=""), append = FALSE,
            row.names = FALSE, col.names = FALSE,quote = FALSE) 


############################################################################################################################################################################
###### 5c. Saving ant group into myrmidon metadata so nurses and foragers can easily be displayed in different colors in Visualization  ####################################
############################################################################################################################################################################

###Now let's assume we want to copy this information into the metadata - DOES NOT WORK YET
##For this we will need to open the myrmidon file in read/write
###First remove object "e" and flush cache
rm(list=c("e"))
gc()
### Second after closing fort-studio, reopen myrmidon file in read/write
e <- fmExperimentOpen("training.myrmidon")
###Third add a new metadata column called group
e$addMetadataColumn("group",fmAntMetadataType$STRING,fmAntStaticString("ant")) 



###here we could potentially have the same discrepancy issue between the order of the antID and the values of the antID - especially if you have accidentally mixed up the columns of the positions$trajectory_summary object.
### so we are going to use a convoluted method to make sure we get this right
###first, we are going to extract the list of ants from the experiment using the function we will use to edit the metadata:
ant_list <- e$ants()$summary
##and as above we are going to create a antID_str column
ant_list$antID_str <-  paste("ant_",ant_list$antID,sep="")

###then let's loop over our positions$trajectory_summary table
for (i in 1:nrow(positions$trajectory_summary)){
  ###get the index correspondence between positions$trajectory_summary and ant_list
  j <- which(ant_list$antID_str == positions$trajectory_summary[i,"antID_str"])

  ###check: we should be able to read the value "group" for that ant - before modification
  e$ants()$objects[[j]]$getValue("group",fmTimeNow()) ### value = ant

  ###now let's modify this value 
  e$ants()$objects[[j]]$setValue(name="group",value=fmAntStaticString(positions$trajectory_summary[i,"group"]),time=fmTimeInf() ) 

  ###check: the value should now have changed
  e$ants()$objects[[j]]$getValue("group",fmTimeNow()) ### value = ant

}

###finally save changes, and reopen in read-only
e$save("/media/bzniks/SAMSUNG/Postprocessing_training/training.myrmidon")
rm(list=c("e"))
gc()

e <- fmExperimentOpenReadOnly("training.myrmidon")

##################################################################################################################################
###### 6. EXAMPLE USE OF TRAJECTORIES: using R libraries to calculate turn angles and home ranges#################################
##################################################################################################################################
####This is just meant as a teaser for what we can do with trajectory analysis in R

####First let's extract a single ant's trajectory
trajectory <-  positions$trajectories[["ant1"]]

###Second - let's convert the relative time (in seconds since start) into absolute times
###For this we need the starting time of that ant's trajectory
###Which we find in object positions$trajectory_summary : positions$trajectory_summary[which(positions$trajectory_summary$antID_str=="ant1"),"start"]
trajectory$time_abs <- positions$trajectory_summary[which(positions$trajectory_summary$antID_str=="ant1"),"start"] + trajectory$time

###Third: create a trajectory object in R
R_traj <- as.ltraj(trajectory[c("x","y")],date=as.POSIXct(trajectory$time_abs ),id="ant1")[[1]]

###This object is a data.frame containing basic information, such as the distance moved at each step, rel.angle the relative turn angle at each step (with negative and positive values depending on whether the ant turns left and right ), abs.angle the direction of each move relative to the x-axis
##Up to you to perform more analyses on the trajectory of each ant depending on your question


###Other example: how to calculate a home range (area in which the ant spends 95% of her time)
trajectory_sp <- SpatialPoints(trajectory[c("x","y")]) ###cretae a SpatialPoints object
trajectory_kde <- kernelUD(trajectory_sp) ###computes a kernel utilisation distribution based on the points provided
trajectory_hr <- getverticeshr(trajectory_kde, 95) ###returns the home range (95%)
trajectory_core <- getverticeshr(trajectory_kde, 50) ###returns the core area (50%)

##################################################################################################################################
###### 7. EXAMPLE USE OF TRAJECTORIES: how to apply a function to all trajectories simultaneously#################################
##################################################################################################################################
###let's assume we want to know the duration of each trajectory and fill in the information into the positions$trajectory_summary file in one go (i.e. without a loop, which is not very effective)

###First we would define our own function:
trajectory_duration                    <- function(trajectory){ return (max(trajectory$time,na.rm=T)-min(trajectory$time,na.rm=T))}
###Second let's apply that function to all trajectories and fill in the results
positions$trajectory_summary$duration <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=trajectory_duration)) ###the match is once again to ensure we will extract the right information for the right ant

###Another example: number of coordinates per trajectory
positions$trajectory_summary$nb_coordinates <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=nrow)) ###the match is once again to ensure we will extract the right information for the right ant

###etc. Generally this will save you time compared to a loop


##################################################################################################################################
###### 8. READING TRAJECTORIES WITH MATCHER#######################################################################################
##################################################################################################################################
###example 1: computes trajectory for a single ant, e.g. ant 5
positions_matcher_1 <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = TRUE,matcher=fmMatcherAntID(5))
View(positions_matcher_1)                     ###shows that only 1 trajectory was output
print(positions_matcher_1$trajectory_summary) ###shows that the trajectory output was that of ant 5

###example 2: computes trajectories of all nurses
positions_matcher_2 <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = TRUE,matcher=fmMatcherAntColumn("group",fmAntStaticString("nurse")))
View(positions_matcher_2) ###returns 14 trajectories (nurses only)

positions_matcher_3 <- fmQueryComputeAntTrajectories(e,start = time_start,end = time_stop,maximumGap = max_gap,computeZones = TRUE,matcher=fmMatcherAntColumn("group",fmAntStaticString("forager")))
View(positions_matcher_3) ###returns 8 trajectories (foragers only)


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
###for example, if we want interactions that involve either body/body or antenna/body interasections:
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
###### 11. BUILDING NETWORK FROM INTERACTIONS #################################
###############################################################################
###to build an aggregated network, we are going to draw an edge between each pair of ants that was seen to interact
###and assign an edge weight to each edge. this can be either the total number of interactions, or the cumulated duration of interactions
###for this it is useful to create a column "number" and "duration" within the interactions data.frame
###where number is always equal to 1 (each interaction counts for 1), and duration = end-start+frame_duration (so that 1-frame interactions don't have a duration of 0)
interactions$number   <- 1

###where can we get the duration of a frame? Either you know the frame rate of your experiment and it will be 1/frame rate, or we can use the results of fmQueryIdentifyFrames
###we ran fmQueryIdentifyFrames earlier and stored the results in Frames_summary 
###We can make sure Frame_summary is sorted in order of increasing time, and then calculate the time difference between each successive frames and takr the mean
Frames_summary <-Frames_summary[order(Frames_summary$time),]
time_differences <- diff(Frames_summary$time)
frame_duration <- mean(time_differences, na.rm=T)

interactions$duration <- interactions$end-interactions$start + frame_duration

###We can now build the network
##1. get the list of nodes involved:
actors <- data.frame(name=as.character(unique(c(interactions$ant1,interactions$ant2))))
###2. build network in first instance 
net <- graph.data.frame(interactions[c("ant1","ant2")],directed=F,vertices=actors)
###3. add edge weights
E(net)$weight <- interactions[,"duration"] ###or : E(net)$weight <- interactions[,"number"]
##this will have created mutliple edges, i.e. if there are several inetractions between the same pair of ats, it will have created as many edges between these two ants as interactions
##so now we want to merge these multiple edges into one per pair, and sum their respective weights
##this is done using function simplify
net <- simplify(net,remove.multiple=TRUE,remove.loop=TRUE,edge.attr.comb="sum")

##Now our network is ready to be plotted or analysed
###follow igraph tutorials for detail

###############################################################################
###### 12. QUERYING METADATA ###################################################
###############################################################################
get_ant_metadata <- function(ant,column,type,time=NULL){
  ant_id <- ant$antID() 
  if (type=="bool"){
    value  <- ant$getValue(column,fmTimeNow())$toBool()
  }else if (type=="string"){
    value  <- ant$getValue(column,fmTimeNow())$toString()
  }else if (type=="double"){
    value  <- ant$getValue(column,fmTimeNow())$toDouble()
  }else if (type=="integer"){
    value  <- ant$getValue(column,fmTimeNow())$toInteger()
  }else{
    print(paste("No method defined for type",type))
    value <- NA
  }
  DF <- data.frame(ant_id=ant_id,value=value);names(DF)=c("Ant",column)
  return(DF)
}
group <- NULL; 
for ( a in e$cAnts()$objects) {
  group <- rbind(group,get_ant_metadata(a,column="group",type="string"))
}
  

###############################################################################
###### 13. READING COLLISIONS ##################################################
###############################################################################
###collisions are for each frame, the list of ants whose shapes intersect one another. Normally not used
collisions <- fmQueryCollideFrames(e, start=time_start, end=time_stop)

###visualise collision object
View(collisions)
###collisions contains 3 objects: frames, positions and collisions
###let's extract them#
collisions_frames    <- collisions$frames      ###data frames giving detail (time, space, width,height) about each frame
collisions_positions <- collisions$positions   ###list containing as many objects as there are frames in collisions_frames. For each one, lists the positions of all detected ants on that frame
                                               ###(those two are the exactly the same as the output from fmQueryIdentifyFrames)
collisions <- collisions$collisions            #### dataframe containing all detected collisions

###let's view the first few lines of collisions
head(collisions) ###columns giving the ID of the two colliding ants, the zone where they are, the types of collisions (all types), and the frames_row_index referring to which frame that collision was detected in (matches the list indices in collisions_positions)
head(collisions)

