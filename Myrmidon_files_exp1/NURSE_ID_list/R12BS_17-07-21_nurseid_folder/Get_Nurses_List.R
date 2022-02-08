rm(list=ls())
gc()

define_ants <- F ###change back to TRUE when you do a new colony

    library("FortMyrmidon") ###latest version is 0.6.1
######## DOCUMENTATION
# https://formicidae-tracker.github.io/studio/docs/latest/api/

# setwd("/media/adriano/DISK2/36h10_bigNest_0_7mlHumidity_test")
setwd("/media/adriano/DISK3/ADRIANO/NURSEID/R12BS_17-07-21_nurseid_folder/")
myrmidon_file <- "R12BS_17-07-21_nurseid.myrmidon"
##########DEFINE ANTS#############
if (define_ants){source("Define_Ant_Identifications.R")}

#########READ TRACKING (READ ONLY)
tracking_data <- fmExperimentOpenReadOnly(myrmidon_file) #if it says it's locked, click Session>Restart R and clear objects from workspace including hidden ones

####define time period to use for defining nurses/foragers
# time_start <- fmTimeParse("2020-11-04T18:49:21Z")
# time_stop <- fmTimeParse("2020-11-05T18:49:21Z")
###above are the times you chose but I recommend the following for the experiment (so it nevers need to be edited):
time_start <- fmTimeCPtrFromAnySEXP(tracking_data$getDataInformations()$end - 24*3600)####last time in tracking minus 48 hours
time_stop  <- fmTimeCPtrFromAnySEXP(tracking_data$getDataInformations()$end) ####last time in tracking
###QUERY 3: fmQueryComputeAntTrajectories()
positions                 <- fmQueryComputeAntTrajectories(tracking_data,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
positions_summaries       <- positions$trajectory_summary
positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
positions_list            <- positions$trajectories

#make sure to define zones correctly. zones (ID = 1 Name = new-zone) and ( ID = 2 Name = new-zone) appear in data as "zone"=0 and 1
###THAT'S INCORRECT - THE FORAGING ZONE HAD NOT BEEN DEFINED PROPERLY AND CORRESPONDS TO ZONE 2! 0 IS FOR ANTS DETECTED OUTSIDE OF ANY ZONE. 
####ALSO. I WILL SHOW YOU HOW TO NAME THE ZONES AS THIS IS IMPORTANT TO FOOL-PROOF THINGS
zones <- tracking_data$cSpaces()$objects[[1]]$cZones()$summary #function to show the Zones present in the Space
foraging_zone <- zones[which(grepl("forag",zones$name)),"ID"] ##fool-proofing - this way we are sure to always select the right zone
print(paste("Foraging zone = zone",foraging_zone))
#for each trajectory, check whether the ant was observed in the foraging zone (positions_list$zone=2) or not. 
# if so the ant is a forager, if not the ant is a nurse
positions_summaries$AntTask <- NA

##before going back and forth between positions_summaries and positions_list:
####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
nrow(positions_summaries)
length(positions_list)
####2/ always make sure that positions_summaries is ordered correctly, using the index column
positions_summaries <- positions_summaries[order(positions_summaries$index),]
###this ensures that the first row in psoitions_summaries corresponds to the first trajectory in positions_list, etc.
for ( ant_index in 1:length(positions_list)) {
  positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]==foraging_zone))
  positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]!=foraging_zone))
  
  
    if ( foraging_zone %in% positions_list[[ant_index]][,"zone"]){
      positions_summaries[ant_index,"AntTask"] <- "forager"
    }else{
      positions_summaries[ant_index,"AntTask"] <- "nurse"
    }
}
#match antID and tagID (live tracking gives tagID). 
IDs <- tracking_data$identificationsAt(fmTimeNow(),FALSE)
IDs <- data.frame(tag_hex_ID=rownames(IDs), IDs["antID"],stringsAsFactors = F)
positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
positions_summaries1 <- positions_summaries  
  
#export nurse list: hexadecimal tagIDs corresponding to nurses in your data frame. 

rm(list=ls()[which(!ls()%in%c("positions_summaries1","myrmidon_file"))]) #close experiment
gc()



tracking_data <- fmExperimentOpenReadOnly(myrmidon_file) #if it says it's locked, click Session>Restart R and clear objects from workspace including hidden ones

####define time period to use for defining nurses/foragers
# time_start <- fmTimeParse("2020-11-04T18:49:21Z")
# time_stop <- fmTimeParse("2020-11-05T18:49:21Z")
###above are the times you chose but I recommend the following for the experiment (so it nevers need to be edited):
time_start <- fmTimeCPtrFromAnySEXP(tracking_data$getDataInformations()$end - 48*3600)####last time in tracking minus 48 hours
time_stop  <- fmTimeCPtrFromAnySEXP(tracking_data$getDataInformations()$end - 24*3600) ####last time in tracking
###QUERY 3: fmQueryComputeAntTrajectories()
positions                 <- fmQueryComputeAntTrajectories(tracking_data,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
positions_summaries       <- positions$trajectory_summary
positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
positions_list            <- positions$trajectories

#make sure to define zones correctly. zones (ID = 1 Name = new-zone) and ( ID = 2 Name = new-zone) appear in data as "zone"=0 and 1
###THAT'S INCORRECT - THE FORAGING ZONE HAD NOT BEEN DEFINED PROPERLY AND CORRESPONDS TO ZONE 2! 0 IS FOR ANTS DETECTED OUTSIDE OF ANY ZONE. 
####ALSO. I WILL SHOW YOU HOW TO NAME THE ZONES AS THIS IS IMPORTANT TO FOOL-PROOF THINGS
zones <- tracking_data$cSpaces()$objects[[1]]$cZones()$summary #function to show the Zones present in the Space
foraging_zone <- zones[which(grepl("forag",zones$name)),"ID"] ##fool-proofing - this way we are sure to always select the right zone
print(paste("Foraging zone = zone",foraging_zone))
#for each trajectory, check whether the ant was observed in the foraging zone (positions_list$zone=2) or not. 
# if so the ant is a forager, if not the ant is a nurse
positions_summaries$AntTask <- NA

##before going back and forth between positions_summaries and positions_list:
####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
nrow(positions_summaries)
length(positions_list)
####2/ always make sure that positions_summaries is ordered correctly, using the index column
positions_summaries <- positions_summaries[order(positions_summaries$index),]
###this ensures that the first row in psoitions_summaries corresponds to the first trajectory in positions_list, etc.
for ( ant_index in 1:length(positions_list)) {
  positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]==foraging_zone))
  positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]!=foraging_zone))
  if ( foraging_zone %in% positions_list[[ant_index]][,"zone"]){
    positions_summaries[ant_index,"AntTask"] <- "forager"
  }else{
    positions_summaries[ant_index,"AntTask"] <- "nurse"
  }
}
#match antID and tagID (live tracking gives tagID). 
IDs <- tracking_data$identificationsAt(fmTimeNow(),FALSE)
IDs <- data.frame(tag_hex_ID=rownames(IDs), IDs["antID"],stringsAsFactors = F)
positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
positions_summaries2 <- positions_summaries


#merge positions_summaries1 and positions_summaries2
names(positions_summaries2)[names(positions_summaries2)%in%c("nb_frames_outside","nb_frames_inside")] <- paste(names(positions_summaries2)[names(positions_summaries2)%in%c("nb_frames_outside","nb_frames_inside")],"_2",sep="")
positions_summaries <- merge(positions_summaries1[c("antID","tag_hex_ID","nb_frames_outside","nb_frames_inside")],
                             positions_summaries2[c("antID","tag_hex_ID","nb_frames_outside_2","nb_frames_inside_2")],
                             all.x=T,all.y=T
                             )
positions_summaries$prop_time_outside <- (positions_summaries$nb_frames_outside+positions_summaries$nb_frames_outside_2)/(positions_summaries$nb_frames_outside+positions_summaries$nb_frames_outside_2+positions_summaries$nb_frames_inside+positions_summaries$nb_frames_inside_2)
positions_summaries[which(positions_summaries$prop_time_outside<=0.01),"AntTask"] <- "nurse"
positions_summaries[which(positions_summaries$prop_time_outside>0.01),"AntTask"] <- "forager"


#export nurse list: hexadecimal tagIDs corresponding to nurses in your data frame. 
nurse_list <- data.frame(tag_hex_id=paste("- ", as.character(positions_summaries[which(positions_summaries$AntTask=="nurse"),"tag_hex_ID"] ),sep=""))


nurse_list
# write.table output syntax to be checked in the tracking room
write.table(nurse_list, file = paste(tracking_data$getDataInformations()[["details"]][["tdd.URI"]],"_nurses_list_1percent.txt",sep=""), append = FALSE,
            row.names = FALSE, col.names = FALSE,quote = FALSE) 

# #---------------Scraps------------------------
# 
# # 1 %in% positions_list[[1]][["zone"]]
# # 1 %in% positions_list[5] #zone is at number 5
# 
# 
# 
# for ( ant in positions_list) {
#   ifelse (0 %in% ant$zone, positions_summaries$AntTask[which(positions_summaries$index==ant)=="nurse"], positions_summaries$AntTask[which(positions_summaries$index==ant)=="bgbgbg"])
#   #ifelse ( 0 %in% ant$zone, positions_summaries$AntTask== "nurse",positions_summaries$AntTask== "forager")
# }
# 
# 
# #--------------------------------------attempt by separating the two functions
# Zone <- lapply(positions_list, function(GetColumn) { GetColumn[,5] }) #isolate needed element
# for (i in 1:length(Zone)){
#   ifelse ( 0 %in% Zone[[i]], print("nurse"), print("forager") )
# } #this sem to work but doesn't attach result to dataframe as needed