
################## GET ZONES usage ###################################
SpaceUse <- function(e, start, end){
  
  #required packages
  require(lubridate)
  require(FortMyrmidon)
  require(mallinfo)
  
  #define zones
  #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
  zones <- e$spaces[[1]]$zones #function to show the Zones present in the Space
  zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
  
  
  foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
  nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
  print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))
  
  #Get complete list of Ants
  AntID_list <- NULL
  for (ant in   1: length(e$ants)) {
    AntID_list <- c(AntID_list,e$ants[[ant]]$ID)}

    #   print(paste(" Time interval from hour",HOUR_start,"to",(HOUR_start-TimeWindow),"until e end",sep = " "))
    positions                 <- fmQueryComputeAntTrajectories(e,start,end,maximumGap = fmHour(24*365),computeZones = TRUE,singleThreaded=FALSE)
    positions_summaries       <- positions$trajectories_summary
    positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
    positions_list            <- positions$trajectories
    
    ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
    if(!(nrow(positions_summaries)==length(positions_list)))stop("number of rows of positions_summaries is not equal to the length of positions_list")
    ####2/ always make sure that positions_summaries is ordered correctly, using the index column
    positions_summaries <- positions_summaries[order(positions_summaries$index),]
    ###this ensures that the first row in positions_summaries corresponds to the first trajectory in positions_list, etc.
    for ( ant_index in 1:length(positions_list)) {
      positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]%in%foraging_zone))
      positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]%in%nest_zone))
      positions_summaries[ant_index,"prop_time_outside"] <- positions_summaries[ant_index,"nb_frames_outside"] /(positions_summaries[ant_index,"nb_frames_outside"] + positions_summaries[ant_index,"nb_frames_inside"])
      positions_summaries[ant_index,"proportion_time_active"] 
      warning("continue from here")
      #is going to be a bit long:
      # get from *_process trajectories all the elements to "2. cut trajectory into bouts of activity vs. inactivity"
      # get sums: sum(bouts$total_dist,na.rm=T)
      positions_summaries[ant_index,"average_bout_speed_pixpersec"]
      positions_summaries[ant_index,"total_distance_travelled_pix"]
      # prop_time_outside = ,
      # proportion_time_active = ,
      # average_bout_speed_pixpersec = ,
      # total_distance_travelled_pix = ,
      }
    #match antID and tagID (live tracking gives tagID). 
    IDs <- e$identificationsAt(fmTimeCreate(offset=fmQueryGetDataInformations(e)$start))
    IDs[sapply(IDs, is.null)] <- NA # assign NA to dead ants
    IDs <- data.frame(tag_hex_ID=unlist(IDs), antID=1:length(IDs),stringsAsFactors = F)
    positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
    
    SpaceUsage <- aggregate(cbind(nb_frames_outside,nb_frames_inside) ~ antID + tag_hex_ID, FUN = sum, na.rm=T,na.action=na.pass,positions_summaries )
   
    # rm(list=ls()[which(!ls()%in%c("SpaceUsage","e","foraging_zone","nest_zone","AntID_list","hour_chunk_start","TimeWindow","loop_N"))]) #close experiment
    # gc()
    # mallinfo::malloc.trim(0L)
    # 
    
  # #add missing ants
  missing_ants <- subset(AntID_list, !(AntID_list %in% SpaceUsage$antID))
  missing_ants_table <- data.frame()
  
  for (MISSING in missing_ants) {
  for (id in length(e$ants[[MISSING]]$identifications)) {
   #print ( e$ants[[MISSING]]$identifications[[id]]$tagValue )
    missing_ants_table <- rbind(missing_ants_table, data.frame(antID=MISSING, tag_hex_ID= e$ants[[MISSING]]$identifications[[id]]$tagValue))
  }}
  
  if (nrow(missing_ants_table)>0) {
  #add empty cols
  missing_ants_table[setdiff(names(SpaceUsage),names(missing_ants_table))] <- NA
  SpaceUsage <- rbind(SpaceUsage, missing_ants_table)
  }
  SpaceUsage <- SpaceUsage[order(SpaceUsage$antID),]
  #rename antId to Tag to comply with Science2018
  colnames(SpaceUsage)[which(colnames(SpaceUsage)=="antID")] <- "Tag"
  SpaceUsage$tag_hex_ID <- NULL
  
  rm(list=ls()[which(!ls()%in%c("SpaceUsage"))]) #close experiment
  gc()
  mallinfo::malloc.trim(0L)
  
  #print("SpaceUse computed")
  ##RETURN OUTPUT
  return(SpaceUsage)

}
