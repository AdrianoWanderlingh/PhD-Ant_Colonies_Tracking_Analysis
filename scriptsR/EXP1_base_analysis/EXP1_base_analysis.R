rm(list=ls())

##################################################################################################
############################# FUNCTIONS ##########################################################
# MOVE THEM TO ANOTHER SCRIPT OT AVOID CONFUSION AND SOURCE THEM


################## GET ANT TASK ###################################
getAntTasks <- function(exp){
  
  ## get 2 12Hours window for the Task calculation
  start <- fmQueryGetDataInformations(exp)$start + 33*3600 ####first time in tracking plus 21 hours, to skip acclimation time + 12 HOURS
  time_start <- fmTimeCreate(offset=start)
  #time_start <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end - 24*3600)####last time in tracking minus 24 hours
  stop  <- fmQueryGetDataInformations(exp)$start + 45*3600 ####pre-tracking period 
  time_stop   <- fmTimeCreate(offset=stop)
  #time_stop  <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end) ####last time in tracking
  ###QUERY 3: fmQueryComputeAntTrajectories()
  positions                 <- fmQueryComputeAntTrajectories(exp,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
  positions_summaries       <- positions$trajectories_summary
  positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
  positions_list            <- positions$trajectories
  
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
  IDs <- exp$identificationsAt(fmTimeNow())
  IDs <- data.frame(tag_hex_ID=rownames(IDs), IDs["antID"],stringsAsFactors = F)
  positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
  #positions_summaries1 <- positions_summaries
  
  positions_summaries1 <- aggregate(positions_summaries[ , 6:7], by = list(positions_summaries$antID), FUN = sum);  colnames(positions_summaries1) [match("Group.1",colnames(positions_summaries1))] <- "antID"
  
  
  #export nurse list: hexadecimal tagIDs corresponding to nurses in your data frame. 
  
  rm(list=ls()[which(!ls()%in%c("positions_summaries1","exp","foraging_zone"))]) #close experiment
  gc()
  
  #exp <- fmExperimentOpen(myrmidon_file) #if it says it's locked, click Session>Restart R and clear objects from workspace including hidden ones
  
  ####define time period to use for defining nurses/foragers
  start <- fmQueryGetDataInformations(exp)$start + 45*3600#### second block, make sure not to overlap with first time block or with exposure time
  time_start <- fmTimeCreate(offset=start)
  #time_start <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end - 24*3600)####last time in tracking minus 24 hours
  stop  <- fmQueryGetDataInformations(exp)$start + 57*3600 ####pre-tracking period 
  time_stop   <- fmTimeCreate(offset=stop)
  
  ###QUERY 3: fmQueryComputeAntTrajectories()
  positions                 <- fmQueryComputeAntTrajectories(exp,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
  positions_summaries       <- positions$trajectories_summary
  positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
  positions_list            <- positions$trajectories
  
  #for each trajectory, check whether the ant was observed in the foraging zone (positions_list$zone=2) or not. 
  # if so the ant is a forager, if not the ant is a nurse
  positions_summaries$AntTask <- NA
  
  ##before going back and forth between positions_summaries and positions_list:
  ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
  nrow(positions_summaries) ==length(positions_list)
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
  IDs <- exp$identificationsAt(fmTimeNow())
  IDs <- data.frame(tag_hex_ID=rownames(IDs), IDs["antID"],stringsAsFactors = F)
  positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
  #positions_summaries2 <- positions_summaries
  positions_summaries2 <- aggregate(positions_summaries[ , 6:7], by = list(positions_summaries$antID), FUN = sum);  colnames(positions_summaries2) [match("Group.1",colnames(positions_summaries2))] <- "antID"
  
  
  #merge positions_summaries1 and positions_summaries2
  names(positions_summaries2)[names(positions_summaries2)%in%c("nb_frames_outside","nb_frames_inside")] <- paste(names(positions_summaries2)[names(positions_summaries2)%in%c("nb_frames_outside","nb_frames_inside")],"_2",sep="")
  positions_summaries <- merge(positions_summaries1[c("antID","nb_frames_outside","nb_frames_inside")], # , "tag_hex_ID"
                               positions_summaries2[c("antID","nb_frames_outside_2","nb_frames_inside_2")], # , "tag_hex_ID"
                               all.x=T,all.y=T
  )
  positions_summaries[is.na(positions_summaries)] <- 0
  
  positions_summaries$prop_time_outside <- (positions_summaries$nb_frames_outside+positions_summaries$nb_frames_outside_2)/(positions_summaries$nb_frames_outside+positions_summaries$nb_frames_outside_2+positions_summaries$nb_frames_inside+positions_summaries$nb_frames_inside_2)
  positions_summaries[which(positions_summaries$prop_time_outside<=0.01),"AntTask"] <- "nurse"
  positions_summaries[which(positions_summaries$prop_time_outside>0.01),"AntTask"] <- "forager"
  
  positions_summaries$AntTask_num <- NA
  positions_summaries[which(positions_summaries$AntTask=="nurse"),"AntTask_num"] <- 1
  positions_summaries[which(positions_summaries$AntTask=="forager"),"AntTask_num"] <- 2
  
  AntTasks <- data.frame(antID=positions_summaries[,"antID"],AntTask= positions_summaries[,"AntTask"],AntTask_num= positions_summaries[,"AntTask_num"])
  
  warning("check that the N of returned ants is equal to length(exp$ants)! \nMaybe that's caused by a dead ant or something that has no coords in the time-span? \nit could be fixed by checking the $ants list and assignign to any missing ant an NA? ")
  return(AntTasks)
}

################## COMPUTE NETWORK ###############################
compute_G <- function(exp, start, end){ # min_cum_duration , link_type, nest_focus, frm_rate
  # convert timestamp of frame into corresponding frame number starting from 1 (with frame#1 at 'start' time)
  # CollideFrames <- fmQueryCollideFrames(exp,start=start,end=end)
  # TimeToFrame <- seq_along(CollideFrames$frames$time)
  Interactions <- fmQueryComputeAntInteractions(exp, start, end, maximumGap=gap, showProgress = TRUE, singleThreaded=FALSE, reportFullTrajectories = F)
  # Interactions$ant_ID1 <- paste("ant_",Interactions$ant1,sep="") ##creates a ID string for each ant: ant1, ant2,...
  # Interactions$ant_ID2 <- paste("ant_",Interactions$ant2,sep="")
  exp.Ants <- exp$ants
  Ant_IDs <- NULL
  # Number of ants
  for (ant in 1:length(exp.Ants)) {Ant_IDs <- c(Ant_IDs, exp.Ants[[ant]]$ID)}
  Max_antID <- max(Ant_IDs)
  # N_ants <- length(exp.Ants)    
  
  # initialise adj-matrix
  adj_matrix <- matrix(0L,nrow=Max_antID, ncol=Max_antID) #np.zeros((N_ants, N_ants))
  rownames(adj_matrix) <- Ant_IDs
  colnames(adj_matrix) <- Ant_IDs
  # Populate network
  for (INT in 1:nrow(Interactions)){
    ANT1 <- Interactions[INT,"ant1"]
    ANT2 <- Interactions[INT,"ant2"]
    
    # Populate adjaciency matrix
    #Choose “link as number of interactions”:
    adj_matrix[ANT1, ANT2] <- adj_matrix[ANT1, ANT2] + 1
    #or “link as total duration of interactions”:
    #  adj_matrix[ant_ID1, ant_ID2] = adj_matrix[ant_ID1, ant_ID2] + int$end – int$start
    
    # if link_type == 'length_inter':
    #   # OPT1
    #   # WEIGHTS: cumulative interaction time
    #   adj_mat[i.IDs[0]-1, i.IDs[1]-1] += (TimeToFrame[fm.Time.ToTimestamp(i.End)] - TimeToFrame[fm.Time.ToTimestamp(i.Start)]) * 1 / frm_rate
    # elif link_type == '#inter':
    #   # OPT2
    #   # WEIGHTS: number of interactions
    #   adj_mat[i.IDs[0]-1, i.IDs[1]-1] += 1
    # else:
    #   raise TypeError('"link_type" not valid')
  }
  
  adj_matrix <- adj_matrix + t(adj_matrix)
  
  # interaction filtering (remove weak connections)
  #adj_mat[adj_mat <  min_cum_duration] = 0
  
  # network build
  G <-  graph_from_adjacency_matrix(adj_matrix)
  
  # # store inverse of weights
  #ADRIANO: look for the function in igraph that works as set_edge_attributes in networkx of python
  # nx.set_edge_attributes(G, 
  #                        {(i,j): 1/adj_mat[j,i] if adj_mat[j,i]>0 else 0 for i in range(len(adj_mat)) for j in range(i)},
  #                        'inv_weight')
  
  return(G)
} # compute_G


##############################################################################

#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

##### LIBRARIES
library(FortMyrmidon) ####R bindings
library(data.table)
library(lubridate)
library(pals)
library(igraph)

#### PARAMETERS
TimeWind        <- 3600 ## in seconds (3600 is an hour)
gap             <- fmSecond(2)


#### FUNCTIONS
#list files recursive up to a certain level (level defined by "n" parameter)
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}

#calculate the distance between successive fixes
# FIX, look at Nathalie version
DISTANCE     <- function(x)  { c(sqrt((x[-nrow(x), "x"] - x[-1, "x"])^2 + (x[-nrow(x), "y"] - x[-1, "y"])^2), NA)}

#### ACCESS FILES
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA", n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]

#### OPEN REPLICATE
# replicate folder
for (REP.n in 1:length(files_list)) {
  # REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = "AntsCreated_AutoOriented_withMetaData")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  #replicate file
  for (REP.FILES in REP.filefolder) {
    # REP.FILES <-  REP.filefolder[2]   #temp
    print(REP.FILES) ##}}
    #open experiment
    exp <- fmExperimentOpen(REP.FILES)
    print("RENAME THE EXP TO \"e\", NOT \"exp\" (FUNCTIONS' VARIABLE NAME) ")
    # exp.Ants <- exp$ants
    
    ## get in R the spaceID / name correspondance
    ##PROBABLY NOT NEEDED. ANYWAY, IT WORKS DIFFERNTLY THAN IN FLORA-bee
    BoxCodes <- exp$spaces[[1]]$zones
    BoxCodes <- data.frame(space=c(BoxCodes[[1]]$name, BoxCodes[[2]]$name), box=c(exp$spaces[[1]]$name), stringsAsFactors = F )
    
    ########### DEFINE ZONES PROPERLY
    zones <- exp$spaces[[1]]$zones #function to show the Zones present in the Space
    zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
    foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
    nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
    print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))

    ########### COMPUTE THE ANT TASKS 
    AntTasks <- getAntTasks(exp)
    
    
    # ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
    # RawContacts_ALL <- NULL; AntMinuteMeans_ALL <- NULL
    # 
    # for (HOUR in seq(from=0, to=45, by=3))  ## increments by 3 hours for 48 hours 
    # {
    #   ## increment the window start & end times by 1 hour
    HOUR <- 0 #TEMP!!! REACTIVAT HOUR LOOP
      From  <- fmQueryGetDataInformations(exp)$end - 48*3600 + (HOUR * TimeWind)
      To    <- fmQueryGetDataInformations(exp)$end - 45*3600 + (HOUR * TimeWind)
      print(paste("Time window, from", From, "to", To))
      start <- fmTimeCreate(offset=From) #end minus 48 hours plus incremental time
      end   <- fmTimeCreate(offset=To) #end minus 45 hours plus incremental time

    #   ## extract 3 hours of contacts
    #   ## WARNING - THIS WILL CRASH IF YOU READ IN TOO MANY HOURS AT ONCE
    #   positions   <- fmQueryComputeAntInteractions(e, start=start, end=end, maximumGap=gap, showProgress = TRUE, singleThreaded=FALSE, reportFullTrajectories = T)#, reportLocalTajectories=FALSE) #reportTrajectories=TRUE
    #   
    #   positions_summaries       <- positions$trajectories_summary
    #   positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
    #   positions_list            <- positions$trajectories
    #   interactions              <- positions$interactions
    #   
    #   # ################### PART 1 : EXTRACTING TRAJECTORY SEGMENTS ##################################
    #   # 
    #   # 
    #   # ## for each ant, vertically stack the trajectory segments from the current HOUR
    #   # AntMinuteMeans_HOUR <- NULL
    #   # for (ID in unique(dd$trajectories_summary$antID))
    #   # {print(paste("stacking trajectory segments for ID", ID))
    #   #   AntIndices <- which(dd$trajectories_summary$antID == ID)     ## find which indices correspond to this ant
    #   #   Traj_Segments <- dd[["trajectories"]]       [AntIndices]  ## each traj segment for ID
    #   #   Traj_Starts   <- dd$trajectories_summary$start [AntIndices]  ## the UNIX start time for each traj segment
    #   #   ## for each trajectory segment,
    #   #   for (SG in 1:length (AntIndices))
    #   #   {
    #   #     Traj_Segments[[SG]]$dt   <- c(diff(Traj_Segments[[SG]]$time),NA) ## the max value will never be any more than 'gap'
    #   #     Traj_Segments[[SG]]$time <- ceiling_date( lubridate::seconds(Traj_Segments[[SG]]$time) + Traj_Starts[SG] , "min")  ## bin using 'round', so 11 min, 10 sec -> 11 min
    #   #     ## calculate the distance between successive fixes
    #   #     Traj_Segments[[SG]]$distance <- DISTANCE(x = Traj_Segments[[SG]])
    #   #     Traj_Segments[[SG]]$speed    <- Traj_Segments[[SG]]$distance / Traj_Segments[[SG]]$dt
    #   #   }
    #   #   ## collapse segments -> one data frame
    #   #   Traj_Segments <-  do.call("rbind", Traj_Segments)
    #   #   ## get mean & sd of speed within each rounded minute
    #   #   AntMinuteMeans          <- aggregate(speed ~ time, FUN=length,        Traj_Segments, na.action=na.pass); names(AntMinuteMeans)[match("speed",names(AntMinuteMeans))] <- "N_observations"
    #   #   AntMinuteMeans$speed    <- aggregate(speed ~ time, FUN=mean, na.rm=T, Traj_Segments, na.action=na.pass)$speed
    #   #   AntMinuteMeans$speed_sd <- aggregate(speed ~ time, FUN=sd,   na.rm=T, Traj_Segments, na.action=na.pass)$speed
    #   #   ## stack the minutely-means
    #   #   AntMinuteMeans_HOUR <- rbind(AntMinuteMeans_HOUR, data.frame(ID, AntMinuteMeans))
    #   # }
    #   # ## stack the hourly segments (faster than just adding each segment to a single object)
    #   # AntMinuteMeans_ALL <- rbind(AntMinuteMeans_ALL, AntMinuteMeans_HOUR)
    #   #
    #   # ## occasional time-series plotting
    #   # if (runif(1)<0.1)
    #   # {
    #   #   par(mfrow=c(2,1), mai=c(0.45,0.45,0.1,0.1), mgp=c(1.3,0.3,0), tcl=-0.2)
    #   #   ColonyACtivity <- aggregate(speed ~ round_date(time,"10 minutes"), FUN=mean, na.rm=T, AntMinuteMeans_ALL)
    #   #   ## Mean activity
    #   #   plot(ColonyACtivity, type="h", xaxs="i", xlab="", yaxt="n", ylab="")
    #   #   ## hacky colony-level heatmap
    #   #   plot(ColonyACtivity, type="n", xaxs="i", xlab="", yaxt="n", ylab="")
    #   #   MeanSpeedPerAnt <- aggregate(speed ~ round_date(time,"10 minutes") + ID, FUN=mean, na.rm=T, AntMinuteMeans_ALL); colnames(MeanSpeedPerAnt)[1] <- "time"; MeanSpeedPerAnt <- xtabs(speed ~ time + ID, MeanSpeedPerAnt)
    #   #   par(new=T)
    #   #   image(x=1:nrow(MeanSpeedPerAnt), y=1:ncol(MeanSpeedPerAnt), z=sqrt(MeanSpeedPerAnt), col=c("white",pals::tol.rainbow()), xaxt="n", yaxt="n", xlab="Time", ylab="Bee"); box()
    #   # }
    # 
    # 
    #   ################### PART 2 : EXTRACTING ANT-TO-ANT INTERACTIONS ##################################
    # 
    #   ##################################################################################################
    #   ################ NOTE: THIS CAN'T REALLY WORK AS FOR EVERY BLOCK OF 3 HOURS THERE WILL BE A DIFFERENT ASSIGNED ID ##
    #   ################ TO OVERCOME THIS, STORE SEPARATELY A DATAFRAME WITH THE nb_frames_outside AND nb_frames_inside AND CALCUULATE ROLE 
    #   ## PER EACH ANT ASSIGNING IT AT RAW_CONTACTS_ALL
    #   
    #   
    #   ########### IMPORTANT!!!!!!!!!!!!! MORE THAN ANT ROLE, GET THE DAMN SPACE IN THE FINAL RAWCONTACTS_ALL!!!!!
    #   # 
    #   # ##before going back and forth between positions_summaries and positions_list:
    #   # ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
    #   # nrow(positions$trajectories_summary)
    #   # length(positions_list)
    #   # ####2/ always make sure that positions_summaries is ordered correctly, using the index column
    #   # positions_summaries <- positions_summaries[order(positions_summaries$index),]
    #   # ###this ensures that the first row in psoitions_summaries corresponds to the first trajectory in positions_list, etc.
    #   # for ( ant_index in 1:length(positions_list)) {
    #   #   positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]==foraging_zone))
    #   #   positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]!=foraging_zone))
    #   #   
    #   #   
    #   #   if ( foraging_zone %in% positions_list[[ant_index]][,"zone"]){
    #   #     positions_summaries[ant_index,"AntTask"] <- "forager"
    #   #   }else{
    #   #     positions_summaries[ant_index,"AntTask"] <- "nurse"
    #   #   }
    #   # }
    #   # 
    #   # positions_summaries$prop_time_outside <- (positions_summaries$nb_frames_outside)/(positions_summaries$nb_frames_outside+positions_summaries$nb_frames_inside)
    #   # positions_summaries[which(positions_summaries$prop_time_outside<=0.01),"AntTask"] <- "nurse"
    #   # positions_summaries[which(positions_summaries$prop_time_outside>0.01),"AntTask"] <- "forager"
    #   # 
    #   # #################################################################################
    #   # ####ASSIGN ANT ROLE BY ADDING COLUMN THAT CHECKS IN POSITIONS_SUMMARY
    #   # #ASSING ROLE: TO BE COPIED IN THE FINAL FILE 
    #   # 
    #   # interactions$ant1_task <- 0
    #   # interactions$ant2_task <- 0
    #   # 
    #   # #loop and assign task
    #   # interactions$ant1_task <- positions_summaries[ant_index,"AntTask"]
    #   # interactions$ant2_task <- positions_summaries[ant_index,"AntTask"]
    # 
    #   ## assign box codes to the data frame
    #   #NOT NEEDED ~~~~~~~~~ IT ASSIGN BOX CODES IF YOU HAVE MULTIPLE TSs
    #   # RawContacts$box <- NA
    #   # RawContacts$box [ RawContacts$space==BoxCodes$space[1]] <- BoxCodes$box[1]
    #   # RawContacts$box [ RawContacts$space==BoxCodes$space[2]] <- BoxCodes$box[2]
    # 
    #   # RawContacts_ALL <- rbind(RawContacts_ALL, data.frame(HOUR, From, To, interactions))
    #   # print(paste("stacked Raw Contacts has",nrow(RawContacts_ALL),"rows"))
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
    #   
    #   
    # }## HOUR

  #COMPUTE NETWORK
  compute_Graph <- compute_G(exp = exp, start = start, end=end)
    
  ######################################################################################################
  #####################################################################################################
  ##### NETWORK PROPERTIES ######################################
      
      ##### ALL TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #Density: The proportion of present edges from all possible edges in the network.
  Density <- edge_density(compute_Graph, loops=F)
  
  ### to fill in!!!! #check Enrico or nath code
  Clustering_coeff  <- transitivity(compute_Graph,type ="barrat", isolates = c("NaN", "zero") )
    
  #Diameter: the longest geodesic distance (length of the shortest path between two nodes) in the network. In igraph, diameter() returns the distance, while get_diameter() returns the nodes along the first found path of that distance.
  #Note that edge weights are used by default, unless set to NA.
  # check enrico code, weights have to be 1/contact duration
  Diameter  <- diameter(compute_Graph, directed=F, weights=NA)
  
  #Network efficiency
  #unclear on what it is in igraph. 
  #Network efficiency: average connection efficiency of all pairs of nodes, where connection efficiency is the reciprocal of the shortest path length between the two nodes
  Net_eff  <- 1/shorteest-path-length-between-nodes
 
  #  Degree centrality:   degree of a vertex is its the number of its adjacent edges.
  Deg_centrality  <- degree(compute_Graph, mode="in") #CHECK for MODE!
  
  
  ### TO FIX AND ACTIVATE ONCE ANT TASKS ARE ASSIGNED!
  Assortativity  <- assortativity_nominal(compute_Graph, types=AntTasks$AntTask_num, directed=F)
  
  ## USE   Louvain (igraph::cluster_louvain) method  for  modularity optimization of weighted networks (see enrico script) 
  ###THIS IS NOT HOW MEMBERSHIP WORKS! (Numeric vector)
  Modularity <- modularity_matrix(compute_Graph, membership=AntTasks$AntTask_num ,weights = NULL, resolution = 1, directed = FALSE)
  
  
  
  ##############################################################################################
  ################### STEAL IDEAS AND STRUCTURE FROM ENRICO!!! BUT NO NEED TO TRANSLATE HIS CODE!######
  
  
  # function to compute netowrk properties
  def G_prop(compute_Graph, exp, start, end, time_win, max_gap): #, name, var, nest_focus, PLOT_HM_check = False
    
    # compute connencted components
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
  
  # Define Giant Component
  GC = G.subgraph(Gcc[0])
  res = 1
  
  # Best partition Louvian Method
  best_partition = nxc.greedy_modularity_communities(G, weight='weight', resolution=res)
  best_partition_CC = nxc.greedy_modularity_communities(GC, weight='weight', resolution=res)
  best_partition_CC_res_09 = nxc.greedy_modularity_communities(GC, weight='weight', resolution=0.9)
  if var != None:
    best_partition_CC_3p = nxc.greedy_modularity_communities(GC, weight='weight', cutoff=3, best_n=3, resolution=res)
  mode_p = mode_communities_dic[var['link_type']][name[11:14]][time_win]
  else:
    best_partition_CC_3p = nxc.greedy_modularity_communities(GC, weight='weight', cutoff=1, best_n=1, resolution=res)
  mode_p = 1
  
  best_partition_CC_mode_p = nxc.greedy_modularity_communities(GC, weight='weight', cutoff=mode_p, best_n=mode_p, resolution=res)
  
  
  # Heatmap partition plotting (optional)
  if PLOT_HM_check:
    
    # Save HM_partition Louvain
    directory = 'plots/HM_partition/unsup_mod/' + myrm_file[7:10]
  if not os.path.exists(directory):
    os.makedirs(directory)
  
  attributes = '_' + var['link_type'] + '_' + str(time_win / 3600) + 'h_' + str(var['h'])
  PLOT_HM_partition(compute_HM_stack(exp, start, end), 
                    best_partition_CC, 
                    directory + '/', 
                    name[7:15] + attributes)
  
  # Save HM_partition 3 partition
  directory = 'plots/HM_partition/sup_mod3p/' + myrm_file[7:10]
  if not os.path.exists(directory):
    os.makedirs(directory)
  
  attributes = '_' + var['link_type'] + '_' + str(time_win / 3600) + 'h_' + str(var['h'])
  PLOT_HM_partition(compute_HM_stack(exp, start, end), 
                    best_partition_CC_3p, 
                    directory + '/', 
                    name[7:15] + attributes)
  
  # Save HM_partition mode parition
  directory = 'plots/HM_partition/sup_modmp/' + myrm_file[7:10]
  if not os.path.exists(directory):
    os.makedirs(directory)
  
  attributes = '_' + var['link_type'] + '_' + str(time_win / 3600) + 'h_' + str(var['h'])
  PLOT_HM_partition(compute_HM_stack(exp, start, end), 
                    best_partition_CC_mode_p, 
                    directory + '/', 
                    name[7:15] + attributes)
  
  return {'rep': int(name[8:10]),
    'exp': name[11:15],
    'start': fm.Time.ToDateTime(start), 
    'time_win': time_win, 
    'nest_focus': nest_focus,
    'max_gap': max_gap,
    'GC': GC.number_of_nodes(),
    'ants': G.number_of_nodes(),
    'cMOD_communities': [len(best_partition_CC[i]) for i in range(len(best_partition_CC))],
    'cmpMOD_communities': [len(best_partition_CC_mode_p[i]) for i in range(len(best_partition_CC_mode_p))],
    'c3pMOD_communities': [len(best_partition_CC_3p[i]) for i in range(len(best_partition_CC_3p))],
    'cMODres09_communities': [len(best_partition_CC_res_09[i]) for i in range(len(best_partition_CC_res_09))],
    'MOD': nxc.modularity(G, best_partition),
    'cMOD': nxc.modularity(GC, best_partition_CC),
    'c3pMOD': nxc.modularity(GC, best_partition_CC_3p),
    'cmpMOD': nxc.modularity(GC, best_partition_CC_mode_p),
    'DEN': nx.density(G), 
    'wDEN': nx.adjacency_matrix(G).sum() / (G.number_of_nodes() * (G.number_of_nodes() - 1) * time_win),  # weighted density = sum all weights /(|V|*(|V|-1)/2 * time_win)
    'DIA': nx.diameter(GC),
    'wDIA': nx.diameter(GC, e=nx.eccentricity(GC, sp=dict(nx.shortest_path_length(GC,weight='inv_weight')))),
    'RAD': nx.radius(GC),
    'wRAD': nx.radius(GC, e=nx.eccentricity(GC, sp=dict(nx.shortest_path_length(GC,weight='inv_weight')))),
    'DEH': np.std([G.degree(n) for n in G.nodes()]),
    'cDEH': np.std([GC.degree(n) for n in GC.nodes()]),
    'wDEH': np.std(nx.adjacency_matrix(G).sum(axis=0)), # strength heterogeneity
    'cwDEH': np.std(nx.adjacency_matrix(GC).sum(axis=0)), # strength heterogeneity
    'CLS': np.mean([c for c in nx.clustering(G, weight='weight').values()]) 
  }
  
  # initialise data-frame with properties
  prop_df = pd.DataFrame(columns=G_prop(nx.star_graph(5),[],fm.Time.Now(),[],1,1,myrm_list[0],var=None, nest_focus=True).keys())
    
  
  
  
  
  ######################################################################################################
  #####################################################################################################
  
  
  
  
  
  
  } }#REP LOOP
    
    
    
    #####################################################################################
    #####################################################################################
    #####################################################################################
    # CONTINUE FROM HERE
    #####################################################################################
    #####################################################################################
    #####################################################################################

    
    
    
    
    ###################### ADD CONDITION HERE: IF THERE IS ONLY 1 SPACE (ALL ANTS ARE IN NEST), SKIP PART AND ASSIGN 100% OF TIME TO NEST)
    
    # ## calculate proportion of time in the foraging box
    # RawContacts_ALL$space_binary <- RawContacts_ALL$space - 1 ## only works if there are exactly 2 cameras
    # Contact_Time_Allocation      <- data.frame(ant=c(RawContacts_ALL$ant1,RawContacts_ALL$ant2), space_binary=c(RawContacts_ALL$space_binary,RawContacts_ALL$space_binary))
    # Box_Time_Allocation          <- aggregate(space_binary ~ ant, FUN=mean, Contact_Time_Allocation) ## mean of c(0,0,1,1,0,0)=1/3
    # 
    # ## count the contacts between each unique pair of ants
    # AggContacts       <- aggregate(end ~ ant1 + ant2, FUN=length, RawContacts_ALL); colnames(AggContacts)[match("end",colnames(AggContacts))] <- "weight"
    # ## count the contacts by the time of day & the box
    # AggContactsByHour <- aggregate(end ~ HOUR + box,  FUN=length, RawContacts_ALL); colnames(AggContactsByHour)[match("end",colnames(AggContactsByHour))] <- "N_contacts"
    # 
    # 
    # ## calculate proportion of time in the foraging box
    # RawContacts_ALL$space_binary <- RawContacts_ALL$space - 1 ## only works if there are exactly 2 cameras
    # Contact_Time_Allocation      <- data.frame(ant=c(RawContacts_ALL$ant1,RawContacts_ALL$ant2), space_binary=c(RawContacts_ALL$space_binary,RawContacts_ALL$space_binary))
    # 
    # ## count the contacts between each unique pair of ants
    # AggContacts       <- aggregate(end ~ ant1 + ant2, FUN=length, RawContacts_ALL); colnames(AggContacts)[match("end",colnames(AggContacts))] <- "weight"
    # ## count the contacts by the time of day & the box
    # AggContactsByHour <- aggregate(end ~ HOUR + box,  FUN=length, RawContacts_ALL); colnames(AggContactsByHour)[match("end",colnames(AggContactsByHour))] <- "N_contacts"
    # 
    # ## create an overall time-aggregated contact network from the counts
    # G <- graph_from_data_frame(d=AggContacts, directed=F); is.weighted(G)
    # 
    # ## calculate the proportion of time spent in the foraging box
    # V(G)$space_binary <- Box_Time_Allocation $ space_binary [match(V(G)$name, Box_Time_Allocation$ant)]
    # 
    # ## node centrality - use to size the nodes
    # V(G)$degree          <- degree(G, mode="all")
    # V(G)$weighted_degree <- graph.strength(G, mode="all")
    # 
    # ## threshold edges - just for plotting
    # G2 <- delete.edges(graph=G, edges=E(G) [E(G)$weight < quantile(E(G)$weight,0.5)] )  ## remove 50% of the weakest edges, for plotting only
    # ## spring-embedded graph layout
    # Layout <- layout_with_fr(graph=G2, weights= log10(E(G2)$weight) ) ##
    # 
    # 
    # 
    # ## PLOTTING
    # par(mfrow=c(2,2), mai=c(0.4,0.4,0.1,0.1), tcl=-0.2, mgp=c(1.3,0.3,0))
    # ## time-series
    # plot  (N_contacts ~ HOUR, AggContactsByHour[AggContactsByHour$box=="Foraging",], type="b", pch=21, bg=1, ylab="N contacts / hour", ylim=range(AggContactsByHour$N_contacts))
    # points(N_contacts ~ HOUR, AggContactsByHour[AggContactsByHour$box=="Nest",],     type="b", pch=21, bg=2)
    # ## distribution of proportion contacts in the foraging box
    # hist(Box_Time_Allocation$space_binary, breaks=25, main="", xlab="Proportion of time in foraging arena", ylab="N bees", col=1)
    # ## degree ~ prop time in forage arena
    # plot(V(G)$weighted_degree ~ V(G)$space_binary, xlab="Proportion of time in foraging arena", ylab="N contacts")
    # ## plot the contact network
    # plot(G, Layout, 
    #      vertex.size= 5 + (8 * V(G)$weighted_degree/max(V(G)$weighted_degree)), 
    #      vertex.color=parula(11)[10 * round(V(G)$space_binary,1)],
    #      vertex.label=NA, 
    #      edge.color=rgb(0.5,0.5,0.5,0.75,maxColorValue =1),
    #      edge.width=2*E(G)$weight/mean(E(G)$weight),
    #      edge.curved=0.2)
    # ## mean of c(0,0,1,1,0,0)=1/3
    # 
    # 

    
  }
}
