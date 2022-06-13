rm(list=ls())

##################################################################################################
############################# FUNCTIONS ##########################################################
# MOVE THEM TO ANOTHER SCRIPT OT AVOID CONFUSION AND SOURCE THEM

################## GET ANT TASK ###################################
AntTasks.ZoneUse <- function(exp,window_shift){
  print("Computing AntTasks based on 24h time-window before exposure")
  #Get complete list of Ants
  AntID_list <- NULL
  for (ant in   1: length(exp$ants)) {
    AntID_list <- c(AntID_list,exp$ants[[ant]]$ID)}
  
  ## get 2 12Hours window for the Task calculation
  ## calcualte the task BEFORE the EXPOSURE
  start <- fmQueryGetDataInformations(exp)$end - 51*3600 - window_shift
  #start <- fmQueryGetDataInformations(exp)$start + 33*3600 ####first time in tracking plus 21 hours, to skip acclimation time + 12 HOURS
  time_start <- fmTimeCreate(offset=start)
  #time_start <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end - 24*3600)####last time in tracking minus 24 hours
  stop <- fmQueryGetDataInformations(exp)$end - 39*3600 - window_shift
  #stop  <- fmQueryGetDataInformations(exp)$start + 45*3600 ####pre-tracking period 
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
  
  
  rm(list=ls()[which(!ls()%in%c("positions_summaries1","exp","foraging_zone","window_shift","AntID_list"))]) #close experiment
  gc()
  
  #exp <- fmExperimentOpen(myrmidon_file) #if it says it's locked, click Session>Restart R and clear objects from workspace including hidden ones
  
  ####define time period to use for defining nurses/foragers
  start <- fmQueryGetDataInformations(exp)$end - 39*3600 - window_shift
  #start <- fmQueryGetDataInformations(exp)$start + 45*3600#### second block, make sure not to overlap with first time block or with exposure time
  time_start <- fmTimeCreate(offset=start)
  #time_start <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end - 24*3600)####last time in tracking minus 24 hours
  stop <- fmQueryGetDataInformations(exp)$end - 27*3600 - window_shift
  #stop  <- fmQueryGetDataInformations(exp)$start + 57*3600 ####pre-tracking period 
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
  
  AntTasks <- data.frame(antID=positions_summaries[,"antID"],AntTask= positions_summaries[,"AntTask"])
  
  print("AntTasks computed")
  
  # #add missing ants as NURSE by DEFAULT
  # missing_ants <- subset(AntID_list, !(AntID_list %in% AntTasks$antID))
  # AntTasks <- rbind(AntTasks, data.frame(antID=missing_ants,AntTask=NA))
  # AntTasks[which(is.na(AntTasks$AntTask)),"AntTask"] <- "nurse"

  AntTasks$AntTask_num <- NA
  AntTasks[which(AntTasks$AntTask=="nurse"),"AntTask_num"] <- 1
  AntTasks[which(AntTasks$AntTask=="forager"),"AntTask_num"] <- 2
  AntTasks <- AntTasks[order(AntTasks$antID),]
  
  rm(list=ls()[which(!ls()%in%c("positions_summaries1","positions_summaries2","exp","foraging_zone","window_shift","AntID_list","AntTasks"))]) #close experiment
  gc()
  
  ################################## ZONE USE
  
  #if (ZoneUsage) {
    print("Computing Zone (nest, foraging area) usage pre-post exposure")
    
    ## get 2 12Hours window for the Task calculation
    ## calcualte the task BEFORE the EXPOSURE
    start <- fmQueryGetDataInformations(exp)$end - 24*3600 - window_shift
    time_start <- fmTimeCreate(offset=start)
    stop <- fmQueryGetDataInformations(exp)$end - 12*3600 - window_shift
    time_stop   <- fmTimeCreate(offset=stop)
    ###QUERY 3: fmQueryComputeAntTrajectories()
    positions                 <- fmQueryComputeAntTrajectories(exp,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
    positions_summaries       <- positions$trajectories_summary
    positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
    positions_list            <- positions$trajectories
    
    
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
    }
    #match antID and tagID (live tracking gives tagID). 
    IDs <- exp$identificationsAt(fmTimeNow())
    IDs <- data.frame(tag_hex_ID=rownames(IDs), IDs["antID"],stringsAsFactors = F)
    positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
    positions_summaries3 <- aggregate(positions_summaries[ , 5:6], by = list(positions_summaries$antID), FUN = sum);  colnames(positions_summaries3) [match("Group.1",colnames(positions_summaries3))] <- "antID"
    
    names(positions_summaries3)[names(positions_summaries3)%in%c("nb_frames_outside","nb_frames_inside")] <- paste(names(positions_summaries3)[names(positions_summaries3)%in%c("nb_frames_outside","nb_frames_inside")],"_3",sep="")
    
    
    rm(list=ls()[which(!ls()%in%c("positions_summaries1","positions_summaries2","positions_summaries3","exp","foraging_zone","window_shift","AntID_list","AntTasks"))]) #close experiment
    gc()
    
    ####define time period to use for defining nurses/foragers
    start <- fmQueryGetDataInformations(exp)$end - 12*3600 - window_shift
    time_start <- fmTimeCreate(offset=start)
    stop <- fmQueryGetDataInformations(exp)$end  - window_shift
    time_stop   <- fmTimeCreate(offset=stop)
    
    ###QUERY 3: fmQueryComputeAntTrajectories()
    positions                 <- fmQueryComputeAntTrajectories(exp,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
    positions_summaries       <- positions$trajectories_summary
    positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
    positions_list            <- positions$trajectories
    
    ##before going back and forth between positions_summaries and positions_list:
    ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
    nrow(positions_summaries) ==length(positions_list)
    ####2/ always make sure that positions_summaries is ordered correctly, using the index column
    positions_summaries <- positions_summaries[order(positions_summaries$index),]
    ###this ensures that the first row in psoitions_summaries corresponds to the first trajectory in positions_list, etc.
    for ( ant_index in 1:length(positions_list)) {
      positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]==foraging_zone))
      positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]!=foraging_zone))
    }
    #match antID and tagID (live tracking gives tagID). 
    IDs <- exp$identificationsAt(fmTimeNow())
    IDs <- data.frame(tag_hex_ID=rownames(IDs), IDs["antID"],stringsAsFactors = F)
    positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
    #positions_summaries2 <- positions_summaries
    positions_summaries4 <- aggregate(positions_summaries[ , 5:6], by = list(positions_summaries$antID), FUN = sum);  colnames(positions_summaries4) [match("Group.1",colnames(positions_summaries4))] <- "antID"
    
    
    #merge positions_summaries1 and positions_summaries2
    names(positions_summaries4)[names(positions_summaries4)%in%c("nb_frames_outside","nb_frames_inside")] <- paste(names(positions_summaries4)[names(positions_summaries4)%in%c("nb_frames_outside","nb_frames_inside")],"_4",sep="")
   
    #put all data frames into list
    psummaries_list <- list(positions_summaries1, positions_summaries2, positions_summaries3,positions_summaries4)      
    
    #merge all data frames together
    positions_summaries <- Reduce(function(x, y) merge(x, y, all=TRUE), psummaries_list)
    positions_summaries[is.na(positions_summaries)] <- 0
    
    positions_summaries$frames_outside_PRE <- positions_summaries$nb_frames_outside + positions_summaries$nb_frames_outside_2
    positions_summaries$frames_outside_POST <- positions_summaries$nb_frames_outside_3 + positions_summaries$nb_frames_outside_4
    
    positions_summaries$frames_inside_PRE <- positions_summaries$nb_frames_inside + positions_summaries$nb_frames_inside_2
    positions_summaries$frames_inside_POST <- positions_summaries$nb_frames_inside_3 + positions_summaries$nb_frames_inside_4

    ZoneUse <- positions_summaries[c("antID","frames_inside_PRE","frames_inside_POST","frames_outside_PRE","frames_outside_POST")]
    
    #merge this output with the AntTasks dataframe
    output_summ_list <- list(AntTasks,ZoneUse)
    
    AntTasks <- Reduce(function(x, y) merge(x, y, all=TRUE), output_summ_list)
    
    ############## MODIFY! KEEP SEPARATED THE PRE AND POST OUTPUTS!!!
    
     #}# ZoneUSe
  
    rm(list=ls()[which(!ls()%in%c("AntTasks"))]) #close experiment
    gc()
    
  ##RETURN OUTPUT
  # warning("Ants that died before the considered time window (pre treatment) will not be assigned a Task by the function. Currently, no task will default to Nurse")
  return(AntTasks)
  
} #, ZoneUsage= TRUE

################## COMPUTE NETWORK ###############################
# it require a "gap" to be defined, should be required by the function
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
  G <-  graph_from_adjacency_matrix(adj_matrix,mode = "undirected")
  actors <- V(G)
  # # store inverse of weights
  #ADRIANO: look for the function in igraph that works as set_edge_attributes in networkx of python
  # nx.set_edge_attributes(G, 
  #                        {(i,j): 1/adj_mat[j,i] if adj_mat[j,i]>0 else 0 for i in range(len(adj_mat)) for j in range(i)},
  #                        'inv_weight')
  
  #### add a column contaning interaction duration in min
  Interactions["duration_min"] <- as.numeric(difftime(Interactions$end, Interactions$start, units = "mins") + 0.125) ###duration in minutes (one frame = 0.125 second)
  ### add edge weights
  E(G)$weight <- Interactions[,"duration_min"]
  ###simplify graph (merge all edges involving the same pair of ants into a single one whose weight = sum of these weights)
  G <- simplify(G,remove.multiple=TRUE,remove.loop=TRUE,edge.attr.comb="sum")
  ##################remove unconnected nodes
  unconnected <-  actors[degree(G)==0]
  G <- G - as.character(unconnected)
  ##################update actor list
  actors <- get.vertex.attribute(G,"name")
  
  # Assign vertex types: "AntTask"
  # include in function eventually
  if (exists("AntTasks")) {
    #create matching of vertices with ants (there could be less ants in a 3-hours window than in the full exp)
    #keep only antTasks corresponding to vertices
    V_AntTasks <- AntTasks[which(AntTasks$antID %in% V(G)$name),]
    G <- set_vertex_attr(G, name="AntTask", index = V(G), value = V_AntTasks$AntTask_num)
  }
  
  return(G)
} # compute_G 

################## COMPUTE NETWORK PROPERTIES ####################
NetProperties <- function(graph){
  
  summary_collective <- NULL
  
  ##### COLLECTIVE NETWORK PROPERTIES ######################################
  #### inherited from Stroeymeyt et al. 2018
  
  ## Assortativity - Task
  #degree of preferential association between workers of the same task group, calculated using Newman’s method
  task_assortativity  <- assortativity_nominal(graph, types= V(graph)$AntTask, directed=F)
  ##Clustering
  clustering <- mean(transitivity(graph,type="barrat",weights=E(graph)$weight,isolates = c("NaN")),na.rm=T)
  ##Degree mean and max
  # Degree centrality:   degree of a vertex is its the number of its adjacent edges.
  degrees         <- degree(graph,mode="all")
  degree_mean     <- mean(degrees,na.rm=T)
  degree_maximum  <- max(degrees,na.rm=T)
  ##Density
  #Density: The proportion of present edges from all possible edges in the network.
  density  <- igraph::edge_density(graph)
  ##Diameter
  #Diameter: the longest geodesic distance (length of the shortest path between two nodes) in the network. In igraph, diameter() returns the distance, while get_diameter() returns the nodes along the first found path of that distance.
  diameter <- igraph::diameter(graph,directed=F,unconnected=TRUE,weights=(1/E(graph)$weight)) ###here use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
  ##Efficiency
  #Network efficiency: average connection efficiency of all pairs of nodes, where connection efficiency is the reciprocal of the shortest path length between the two nodes
  net_dist                    <- shortest.paths(graph, weights=1/E(graph)$weight, mode="all") ##again use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
  net_dist[net_dist==0]       <- NA ##remove distances to self
  efficiency                  <- 1/net_dist ##transform each distance into an efficiency
  efficiency <- (1/((vcount(graph)*(vcount(graph)-1))))*(sum(efficiency,na.rm=TRUE))
  ## Modularity
  communities             <- cluster_louvain(graph, weights = E(graph)$weight)
  community_membership    <- communities$membership
  modularity              <- modularity(graph,community_membership,weights=E(graph)$weight)
  
  
  
  #FINAL OUTPUT #DATAFRAME with the network properties per each 3 hours timeslot
  summary_collective <- rbind(summary_collective,data.frame(randy=REP.FILES,colony=COLONY,colony_size=COLONY_SIZE,treatment=TREATMENT,period=PERIOD,time_hours=HOUR, From, To,#time_of_day=time_of_day,
                                                            task_assortativity=task_assortativity,
                                                            clustering=clustering,
                                                            degree_mean=degree_mean,
                                                            degree_maximum=degree_maximum,
                                                            density=density,
                                                            diameter=diameter,
                                                            efficiency=efficiency,
                                                            modularity=modularity,stringsAsFactors = F))
  
  return(summary_collective)
  
  ########### EXPANSION ##########################
  ####Part 2: individual network properties ####
  # look at line 218 on from /13_network_analysis.R
  # (path length to queen, etc...)
  
}


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
gap             <- fmSecond(10)
## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED! 
window_shift <- 60*10 #approx N of minutes that where given at the end as leeway, can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition

## some initialization
Period_dataframe <- NULL #checking time correspondances

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


#### ACCESS FILES
WORKDIR <- "/media/cf19810/DISK4/ADRIANO"
DATADIR <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/")

#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]

#start fresh
Network_properties            <- data.frame()

###initialise general output folder
###remove folder if already exists to make sure we don't mix things up
if (file.exists(file.path(DATADIR, "NetworkAnalysis_outcomes"))){
  unlink(file.path(DATADIR, "NetworkAnalysis_outcomes"),recursive=T)
} 
###create folder
dir.create(file.path(DATADIR, "NetworkAnalysis_outcomes"),recursive = T)
###define name of general output table containing Network_properties
output_name <- file.path(DATADIR, "NetworkAnalysis_outcomes","Network_properties.txt")
AntTasks_SpaceUse <- file.path(DATADIR, "NetworkAnalysis_outcomes","AntTasks_SpaceUse.txt")

##### RUNNING TIME
loop_start_time <- Sys.time()

####define to_keep variables to keep clearing memory between runs
to_keep <- c(ls(),c("to_keep"))

#little improvement to do
print("RENAME THE EXP TO \"e\", NOT \"exp\" (FUNCTIONS' VARIABLE NAME) ")

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
    # exp.Ants <- exp$ants
    exp_end <- fmQueryGetDataInformations(exp)$end - window_shift
    
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

    ########### COMPUTE THE ANT TASKS (24h before exposure)  and the zone use (pre-post exposure)
    AntTasks <- AntTasks.ZoneUse(exp=exp,window_shift=window_shift)
    
    ########## GET EXPOSED ANTS 
    exp.Ants <- exp$ants
    AntTasks$Exposed <- "no"
    # check who's the queen
    for (ant in exp.Ants){
      individual  <- ant$ID
      #print(ant)
        if (TRUE %in% ant$getValues("Exposed")[,"values"]) { AntTasks[individual,"Exposed"] <- "exposed" }
    }

    #base file info
    #TREATMENT
    COLONY <- sub("\\_.*", "", basename(REP.FILES))
    TREATMENT <- substr(COLONY,(nchar(COLONY)+1)-2,nchar(COLONY))
    #
    COLONY_SIZE <- nrow(AntTasks)

    # ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
    for (HOUR in seq(from=0, to=48, by=3)){  ## increments by 3 hours for 48 hours
      
    #   ## increment the window start & end times by 1 hour
    #HOUR <- 0 #TEMP!!! REACTIVAT HOUR LOOP
      From  <- fmQueryGetDataInformations(exp)$end - 51*3600 + (HOUR * TimeWind) - window_shift
      To    <- fmQueryGetDataInformations(exp)$end - 48*3600 + (HOUR * TimeWind) - window_shift
      print(paste0("computing hour ",HOUR))
      print(paste("Time window, from", From, "to", To))
      start <- fmTimeCreate(offset=From) #end minus 48 hours plus incremental time
      end   <- fmTimeCreate(offset=To) #end minus 45 hours plus incremental time
      
      #base file information
      #PERIOD
      TimeDiff <- difftime(exp_end, To, units = "hours")
      if(TimeDiff < 24){ PERIOD <- "POST"
      }else if ( TimeDiff >= 27 & TimeDiff < 51) { PERIOD <- "PRE"  } else{ PERIOD <- "EXPOSURE_GAP"}
  
      Period_dt<- data.frame(From, To, PERIOD)
      Period_dataframe <- rbind(Period_dataframe, Period_dt)
      
    # }      
    # 
    # table(Period_dataframe$PERIOD) # shall be equal!
    
    # RUN FOR PRE AND POST (skip the 3h exposure gap)
    if (!Period_dt$PERIOD=="EXPOSURE_GAP") {

  #COMPUTE NETWORK
  Graph <- compute_G(exp = exp, start = start, end=end)

  # COMPUTE NETWORK PROPERTIES
  Network_prop_hour <- NetProperties(graph=Graph)
  
  Network_properties <- rbind(Network_properties,Network_prop_hour)
  
  } }#REP LOOP
    
    
    
    #keep relevant exp info
    AntTasks <- data.frame(randy=REP.FILES,colony=COLONY,colony_size=COLONY_SIZE,treatment=TREATMENT, AntTasks)
    
    ########################################
    ### prop of time Exposed ants spend inside the nest, pre-post exposure
    AntTasks$prop_inside_PRE <- AntTasks$frames_inside_PRE /(AntTasks$frames_inside_PRE + AntTasks$frames_outside_PRE)
    AntTasks$prop_inside_POST <- AntTasks$frames_inside_POST /(AntTasks$frames_inside_POST + AntTasks$frames_outside_POST)
    AntTasks$delta_time_inside <- AntTasks$prop_inside_POST - AntTasks$prop_inside_PRE
    
    ########################################
    ##### SAVE FILES IN FOLDER #############

    ## Network properties save
    if (file.exists(output_name)){
      write.table(Network_properties,file=output_name,append=T,col.names=F,row.names=F,quote=T)
    }else{
      write.table(Network_properties,file=output_name,append=F,col.names=T,row.names=F,quote=T)
    }
    
    ## AntTasks save
    if (file.exists(AntTasks_SpaceUse)){
      write.table(AntTasks,file=AntTasks_SpaceUse,append=T,col.names=F,row.names=F,quote=T)
    }else{
      write.table(AntTasks,file=AntTasks_SpaceUse,append=F,col.names=T,row.names=F,quote=T)
    }
    
    
    #cleaning
    rm(   list   =  ls()[which(!ls()%in%to_keep)]    )
    gc()
    
    
    ######################################################
    ### PLOTTING (1 col..., should be all)
    #plot(AntTasks$delta_time_inside, col=as.factor(AntTasks$Exposed))
    
    
    ####################################################################################
    ### plot MEAN DELTA PER ANT SEPARATED BETWEEN NON EXPOSED AND EXPOSED, WITH COL. SIZE COMPARISON
    ###################################################################################
    
    
    ### SEE NATHALIE WORK FOR PLOTTING INSPIRATION
    
    
    
    
    
    
    
    
    
    
    ############## RANDOM PLOTTING COPIED FROM BEE_FLORA
    
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
    
  }
}


loop_end_time <- Sys.time()
print (paste("loop took ",as.numeric(difftime(loop_end_time, loop_start_time, units = "mins"))," minutes to complete"))

