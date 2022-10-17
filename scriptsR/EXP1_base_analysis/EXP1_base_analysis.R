rm(list=ls())

##################################################################################################
############################# FUNCTIONS ##########################################################
# MOVE THEM TO ANOTHER SCRIPT OT AVOID CONFUSION AND SOURCE THEM

################## GET ANT TASK ###################################

## IMPORTANT CONSIDERATIONS:

# - DATA HAS TO BE IN CHUNKS OF 3 HOURS!!
# - time_hours ranges from -30 to 21, with 0 being the first hour post exposure and -3 missing as it is the EXPOSURE_TIME
# - the computation could be sped-up only doing calculations for nurses? BUT it is USEFUL to keep the comparison between treated and non treated


########################################################################################################
################### IMPORTANT THINGS TO MODIFY #########################################################
# ADD THE FOLLOWING VARS:
# period ("after","before")
# SHOULD BE ASSIGNED IN THE ANT_TASK FUNCTION

# time_hours
# FIND A WAY TO ASSIGN TIMING TO ANT_TASK AND UPDATE NETWORK_ANALYSIS TIMING TO MATCH THE Time_dictionary

# time_of_day
# SIMPLY USE THE Time_dictionary AS GUIDE -> EVERY time_hours HAS A SPECIFIC time_of_day BUT THE OPPOSITE IS NOT TRUE (MULTIPLE DAYS PRESENT, REPETED time_of_day VALUES)

# CHECK Plot_Grooming_Pre-Post.R to see how these things are used there!


################## GET ANT TASK ###################################
# AntTasks.ZoneUse <- function(exp,window_shift){
#   print("Computing AntTasks based on 48h time-window before exposure")
#   #Get complete list of Ants
#   AntID_list <- NULL
#   for (ant in   1: length(exp$ants)) {
#     AntID_list <- c(AntID_list,exp$ants[[ant]]$ID)}
#   
#   
#   hour_chunk_start <- c(75,63,51,39)
#   positions_summaries_list <- list()
#   loop_N <- 0
#   
#   for (HOUR_start in hour_chunk_start) {
#     loop_N <- loop_N + 1 
#     ## get 2 12Hours window for the Task calculation
#     ## calcualte the task BEFORE the EXPOSURE
#     start <- fmQueryGetDataInformations(exp)$end - HOUR_start*3600 - window_shift
#     #start <- fmQueryGetDataInformations(exp)$start + 33*3600 ####first time in tracking plus 21 hours, to skip acclimation time + 12 HOURS
#     time_start <- fmTimeCreate(offset=start)
#     #time_start <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end - 24*3600)####last time in tracking minus 24 hours
#     stop <- fmQueryGetDataInformations(exp)$end - (HOUR_start-12)*3600 - window_shift
#     #stop  <- fmQueryGetDataInformations(exp)$start + 45*3600 ####pre-tracking period 
#     time_stop   <- fmTimeCreate(offset=stop)
#     #time_stop  <- fmTimeCPtrFromAnySEXP(exp$getDataInformations()$end) ####last time in tracking
#     ###QUERY 3: fmQueryComputeAntTrajectories()
#     positions                 <- fmQueryComputeAntTrajectories(exp,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
#     positions_summaries       <- positions$trajectories_summary
#     positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
#     positions_list            <- positions$trajectories
#     
#     #for each trajectory, check whether the ant was observed in the foraging zone (positions_list$zone=2) or not. 
#     # if so the ant is a forager, if not the ant is a nurse
#     positions_summaries$AntTask <- NA
#     
#     ##before going back and forth between positions_summaries and positions_list:
#     ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
#     nrow(positions_summaries)==length(positions_list)
#     ####2/ always make sure that positions_summaries is ordered correctly, using the index column
#     positions_summaries <- positions_summaries[order(positions_summaries$index),]
#     ###this ensures that the first row in positions_summaries corresponds to the first trajectory in positions_list, etc.
#     for ( ant_index in 1:length(positions_list)) {
#       positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]==foraging_zone))
#       positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]!=foraging_zone))
#       
#       if ( foraging_zone %in% positions_list[[ant_index]][,"zone"]){
#         positions_summaries[ant_index,"AntTask"] <- "forager"
#       }else{
#         positions_summaries[ant_index,"AntTask"] <- "nurse"
#       }
#     }
#     #match antID and tagID (live tracking gives tagID). 
#     IDs <- exp$identificationsAt(fmTimeNow()) #this skips dead ants
#     IDs[sapply(IDs, is.null)] <- NA # assign NA to dead ants
#     IDs <- data.frame(tag_hex_ID=unlist(IDs), antID=1:length(IDs),stringsAsFactors = F)
#     positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
#     #positions_summaries1 <- positions_summaries
#     #positions_summaries_LOOP <- aggregate(cbind(positions_summaries$nb_frames_outside,positions_summaries$nb_frames_inside), by = list(positions_summaries$antID), FUN = sum);  colnames(positions_summaries_LOOP) [match("Group.1",colnames(positions_summaries_LOOP))] <- "antID"
#     positions_summaries_LOOP <- aggregate(cbind(nb_frames_outside,nb_frames_inside) ~ antID + tag_hex_ID, FUN = sum, na.rm=T,na.action=na.pass,positions_summaries )
#     #add names that help in merging later on
#     colnames(positions_summaries_LOOP) [match("nb_frames_outside",colnames(positions_summaries_LOOP))] <- paste0("nb_frames_outside",loop_N)
#     colnames(positions_summaries_LOOP) [match("nb_frames_inside",colnames(positions_summaries_LOOP))] <- paste0("nb_frames_inside",loop_N)
#     
#     
#     #positions_summaries_list <- append(positions_summaries_list, positions_summaries_LOOP)
#     positions_summaries_list[[loop_N]] <-  positions_summaries_LOOP
#     
#     
#     rm(list=ls()[which(!ls()%in%c("positions_summaries_list","exp","foraging_zone","window_shift","AntID_list","hour_chunk_start","loop_N"))]) #close experiment
#     gc()
#     mallinfo::malloc.trim(0L)
#     
#   }
#   
#   #merge all data frames together
# 
#   # # non-good looking recursive merging
#   # positions_summaries_mergA <- merge(positions_summaries_list[[1]][c("antID","nb_frames_outside1","nb_frames_inside1")], # , "tag_hex_ID"
#   #                                    positions_summaries_list[[2]][c("antID","nb_frames_outside2","nb_frames_inside2")], # , "tag_hex_ID"
#   #                                    all.x=T,all.y=T)
#   # 
#   # positions_summaries_mergB <-  merge(positions_summaries_list[[3]][c("antID","nb_frames_outside3","nb_frames_inside3")], # , "tag_hex_ID"
#   #                                     positions_summaries_list[[4]][c("antID","nb_frames_outside4","nb_frames_inside4")], # , "tag_hex_ID"
#   #                                     all.x=T,all.y=T)
#   
#   positions_summaries <- Reduce(function(x, y) merge(x, y, all=TRUE), positions_summaries_list)
#   
#   positions_summaries <- as.data.frame(sapply(positions_summaries,as.numeric))
#   
#   positions_summaries[is.na(positions_summaries)] <- 0
#   
#   #positions_summaries$prop_time_outside <- (positions_summaries$nb_frames_outside1+positions_summaries$nb_frames_outside2)/(positions_summaries$nb_frames_outside1+positions_summaries$nb_frames_outside2+positions_summaries$nb_frames_inside1+positions_summaries$nb_frames_inside2)
#   
#   #sum inside & outside
#   positions_SUMS<- data.frame(antID=positions_summaries$antID,tag_hex_ID=positions_summaries$tag_hex_ID, outside=rowSums(positions_summaries[, grep("outside", colnames(positions_summaries))]),
#              inside=rowSums(positions_summaries[, grep("inside", colnames(positions_summaries))]))
#   
#   positions_SUMS$prop_time_outside <- positions_SUMS$outside/(positions_SUMS$outside+positions_SUMS$inside)
#   
#   positions_SUMS[which(positions_SUMS$prop_time_outside<=0.01),"AntTask"] <- "nurse"
#   positions_SUMS[which(positions_SUMS$prop_time_outside>0.01),"AntTask"] <- "forager"
#   
#   AntTasks <- data.frame(antID=positions_SUMS[,"antID"],tag_hex_ID=positions_SUMS[,"tag_hex_ID"],AntTask= positions_SUMS[,"AntTask"])
#   
#   print("AntTasks computed")
# 
#   # #add missing ants as NURSE by DEFAULT
#   # missing_ants <- subset(AntID_list, !(AntID_list %in% AntTasks$antID))
#   # AntTasks <- rbind(AntTasks, data.frame(antID=missing_ants,AntTask=NA))
#   # AntTasks[which(is.na(AntTasks$AntTask)),"AntTask"] <- "nurse"
#   
#   AntTasks$AntTask_num <- NA
#   AntTasks[which(AntTasks$AntTask=="nurse"),"AntTask_num"] <- 1
#   AntTasks[which(AntTasks$AntTask=="forager"),"AntTask_num"] <- 2
#   AntTasks <- AntTasks[order(AntTasks$antID),]
#   
#   rm(list=ls()[which(!ls()%in%c("positions_summaries_list","positions_SUMS","exp","foraging_zone","window_shift","AntID_list","AntTasks"))]) #close experiment
#   gc()
#   mallinfo::malloc.trim(0L)
# 
#   #  POSITION SUMMARIES LAST TWO ELEMENTS (IN CASE OF BLOCKS OF 12H) USED FOR ZONE USAGE
#   
#   ################################## ZONE USE
#   
#   #if (ZoneUsage) {
#   print("Computing Zone (nest, foraging area) usage pre-post exposure")
#   
#   ## get 2 12Hours window for the Task calculation
#   ## calcualte the task AFTER the EXPOSURE
#   start <- fmQueryGetDataInformations(exp)$end - 24*3600 - window_shift
#   time_start <- fmTimeCreate(offset=start)
#   stop <- fmQueryGetDataInformations(exp)$end - 12*3600 - window_shift
#   time_stop   <- fmTimeCreate(offset=stop)
#   ###QUERY 3: fmQueryComputeAntTrajectories()
#   positions                 <- fmQueryComputeAntTrajectories(exp,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
#   positions_summaries       <- positions$trajectories_summary
#   positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
#   positions_list            <- positions$trajectories
#   
#   
#   ##before going back and forth between positions_summaries and positions_list:
#   ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
#   nrow(positions_summaries)
#   length(positions_list)
#   ####2/ always make sure that positions_summaries is ordered correctly, using the index column
#   positions_summaries <- positions_summaries[order(positions_summaries$index),]
#   ###this ensures that the first row in positions_summaries corresponds to the first trajectory in positions_list, etc.
#   for ( ant_index in 1:length(positions_list)) {
#     positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]==foraging_zone))
#     positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]!=foraging_zone))
#   }
# 
#   #match antID and tagID (live tracking gives tagID). 
#   IDs <- exp$identificationsAt(fmTimeNow()) #this skips dead ants
#   IDs[sapply(IDs, is.null)] <- NA # assign NA to dead ants
#   IDs <- data.frame(tag_hex_ID=unlist(IDs), antID=1:length(IDs),stringsAsFactors = F)
#   positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
#   #positions_summaries1 <- positions_summaries
#   positions_summaries1post <- aggregate(cbind(nb_frames_outside,nb_frames_inside) ~ antID + tag_hex_ID, FUN = sum, na.rm=T,na.action=na.pass,positions_summaries )
#   #add names that help in merging later on
#   names(positions_summaries1post)[names(positions_summaries1post)%in%c("nb_frames_outside","nb_frames_inside")] <- paste(names(positions_summaries1post)[names(positions_summaries1post)%in%c("nb_frames_outside","nb_frames_inside")],"1post",sep="")
#   
#   
#   rm(list=ls()[which(!ls()%in%c("positions_summaries_list","positions_SUMS","positions_summaries1post","exp","foraging_zone","window_shift","AntID_list","AntTasks"))]) #close experiment
#   gc()
#   mallinfo::malloc.trim(0L)
#   
#   ####define time period to use for defining nurses/foragers
#   start <- fmQueryGetDataInformations(exp)$end - 12*3600 - window_shift
#   time_start <- fmTimeCreate(offset=start)
#   stop <- fmQueryGetDataInformations(exp)$end  - window_shift
#   time_stop   <- fmTimeCreate(offset=stop)
#   
#   ###QUERY 3: fmQueryComputeAntTrajectories()
#   positions                 <- fmQueryComputeAntTrajectories(exp,start = time_start,end = time_stop,maximumGap = fmHour(24*365),computeZones = TRUE)
#   positions_summaries       <- positions$trajectories_summary
#   positions_summaries       <- data.frame(index=1:nrow(positions_summaries),positions_summaries,stringsAsFactors = F)
#   positions_list            <- positions$trajectories
#   
#   ##before going back and forth between positions_summaries and positions_list:
#   ####1/ always check that the number of rows of positions_summaries is equal to the length of positions_list
#   nrow(positions_summaries) ==length(positions_list)
#   ####2/ always make sure that positions_summaries is ordered correctly, using the index column
#   positions_summaries <- positions_summaries[order(positions_summaries$index),]
#   ###this ensures that the first row in psoitions_summaries corresponds to the first trajectory in positions_list, etc.
#   for ( ant_index in 1:length(positions_list)) {
#     positions_summaries[ant_index,"nb_frames_outside"]  <- length(which(positions_list[[ant_index]][,"zone"]==foraging_zone))
#     positions_summaries[ant_index,"nb_frames_inside"] <- length(which(positions_list[[ant_index]][,"zone"]!=foraging_zone))
#   }
# 
#   #match antID and tagID (live tracking gives tagID).
#   IDs <- exp$identificationsAt(fmTimeNow()) #this skips dead ants
#   IDs[sapply(IDs, is.null)] <- NA # assign NA to dead ants
#   IDs <- data.frame(tag_hex_ID=unlist(IDs), antID=1:length(IDs),stringsAsFactors = F)
#   positions_summaries$tag_hex_ID <- IDs[ match(positions_summaries$antID,IDs$antID)     , "tag_hex_ID"]
#   #positions_summaries1 <- positions_summaries
#   positions_summaries2post <- aggregate(cbind(nb_frames_outside,nb_frames_inside) ~ antID + tag_hex_ID, FUN = sum, na.rm=T,na.action=na.pass,positions_summaries )
#   #add names that help in merging later on
#   names(positions_summaries2post)[names(positions_summaries2post)%in%c("nb_frames_outside","nb_frames_inside")] <- paste(names(positions_summaries2post)[names(positions_summaries2post)%in%c("nb_frames_outside","nb_frames_inside")],"2post",sep="")
#   
#   #merge positions_summaries1 and positions_summaries2
#   #put all data frames into list
#   psummaries_list <- list(positions_summaries_list[[3]], positions_summaries_list[[4]], positions_summaries1post,positions_summaries2post)      
#   
#   #merge all data frames together
#   positions_summaries <- Reduce(function(x, y) merge(x, y, all=TRUE), psummaries_list)
#   positions_summaries[is.na(positions_summaries)] <- 0
#   
#   positions_summaries$frames_outside_24hPRE <- positions_summaries$nb_frames_outside3 + positions_summaries$nb_frames_outside4
#   positions_summaries$frames_outside_24hPOST <- positions_summaries$nb_frames_outside1post + positions_summaries$nb_frames_outside2post
#   
#   positions_summaries$frames_inside_24hPRE <- positions_summaries$nb_frames_inside3 + positions_summaries$nb_frames_inside4
#   positions_summaries$frames_inside_24hPOST <- positions_summaries$nb_frames_inside1post + positions_summaries$nb_frames_inside2post
#   
#   ZoneUse <- positions_summaries[c("antID","frames_inside_24hPRE","frames_inside_24hPOST","frames_outside_24hPRE","frames_outside_24hPOST")]
#   
#   #merge this output with the AntTasks dataframe
#   output_summ_list <- list(AntTasks,ZoneUse)
#   
#   AntTasks <- Reduce(function(x, y) merge(x, y, all=TRUE), output_summ_list)
#   
#   ############## MODIFY! KEEP SEPARATED THE PRE AND POST OUTPUTS!!!
#   
#   #}# ZoneUSe
#   
#   rm(list=ls()[which(!ls()%in%c("AntTasks"))]) #close experiment
#   gc()
#   mallinfo::malloc.trim(0L)
#   
#   ##RETURN OUTPUT
#   # warning("Ants that died before the considered time window (pre treatment) will not be assigned a Task by the function. Currently, no task will default to Nurse")
#   return(AntTasks)
#   
# } #, ZoneUsage= TRUE

################## COMPUTE NETWORK ###############################
# it requires a "gap" to be defined, should be required by the function
compute_G <- function(exp, start, end){ # min_cum_duration , link_type, nest_focus, frm_rate
  # convert timestamp of frame into corresponding frame number starting from 1 (with frame#1 at 'start' time)
  # CollideFrames <- fmQueryCollideFrames(exp,start=start,end=end)
  # TimeToFrame <- seq_along(CollideFrames$frames$time)
  Interactions <- fmQueryComputeAntInteractions(exp, start, end, maximumGap=gap, singleThreaded=FALSE, reportFullTrajectories = F)
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
#install.packages("mallinfo", repos = "http://www.rforge.net/")
library(mallinfo)
library(reshape)

#### PARAMETERS
TimeWind        <- 3600 ## in seconds (3600 is an hour)
gap             <- fmSecond(10)
## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED! 
window_shift <- 60*15 #approx N of minutes that where given at the end as leeway, minutes can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition

## some initialization
Period_dataframe <- NULL #checking time correspondances
#start fresh
Network_properties            <- data.frame()

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
WORKDIR   <- "/media/cf19810/DISK4/ADRIANO"
DATADIR   <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/")
SCRIPTDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP1_analysis scripts"
metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021_2022-10-12.txt",sep=""),header=T,stringsAsFactors = F, sep=",")


###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"AntTasks_v082.R",sep="/"))



#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]


###initialise general output folder
###remove folder if already exists to make sure we don't mix things up
if (file.exists(file.path(DATADIR, "NetworkAnalysis_outcomes"))){
  unlink(file.path(DATADIR, "NetworkAnalysis_outcomes"),recursive=T)
} 
###create folder
dir.create(file.path(DATADIR, "NetworkAnalysis_outcomes"),recursive = T)
###define name of general output table containing Network_properties
output_name <- file.path(DATADIR, "NetworkAnalysis_outcomes","Network_properties.txt") # (saved INSIDE the Network_analysis folder)
AntTasks_SpaceUse <- file.path(DATADIR,"AntTasks_SpaceUse_july2022.txt") # (saved OUTSIDE the Network_analysis folder)

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
  REP.files       <- list.files(REP.folder, pattern = "AntsCreated_AutoOriented_withMetaData_NS_NS_q")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  #replicate file
  for (REP.FILES in REP.filefolder) {
    # REP.FILES <-  REP.filefolder[2]   #temp
    print(REP.FILES) ##}}
    #open experiment
    exp <- fmExperimentOpen(REP.FILES)
    # exp.Ants <- exp$ants
    print(paste0("Processing ",basename(REP.FILES)))
    exp_end <- fmQueryGetDataInformations(exp)$end - window_shift
    
    ## get in R the spaceID / name correspondance
    ##PROBABLY NOT NEEDED. ANYWAY, IT WORKS DIFFERNTLY THAN IN FLORA-bee
    BoxCodes <- exp$spaces[[1]]$zones
    BoxCodes <- data.frame(space=c(BoxCodes[[1]]$name, BoxCodes[[2]]$name), box=c(exp$spaces[[1]]$name), stringsAsFactors = F )
    
    # ########### DEFINE ZONES PROPERLY - DEFINED INSIDE THE ANT TASKS FUNCTION, SHOULD DO SAME FOR NEST USE
    # #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
    # zones <- exp$spaces[[1]]$zones #function to show the Zones present in the Space
    # zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
    # foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
    # nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
    # print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))

    
    
    ########## GET EXPOSED ANTS from METADATA
    exp.Ants <- exp$ants
    metadata$Exposed <- "no"
    for (ant in exp.Ants){
      individual  <- ant$ID
      #print(ant)
        if (TRUE %in% ant$getValues("Exposed")[,"values"]) { metadata[individual,"Exposed"] <- "exposed" }
    }
    
    ########## GET QUEENS 
    metadata$IsQueen <- "no"
    for (ant in exp.Ants){
      individual  <- ant$ID
      #print(ant)
      if (TRUE %in% ant$getValues("IsQueen")[,"values"]) { metadata[individual,"IsQueen"] <- "queen" }
    }

    #base file info
    #TREATMENT
    COLONY <- sub("\\_.*", "", basename(REP.FILES))
    TREATMENT <- substr(COLONY,(nchar(COLONY)+1)-2,nchar(COLONY))
    #
    COLONY_SIZE <- unique(metadata[which(metadata$REP_treat == COLONY),"colony_size"])

    # ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
    print(paste0("Compute 3-hours analysis"))
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
  print(paste0("Computing networks and properties"))
  Graph <- compute_G(exp = exp, start = start, end=end)

  # COMPUTE NETWORK PROPERTIES
  Network_prop_hour <- NetProperties(graph=Graph)
  
  Network_properties <- rbind(Network_properties,Network_prop_hour)
  
  
  #####################################################################################################
  
  ######TO FIX!!!!!
  
  # SpaceUsage should produce as output a 3h blocks list of the ants, run by run of the  HOUR loop. 
  # similarly to Network_prop_hour (but with many rows) there should be a stacking, keeping the "hour" label and then saving the output outside per each replicate!
  # RE RUN SPACE USAGE TO MAKE SURE IT PRODUCES AS OUTPUT A ROW PER ANT PER HOURBLOCK!
  
  #keep relevant exp info for SpaceUsage
  SpaceUsage <- data.frame(randy=REP.FILES,colony=COLONY,colony_size=COLONY_SIZE,treatment=TREATMENT,period=PERIOD,time_hours=HOUR, From, To, SpaceUsage)
  
  #####################################################################################################
  
  } }#REP LOOP
    

    ########################################
    ### prop of time Exposed ants spend inside the nest, pre-post exposure
    # AntTasks$prop_inside_24hPRE <- AntTasks$frames_inside_24hPRE /(AntTasks$frames_inside_24hPRE + AntTasks$frames_outside_24hPRE)
    # AntTasks$prop_inside_24hPOST <- AntTasks$frames_inside_24hPOST /(AntTasks$frames_inside_24hPOST + AntTasks$frames_outside_24hPOST)
    # AntTasks$delta_time_inside <- AntTasks$prop_inside_24hPOST - AntTasks$prop_inside_24hPRE
    # 
    ########################################
    ##### SAVE FILES IN FOLDER #############

    ## Network properties save (saved INSIDE the Network_analysis folder)
    if (file.exists(output_name)){
      write.table(Network_properties,file=output_name,append=T,col.names=F,row.names=F,quote=T,sep=",")
    }else{
      write.table(Network_properties,file=output_name,append=F,col.names=T,row.names=F,quote=T,sep=",")
    }
    
    
    ###########################################################################################
    ###########################################################################################
    #### IMPORTANT FOR SPACEUSAGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # ## AntTasks save (saved OUTSIDE the Network_analysis folder)
    # if (file.exists(AntTasks_SpaceUse)){
    #   write.table(AntTasks,file=AntTasks_SpaceUse,append=T,col.names=F,row.names=F,quote=T,sep=",")
    # }else{
    #   write.table(AntTasks,file=AntTasks_SpaceUse,append=F,col.names=T,row.names=F,quote=T,sep=",")
    # }
    #
    
    
    #start fresh
    Network_properties            <- data.frame()
    
    #cleaning
    rm(   list   =  ls()[which(!ls()%in%to_keep)]    )
    gc()
    mallinfo::malloc.trim(0L)
    
    
    ######################################################
    ### PLOTTING (1 col..., should be all)
    #plot(AntTasks$delta_time_inside, col=as.factor(AntTasks$Exposed))
    
    ####################################################################################
    ### plot MEAN DELTA PER ANT SEPARATED BETWEEN NON EXPOSED AND EXPOSED, WITH COL. SIZE COMPARISON
    ###################################################################################
    

     
  }
}


loop_end_time <- Sys.time()
print (paste("loop took ",as.numeric(difftime(loop_end_time, loop_start_time, units = "mins"))," minutes to complete"))


#######
## SECTION OF MODS TO BE APPLIED!!

##############################################################################################################
##############################################################################################################
###### IMPORTANT #############################################################################################
## FOR "period" , "time_hours" AND "time_of_day"TO BE ASSIGNED WE NEED TO CALCULATE THINGS ON 3H CHUNKS (NOT 12!)
###
## THESE CALCULATIONS HAVE TO BE MADE INSIDE THE SCRIPT, BEFORE THE FILE IS SAVED AT EVERY LOOP ITERATION ####
##############################################################################################################
##############################################################################################################
##############################################################################################################


#reference
individual_behavioural_data <- read.table("/home/cf19810/Documents/TEMP/individual_behavioural_data.txt",header=T,stringsAsFactors = F)

Time_dictionary <- unique(individual_behavioural_data[c("time_hours","time_of_day","period")])

#treatment
AntTasks$treatment_new <- NA
for (ROW in 1:nrow(AntTasks)) {
  if(AntTasks[ROW,"treatment"]  %in% c("BP","SP")){  AntTasks[ROW,"treatment_new"] <- "pathogen"
  }else if ( AntTasks[ROW,"treatment"] %in% c("BS","SS")) { AntTasks[ROW,"treatment_new"] <- "control" } else{ print("ERROR")}
  
}

#status
AntTasks$status <- NA
for (ROW in 1:nrow(AntTasks)) {
  if(AntTasks[ROW,"treatment"]  %in% c("BP","BS")){  AntTasks[ROW,"status"] <- "large"
  }else if ( AntTasks[ROW,"treatment"] %in% c("SP","SS")) { AntTasks[ROW,"status"] <- "small" } else{ print("ERROR")}
  
}

#age
AntTasks$age <- NA
for (ROW in 1:nrow(AntTasks)) {
  if(AntTasks[ROW,"treatment"]  %in% c("BP","BS")){  AntTasks[ROW,"status"] <- "old"
  }else if ( AntTasks[ROW,"treatment"] %in% c("SP","SS")) { AntTasks[ROW,"status"] <- "young" } else{ print("ERROR")}
  
}


# period ("after","before")
# SHOULD BE ASSIGNED IN THE ANT_TASK FUNCTION

# time_hours
# FIND A WAY TO ASSIGN TIMING TO ANT_TASK AND UPDATE NETWORK_ANALYSIS TIMING TO MATCH THE Time_dictionary

# time_of_day
# SIMPLY USE THE Time_dictionary AS GUIDE -> EVERY time_hours HAS A SPECIFIC time_of_day BUT THE OPPOSITE IS NOT TRUE (MULTIPLE DAYS PRESENT, REPETED time_of_day VALUES)



# COLUMNS: 
# colony    colony_size x
# treatment (two values: pathogen and control)  x
# tag    (only include treated workers) x
# age    x
# status    ( replace the content of the "status" column with either "large" or "small" (rather than "treated" or "untreated") ) x
# period    ( pre/post chunks corresponding to the same time of day should have the same value in column "time_of_day") 
# time_hours    
# time_of_day
# 
# EXTRA FROM ME:
#   prop_time_outside
# Remaing cols

# OPTIONAL:
# proportion_time_active    
# average_bout_speed_pixpersec    
# total_distance_travelled_pix



## time_dictionary
# time_hours time_of_day period
# 1              0          12  after
# 69             3          15  after
# 137            6          18  after
# 205            9          21  after
# 273           12           0  after
# 341           15           3  after
# 409           18           6  after
# 477           21           9  after
# 545          -30           6 before
# 613          -27           9 before
# 681          -24          12 before
# 749          -21          15 before
# 817          -18          18 before
# 885          -15          21 before
# 953          -12           0 before
# 1021          -9           3 before
# 11385         -6           6 before
