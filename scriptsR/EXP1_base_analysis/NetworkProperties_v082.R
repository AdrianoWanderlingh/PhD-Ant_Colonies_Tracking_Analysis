##################################################################
################## COMPUTE NETWORK ###############################
##################################################################

# it requires a "gap" to be defined, should be required by the function
#Gap is 10 seconds as these are interactions (he maximum gap in tracking before cutting the trajectory or interactions in two different object.)
compute_G <- function(exp, start, end, gap){ # min_cum_duration , link_type, nest_focus, frm_rate
  print(paste0("Computing networks"))
  
  
  #required packages
  require(igraph)
  require(FortMyrmidon)
  
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
  
  ################# assign vertex attributes
  #subset metadata
  METADATA <- subset(metadata,metadata$REP_treat==REP_TREAT)
  #create matching of vertices with ants (there could be less ants in a 3-hours window than in the full exp)

  ############# Assign vertex types: "AntTask"
  if ("AntTask" %in% colnames(metadata)) {
    METADATA$AntTask_num <- NA
    METADATA[which(METADATA$AntTask=="nurse"),"AntTask_num"] <- 1
    METADATA[which(METADATA$AntTask=="forager"),"AntTask_num"] <- 2
    
    #keep only ants corresponding to vertices
    V_METADATA <- METADATA[which(METADATA$antID %in% V(G)$name),]
 
    G <- set_vertex_attr(G, name="AntTask", index = V(G), value = V_METADATA$AntTask_num)
    #REMOVE nodes without ANTTASK (missing ants from observation as likely died at early stage of experiment - not appearing in the 48h pre-exposure)
    G <- delete_vertices(G, is.na(V(G)$AntTask))
  }else{ 
  warning("No metadata provided, impossible to assign Vertex types to Graph")}
  
  #############  Assign vertex types: "Exposed"
  G <- set_vertex_attr(G, name="Exposed", index = V(G), value = V_METADATA$Exposed)
  
  #############  Assign vertex types: "IsQueen"
  G <- set_vertex_attr(G, name="IsQueen", index = V(G), value = V_METADATA$IsQueen)
  
  #############  Assign vertex types: "antID"
  G <- set_vertex_attr(G, name="antID", index = V(G), value = V_METADATA$antID)
  
  return(G)
} # compute_G 

##################################################################
################## COMPUTE NETWORK PROPERTIES ####################
##################################################################

NetProperties <- function(graph){
  print(paste0("Computing network properties"))
  
  #required packages
  require(igraph)
  require(FortMyrmidon)
  
  
  summary_collective <- NULL
  summary_individual <- NULL
  
  
  ##### COLLECTIVE NETWORK PROPERTIES ######################################
  #### inherited from Stroeymeyt et al. 2018
  
  ## Assortativity - Task
  #degree of preferential association between workers of the same task group, calculated using Newman’s method
  #if (!is.n(V(graph)$AntTask)) {
  task_assortativity  <- assortativity_nominal(graph, types= V(graph)$AntTask, directed=F)#}else{
    #warning("No  Task Vertex information, impossible to calculate task_assortativity")}
  
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
  summary_collective <- rbind(summary_collective,data.frame(
                                                            task_assortativity=task_assortativity,
                                                            clustering=clustering,
                                                            degree_mean=degree_mean,
                                                            degree_maximum=degree_maximum,
                                                            density=density,
                                                            diameter=diameter,
                                                            efficiency=efficiency,
                                                            modularity=modularity,stringsAsFactors = F))
  
  ##### INDIVIDUAL NETWORK PROPERTIES ######################################
  #### inherited from Stroeymeyt et al. 2018
  ####Part 2: individual network properties ####

  ####Part 2: individual network properties ####
  ###prepare table
  individual <- data.frame(antID=V(graph)$antID,
                           IsQueen=V(graph)$IsQueen,
                           Exposed=V(graph)$Exposed,
                           AntTask=V(graph)$AntTask,
                           degree=NA,
                           aggregated_distance_to_queen=NA#,
                           #mean_aggregated_distance_to_treated=NA,
                           #same_community_as_queen=NA
                           )
  ##degree
  individual[match(names(degrees),individual$antID),"degree"] <- degrees
  # ##same community as queen
  # queen_comm <- community_membership[which(V(net)$name==queenid)]
  # community_membership <- community_membership==queen_comm
  # individual[match(V(net)$name,individual$tag),"same_community_as_queen"] <- community_membership
  ##path length to queen
    path_length_to_queen <- t(shortest.paths(graph,v=V(graph),to=V(graph)[which(V(graph)$IsQueen==TRUE)] ,weights=1/E(graph)$weight))
    individual[match(colnames(path_length_to_queen),individual$antID),"aggregated_distance_to_queen"] <- as.numeric(path_length_to_queen )
  ########Mean path length to treated; aggregated_network
    path_length_to_treated                             <- as.data.frame(as.matrix(shortest.paths(graph,v=V(graph),to=V(graph)[which(V(graph)$Exposed==TRUE)] ,weights=1/E(graph)$weight)))
    path_length_to_treated["mean_distance_to_treated"] <- NA
    path_length_to_treated$mean_distance_to_treated    <- as.numeric(rowMeans(path_length_to_treated,na.rm=T))
    individual[match(rownames(path_length_to_treated),individual$antID),"mean_aggregated_distance_to_treated"] <- path_length_to_treated[,"mean_distance_to_treated"]
  ###Add data to main data table
  summary_individual <- rbind(summary_individual,individual)
  
  ####################
  #return a list object.
  
  summaries_network <- list(summary_collective=summary_collective, summary_individual=summary_individual)
  return(summaries_network)
}