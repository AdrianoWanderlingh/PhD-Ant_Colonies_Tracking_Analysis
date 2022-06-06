rm(list=ls())

#"https://formicidae-tracker.github.io/myrmidon/latest/index.html"

##### LIBRARIES
library(FortMyrmidon) ####R bindings
#library(Rcpp)
library(data.table)
library(lubridate)
library(pals)

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
DISTANCE     <- function(x)  { c(sqrt((x[-nrow(x), "x"] - x[-1, "x"])^2 + (x[-nrow(x), "y"] - x[-1, "y"])^2), NA)}


#### ACCESS FILES
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA", n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]

#### OPEN REPLICATE
# replicate folder
for (REP.n in 1:length(files_list)) {
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = ".myrmidon")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  #replicate file
  for (REP.FILES in REP.filefolder) {
    print(REP.FILES)
    #open experiment
    e <- fmExperimentOpen(REP.FILES)
    
    ## use the e$cSpaces() method to get in R the spaceID / name correspondance
    BoxCodes <- paste0(capture.output(e$cSpaces()),  sep=" ", collapse=" ")
    BoxCodes <- gsub("fmCSpace","",BoxCodes); BoxCodes <- gsub("ByID","",BoxCodes)
    BoxCodes <- strsplit(BoxCodes, "->")[[1]]
    BoxCodes <- unlist(strsplit(BoxCodes," ID = "))
    BoxCodes <- unlist(BoxCodes [ grep("Name ",BoxCodes)])
    BoxCodes <- unlist(strsplit(BoxCodes, "\\)"))
    BoxCodes <- unlist(BoxCodes [ grep("Name ",BoxCodes)])
    BoxCodes <- gsub("Name =","",BoxCodes)
    BoxCodes <- strsplit(BoxCodes," ")
    BoxCodes <- data.frame(space=c(BoxCodes[[1]][1], BoxCodes[[2]][1]), box=c(BoxCodes[[1]][3], BoxCodes[[2]][3]), stringsAsFactors = F )

    ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
    RawContacts_ALL <- NULL; AntMinuteMeans_ALL <- NULL
    
    for (HOUR in seq(from=0, to=48, by=3))  ## increments by 3 hours for 48 hours 
    {
      ## increment the window start & end times by 1 hour
      From  <- fmQueryGetDataInformations(e)$end - 48*3600 + (HOUR * TimeWind)
      To    <- fmQueryGetDataInformations(e)$end - 45*3600 + (HOUR * TimeWind)
      print(paste("Time window, from", From, "to", To))
      start <- fmTimeCreate(offset=From) #end minus 48 hours plus incremental time
      end   <- fmTimeCreate(offset=To) #end minus 45 hours plus incremental time
      
      ## extract 3 hours of contacts
      ## WARNING - THIS WILL CRASH IF YOU READ IN TOO MANY HOURS AT ONCE
      dd   <- fmQueryComputeAntInteractions(e, start=start, end=end, maximumGap=gap, showProgress = TRUE, singleThreaded=FALSE, reportFullTrajectories = T)#, reportLocalTajectories=FALSE) #reportTrajectories=TRUE
      
      # ################### PART 1 : EXTRACTING TRAJECTORY SEGMENTS ##################################
      # 
      # 
      # ## for each ant, vertically stack the trajectory segments from the current HOUR
      # AntMinuteMeans_HOUR <- NULL
      # for (ID in unique(dd$trajectories_summary$antID))
      # {print(paste("stacking trajectory segments for ID", ID))
      #   AntIndices <- which(dd$trajectories_summary$antID == ID)     ## find which indices correspond to this ant
      #   Traj_Segments <- dd[["trajectories"]]       [AntIndices]  ## each traj segment for ID
      #   Traj_Starts   <- dd$trajectories_summary$start [AntIndices]  ## the UNIX start time for each traj segment
      #   ## for each trajectory segment,
      #   for (SG in 1:length (AntIndices))
      #   {
      #     Traj_Segments[[SG]]$dt   <- c(diff(Traj_Segments[[SG]]$time),NA) ## the max value will never be any more than 'gap'
      #     Traj_Segments[[SG]]$time <- ceiling_date( lubridate::seconds(Traj_Segments[[SG]]$time) + Traj_Starts[SG] , "min")  ## bin using 'round', so 11 min, 10 sec -> 11 min
      #     ## calculate the distance between successive fixes
      #     Traj_Segments[[SG]]$distance <- DISTANCE(x = Traj_Segments[[SG]])
      #     Traj_Segments[[SG]]$speed    <- Traj_Segments[[SG]]$distance / Traj_Segments[[SG]]$dt
      #   }
      #   ## collapse segments -> one data frame
      #   Traj_Segments <-  do.call("rbind", Traj_Segments)
      #   ## get mean & sd of speed within each rounded minute
      #   AntMinuteMeans          <- aggregate(speed ~ time, FUN=length,        Traj_Segments, na.action=na.pass); names(AntMinuteMeans)[match("speed",names(AntMinuteMeans))] <- "N_observations"
      #   AntMinuteMeans$speed    <- aggregate(speed ~ time, FUN=mean, na.rm=T, Traj_Segments, na.action=na.pass)$speed
      #   AntMinuteMeans$speed_sd <- aggregate(speed ~ time, FUN=sd,   na.rm=T, Traj_Segments, na.action=na.pass)$speed
      #   ## stack the minutely-means
      #   AntMinuteMeans_HOUR <- rbind(AntMinuteMeans_HOUR, data.frame(ID, AntMinuteMeans))
      # }
      # ## stack the hourly segments (faster than just adding each segment to a single object)
      # AntMinuteMeans_ALL <- rbind(AntMinuteMeans_ALL, AntMinuteMeans_HOUR)
      #
      # ## occasional time-series plotting
      # if (runif(1)<0.1)
      # {
      #   par(mfrow=c(2,1), mai=c(0.45,0.45,0.1,0.1), mgp=c(1.3,0.3,0), tcl=-0.2)
      #   ColonyACtivity <- aggregate(speed ~ round_date(time,"10 minutes"), FUN=mean, na.rm=T, AntMinuteMeans_ALL)
      #   ## Mean activity
      #   plot(ColonyACtivity, type="h", xaxs="i", xlab="", yaxt="n", ylab="")
      #   ## hacky colony-level heatmap
      #   plot(ColonyACtivity, type="n", xaxs="i", xlab="", yaxt="n", ylab="")
      #   MeanSpeedPerAnt <- aggregate(speed ~ round_date(time,"10 minutes") + ID, FUN=mean, na.rm=T, AntMinuteMeans_ALL); colnames(MeanSpeedPerAnt)[1] <- "time"; MeanSpeedPerAnt <- xtabs(speed ~ time + ID, MeanSpeedPerAnt)
      #   par(new=T)
      #   image(x=1:nrow(MeanSpeedPerAnt), y=1:ncol(MeanSpeedPerAnt), z=sqrt(MeanSpeedPerAnt), col=c("white",pals::tol.rainbow()), xaxt="n", yaxt="n", xlab="Time", ylab="Bee"); box()
      # }


      ################### PART 2 : EXTRACTING ANT-TO-ANT INTERACTIONS ##################################

      ## vertically stack the hourly aggregated contacts
      RawContacts <- dd$interactions

      # ## assign box codes to the data frame
      # RawContacts$box <- NA
      # RawContacts$box [ RawContacts$space==BoxCodes$space[1]] <- BoxCodes$box[1]
      # RawContacts$box [ RawContacts$space==BoxCodes$space[2]] <- BoxCodes$box[2]

      RawContacts_ALL <- rbind(RawContacts_ALL, data.frame(HOUR, From, To, RawContacts))
      print(paste("stacked Raw Contacts has",nrow(RawContacts_ALL),"rows"))
    }## HOUR

    #####################################################################################
    #####################################################################################
    #####################################################################################
    # CONTINUE FROM HERE
    #####################################################################################
    #####################################################################################
    #####################################################################################
    
    
    ## calculate proportion of time in the foraging box
    RawContacts_ALL$space_binary <- RawContacts_ALL$space - 1 ## only works if there are exactly 2 cameras
    Contact_Time_Allocation      <- data.frame(ant=c(RawContacts_ALL$ant1,RawContacts_ALL$ant2), space_binary=c(RawContacts_ALL$space_binary,RawContacts_ALL$space_binary))
    Box_Time_Allocation          <- aggregate(space_binary ~ ant, FUN=mean, Contact_Time_Allocation) ## mean of c(0,0,1,1,0,0)=1/3
    
    ## count the contacts between each unique pair of ants
    AggContacts       <- aggregate(end ~ ant1 + ant2, FUN=length, RawContacts_ALL); colnames(AggContacts)[match("end",colnames(AggContacts))] <- "weight"
    ## count the contacts by the time of day & the box
    AggContactsByHour <- aggregate(end ~ HOUR + box,  FUN=length, RawContacts_ALL); colnames(AggContactsByHour)[match("end",colnames(AggContactsByHour))] <- "N_contacts"
    
    
    ## calculate proportion of time in the foraging box
    RawContacts_ALL$space_binary <- RawContacts_ALL$space - 1 ## only works if there are exactly 2 cameras
    Contact_Time_Allocation      <- data.frame(ant=c(RawContacts_ALL$ant1,RawContacts_ALL$ant2), space_binary=c(RawContacts_ALL$space_binary,RawContacts_ALL$space_binary))
    
    ## count the contacts between each unique pair of ants
    AggContacts       <- aggregate(end ~ ant1 + ant2, FUN=length, RawContacts_ALL); colnames(AggContacts)[match("end",colnames(AggContacts))] <- "weight"
    ## count the contacts by the time of day & the box
    AggContactsByHour <- aggregate(end ~ HOUR + box,  FUN=length, RawContacts_ALL); colnames(AggContactsByHour)[match("end",colnames(AggContactsByHour))] <- "N_contacts"
    
    ## create an overall time-aggregated contact network from the counts
    G <- graph_from_data_frame(d=AggContacts, directed=F); is.weighted(G)
    
    ## calculate the proportion of time spent in the foraging box
    V(G)$space_binary <- Box_Time_Allocation $ space_binary [match(V(G)$name, Box_Time_Allocation$ant)]
    
    ## node centrality - use to size the nodes
    V(G)$degree          <- degree(G, mode="all")
    V(G)$weighted_degree <- graph.strength(G, mode="all")
    
    ## threshold edges - just for plotting
    G2 <- delete.edges(graph=G, edges=E(G) [E(G)$weight < quantile(E(G)$weight,0.5)] )  ## remove 50% of the weakest edges, for plotting only
    ## spring-embedded graph layout
    Layout <- layout_with_fr(graph=G2, weights= log10(E(G2)$weight) ) ##
    
    
    
    ## PLOTTING
    par(mfrow=c(2,2), mai=c(0.4,0.4,0.1,0.1), tcl=-0.2, mgp=c(1.3,0.3,0))
    ## time-series
    plot  (N_contacts ~ HOUR, AggContactsByHour[AggContactsByHour$box=="Foraging",], type="b", pch=21, bg=1, ylab="N contacts / hour", ylim=range(AggContactsByHour$N_contacts))
    points(N_contacts ~ HOUR, AggContactsByHour[AggContactsByHour$box=="Nest",],     type="b", pch=21, bg=2)
    ## distribution of proportion contacts in the foraging box
    hist(Box_Time_Allocation$space_binary, breaks=25, main="", xlab="Proportion of time in foraging arena", ylab="N bees", col=1)
    ## degree ~ prop time in forage arena
    plot(V(G)$weighted_degree ~ V(G)$space_binary, xlab="Proportion of time in foraging arena", ylab="N contacts")
    ## plot the contact network
    plot(G, Layout, 
         vertex.size= 5 + (8 * V(G)$weighted_degree/max(V(G)$weighted_degree)), 
         vertex.color=parula(11)[10 * round(V(G)$space_binary,1)],
         vertex.label=NA, 
         edge.color=rgb(0.5,0.5,0.5,0.75,maxColorValue =1),
         edge.width=2*E(G)$weight/mean(E(G)$weight),
         edge.curved=0.2)
    ## mean of c(0,0,1,1,0,0)=1/3
    
    

    
  }
}
