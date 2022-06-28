rm(list=ls())
library(Rcpp); library(FortMyrmidon); library(igraph); library(pals); library(lubridate)

## Parameters
TimeWind        <- 3600 ## in seconds
StartTime       <- "2020-05-15 01:00:12"  ## manually define the start time - surely there's a better way of doing this (get the first frame in the data..?)
StartTime_POSIX <- as.POSIXct(StartTime, tz = "CET") ## IS THIE TIME-ZONE CORECT..?
gap             <- fmSecond(2)

## FUNCTIONS
DISTANCE     <- function(x)  { c(sqrt((x[-nrow(x), "x"] - x[-1, "x"])^2 + (x[-nrow(x), "y"] - x[-1, "y"])^2), NA)}

## Load the myrmidon files & define the time window
e     <- fmExperimentOpenReadOnly("/media/tom/HDD/HoneyBee/HoneyBee_Flora-fixed.myrmidon")  

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

### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts hour-by-hour, stacking vertically
RawContacts_ALL <- NULL; AntMinuteMeans_ALL <- NULL

for (HOUR in 1:(5*24))  ## assumes the tracking lasted for 5 full days
  {
  ## increment the window start & end times by 1 hour
  From <- StartTime_POSIX + (HOUR * TimeWind)  
  To   <- From + TimeWind
  print(paste("Time window, from", From, "to", To))
  ## fmTimeParse wants the time in the format, "2020-05-15T01:00:12.000Z", so make it so ...
  From_Formatted <- paste0(gsub(" ","T",as.character(From)),".000Z",sep="")
  To_Formatted   <- paste0(gsub(" ","T",as.character(To)),  ".000Z",sep="")
  ## Define the time window here. 
  start <- fmTimeParse(From_Formatted)  
  end   <- fmTimeParse(To_Formatted)    
  
  ## extract 1 hour of contacts
  ## WARNING - THIS WILL CRASH IF YOU READ IN "TOO MUCH 24 HOURS OF DATA AT ONCE!!!
  dd   <- fmQueryComputeAntInteractions(e, start=start, end=end, maximumGap=gap, showProgress = TRUE, singleThreaded=FALSE, reportGlobalTrajectories=TRUE)#, reportLocalTajectories=FALSE) #reportTrajectories=TRUE
  
  ################### PART 1 : EXTRACTING TRAJECTORY SEGMENTS ##################################
  
  ## for each ant, vertically stack the trajectory segments from the current HOUR
  AntMinuteMeans_HOUR <- NULL
  for (ID in unique(dd$summaryTrajectory$antID))
    {print(paste("stacking trajectory segments for ID", ID))
    AntIndices <- which(dd$summaryTrajectory$antID == ID)     ## find which indices correspond to this ant
    Traj_Segments <- dd[["trajectories"]]       [AntIndices]  ## each traj segment for ID
    Traj_Starts   <- dd$summaryTrajectory$start [AntIndices]  ## the UNIX start time for each traj segment
    ## for each trajectory segment, 
    for (SG in 1:length (AntIndices))
      {
      Traj_Segments[[SG]]$dt   <- c(diff(Traj_Segments[[SG]]$time),NA) ## the max value will never be any more than 'gap'
      Traj_Segments[[SG]]$time <- ceiling_date( Traj_Segments[[SG]]$time + Traj_Starts[SG] , "min")  ## bin using 'round', so 11 min, 10 sec -> 11 min
      ## calculate the distance between successive fixes
      Traj_Segments[[SG]]$distance <- DISTANCE(x = Traj_Segments[[SG]])
      Traj_Segments[[SG]]$speed    <- Traj_Segments[[SG]]$distance / Traj_Segments[[SG]]$dt
      }  
    ## collapse segments -> one data frame
    Traj_Segments <-  do.call("rbind", Traj_Segments)
    ## get mean & sd of speed within each rounded minute
    AntMinuteMeans          <- aggregate(speed ~ time, FUN=length,        Traj_Segments, na.action=na.pass); names(AntMinuteMeans)[match("speed",names(AntMinuteMeans))] <- "N_observations"
    AntMinuteMeans$speed    <- aggregate(speed ~ time, FUN=mean, na.rm=T, Traj_Segments, na.action=na.pass)$speed
    AntMinuteMeans$speed_sd <- aggregate(speed ~ time, FUN=sd,   na.rm=T, Traj_Segments, na.action=na.pass)$speed
    ## stack the minutely-means
    AntMinuteMeans_HOUR <- rbind(AntMinuteMeans_HOUR, data.frame(ID, AntMinuteMeans))
    }
  ## stack the hourly segments (faster than just adding each segment to a single object)
  AntMinuteMeans_ALL <- rbind(AntMinuteMeans_ALL, AntMinuteMeans_HOUR)
  
  ## occasional time-series plotting
  if (runif(1)<0.1)
    {
    par(mfrow=c(2,1), mai=c(0.45,0.45,0.1,0.1), mgp=c(1.3,0.3,0), tcl=-0.2)
    ColonyACtivity <- aggregate(speed ~ round_date(time,"10 minutes"), FUN=mean, na.rm=T, AntMinuteMeans_ALL)
    ## Mean activity
    plot(ColonyACtivity, type="h", xaxs="i", xlab="", yaxt="n", ylab="")
    ## hacky colony-level heatmap
    plot(ColonyACtivity, type="n", xaxs="i", xlab="", yaxt="n", ylab="")
    MeanSpeedPerAnt <- aggregate(speed ~ round_date(time,"10 minutes") + ID, FUN=mean, na.rm=T, AntMinuteMeans_ALL); colnames(MeanSpeedPerAnt)[1] <- "time"; MeanSpeedPerAnt <- xtabs(speed ~ time + ID, MeanSpeedPerAnt)
    par(new=T)
    image(x=1:nrow(MeanSpeedPerAnt), y=1:ncol(MeanSpeedPerAnt), z=sqrt(MeanSpeedPerAnt), col=c("white",pals::tol.rainbow()), xaxt="n", yaxt="n", xlab="Time", ylab="Bee"); box()
    }
  
  
  ################### PART 2 : EXTRACTING ANT-TO-ANT INTERACTIONS ##################################
  
  ## vertically stack the hourly aggregated contacts
  RawContacts <- dd[["interactions"]]
  
  ## assign box codes to the data frame
  RawContacts$box <- NA
  RawContacts$box [ RawContacts$space==BoxCodes$space[1]] <- BoxCodes$box[1]
  RawContacts$box [ RawContacts$space==BoxCodes$space[2]] <- BoxCodes$box[2]
  
  RawContacts_ALL <- rbind(RawContacts_ALL, data.frame(HOUR, From, To, RawContacts))
  print(paste("stacked Raw Contacts has",nrow(RawContacts_ALL),"rows"))
  }## HOUR



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



