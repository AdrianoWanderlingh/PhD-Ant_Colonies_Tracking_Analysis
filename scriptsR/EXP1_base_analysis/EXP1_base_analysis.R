rm(list=ls())

# NOTE ON DIFFERNCES FROM SCIENCE 2018 PAPER:
# period is ("pre","post"), not ("after","before")
# treatment is ("Pathogen","Sham"), not ("pathogen","sham")
# status is linked to size of the colony, not ("treated","untreated") ant

# missing:
# - prop_time_outside (not calculated but data is there from SpaceUse)
# - proportion_time_active 
# - average_bout_speed_pixpersec
# total_distance_travelled_pix

# ------------------------------------
#DONE:
# COLUMNS: 
# colony    colony_size x
# treatment (two values: pathogen and control)  x
# tag    (only include treated workers) x
# age    x
# status    ( replace the content of the "status" column with either "large" or "small" (rather than "treated" or "untreated") ) x
# period    ( pre/post chunks corresponding to the same time of day should have the same value in column "time_of_day") x
# time_hours    x
# time_of_day x
# 


#reference
individual_behavioural_data <- read.table("/home/cf19810/Documents/TEMP/individual_behavioural_data.txt",header=T,stringsAsFactors = F)
Time_dictionary_N <- unique(individual_behavioural_data[c("time_hours","time_of_day","period")])

########################################################################################################
################### IMPORTANT THINGS TO MODIFY #########################################################

# UNIFORM Plot_Grooming_Pre-Post.R to the output of this script and the sourced functions
# AGGREGATE SPACE USE FUNCTION FOR tag_hex_ID (we want 1 row per ant!)

# ########### DEFINE ZONES PROPERLY - DEFINED INSIDE THE ANT TASKS FUNCTION, SHOULD DO SAME FOR NEST USE
# #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
# zones <- exp$spaces[[1]]$zones #function to show the Zones present in the Space
# zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
# foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
# nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
# print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))


#little improvement to do
print("RENAME THE EXP TO \"e\", NOT \"exp\" (FUNCTIONS' VARIABLE NAME) ")


##############################################################################
######
#######
##########
#############
##########
#######
######
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
library(dplyr)

#### PARAMETERS
TimeWind        <- 3600 ## in seconds (3600 is an hour)
gap             <- fmSecond(10) # for ComputeGraph function
Time_dictionary <- data.frame(time_hours= -36:35, time_of_day= rep(0:23,3)) #time_hours= time since exposure (negatives pre-exposire, positive post exposure), time_of_day=time hour slot of the day
Time_dictionary <- Time_dictionary[which(Time_dictionary$time_hours<=21 & Time_dictionary$time_hours>=-27),]

## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED! 
window_shift <- 60*15 #approx N of minutes that where given at the end as leeway, minutes can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition

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
metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021_2022-10-20.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"SpaceUse_v082.R",sep="/")) # SPACE_USE
source(paste(SCRIPTDIR,"NetworkProperties_v082.R",sep="/")) # NET_properties
#list.functions.in.file(paste(SCRIPTDIR,"NetworkProperties_v082.R",sep="/"), alphabetic = TRUE)


#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]

###initialise general output folder
###remove folder if already exists to make sure we don't mix things up
# if (file.exists(file.path(DATADIR, "NetworkAnalysis_outcomes"))){
#   unlink(file.path(DATADIR, "NetworkAnalysis_outcomes"),recursive=T)
# } 
###create folder
dir.create(file.path(DATADIR, paste0("Exp1_Results_", Sys.Date())),recursive = T)
###define name of general output files
NET_properties <- file.path(DATADIR, paste0("Exp1_Results_", Sys.Date()),"Network_properties.txt") # (saved INSIDE the folder)
SPACE_USE <- file.path(DATADIR, paste0("Exp1_Results_", Sys.Date()),"Space_Usage.txt") # (saved INSIDE the folder)


##### RUNNING TIME
loop_start_time <- Sys.time()

####define to_keep variables to keep clearing memory between runs
to_keep <- c(ls(),c("to_keep"))

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
    
    ## some initialization
    Period_dataframe <- NULL #checking time correspondances
    #start fresh
    Network_properties            <- data.frame()
    SpaceUsage                    <- data.frame()
    
    
    print(REP.FILES) ##}}
    #open experiment
    exp <- fmExperimentOpen(REP.FILES)
    # exp.Ants <- exp$ants
    print(paste0("Processing ",basename(REP.FILES)))
    exp_end <- fmQueryGetDataInformations(exp)$end - window_shift
    
    ## get in R the spaceID / name correspondance
    ##PROBABLY NOT NEEDED. ANYWAY, IT WORKS DIFFERNTLY THAN IN FLORA-bee
    # BoxCodes <- exp$spaces[[1]]$zones
    # BoxCodes <- data.frame(space=c(BoxCodes[[1]]$name, BoxCodes[[2]]$name), box=c(exp$spaces[[1]]$name), stringsAsFactors = F )
    
    
    #base file info
    REP_TREAT <- sub("\\_.*", "", basename(REP.FILES))
    #SIZE_TREAT <- substr(REP_TREAT,(nchar(REP_TREAT)+1)-2,nchar(REP_TREAT))
    #
    COLONY_SIZE <- unique(metadata[which(metadata$REP_treat == REP_TREAT),"colony_size"])

    # ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
    print(paste0("Compute 3-hours analysis"))
    for (TIME_HOURS in Time_dictionary$time_hours[seq(1, length(Time_dictionary$time_hours), 3)] ){  ## increments by 3 hours for 48 hours
      # HOUR <- seq(from=0, to=48, by=3)
      # From  <- fmQueryGetDataInformations(exp)$end - 51*TimeWind + (HOUR * TimeWind) - window_shift
      # To    <- fmQueryGetDataInformations(exp)$end - 48*TimeWind + (HOUR * TimeWind) - window_shift
      # 
      From  <- fmQueryGetDataInformations(exp)$end + (TIME_HOURS - 24)*TimeWind - window_shift
      To    <- fmQueryGetDataInformations(exp)$end + (TIME_HOURS - 21)*TimeWind - window_shift
      #############################
      
      print(paste0("computing hour ",TIME_HOURS))
      print(paste("Time window, from", From, "to", To))
      start <- fmTimeCreate(offset=From) #end minus 48 hours plus incremental time
      end   <- fmTimeCreate(offset=To) #end minus 45 hours plus incremental time
      
      #base file information
      #PERIOD
      TimeDiff <- difftime(exp_end, To, units = "hours")
      if(TimeDiff < 24){ PERIOD <- "post"
      }else if ( TimeDiff >= 27 & TimeDiff < 51) { PERIOD <- "pre"  } else{ PERIOD <- "EXPOSURE_GAP"}
      Period_dt<- data.frame(From, To, PERIOD)
      Period_dataframe <- rbind(Period_dataframe, Period_dt)
      #TIME_OF_DAY
      TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS),"time_of_day"]
    #}      
    # 
    # table(Period_dataframe$PERIOD) # shall be equal!
    
    # RUN FOR PRE AND POST (skip the 3h exposure gap)
    if (!Period_dt$PERIOD=="EXPOSURE_GAP") {

  ############ NETWORK PROPERTIES: INDIVIDUAL & COLONY LEVEL ##########################    
      
  # COMPUTE NETWORK
  Graph               <- compute_G(exp = exp, start = start, end=end)
  # COMPUTE NETWORK PROPERTIES
  Network_prop_hour   <- NetProperties(graph=Graph)
  #Add metadata info
  Network_prop_hour <- cbind(data.frame(randy=REP.FILES,REP_treat=REP_TREAT,colony_size=COLONY_SIZE,period=PERIOD,time_hours=TIME_HOURS, time_of_day=TIME_OF_DAY,From, To,
                                        Network_prop_hour,
                                        stringsAsFactors = F))
  #stack
  Network_properties  <- rbind(Network_properties,Network_prop_hour)
  
  ############ INDIVIDUAL LEVEL ##########################    
  # SPACE USAGE
  SpaceUsage_hour     <- SpaceUse(exp = exp, start = start, end=end)
  
  #Add metadata info
  SpaceUsage_hour   <- cbind(data.frame(randy=REP.FILES,REP_treat=REP_TREAT,colony_size=COLONY_SIZE,period=PERIOD,time_hours=TIME_HOURS, time_of_day=TIME_OF_DAY,From, To,
                                        SpaceUsage_hour,
                                        stringsAsFactors = F))
  
  #stack
  SpaceUsage          <- rbind(SpaceUsage,SpaceUsage_hour)

  #####################################################################################################
  
  } }#REP LOOP
    
    #### ADD EXTRA REPLICATE INFORMATIONS

    #add status (large, small) info to Network_properties
    Network_properties <- dplyr::left_join(Network_properties,metadata[c("size_treat","status","treatment","REP_treat")], by = "REP_treat")         # Apply left_join dplyr function
    
    #add task, exposure and status (large, small)  info to SpaceUsage
    SpaceUsage <- dplyr::left_join(SpaceUsage,metadata[c("size_treat","status","treatment","REP_treat","antID","Exposed","AntTask")], by = c("REP_treat","antID"))         # Apply left_join dplyr function

    
    ########################################
    ##### SAVE FILES IN FOLDER #############

    ## Network properties save (saved INSIDE the Network_analysis folder)
    if (file.exists(NET_properties)){
      write.table(Network_properties,file=NET_properties,append=T,col.names=F,row.names=F,quote=T,sep=",")
    }else{
      write.table(Network_properties,file=NET_properties,append=F,col.names=T,row.names=F,quote=T,sep=",")
    }
    
    ## Space Use save (saved INSIDE the Network_analysis folder)
    if (file.exists(SPACE_USE)){
      write.table(SpaceUsage,file=SPACE_USE,append=T,col.names=F,row.names=F,quote=T,sep=",")
    }else{
      write.table(SpaceUsage,file=SPACE_USE,append=F,col.names=T,row.names=F,quote=T,sep=",")
    }
    
    
    
    #start fresh
    Network_properties            <- data.frame()
    SpaceUsage                    <- data.frame()
    
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










################################################################################################################
########## THE INFORMATION ON METADATA IS USELESS AND MISLEADING, REMOVE ######################################
###############################################################################################################

# ########## GET EXPOSED ANTS from METADATA
# exp.Ants <- exp$ants
# metadata$Exposed <- "no"
# for (ant in exp.Ants){
#   individual  <- ant$ID
#   #print(ant)
#     if (TRUE %in% ant$getValues("Exposed")[,"values"]) { metadata[individual,"Exposed"] <- "exposed" }
# }
# 
# ########## GET QUEENS 
# metadata$IsQueen <- "no"
# for (ant in exp.Ants){
#   individual  <- ant$ID
#   #print(ant)
#   if (TRUE %in% ant$getValues("IsQueen")[,"values"]) { metadata[individual,"IsQueen"] <- "queen" }
# }
