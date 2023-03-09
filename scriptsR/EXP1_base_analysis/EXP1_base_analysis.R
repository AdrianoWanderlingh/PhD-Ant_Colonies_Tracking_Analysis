rm(list = ls())
gc()
mallinfo::malloc.trim(0L)

# NOTE ON DIFFERENCES FROM SCIENCE 2018 PAPER:
# period is ("pre","post"), not ("after","before")
# treatment is ("Pathogen","Sham"), not ("pathogen","sham")
# status is linked to size of the colony, not ("treated","untreated") ant

# to do in MARCH 2023
# - proportion_time_active
# - average_bout_speed_pixpersec
# total_distance_travelled_pix

# ------------------------------------
# DONE in NETWORK OUTPUT:
# COLUMNS:
# colony    colony_size x
# treatment (two values: pathogen and control)  x
# tag    (only include treated workers) x
# age    x
# status    ( replace the content of the "status" column with either "large" or "small" (rather than "treated" or "untreated") ) x
# period    ( pre/post chunks corresponding to the same time of day should have the same value in column "time_of_day") x
# time_hours    x
# time_of_day x

##############################################################################
################### IMPORTANT THINGS TO MODIFY ###############################

# UNIFORM Plot_Grooming_Pre-Post.R to the output of this script and the sourced functions
# AGGREGATE SPACE USE FUNCTION FOR tag_hex_ID (we want 1 row per ant!)

# ########### DEFINE ZONES PROPERLY - DEFINED INSIDE THE ANT TASKS FUNCTION, SHOULD DO SAME FOR NEST USE
# #### EXPAND THIS TO EXTRACT ALL ZONES (WATER, SUGAR, ETC)
# zones <- e$spaces[[1]]$zones #function to show the Zones present in the Space
# zones_tab <- data.frame(ID =c(zones[[1]]$ID, zones[[2]]$ID), name=c(zones[[1]]$name, zones[[2]]$name))
# foraging_zone <- zones_tab[which(grepl("forag",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
# nest_zone <- zones_tab[which(grepl("nest",zones_tab$name)),"ID"]##fool-proofing - this way we are sure to always select the right zone
# print(paste("Foraging zone = zone",foraging_zone, "& Nest zone = zone",nest_zone))

# change window_shift to actual date? As per Grooming (new exact data info to be gathered)
# THE EXACT DATE IS IN /home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data/Grooming_Classifier_CrossVal_Adriano2022_RETURN_EXP_TIME_ZULU.csv

##############################################################################
##############################################################################
# "https://formicidae-tracker.github.io/myrmidon/latest/index.html"

##### LIBRARIES
library(data.table)
library(lubridate)
library(pals)
library(dplyr)
# library(reshape)
library(adehabitatHR)
library(adehabitatLT)
library(changepoint)
library(e1071)
library(igraph)
library(gtools)
library(Rcpp)
library(survival)

# starting params
# USER <- "2A13_Office" # Nath_office 
USER <- "2A13_Office" # Nath_office 

if (USER == "2A13_Office") {
  usr <- "cf19810"
} else {
  usr <- "bzniks"
}
  
SAVEDIR <- paste("/media/",usr,"/DISK4/EXP1_base_analysis/EXP_summary_data",sep="")
#where to save the interactions in the same structure as for Science 2018
#I copied the folder structure from Nathalie, ideally one should create it via the terminal:
# use the file /media/cf19810/DISK4/Lasius-Bristol_pathogen_experiment/dirs.txt
# and in the terminal use: xargs mkdir -p < dirs.txt
INTDIR <- paste("/media/",usr,"/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment/intermediary_analysis_steps",sep="")
BEHDIR <- paste("/media/",usr,"/DISK4/Lasius-Bristol_pathogen_experiment/main_experiment/processed_data/individual_behaviour",sep="")

  WORKDIR <- paste("/media",usr,"DISK4/ADRIANO",sep="/")
  DATADIR <- paste(WORKDIR, "EXPERIMENT_DATA", sep = "/")
  SCRIPTDIR <- paste("/media",usr,"DISK4/EXP1_base_analysis/EXP1_analysis scripts",sep="/") # "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP1_analysis scripts"
  metadata <- read.table(paste(DATADIR, "/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
  BEH_FUNCTIONS <-  paste("/media",usr,"DISK4/Ants_behaviour_analysis/ScriptsR",sep="/")
  BODYLENGTH_FILE <- paste(BEH_FUNCTIONS,"Mean_ant_length_per_TrackingSystem.txt", sep = "/")
  metadata_info         <- read.csv(paste(DATADIR,"Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU_All_colonies.csv",sep = "/"), sep = ",")

##metadata
  metadata$status_ant <- NA
  metadata$status_ant <-ifelse(metadata$Expose==TRUE,"treated","untreated")
  
# #### ACCESS FILES
# if (USER == "Adriano") {
#   WORKDIR <- "/media/cf19810/DISK4/ADRIANO"
#   DATADIR <- paste(WORKDIR, "EXPERIMENT_DATA", sep = "/")
#   SCRIPTDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP1_analysis scripts"
#   metadata <- read.table(paste(DATADIR, "/Metadata_Exp1_2021_2022-10-20.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
#   warning("METADATA MISSING FOR R9SS AND R9BS")
#   BEH_FUNCTIONS <- "/home/cf19810/Documents/scriptsR/Ants_behaviour_analysis/ScriptsR"
#   BODYLENGTH_FILE <- paste(BEH_FUNCTIONS,"Mean_ant_length_per_TrackingSystem.txt", sep = "/")
#   metadata_info         <- read.csv(paste(DATADIR,"Grooming_Classifier_CrossVal_RETURN_EXP_TIME_ZULU_All_colonies.csv",sep = "/"), sep = ",")
#   
# }

### source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR, "SpaceUse_v082.R", sep = "/")) # SPACE_USE
source(paste(SCRIPTDIR, "NetworkProperties_v082.R", sep = "/")) # NET_properties collective + individual
suppressMessages(source(paste(BEH_FUNCTIONS, "BEH_libraries.R", sep = "/")))
source(paste(BEH_FUNCTIONS, "interaction_detection.R", sep = "/"))
source(paste(BEH_FUNCTIONS, "trajectory_extraction.R", sep = "/"))
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # required by merge_interactions.cpp
sourceCpp(paste(BEH_FUNCTIONS, "merge_interactions.cpp", sep = "/"))
# list.functions.in.file(paste(SCRIPTDIR,"NetworkProperties_v082.R",sep="/"), alphabetic = TRUE)

#### FUNCTIONS
# list files recursive up to a certain level (level defined by "n" parameter)
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n - 1)
    c(res, add)
  } else {
    res
  }
}

# body_lenghts
get_body_lengths <- function(e, all_body_lengths) {
  ### check how many spaces there are
  space_list <- e$spaces
  body_lengths <- NULL
  for (space_ID in 1:length(space_list)) {
    body_lengths <- rbind(body_lengths, data.frame(
      space_ID = space_list[[space_ID]]$ID,
      tracking_system = space_list[[space_ID]]$name,
      body_length_px = all_body_lengths[which(all_body_lengths$TS == space_list[[space_ID]]$name), "mean_worker_length_px"],
      body_length_mm = all_body_lengths[which(all_body_lengths$TS == space_list[[space_ID]]$name), "mean_worker_length_mm"]
    ))
  }
  return(body_lengths)
}


#### PARAMETERS
all_body_lengths <- read.table(BODYLENGTH_FILE, header = T, stringsAsFactors = F, sep = ",") ### body length information
TimeWind <- 3600 ## in seconds (3600 is an hour)
N_DECIMALS <- 3
MAX_INTERACTION_GAP <- 10 # fmSecond(10) # the maximum gap in tracking before cutting the interactions in two different object.  Used in ComputeGraph function
MAX_DIST_GAP_MM <- 0.6
FRAME_RATE <- 8
interactions_of_interest <- list(c("head", "body")) ### if interested in more than one type, list them as follows: interactions_of_interest <- c(list(c("head","body")),list(c("head","head")))

Time_dictionary <- data.frame(time_hours = -36:35, time_of_day = rep(0:23, 3)) # time_hours= time since exposure (negatives pre-exposire, positive post exposure), time_of_day=time hour slot of the day
Time_dictionary <- Time_dictionary[which(Time_dictionary$time_hours <= 21 & Time_dictionary$time_hours >= -27), ]
Time_dictionary$period <- ifelse(Time_dictionary$time_hours<0, "pre", "post")

#OUTPUT_FOLDER <-  "Exp1_Results_2023" # paste0("Exp1_Results_", Sys.Date()) #OUT OF DATE?
SPACE_USE     <-  file.path(BEHDIR,"pre_vs_post_treatment","individual_behavioural_data.txt")

#### FLAGS
RUN_INTERACT     <- FALSE
RUN_SPACEUSE     <- TRUE
RUN_NETWORKS     <- FALSE
warning(paste("RUN_INTERACT is set to:",RUN_INTERACT,
              "\nRUN_SPACEUSE is set to:",RUN_SPACEUSE,
              "\nRUN_NETWORKS is set to:",RUN_NETWORKS,sep="\t"))

## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED!
#window_shift <- 60 * 15 # approx N of minutes that where given at the end as leeway, minutes can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition

# list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
# select REP folders
files_list <- files_list[grep("REP", files_list)]

### initialise general output folder
### remove folder if already exists to make sure we don't mix things up
# if (file.exists(file.path(DATADIR, "NetworkAnalysis_outcomes"))){
#   unlink(file.path(DATADIR, "NetworkAnalysis_outcomes"),recursive=T)
# }

# ### create folder
# dir.create(file.path(SAVEDIR, OUTPUT_FOLDER), recursive = T) 
# ### define name of general output files
# NET_properties_collective <- file.path(SAVEDIR, OUTPUT_FOLDER, "NetworkProp_collective.txt") # (saved INSIDE the folder)
# NET_properties_individual <- file.path(SAVEDIR, OUTPUT_FOLDER, "NetworkProp_individual.txt") # (saved INSIDE the folder)

##### RUNNING TIME
loop_start_time <- Sys.time()


###create object that will contain time limits for annotations_all (contains all annotations for all behaviour so true time limit for period analysed)
# time_window_all        <- merge(aggregate(T_start_UNIX ~ PERIOD + REPLICATE, FUN=min, data=annotations_all),aggregate(T_stop_UNIX  ~ PERIOD + REPLICATE, FUN=max, data=annotations_all))
# names(time_window_all) <- c("PERIOD","REPLICATE","time_start","time_stop")
time_window_all        <- data.frame(
  PERIOD = "all"
  ,REPLICATE = metadata_info$REP_treat
  ,return_time = metadata_info$ReturnExposed_time #time_start
)
# 
time_window_all$return_time <- as.POSIXct(time_window_all$return_time, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
time_window_all$time_start <- time_window_all$return_time - 27*60*60 #(24h + 3h gap)
time_window_all$time_stop  <- time_window_all$return_time + 24*60*60
# 

#### define to_keep variables to keep clearing memory between runs
to_keep <- c(ls(), c("to_keep"))

#### OPEN REPLICATE
# replicate folder
for (REP.n in 1:length(files_list)) {
  # REP.n <- 1    #temp
  REP.folder <- files_list[REP.n]
  REP.files <- list.files(REP.folder, pattern = "CapsuleDef2018_q") #AntsCreated_AutoOriented_withMetaData_NS_NS_q_
  ##
  ### select only for the science2018 capsules, the name to grep for is in NOTION
  ##
  #REP.files <- REP.files[!grepl("CapsuleDef", REP.files)]
  REP.filefolder <- paste(REP.folder, REP.files, sep = "/")

  # replicate file
  for (REP.FILES in REP.filefolder) {
    # REP.FILES <-  REP.filefolder[2]   #temp

    ## some initialization
    Period_dataframe <- NULL # checking time correspondances
    # start fresh
    NetworkProp_collective <- data.frame()
    NetworkProp_collective_hour <- data.frame()
    NetworkProp_individual <- data.frame()
    NetworkProp_individual_hour <- data.frame()
    
    #Interactions_total <- data.frame()
    Interactions_REP_TREAT <- data.frame()
    #SpaceUsage <- data.frame()

    cat(paste("########################################\n",basename(REP.FILES)),sep="") ## }}
    # open experiment
    e <- fmExperimentOpen(REP.FILES)
    # e.Ants <- e$ants
    print(paste0("Processing ", basename(REP.FILES)))

    ### get body length
    body_lengths <- get_body_lengths(e, all_body_lengths)

    # base file info
    # Split the string into a vector based on "_"
    input_vector <- strsplit(basename(REP.FILES), "_")[[1]]
    # Extract the first 2 elements of the vector
    # Use grep to find the element containing "R"
    REP_TREAT <- input_vector[grep("R", input_vector[1:2])]
    # REP_TREAT <- sub("\\_.*", "", )
    # SIZE_TREAT <- substr(REP_TREAT,(nchar(REP_TREAT)+1)-2,nchar(REP_TREAT))
  
  #   
  #   REPS <- c(REPS, REP_TREAT) #
  # 
  # }}

    
    #####
    # end time is return time + 24h
    exp_end <- time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"return_time"] + 60*60*24 # 24h
    #exp_end <- fmQueryGetDataInformations(e)$end - window_shift
    
#     if (RUN_NETWORKS) {
#     COLONY_SIZE <- unique(metadata[which(metadata$REP_treat == REP_TREAT), "colony_size"])
# }
# 
#     ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
#     print(paste0("Compute 3-hours analysis"))
#     # TIME_HOURS zero is the moment of exposed ants return
#     for (TIME_HOURS in Time_dictionary$time_hours[seq(1, length(Time_dictionary$time_hours), 3)]) { ## increments by 3 hours for 48 hours
#       # HOUR <- seq(from=0, to=48, by=3)
#       # From  <- fmQueryGetDataInformations(e)$end - 51*TimeWind + (HOUR * TimeWind) - window_shift
#       # To    <- fmQueryGetDataInformations(e)$end - 48*TimeWind + (HOUR * TimeWind) - window_shift
#       #
#       From <- fmQueryGetDataInformations(e)$end + (TIME_HOURS - 24) * TimeWind - window_shift
#       To <- fmQueryGetDataInformations(e)$end + (TIME_HOURS - 21) * TimeWind - window_shift
#       #############################
# 
#       print(paste0("computing hour ", TIME_HOURS))
#       print(paste("Time window, from", From, "to", To))
#       time_start <- fmTimeCreate(offset = From) # time_stop minus 48 hours plus incremental time
#       time_stop <- fmTimeCreate(offset = To) # time_stop minus 45 hours plus incremental time
# 
#       # base file information
#       # PERIOD
#       TimeDiff <- difftime(exp_end, To, units = "hours")
#       if (TimeDiff < 24) {
#         PERIOD <- "post"
#       } else if (TimeDiff >= 27 & TimeDiff < 51) {
#         PERIOD <- "pre"
#       } else {
#         PERIOD <- "EXPOSURE_GAP"
#       }
#       Period_dt <- data.frame(From, To, PERIOD)
#       Period_dataframe <- rbind(Period_dataframe, Period_dt)
#       # TIME_OF_DAY
#       TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS), "time_of_day"]
#       # }
#       #
#       # table(Period_dataframe$PERIOD) # shall be equal!
# 
#       if (RUN_SPACEUSE) {
# 
#       # RUN FOR PRE AND POST (skip the 3h exposure gap)
#       if (!Period_dt$PERIOD == "EXPOSURE_GAP") {
# 
# 
#         # ############ NETWORK PROPERTIES: INDIVIDUAL & COLONY LEVEL ##########################
#         #
#         # # COMPUTE NETWORK
#         # Graph <- compute_G(e = e, start = time_start, end = time_stop, gap = MAX_INTERACTION_GAP)
#         # # COMPUTE NETWORK PROPERTIES (collective and individual)
#         # Network_summ_prop_hour <- NetProperties(graph = Graph)
#         #
#         # ### Collective
#         # # Add metadata info
#         # NetworkProp_collective_hour <- cbind(data.frame(
#         #   randy = REP.FILES, REP_treat = REP_TREAT, colony_size = COLONY_SIZE, period = PERIOD, time_hours = TIME_HOURS, time_of_day = TIME_OF_DAY, From, To,
#         #   Network_summ_prop_hour$summary_collective,
#         #   stringsAsFactors = F
#         # ))
#         # # stack
#         # NetworkProp_collective <- rbind(NetworkProp_collective, NetworkProp_collective_hour)
#         #
#         # ### Individual
#         # # Add metadata info
#         # NetworkProp_individual_hour <- cbind(data.frame(
#         #   randy = REP.FILES, REP_treat = REP_TREAT, colony_size = COLONY_SIZE, period = PERIOD, time_hours = TIME_HOURS, time_of_day = TIME_OF_DAY, From, To,
#         #   Network_summ_prop_hour$summary_individual,
#         #   stringsAsFactors = F
#         # ))
#         # # stack
#         # NetworkProp_individual <- rbind(NetworkProp_individual, NetworkProp_individual_hour)
# 
#       }
#         }
#       } # REP LOOP
    

    ### PERFORM THE FULL INTERACTION AND NETWORK ANALYSIS OUTSIDE OF THE HOURLY LOOP #########################
    
    # # # Extract the minimum "From" time and maximum "To" time for each PERIOD
    Period_windows <- data.frame(Period = c("pre","post")
                                 , From = c(time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"time_start"], #skip 3h gap
                                            time_window_all[which(time_window_all$REPLICATE==REP_TREAT),"return_time"]
                                 ))
    Period_windows$To <-   Period_windows$From + 24*60*60
    
    
    for (PERIOD in c("pre","post")) {
    
      #segment the time dictionary
      Time_dictionary_PERIOD <- Time_dictionary[which(Time_dictionary$period==PERIOD),]
      #conform naming to science2018
      # colony code
      REP_NUM           <- substring(REP_TREAT, 2, nchar(REP_TREAT))
      colony            <- paste("colony", paste(rep(0,4-nchar(REP_NUM)),collapse=""), REP_NUM,sep="")
      # colony_status
      status_char       <- substr(REP_TREAT, nchar(REP_TREAT), nchar(REP_TREAT))
      colony_status     <- ifelse(status_char == "P", "pathogen", ifelse(status_char == "S", "control", NA))
      # treatment code
      period_code       <- ifelse(PERIOD == "pre", "PreTreatment", ifelse(PERIOD == "post", "PostTreatment", NA))
      size_char         <- substr(REP_TREAT, nchar(REP_TREAT)-1, nchar(REP_TREAT)-1)
      size_status       <- ifelse(size_char == "S", "small", ifelse(size_char == "B", "big", NA))
      treatment_code    <- paste(colony_status,size_status,sep=".")
      
      # Select the full PERIOD (24h)
      cat(paste("#######","period:", PERIOD, sep= " "))
      #time_start <- as.numeric(Period_windows[which(Period_windows$Period==PERIOD),"From"]) # time_stop minus 48 hours plus incremental time
      time_start <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"From"]) # time_stop minus 48 hours plus incremental time
      time_stop <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"To"]) # time_stop minus 45 hours plus incremental time
      ##TEMP TIME STOP
      # warning("using tiny time window for testing")
      # time_stop <- fmTimeCreate(offset = Period_windows[which(Period_windows$Period==PERIOD),"From"] + 2*60)
      
      
      if (RUN_INTERACT) {
      # DESTINATION FOLDER
      INTERACTIONS_FULL   <-  file.path(INTDIR,"full_interaction_lists",period_code,"observed", 
                                 paste(colony,treatment_code,period_code,"interactions.txt",sep="_"))

    ############ GET INTERACTIONS ##########################
    if(file.exists(INTERACTIONS_FULL)){
      warning(paste(REP_TREAT,"|",colony,treatment_code,period_code,"| already present, SKIP >>", sep=" "))
    }
      
    if(!file.exists(INTERACTIONS_FULL)){
    
      # INTERACTIONS IN THIS FUNCTION ARE CALCULATED ACCORDING TO STROEYMEYT ET AL, SCIENCE 2018
      Interactions <- compute_Interactions(e = e, start = time_start, end = time_stop, max_time_gap = MAX_INTERACTION_GAP)
      
      # Remove interactions involving dead ants
      dead_by_REP <- metadata[which(metadata$IsAlive==FALSE & metadata$REP_treat==REP_TREAT),"antID"]
      Interactions <- Interactions[!(Interactions$ant1 %in% dead_by_REP | Interactions$ant2 %in% dead_by_REP), ]
      
      # create time vars
      Interactions$time_hours   <- NA
      Interactions$time_of_day  <- NA
      # TIME_HOURS zero is the moment of exposed ants return
      warning("interaction binning loop should be fixed as done later for the SpaceUse, see TIME_HOURS ")
      for (TIME_HOURS in Time_dictionary$time_hours[seq(1, length(Time_dictionary$time_hours), 3)]) { ## increments by 3 hours for 48 hours
        
        # TIME_OF_DAY
        TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS), "time_of_day"]
        #time windows
        From_TIME_HOURS <- as.numeric((time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind),N_DECIMALS)
        To_TIME_HOURS <- as.numeric((time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 21) * TimeWind),N_DECIMALS)
        
        #assing time labels
        Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"time_hours"]   <- TIME_HOURS
        Interactions[which(Interactions$Starttime >= From_TIME_HOURS & Interactions$Starttime <= To_TIME_HOURS),"time_of_day"]  <- TIME_OF_DAY
      }
      
      # Add metadata info
      Interactions <- cbind(data.frame(
        period = PERIOD,Interactions, REP_treat = REP_TREAT, colony=colony, treatment=treatment_code,
        stringsAsFactors = F
      ))
      
      # extras not present in Science files: (added at the end of the file output)
      #"REP_treat","period","ant1.zones","ant2.zones","duration"
      Interactions <- Interactions[,c("Tag1","Tag2","Startframe","Stopframe","Starttime","Stoptime","Box","Xcoor1","Ycoor1","Angle1","Xcoor2","Ycoor2","Angle2","Direction","Detections","time_hours","time_of_day","colony","treatment","REP_treat","period","ant1.zones","ant2.zones","duration")]
      
      ## Interactions save (saved INSIDE the Network_analysis folder)
      #if (file.exists(INTERACTIONS_FULL)) {
      #  write.table(Interactions, file = INTERACTIONS_FULL, append = T, col.names = F, row.names = F, quote = F, sep = ",")
      #} else {
        write.table(Interactions, file = INTERACTIONS_FULL, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
      #}
      
      ### split output into bins
        ### HOURLY LOOP; loading 24 hours of data requires >32 GB RAM & causes R to crash, so we must load the contacts in 3h blocks, stacking vertically
        print(paste0("split files into 3-hours bins"))
        
        for (TIME_HOURS in unique(Interactions$time_hours)) {
          #labels for subsets
          TH <- paste0("TH",unique(Interactions[which(Interactions$time_hours==TIME_HOURS),"time_hours"]))
          TD <- paste0("TD",unique(Interactions[which(Interactions$time_hours==TIME_HOURS),"time_of_day"]))
          print(paste("TIME_HOURS", TH,"TIME_OF_DAY", TD,"PERIOD",PERIOD,sep=" "))
          
          INTERACTIONS_BINNED <-  file.path(INTDIR,"binned_interaction_lists",period_code,"observed", 
                                            paste(colony,treatment_code,period_code,TH,TD,"interactions.txt",sep="_"))
                    #save object by TH and TD
          write.table(Interactions[which(Interactions$time_hours==TIME_HOURS),], file = INTERACTIONS_BINNED, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
          }
    }
      }# RUN_INTERACT
      
      if (file.exists(SPACE_USE)) {
        print("file exists")
        spaceUse_done <- read.table(SPACE_USE, header = T, stringsAsFactors = F, sep = "\t")
        spaceUse_done <- spaceUse_done[which(spaceUse_done$colony==colony),"time_hours"]
      }else{
        spaceUse_done <- c() # placeholder
      }
        
      
      if (RUN_SPACEUSE) {
        print("Computing SpaceUse based on 24h time-window pre AND post exposure")
        warning("the Space Use Script can be substantially sped up by computing trajectories beforehand, then cutting the result by 3h chunks and performing operations on them.
                \n -This is because trajectories stitching is slow so may be better to perform it once. see inspiration from 8_process_trajectory_files.R
                \n -This will require saving the summary stats per ant in a safe place (possibly the trajectories summary as now but with clearing of the vars at each time cycle)
                \n -parallelisation of ant computation greatly improve speed. Before parallelisation of SpaceUse calculations: 2.1 sec per ant per 3h chunk (70h), after 0.21 sec with 10 cores. this includes only the summary stats computations (traj stitching excluded)")
        SpaceUse_loop_start_time <- Sys.time()
        # TIME_HOURS zero is the moment of exposed ants return
        for (TIME_HOURS in Time_dictionary_PERIOD$time_hours[seq(1, length(Time_dictionary_PERIOD$time_hours), 3)]) { ## increments by 3 hours for 48 hours
          
          
          # TIME_OF_DAY
          TIME_OF_DAY <- Time_dictionary[which(Time_dictionary$time_hours == TIME_HOURS), "time_of_day"]
          print(paste("TIME_HOURS", TIME_HOURS,"TIME_OF_DAY", TIME_OF_DAY,"PERIOD",PERIOD,sep=" "))
          
            # has the time_hour been computed already? 
            if (TIME_HOURS %in% spaceUse_done) {
              print("DONE, SKIP TO NEXT 3-HOURS BLOCK >>")
                      }else{
          
          #time windows
          time_start_h <- fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind)) # time_stop minus 48 hours plus incremental time
          time_stop_h  <- fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 21) * TimeWind)) # time_stop minus 45 hours plus incremental time
          ##TEMP TIME STOP
          # warning("using tiny time window for testing")
          #time_stop_h <-  fmTimeCreate(offset = (time_window_all[time_window_all$REPLICATE==REP_TREAT,"time_stop"] + (TIME_HOURS - 24) * TimeWind + 40*60) ) 
          
         
          ############ SPACE USE ##########################
          # SPACE USAGE
          SpaceUsage <- SpaceUse(e = e, start = time_start_h, end = time_stop_h)
          
          SpaceUsage <- dplyr::left_join(SpaceUsage, metadata[which(metadata$REP_treat==REP_TREAT),c("status_ant", "antID","IsAlive")], by = "antID") # Apply left_join dplyr function
          colnames(SpaceUsage)[which(colnames(SpaceUsage)=="antID")] <- "tag"
          colnames(SpaceUsage)[which(colnames(SpaceUsage)=="status_ant")] <- "status"
          
          #remove dead ants
          SpaceUsage <- SpaceUsage[which(SpaceUsage$IsAlive==T),]
          SpaceUsage$IsAlive <- NULL

          # Add metadata info
          SpaceUsage <- cbind(data.frame(
            colony = colony,
            colony_size = unique(metadata[which(metadata$REP_treat==REP_TREAT),"colony_size"]),
            treatment = treatment_code,
            age = NA,
            period = PERIOD,
            time_hours = TIME_HOURS,
            time_of_day = TIME_OF_DAY,
            SpaceUsage,
            REP_treat = REP_TREAT,
            stringsAsFactors = F
          ))

            ## Space Use save (saved INSIDE the Network_analysis folder)
            if (file.exists(SPACE_USE)) {
              write.table(SpaceUsage, file = SPACE_USE, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
            } else {
              write.table(SpaceUsage, file = SPACE_USE, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
            }
          } # if present already
        } # hour loop
        SpaceUse_loop_end_time <- Sys.time()
        print(paste("spaceUse 3h chunk took ", as.numeric(difftime(SpaceUse_loop_end_time, SpaceUse_loop_start_time, units = "secs")), " sec to complete"))
        }
      } # PERIOD
    }
    
  
    #### ADD EXTRA REPLICATE INFORMATIONS
  if (RUN_NETWORKS) {
    # add status (large, small) info to NetworkProp_collective
    NetworkProp_collective <- dplyr::left_join(NetworkProp_collective, unique(metadata[c("size_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function

    # add status (large, small) info to NetworkProp_individual
    NetworkProp_individual <- dplyr::left_join(NetworkProp_individual, unique(metadata[c("size_treat", "status", "treatment", "REP_treat")]), by = "REP_treat") # Apply left_join dplyr function

    ########################################
    ##### SAVE FILES IN FOLDER #############
  
    ## Network properties Collective save (saved INSIDE the Network_analysis folder)
    if (file.exists(NET_properties_collective)) {
      write.table(NetworkProp_collective, file = NET_properties_collective, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
    } else {
      write.table(NetworkProp_collective, file = NET_properties_collective, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
    }

    ## Network properties Individual save (saved INSIDE the Network_analysis folder)
    if (file.exists(NET_properties_individual)) {
      write.table(NetworkProp_individual, file = NET_properties_individual, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
    } else {
      write.table(NetworkProp_individual, file = NET_properties_individual, append = F, col.names = T, row.names = F, quote = F, sep = "\t")
    }
}

    # start fresh
    NetworkProp_collective <- data.frame()
    NetworkProp_individual <- data.frame()
    Interactions           <- data.frame()
    SpaceUsage             <- data.frame()

    # cleaning
    rm(list = ls()[which(!ls() %in% to_keep)])
    gc()
    mallinfo::malloc.trim(0L)
  }

loop_end_time <- Sys.time()
print(paste("loop took ", as.numeric(difftime(loop_end_time, loop_start_time, units = "mins")), " minutes to complete"))
