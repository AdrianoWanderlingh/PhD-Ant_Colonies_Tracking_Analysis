# Metadata extraction file

####################################################################################
#### THIS SCRIPT CONTAINS:
#### METADATA EXTRACTION FROM MYRMIDON FILES
####################################################################################

#### LIBRARIES
library(FortMyrmidon) ####R bindings

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

############# GET METADATA FUNCTION #########################
Get_metadata <- function() {
  print(REP.FILES) ##}}
  #open experiment
  e <- fmExperimentOpen(REP.FILES)
  exp.Ants  <- e$ants
  exp_end   <- fmQueryGetDataInformations(e)$end
  exp_start <- fmQueryGetDataInformations(e)$start
  
  ########### COMPUTE THE ANT TASKS (48h before exposure)
  #make sure that the file does not exist already
  print(paste0("Computing Ant Tasks and Zone Uses"))
  #re-source here as it keeps getting deleted
  source(paste(SCRIPTDIR,"AntTasks_v082.R",sep="/"))
  AntTasks <- AntTasks(e=e)
  
  ############# CREATE BASE FILE
  metadata <- NULL
  
  for (ant in exp.Ants){
    for (id in ant$identifications){
      metadata <- rbind(metadata,data.frame(REP_treat        = REP_treat,
                                            size_treat       = substr(REP_treat,(nchar(REP_treat)+1)-2,nchar(REP_treat)),
                                            status           = NA,
                                            antID            = ant$ID,
                                            tagIDdecimal     = id$tagValue,
                                            identifStart     = capture.output(print(id$start)), 
                                            identifEnd       = capture.output(print(id$end)), 
                                            AntTask1perc     = AntTasks[which(AntTasks$antID==ant$ID),"AntTask"],
                                            prop_time_outside= round(AntTasks[which(AntTasks$antID==ant$ID),"prop_time_outside"],3),
                                            box              = e$spaces[[1]]$name, # tracking system
                                            #AntTask_num          = AntTasks[which(AntTasks$antID==ant$ID),"AntTask_num"],
                                            stringsAsFactors = F))
      }
    }
  
  metadata[grepl("∞", metadata$identifEnd),"identifEnd"] <- NA
  metadata[grepl("∞", metadata$identifStart),"identifStart"] <- NA
  
  #transform time in GMT
  metadata$identifEnd <- as.POSIXct( metadata$identifEnd , format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
  
  #colony-wide metadata
  metadata$colony_size <- length(unique(metadata$antID))
  
  metadata$status[metadata$size_treat=="BS"] <- "Big Sham"
  metadata$status[metadata$size_treat=="BP"] <- "Big Pathogen"
  metadata$status[metadata$size_treat=="SS"] <- "Small Sham"
  metadata$status[metadata$size_treat=="SP"] <- "Small Pathogen" 
  
  #split in two cols output
  treat_colz <-     data.frame(do.call('rbind',strsplit(metadata$status,' ',fixed=TRUE)))
  metadata$status <- treat_colz$X1
  metadata$treatment <- treat_colz$X2
  
  #empty metadata columns
  metadata$Comment     <- NA
  metadata$Exposed     <- NA 
  metadata$IsAlive     <- NA 
  metadata$IsQueen     <- NA 
  metadata$surviv_time <- NA # time of survival: time of death or end of exp
  
  ############# METADATA KEYS AND VALUES #########################
  list_keys <- list()
  #assign metadata keys
  for (KEY in   1:length(e$metaDataKeys)) {
    key <- names(e$metaDataKeys[KEY])
    #defaultvalue <- unname(e$metaDataKeys[KEY][[1]])
    list_keys <- c(list_keys,key)
  }
  
  ## assign metadata values
  for (ant in exp.Ants){
    individual  <- ant$ID
    #extract metadata info per key
    for (METADATA_KEY in list_keys) {
      #for more than one row, always last row will be picked (the relevant one if there is a timed change or the default one if there is not)
      for(ROW in 1:nrow(ant$getValues(METADATA_KEY))) {
        #assign metadata_key value when ID corresponds
        metadata[which(metadata$antID==ant$ID),METADATA_KEY] <- ant$getValues(METADATA_KEY)[ROW,"values"]
        #if the ant died
        if (METADATA_KEY=="IsAlive") {
          if (ant$getValues("IsAlive")[ROW,"values"]==FALSE) {
            metadata[which(metadata$antID==ant$ID),"surviv_time"] <- ant$getValues("IsAlive")[ROW,"times"]
            # if didn't die    
          }else if (ant$getValues("IsAlive")[ROW,"values"]==TRUE) {
            metadata[which(metadata$antID==ant$ID),"surviv_time"] <- as.POSIXct( exp_end,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
          }} # IsAlive check
      } # ROW
    } #METADATA_KEY
  }#ant
  metadata$surviv_time <- as.POSIXct(metadata$surviv_time, origin="1970-01-01", tz="GMT")
  
  # relabel Q as queen, not nurse
  metadata[which(metadata$IsQueen==TRUE),"AntTask1perc"] <- "queen"
  
  # N_exposed in the colony
  metadata$N_exposed <- sum(metadata$Exposed)
  #start and end times
  metadata$ExpStart  <-  exp_start
  attr(metadata$ExpStart,"tzone") <- "GMT"
  metadata$ExpEnd    <-  exp_end
  attr(metadata$ExpEnd,"tzone") <- "GMT"
  
  # save
  if (file.exists(Metadata_Exp1)){
    write.table(metadata,file=Metadata_Exp1,append=T,col.names=F,row.names=F,quote=T,sep=",")
  }else{
    write.table(metadata,file=Metadata_Exp1,append=F,col.names=T,row.names=F,quote=T,sep=",")
  }
  gc()
  #remove experiment and all
  #clean up
  rm(list=ls()[which(!ls()%in%c("REP.folder","REP.files","REP.filefolder","files_list","Metadata_Exp1","SCRIPTDIR","WORKDIR","DATADIR","list.dirs.depth.n","AntTasks"))]) #close experiment
}# get metadata



############# ACCESS FILES #########################
USER <- "supercompAdriano"

if (USER=="Adriano") {
WORKDIR <- "/media/cf19810/DISK4/ADRIANO" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
DATADIR <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
SCRIPTDIR <- "/media/cf19810/DISK4/EXP1_base_analysis/EXP1_analysis scripts"
}
if (USER=="supercompAdriano") {
  WORKDIR <- "/media/cf19810/DISK4/ADRIANO" # "/media/cf19810/Seagate Portable Drive/ADRIANO"
  DATADIR <- paste(WORKDIR,"EXPERIMENT_DATA",sep="/") # paste(WORKDIR,"EXPERIMENT_DATA_EXTRAPOLATED",sep="/")
  SCRIPTDIR <- "/media/cf19810/DISK4/EXP1_base_analysis/EXP1_analysis scripts"
}

###source function scripts
print("Loading functions and libraries...")
source(paste(SCRIPTDIR,"AntTasks_v082.R",sep="/"))

###define name of general output table
Metadata_Exp1 <- file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt") 

# #inferred$return_time <- NA
# inferred$end_time <- NA
# inferred$exposed <- "no"
# inferred$dead <- "no"
#inferred$time_since_start <- NA

### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n(DATADIR, n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]

# replicate folder
for (REP.n in 1:length(files_list)) {
  # REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = "CapsuleDef2018_q.myr")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  #replicate file
  for (REP.FILES in REP.filefolder) {
    #get substring in variable until R9BS_
    # REP.FILES <- REP.filefolder[2]
    REP_treat <- sub("\\_.*", "", basename(REP.FILES))
    #treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))
    
    # #check if the colony has already been recorded in the metadata{}
    # if(file.exists(file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt"))){
    # metadata_present <- read.table(file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt"),header=T,stringsAsFactors = F, sep=",")
    # if (REP_treat %in% unique(metadata_present$REP_treat) ) {
    #   print(paste0(REP_treat," already present in metadata, skip"))
    # }
    # } else {
      
      #check if the metadata file exists
      if(file.exists(file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt"))) {
        metadata_present <- read.table(file.path(DATADIR,"Metadata_Exp1_2021_2023-02-27.txt"),header=T,stringsAsFactors = F, sep=",")
        #check if the colony has already been recorded in the metadata{}
        if (REP_treat %in% unique(metadata_present$REP_treat)) {
          print(paste0(REP_treat," already present in metadata, skip"))
        } else {
          Get_metadata() # the code that you want to execute if REP_treat is not present in metadata_present$REP_treat
        }
      } else {
        Get_metadata() # the code that you want to execute if the file does not exist
      }
    } # REP folders
  } # REP by REP
