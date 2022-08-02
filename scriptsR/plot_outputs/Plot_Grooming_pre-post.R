####################################################################################
#### THIS SCRIPT CONTAINS:
#### DATA MANIPULATION TO PLOT INFERRED GROOMING
#### extra: SOME BASE PLOTS USED IN THE IUSSI SAN DIEGO 2022 PRESENTATION (MOVE THEM?)
####################################################################################

gc()
mallinfo::malloc.trim(0L)


library(FortMyrmidon)
library("ggplot2")
library(lubridate)
library(plotrix)
library(scales)
library(car)
library(lme4)
library(Hmisc)
library("viridis")
library(stringr)
library(dplyr)

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


### DIRECTORIES
WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data"

### EXTRAS
Presentation_plots             <- FALSE #script for single interactions heatmaps

###########################################################################################################
###### THIS PART OF THE SCRIPT COLLECTS INFO FROM THE MYRMIDON FILES AND IT IS TURNED OFF AS THE
###### PRELIMINARY OUTPUT HAS BEEN SAVED ON A SPECIFIC THIRD FILE inferred_groomings_ALL_withCommonStart
###########################################################################################################

if (!file.exists(paste0(WORKDIR,"/Data/inferred_groomings_ALL_withCommonStart.txt"))) {

## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED!
window_shift <- 60*20 #approx N of minutes that where given at the end as leeway, minutes can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition
#window_shift_UNIX <- as.POSIXct(window_shift,origin = "1970-01-01",tz = "GMT")

# HAND-COLLECTED RETURN TIMES (MIN) - MORE PRECISE THAN THE SHIFT
ReturnTime_mins <- read.table(paste(WORKDIR,"/Data/ReturnTime_mins.csv",sep=""),header=T,stringsAsFactors = F, sep=",")
ReturnTime_mins$Return.hour <- paste0(ReturnTime_mins$Return.hour,":00")
ReturnTime_mins$REP_treat <- paste0("R",ReturnTime_mins$REP_treat )
ReturnTime_mins <- ReturnTime_mins[,c("REP_treat","Return.hour")]

### DIRECTORIES
REP1to7_SmallOnly <- read.table(paste(WORKDIR,"/Data/inferred_groomings_REP1to7_SmallOnly.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
REP8to14_SmallOnly <- read.table(paste(WORKDIR,"/Data/inferred_groomings_REP8to14_SmallOnly.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

#List all text files in the working directory
filenames <- list.files(DATADIR,pattern = '\\whole_experiment.txt$')
filesdir <- paste(DATADIR,filenames, sep="/")
#Read every text file with header, skipping the 1st row.
#Keep only the 5th column after reading the data.
result <- lapply(filesdir, function(x) read.table(x,header=T,stringsAsFactors = F, sep=","))
inferred_LargeCols <- do.call(rbind, result)

inferred <- rbind(REP1to7_SmallOnly,REP8to14_SmallOnly,inferred_LargeCols)

inferred$T_start_UNIX <- as.POSIXct(inferred$T_start_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
inferred$T_stop_UNIX <- as.POSIXct(inferred$T_stop_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )

#base file info
#TREATMENT
inferred$TREATMENT <- substr(inferred$REPLICATE,(nchar(inferred$REPLICATE)+1)-2,nchar(inferred$REPLICATE))
#get rep-treat to match with exp filenames
inferred$REP_treat <- sub(".*\\/", "", inferred$REPLICATE)

#inferred$return_time <- NA
inferred$end_time <- NA
inferred$exposed <- "no"
inferred$dead <- "no"
#inferred$time_since_start <- NA

### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!
#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n("/media/cf19810/DISK3/ADRIANO/EXPERIMENT_DATA", n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]

Reps_N_exposed <- data.frame(REP_treat= unique(inferred$REP_treat), N_ants = NA)

# replicate folder
for (REP.n in 1:length(files_list)) {
  # REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = "CapDef3.myr")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")

  #replicate file
  for (REP.FILES in REP.filefolder) {
    #get substring in variable until R9BS_
    # REP.FILES <- REP.filefolder[1]
    REP_treat <- sub("\\_.*", "", basename(REP.FILES))
    #treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))

    # REP.FILES <-  REP.filefolder[1]   #temp
    print(REP.FILES) ##}}
    #open experiment
    exp <- fmExperimentOpen(REP.FILES)
    # exp.Ants <- exp$ants
    exp_end <- fmQueryGetDataInformations(exp)$end

    ########## GET EXPOSED ANTS # AW 17June2022
    e.Ants <- exp$ants
    Exposed_list <- vector()
    if (REP_treat %in% Reps_N_exposed$REP_treat) {
    for (ant in e.Ants){
      #exclude dead ants
      if (FALSE %in% ant$getValues("IsAlive")[,"values"] ) { #if FALSE doesn't appear
        # assign DEAD flag in inferred
        inferred[which(inferred$REP_treat==REP_treat & inferred$Rec_Name==paste0("ant_",ant$ID)),"dead"] <- "dead"
      }#exclude dead
        
      if (TRUE %in% ant$getValues("Exposed")[,"values"]) {
        exposed <-ant$ID
        Exposed_list <- c(Exposed_list, exposed) 
        # assign EXPOSED flag in inferred
        inferred[which(inferred$REP_treat==REP_treat & inferred$Rec_Name==paste0("ant_",ant$ID)),"exposed"] <- "exposed"
      } # list exposed
    }
    explist<- data.frame(Exposed_list,stringsAsFactors = F)
    #Collapse
    explist2 <- data.frame(val=paste0(explist$Exposed_list,collapse = ', '),stringsAsFactors = F)

      Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP_treat) ,"N_ants"] <- explist2
    }

    for (ROW in 1:nrow(inferred)) {
      if (inferred[ROW,"REP_treat"]==REP_treat) {
        #add end time minus
        inferred[ROW,"end_time"] <-  as.POSIXct(exp_end ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
        # inferred[ROW,"return_time"] <-  as.POSIXct( exp_end  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
        #assign a new common time for analyses
        #inferred[ROW,"time_since_start"] <-
      }
    }
    rm(list=(c("exp"))) #remove experiment
  }}

#subtract end window
inferred$end_time <- inferred$end_time - window_shift
## force unix_time
inferred$end_time <- as.POSIXct(inferred$end_time,  origin="1970-01-01", tz="GMT" )
inferred$return_time <- inferred$end_time - 24*3600

## WINDOW SHIFT GETS US THE APPROX TIME, ASSIGN THE FINAL RETURN TIME
##substitute time with hand-collected return time for increased precision
#paste the day of inferred$return_time with the time of ReturnTime_mins
inferred$Date_only <- as.character(as.Date(inferred$return_time, format = "%Y-%m-%d"))
inferred <- left_join(inferred, ReturnTime_mins, by = "REP_treat")



inferred$return_time  <- paste(inferred$Date_only,inferred$Return.hour,"GMT")
inferred$return_time  <- as.POSIXct(inferred$return_time,  origin="1970-01-01", tz="GMT" )
#the hand-collected return time is the time at which the return operation is complete. As sometimes it could have taken a few minutes, subtract 5 mins to avoid loosing grooming events
inferred$return_time  <- inferred$return_time - (60*10) #mins

## add a common time to all rows
# Negative before exposure, positive after
inferred$time_stop_since_treat <- as.numeric(difftime(inferred$T_stop_UNIX,inferred$return_time, units = "secs"))

inferred$PERIOD <- NA
# Assign Pre and Post labels
# The time between the nurses sampling from the nest and the return time is labeled as EXPOSURE_GAP (3h window) 
for (ROW in 1:nrow(inferred)) {
  #since these are times in seconds, it is ok to  have strict ">" instead of ">="
if(inferred[ROW,"time_stop_since_treat"] > 0){ inferred[ROW,"PERIOD"] <- "post"
}else if ( inferred[ROW,"time_stop_since_treat"] < (-3*3600) & inferred[ROW,"time_stop_since_treat"] > (-27*3600)) {  inferred[ROW,"PERIOD"] <- "pre"  } else{  inferred[ROW,"PERIOD"] <- "EXPOSURE_GAP"}
}

### Given that this operation requires the HD to be plugged in, save the output
write.table(inferred,file=paste(WORKDIR,"/Data/inferred_groomings_ALL_withCommonStart.txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")
write.table(Reps_N_exposed,file=paste(WORKDIR,"/Data/N_ants_exposed_xREP.txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")
}  
  
#######################################################
#start from here

###### AGGREGATE ALL VALUES FOR PRE.POST  #####
inferred <- read.table(paste(WORKDIR,"/Data/inferred_groomings_ALL_withCommonStart.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
Reps_N_exposed <- read.table(paste(WORKDIR,"/Data/N_ants_exposed_xREP.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
# Rename by name
inferred$TREATMENT <- as.factor(inferred$TREATMENT)
levels(inferred$TREATMENT)[levels(inferred$TREATMENT)=="BS"] <- "Big Sham"
levels(inferred$TREATMENT)[levels(inferred$TREATMENT)=="BP"] <- "Big Pathogen"
levels(inferred$TREATMENT)[levels(inferred$TREATMENT)=="SS"] <- "Small Sham"
levels(inferred$TREATMENT)[levels(inferred$TREATMENT)=="SP"] <- "Small Pathogen"

#SELECT RELEVANT ROWS
#Remove non-exposed reveivers and dead ants
inferred <- inferred[which(inferred$exposed=="exposed"),]
inferred <- inferred[which(inferred$dead=="no"),]
# #REMOVE EXP_GAP data
inferred <- inferred[which(inferred$PERIOD!="EXPOSURE_GAP"),]
#add count column
inferred$Count_byAnt <- 1

#N of reps (some missing!)
table(str_sub( unique(inferred$REP_treat),-2,-1))


## count the number of observations by ant and get mean
inferred_count_summary    <- aggregate(Count_byAnt ~ PERIOD + TREATMENT  + REP_treat + Rec_Name, FUN=length, na.action=na.pass, inferred)
## calculate mean durations for REP
inferred_dur_summary      <- aggregate(duration ~ PERIOD + TREATMENT  + REP_treat + Rec_Name, FUN=mean, na.rm=T, na.action=na.pass, inferred)
#sum by ant (SUM DUR MAY BE MORE INFORMATIVE AS THE GROOMING DETECTION MAY RESULT FAGMENTED)
inferred_SUM              <- aggregate(duration ~ PERIOD + TREATMENT  + REP_treat + Rec_Name, FUN=sum, na.rm=T, na.action=na.pass, inferred); names(inferred_SUM)[names(inferred_SUM) == 'duration'] <- 'SUM_duration'

#MERGE - REP_treatments
inferred_count_summ1 <- list(inferred_count_summary,inferred_dur_summary,inferred_SUM)
inferred_count_summ1 <- Reduce(function(x, y) merge(x, y, all=TRUE), inferred_count_summ1)

## merge counts & durations carefully
#inferred_count_summ1    <-  plyr::join(x=inferred_count_summary, y=inferred_dur_summary, type = "full", match = "all")

#calculate mean by REP
inferred_count_summ1    <- aggregate(cbind(Count_byAnt,duration,SUM_duration) ~ PERIOD + TREATMENT  + REP_treat, FUN=mean, na.action=na.pass, inferred_count_summ1)

## create a data frame with all combinations of the conditioning variables
all_combos1 <- expand.grid ( PERIOD=unique(inferred$PERIOD), REP_treat=unique(inferred$REP_treat))

## add the missing cases
Counts_by_Behaviour_AllCombos1 <- plyr::join (x = inferred_count_summ1 , y=all_combos1, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            

## replace the NAs with 0 counts   
Counts_by_Behaviour_AllCombos1$Count_byAnt[which(is.na(Counts_by_Behaviour_AllCombos1$Count_byAnt))] <- 0

#add N events and SD to calculate the Std Err of the mean differences
infer_Nevents  <- aggregate(Count_byAnt ~ PERIOD  + TREATMENT,                 FUN=length,      na.action=na.pass, Counts_by_Behaviour_AllCombos1); colnames(infer_Nevents) [match("Count_byAnt",colnames(infer_Nevents))] <- "N_Count_REP"
infer_SD  <- aggregate(cbind(Count_byAnt,duration) ~ PERIOD  + TREATMENT,                 FUN=sd,      na.rm=T, na.action=na.pass, Counts_by_Behaviour_AllCombos1); colnames(infer_SD) [match("Count_byAnt",colnames(infer_SD))] <- "SD_Count_byAnt" ; colnames(infer_SD) [match("duration",colnames(infer_SD))] <- "SD_duration"
## finally, get the mean & S.E. for each behav before/after  for barplots
infer_MEAN  <- aggregate(cbind(Count_byAnt,duration,SUM_duration)  ~ PERIOD  + TREATMENT,    FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1) 
infer_SE    <- aggregate(cbind(Count_byAnt,duration,SUM_duration) ~ PERIOD  + TREATMENT,     FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1) 
#rename cols
colnms <- c("Count_byAnt","duration","SUM_duration")
colnames(infer_MEAN)[match(colnms, colnames(infer_MEAN))] <- paste0("Mean_",colnms)
colnames(infer_SE)[match(colnms, colnames(infer_SE))] <- paste0("SE_",colnms)

#MERGE everything
infer_full <- list(infer_Nevents,infer_SD,infer_MEAN,infer_SE)
infer_full <- Reduce(function(x, y) merge(x, y, all=TRUE), infer_full)
infer_full[is.na(infer_full)] <- 0
#reorder levels to fix strange behaviour in plotting
infer_full$TREATMENT <- factor(infer_full$TREATMENT, levels = c("Big Pathogen","Big Sham","Small Pathogen","Small Sham"))

#----------------------------------
#ONGOING
# https://rpubs.com/brouwern/SEdiff2means

## diff after relative to before
#StdError_diff_between_means

#Formula for POOLED standard deviation (to calculate the SE of the difference between means)
## Note the formulas squares SD to get variance
var.pooled <- function(N1,N2,SD1,SD2){
  (N1*SD1^2 + N2*SD2^2)/(N1+N2)
}
# Standard error of difference
## Note that this uses sample size, NOT degrees of freedom (N)
SE.diff <- function(var.pool, n1,n2){
  sqrt(var.pool*(1/n1 + 1/n2))
}

#Apply function
## Note: uses sample size
#se.dif <- SE.diff(var.pool,n1 = N[1],n2 = N[2]) 


##### STD ERR SHOULD NOT BE TREATED LIKE THIS!!!!
#use the functions above to calculate the pooled SD and the N samples! 
infer_full_DELTA <- infer_full %>%
  group_by(TREATMENT) %>%
  dplyr::summarise(Mean_Count_PostPre = Mean_Count_byAnt[match("post", PERIOD)] - Mean_Count_byAnt[match("pre", PERIOD)],
                   Mean_dur_PostPre = Mean_duration[match("post", PERIOD)] - Mean_duration[match("pre", PERIOD)],
                   SE_Count_PostPre = SE_Count_byAnt[match("post", PERIOD)] - SE_Count_byAnt[match("pre", PERIOD)],
                   SE_dur_PostPre = SE_duration[match("post", PERIOD)] - SE_duration[match("pre", PERIOD)],
                   SUM.dur_PostPre = SUM_duration[match("post", PERIOD)] - SUM_duration[match("pre", PERIOD)],
  )
#---------------------------------------



###### AGGREGATE TIME BINS VALUES  #####

#time bins for plotting
time.break <- c("h4","hour") # ,"10min"     "h24",
tokeep <- NA

for (TIME in time.break) {
  #clean all
  rm(list=setdiff(ls(),c("inferred","Reps_N_exposed","infer_full","WORKDIR","DATADIR","TIME","time.break","tokeep",tokeep)))
  #create local copy
  inferred_bin <- inferred
  
  #add time.breaks to dataframe
  if (TIME=="hour") {
    inferred_bin$timespan <- inferred_bin$time_stop_since_treat/3600
  }else if(TIME=="10min"){ inferred_bin$timespan <- inferred_bin$time_stop_since_treat/600
  } else if (TIME=="h4") { 
    #create interval of 4h
    inferred_bin$timespan <- inferred_bin$time_stop_since_treat/14400
    }
  
  inferred_bin$timespan <- round(inferred_bin$timespan,0)
  
  if (TIME=="hour") {
  ## TEMP: assign a time of the day!
  Time_dictionary <- data.frame(timespan= -36:35, time_of_day= rep(0:23,3))
  inferred_bin <- left_join(inferred_bin, Time_dictionary, by = "timespan")
  } else   if (TIME=="h4") {
    Time_dictionary1 <- data.frame(timespan= c(-7:6), time_of_day= c(8, 12, 16, 20, 0, 4, 8,  12,  16,  20,  0,  4,  8, 12))
    inferred_bin <- left_join(inferred_bin, Time_dictionary1, by = "timespan")
  }
  
  
  ## count the number of observations by ant and get mean
  inferred_count_bin_summary    <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + time_of_day + timespan + REP_treat + Rec_Name, FUN=length, na.action=na.pass, inferred_bin)
  ## calculate mean durations for REP
  inferred_dur_bin_summary      <- aggregate(duration ~ PERIOD + TREATMENT + time_of_day + timespan + REP_treat + Rec_Name, FUN=mean, na.rm=T, na.action=na.pass, inferred_bin)
  #sum by ant (SUM DUR MAY BE MORE INFORMATIVE AS THE GROOMING DETECTION MAY RESULT FAGMENTED)
  inferred_bin_SUM              <- aggregate(duration ~ PERIOD + TREATMENT + time_of_day + timespan + REP_treat + Rec_Name, FUN=sum, na.rm=T, na.action=na.pass, inferred_bin); names(inferred_bin_SUM)[names(inferred_bin_SUM) == 'duration'] <- 'SUM_duration'
  
  #MERGE - REP_treatments
  inferred_count_bin_summ1 <- list(inferred_count_bin_summary,inferred_dur_bin_summary,inferred_bin_SUM)
  inferred_count_bin_summ1 <- Reduce(function(x, y) merge(x, y, all=TRUE), inferred_count_bin_summ1)
  
  ## merge counts & durations carefully
  #inferred_count_bin_summ1    <-  plyr::join(x=inferred_count_bin_summary, y=inferred_dur_bin_summary, type = "full", match = "all")
  
  #calculate mean by REP
  inferred_count_bin_summ1    <- aggregate(cbind(Count_byAnt,duration,SUM_duration) ~ PERIOD + TREATMENT + time_of_day + timespan + REP_treat, FUN=mean, na.action=na.pass, inferred_count_bin_summ1)

  ## create a data frame with all combinations of the conditioning variables
  all_combos1 <- expand.grid (TREATMENT=unique(inferred_bin$TREATMENT) , timespan= unique(inferred_bin$timespan)) #), PERIOD=unique(inferred_bin$PERIOD)
  
  ## add the missing cases
  Counts_by_Behaviour_AllCombos1 <- plyr::join (x = inferred_count_bin_summ1 , y=all_combos1, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            

  ## replace the NAs with 0 counts   
  Counts_by_Behaviour_AllCombos1$Count_byAnt[which(is.na(Counts_by_Behaviour_AllCombos1$Count_byAnt))] <- 0
  
  #add N events and SD to calculate the Std Err of the mean differences
  infer_bin_Nevents  <- aggregate(Count_byAnt ~ PERIOD + time_of_day + timespan + TREATMENT,                 FUN=length,      na.action=na.pass, Counts_by_Behaviour_AllCombos1); colnames(infer_bin_Nevents) [match("Count_byAnt",colnames(infer_bin_Nevents))] <- "N_Count_REP"
  infer_bin_SD  <- aggregate(cbind(Count_byAnt,duration) ~ PERIOD + time_of_day + timespan + TREATMENT,                 FUN=sd,      na.rm=T, na.action=na.pass, Counts_by_Behaviour_AllCombos1); colnames(infer_bin_SD) [match("Count_byAnt",colnames(infer_bin_SD))] <- "SD_Count_byAnt" ; colnames(infer_bin_SD) [match("duration",colnames(infer_bin_SD))] <- "SD_duration"
  ## finally, get the mean & S.E. for each behav before/after  for barplots
  infer_bin_MEAN  <- aggregate(cbind(Count_byAnt,duration,SUM_duration)  ~ PERIOD + time_of_day + timespan + TREATMENT,    FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1) 
  infer_bin_SE    <- aggregate(cbind(Count_byAnt,duration,SUM_duration) ~ PERIOD + time_of_day + timespan + TREATMENT,     FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1) 
  #rename cols
  colnms <- c("Count_byAnt","duration","SUM_duration")
  colnames(infer_bin_MEAN)[match(colnms, colnames(infer_bin_MEAN))] <- paste0("Mean_",colnms)
  colnames(infer_bin_SE)[match(colnms, colnames(infer_bin_SE))] <- paste0("SE_",colnms)
  
  #MERGE everything
  infer_bin_full <- list(infer_bin_Nevents,infer_bin_SD,infer_bin_MEAN,infer_bin_SE)
  infer_bin_full <- Reduce(function(x, y) merge(x, y, all=TRUE), infer_bin_full)
  infer_bin_full[is.na(infer_bin_full)] <- 0
  #reorder levels to fix strange behaviour in plotting
  infer_bin_full$TREATMENT <- factor(infer_bin_full$TREATMENT, levels = c("Big Pathogen","Big Sham","Small Pathogen","Small Sham"))
  
  #save output
  if (TIME=="h4") {
    infer_bin_h4 <- infer_bin_full
    tokeep <- c("tokeep","infer_bin_h4")
    
  } else if(TIME=="hour") { 
    #add break for line plots
    GAP <- expand.grid(PERIOD= "pre", timespan=c(-2,-1),time_of_day=c(10,11), TREATMENT=c("Small Pathogen","Small Sham","Big Pathogen","Big Sham"), 
                       N_Count_REP = NA, SD_Count_byAnt = NA, SD_duration = NA, Mean_Count_byAnt=NA, Mean_duration=NA, Mean_SUM_duration = NA, SE_Count_byAnt=NA, SE_duration=NA, SE_SUM_duration=NA )
    infer_bin_1h <- rbind(infer_bin_full,GAP)
    tokeep <- c("tokeep","infer_bin_1h")
    
  }else if(TIME=="10min") {
    #add break  for line plots
    #redo proper breaks, it should be a repeat between -12 and -6 I guess
    #likely not needed anyway as it will be a barplot
    GAP <- expand.grid(PERIOD= "pre", timespan=c(-2*6,-1*6),time_of_day=NA, TREATMENT=c("Small Pathogen","Small Sham","Big Pathogen","Big Sham"), 
                       N_Count_REP = NA, SD_Count_byAnt = NA, SD_duration = NA, Mean_Count_byAnt=NA, Mean_duration=NA, Mean_SUM_duration = NA, SE_Count_byAnt=NA, SE_duration=NA, SE_SUM_duration=NA )
    infer_bin_10min<- rbind(infer_bin_full,GAP)
    tokeep <- c("tokeep","infer_bin_10min")
  }
  
  
} #TIME

#reorder levels
infer_bin_1h$PERIOD <- factor(infer_bin_1h$PERIOD, levels = c("pre","post"))
infer_bin_h4$PERIOD <- factor(infer_bin_h4$PERIOD, levels = c("pre","post"))

#TRIMMED DATA
infer_bin_1h_trim <- infer_bin_1h[which(infer_bin_1h$time_of_day>=12 & infer_bin_1h$time_of_day <=16),]
#remove anything after timespan = 4 to exclude next day!
infer_bin_1h_trim <- infer_bin_1h_trim[which(infer_bin_1h_trim$timespan <=4),]

#TRIMMED DATA
infer_bin_h4_trim <- infer_bin_h4[which(infer_bin_h4$time_of_day==12),]
#remove anything after timespan = 4 to exclude next day!
infer_bin_h4_trim <- infer_bin_h4_trim[which(infer_bin_h4_trim$timespan <=1),]


############################################
# + pathogen−induced changes relative to sham−induced changes



##correct and copy from infer_full
# ##### STD ERR SHOULD NOT BE TREATED LIKE THIS!!!!
# infer_bin_1h_trim_DELTA <- infer_bin_1h_trim %>%
#   group_by(TREATMENT, time_of_day) %>%
#   dplyr::summarise(Mean_Count_PostPre = Mean_Count_byAnt[match("post", PERIOD)] - Mean_Count_byAnt[match("pre", PERIOD)],
#                    Mean_dur_PostPre = Mean_duration[match("post", PERIOD)] - Mean_duration[match("pre", PERIOD)],
#                    SE_Count_PostPre = SE_Count_byAnt[match("post", PERIOD)] - SE_Count_byAnt[match("pre", PERIOD)],
#                    SE_dur_PostPre = SE_duration[match("post", PERIOD)] - SE_duration[match("pre", PERIOD)],
#                    SUM.dur_PostPre = SUM_duration[match("post", PERIOD)] - SUM_duration[match("pre", PERIOD)],
#                    )






############# PLOTS ################
#style

STYLE <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw(),
              scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
              )

STYLE_continous <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw()
)


STYLE_NOVIR <- list(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw(),
              scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)


### GROOMING LOCATION
#mean by ant
Groom_location <- aggregate(REP_treat ~ TREATMENT + PERIOD + ant1.zones + Rec_Name , FUN=length, na.action=na.pass, inferred);  colnames(Groom_location) [match("REP_treat",colnames(Groom_location))] <- "Count_byAnt"
Groom_location$ant1.zones <- str_replace(Groom_location$ant1.zones,"1","Nest Area")
Groom_location$ant1.zones <- str_replace(Groom_location$ant1.zones,"2","Foraging Area")
#mean by nest
Groom_location_MEAN    <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + ant1.zones , FUN=mean, na.rm=T, na.action=na.pass, Groom_location) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count_byAnt"
Groom_location_SE    <- aggregate(Count_byAnt ~ PERIOD + TREATMENT + ant1.zones , FUN=std.error, na.rm=T, na.action=na.pass, Groom_location) ; colnames(Groom_location_SE) [match("Count_byAnt",colnames(Groom_location_SE))] <- "SD_Count"
##JOIN
Groom_location    <-  plyr::join(x=Groom_location_MEAN, y=Groom_location_SE, type = "full", match = "all")
### Grooming Location
ggplot(Groom_location, aes(x=TREATMENT, y=Count_byAnt, fill= PERIOD))+
  #geom_bar(stat='identity',position = position_dodge())+
  geom_errorbar( aes(x=TREATMENT,ymin=Count_byAnt-SD_Count, ymax=Count_byAnt+SD_Count),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~ant1.zones) +
  labs(title= "Grooming Location",subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate")), y = "Mean Freq by ant")+
  STYLE_NOVIR




############################################################
############################################################
############################################################
#####  BARPLOTS FOR PRE-POST ##

## MEAN FREQUENCY
ggplot(infer_full, aes(x=PERIOD, y=Mean_Count_byAnt, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_Count_byAnt-SE_Count_byAnt, ymax=Mean_Count_byAnt+SE_Count_byAnt),position=position_dodge2(width=0.8, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.8, preserve = "single")) +
  #facet_wrap(~ PERIOD) + #, labeller = as_labeller(time_of_day,text.add)
  STYLE +
  labs(title= "Frequency of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate (full period) ")), y = "Mean Freq by ant")

## MEAN DURATION
ggplot(infer_full, aes(x=PERIOD, y=Mean_duration, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_duration-SE_duration, ymax=Mean_duration+SE_duration),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  #facet_wrap(~ time_of_day) + #, labeller = as_labeller(time_of_day,text.add)
  STYLE +
  labs(title= "Duration of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate (full period)")), y = "Mean duration (s) by ant")

## TOTAL DURATION
ggplot(infer_full, aes(x=PERIOD, y=Mean_SUM_duration, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_SUM_duration-SE_SUM_duration, ymax=Mean_SUM_duration+SE_SUM_duration),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  #facet_wrap(~ time_of_day) + #, labeller = as_labeller(time_of_day,text.add)
  STYLE +
  labs(title= "Duration of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Total", phantom(.)%+-%phantom(.), "SE by replicate (full period)")), y = "Total duration (s) by ant")

############################################################
############################################################
############################################################

############################################################
#####  BARPLOTS FOR 4H BINS ##

## MEAN FREQUENCY

ggplot(infer_bin_h4_trim, aes(x=PERIOD, y=Mean_Count_byAnt, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_Count_byAnt-SE_Count_byAnt, ymax=Mean_Count_byAnt+SE_Count_byAnt),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  STYLE +
  labs(title= "Frequency of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate in 4h BLOCK (since 12:00)")), y = "Mean Freq by ant")


## MEAN DURATION
ggplot(infer_bin_h4_trim, aes(x=PERIOD, y=Mean_duration, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_duration-SE_duration, ymax=Mean_duration+SE_duration),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  STYLE +
  labs(title= "Duration of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate in 4h BLOCK (since 12:00)")), y = "Mean duration (s) by ant")


## TOTAL DURATION
ggplot(infer_bin_h4_trim, aes(x=PERIOD, y=Mean_SUM_duration, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_SUM_duration-SE_SUM_duration, ymax=Mean_SUM_duration+SE_SUM_duration),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  STYLE +
  #geom_text("",line= 5) +
  labs(title= "Duration of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Total", phantom(.)%+-%phantom(.), "SD by replicate in 4h BLOCK (since 12:00)")), y = "Total duration (s) by ant")



############################################################
#####  BARPLOTS FOR 1H BINS ##

## MEAN FREQUENCY
#make clear that those are hours
text.add <-":00"
infer_bin_1h_trim$time_of_day <- paste0(infer_bin_1h_trim$time_of_day,text.add)

ggplot(infer_bin_1h_trim, aes(x=PERIOD, y=Mean_Count_byAnt, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_Count_byAnt-SE_Count_byAnt, ymax=Mean_Count_byAnt+SE_Count_byAnt),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~ time_of_day) + #, labeller = as_labeller(time_of_day,text.add)
  STYLE +
  labs(title= "Frequency of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate per HOUR")), y = "Mean Freq by ant")


## MEAN DURATION
ggplot(infer_bin_1h_trim, aes(x=PERIOD, y=Mean_duration, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_duration-SE_duration, ymax=Mean_duration+SE_duration),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~ time_of_day) + #, labeller = as_labeller(time_of_day,text.add)
  STYLE +
  labs(title= "Duration of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate per HOUR")), y = "Mean duration (s) by ant")


## TOTAL DURATION
ggplot(infer_bin_1h_trim, aes(x=PERIOD, y=Mean_SUM_duration, fill= TREATMENT))+
  geom_errorbar( aes(x=PERIOD,ymin=Mean_SUM_duration-SE_SUM_duration, ymax=Mean_SUM_duration+SE_SUM_duration),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~ time_of_day) + #, labeller = as_labeller(time_of_day,text.add)
  STYLE +
  #geom_text("",line= 5) +
  labs(title= "Duration of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Total", phantom(.)%+-%phantom(.), "SD by individual per HOUR")), y = "Total duration (s) by ant")



############################################################
#####  BARPLOTS FOR 1H BINS DELTA POST-PRE ##


ggplot(infer_bin_1h_trim_DELTA, aes(x=TREATMENT, y=Mean_Count_PostPre, fill= TREATMENT))+
  geom_errorbar( aes(x=TREATMENT,ymin=Mean_Count_PostPre-SE_Count_PostPre, ymax=Mean_Count_PostPre+SE_Count_PostPre),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~ time_of_day) + #, labeller = as_labeller(time_of_day,text.add)
  STYLE +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title= "Frequency of Grooming in Sham and Pathogen treated colonies",subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate per HOUR")), y = "delta mean after relative to before")
  

############################################################



##### LINE PLOTS FOR 1H BINS ##

#FREQUENCY
ggplot(infer_bin_1h,
       aes(x = timespan, y = Mean_Count_byAnt,group = TREATMENT,color = TREATMENT)) +
  geom_vline(xintercept = 0,color = "red")+
  geom_point(size=1) +
  geom_smooth(data=subset(infer_bin_1h, PERIOD=="pre"), method = "lm") + #, formula = y ~ x + I(x^2)
  geom_smooth(data=subset(infer_bin_1h, PERIOD=="post"), method = "lm") + #, formula = y ~ x + I(x^2)
  STYLE_continous +
  labs(title = "Frequency of Grooming in Sham and Pathogen treated colonies",
       subtitle = "mean by ant and mean by rep",
       x = "time from treatment", y = "Freq by Hour by Ant") #+
  #facet_wrap(~ PERIOD) #, labeller = as_labeller(time_of_day,text.add)

#MEAN DURATION
ggplot(infer_bin_1h,
       aes(x = timespan, y = Mean_duration,group = TREATMENT,color = TREATMENT)) +
  geom_vline(xintercept = 0,color = "red")+
  geom_point(size=1) +
  geom_smooth(data=subset(infer_bin_1h, PERIOD=="pre"), method = "lm") + #, formula = y ~ x + I(x^2)
  geom_smooth(data=subset(infer_bin_1h, PERIOD=="post"), method = "lm") + #, formula = y ~ x + I(x^2)
  STYLE_continous +
  labs(title = "Duration of Grooming in Sham and Pathogen treated colonies",
       subtitle = "mean by ant and mean by rep",
       x = "time from treatment", y = "Mean duration by Hour") #+
#facet_wrap(~ PERIOD) #, labeller = as_labeller(time_of_day,text.add)

#SUM DURATION
ggplot(infer_bin_1h,
       aes(x = timespan, y = Mean_SUM_duration,group = TREATMENT,color = TREATMENT)) +
  geom_vline(xintercept = 0,color = "red")+
  geom_point(size=1) +
  geom_smooth(data=subset(infer_bin_1h, PERIOD=="pre"), method = "lm") + #, formula = y ~ x + I(x^2)
  geom_smooth(data=subset(infer_bin_1h, PERIOD=="post"), method = "lm") + #, formula = y ~ x + I(x^2)
  STYLE_continous +
  labs(title = "Duration of Grooming in Sham and Pathogen treated colonies",
       subtitle = "sum by ant and mean by rep",
       x = "time from treatment", y = "Total duration (sec) by Hour") #+
#facet_wrap(~ PERIOD) #, labeller = as_labeller(time_of_day,text.add)



###COUNT N OF EXPOSED RECEIVERS
Reps_N_exposed$N_received <- NA
for (REP.TREAT in unique(Reps_N_exposed$REP_treat) ) {
  AntList <- as.numeric(unlist(strsplit(Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP.TREAT),"N_ants"],",")))
  GroomingRec <- inferred[which(inferred$REP_treat==REP.TREAT),"Rec_Name"]
  GroomingRec <- unique(gsub("ant_", "", GroomingRec))
  # How many exposed ants received grooming
  Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP.TREAT),"N_received"] <- length(intersect(unique(GroomingRec),AntList))
}


### N of exposed ants in the colony
Reps_N_exposed$treat <- str_sub(Reps_N_exposed$REP_treat,-2,-1)
Mean_ants_exp <- aggregate(N_received ~ treat, FUN=mean, na.rm=T, na.action=na.pass, Reps_N_exposed)
SD_ants_exp <- aggregate(N_received ~ treat, FUN=sd, na.rm=T, na.action=na.pass, Reps_N_exposed); colnames(SD_ants_exp) [match("N_received",colnames(SD_ants_exp))] <- "SD_received"
N.ants.exposed    <-  plyr::join(x=Mean_ants_exp, y=SD_ants_exp, type = "full", match = "all")
data.frame(treat=N.ants.exposed$treat, N_exposed=sprintf("%.2f \U00B1 %.2f",N.ants.exposed$N_received,N.ants.exposed$SD_received))





##########################################################################################
### EXTRA UNRELATED PLOTS FOR PRESENTATIONS ##############################################

if (Presentation_plots) {

#######################################################################
###### PLOTTING A SINGLE GROOMING INTERACTION #########################

# Load the relevant libraries - do this every time
library(plyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library("scales")
library("viridis")
library(RSBID) #samp;ing

Grooming_int <- read.table("/home/cf19810/Documents/presentation_plot_data/Act28_Rec29_frames12357-12586-ROW15_uniqueID115_Grooming",header=T,stringsAsFactors = F, sep=",")
Grooming_int <- Grooming_int %>% dplyr::select(-contains(c(".x",".y",".angle",".zone",".time_sec","tim-interval","dt_FRAME","time_interval")))

for (variable in names(Grooming_int)){
  if (variable!="frame") {
    Grooming_int[variable]<- as.numeric(rescale(Grooming_int[,variable]))
  }
}

Grooming_int <- Grooming_int[180:229,]
Grooming_int$frame <- 1:length(Grooming_int$frame)

#long format
Grooming_int_long <- melt(setDT(Grooming_int), id.vars = c("frame"), variable.name = "Vars")

## heatmap
ggplot(Grooming_int_long, aes(frame, Vars)) + geom_tile(aes(fill = value),colour = "white", na.rm = T) +
  scale_fill_viridis(option = "B") +  
  #guides(fill=guide_legend(title="Total Incidents")) +
  theme_bw() + theme_minimal() + 
  labs(title = "",
       x = "frame", y = "variable") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",plot.title.position = "plot")

#################################################
#### PLOTTING ALL HIT MISS DATA FOR A SUBSET ####

Train_HitMiss <- read.table("/home/cf19810/Documents/presentation_plot_data/Training_R9SP_HitMiss_Auto_Grooming.txt",header=T,stringsAsFactors = F, sep=",")
Train_HitMiss <- Train_HitMiss %>% select(-contains(c("REPLICATE","PERIOD","BEH","ROW","Name","unique_interaction_id","ant","frame","duration_sec","pair","disagreement")))

### get labels of Relief selected var
RELIEF_selected <- c("prop_time_undetected_REC","sum_moved_distance_px_ACT","prop_time_undetected_ACT","transfmean_speed_pxpersec_REC","transfmean_abs_Body_Rotation_REC","transfSD_abs_jerk_PxPerSec3_REC" ,"stDev_turnAngle_REC","transfmean_abs_accel_pxpersec2_REC","transfmean_abs_jerk_PxPerSec3_REC" ,"stDev_Body_Rotation_REC","transfmean_abs_ang_Velocity_Body_REC","stDev_ang_Velocity_Body_ACT","transfmean_abs_ang_Velocity_Body_ACT","transfmean_abs_Body_Rotation_ACT","transfSD_abs_Body_Rotation_REC","stDev_Body_Rotation_ACT","transfSD_abs_ang_Velocity_Body_ACT","chull_area_REC","sum_moved_distance_px_REC","transfSD_abs_Body_Rotation_ACT","transfmean_speed_pxpersec_ACT","chull_area_ACT","transfSD_abs_ang_Velocity_Movement_REC","transfmean_abs_accel_pxpersec2_ACT", "stDev_turnAngle_ACT","root_mean_square_deviation_px_ACT","root_mean_square_deviation_px_REC" )

#scaling for plotting
for (variable in names(Train_HitMiss)){
  if (variable!="Hit") {
    Train_HitMiss[variable]<- as.numeric(rescale(Train_HitMiss[,variable]))
  }
}

Train_HitMiss <- Train_HitMiss[complete.cases(Train_HitMiss), ]
Train_HitMiss$Hit <- as.factor(Train_HitMiss$Hit)
Train_HitMiss_sampled   <- SBC(Train_HitMiss,   "Hit") #Under-Sampling Based on Clustering (SBC)

#add colour highlight
a <- ifelse(names(Train_HitMiss) %in% RELIEF_selected, "black", "gray")
Alpha <- ifelse(names(Train_HitMiss) %in% RELIEF_selected, "1", "0")

#Train_HitMiss_sampled$unique_interaction_id <- 1:length(Train_HitMiss_sampled$unique_interaction_id)

#assign numeration by hit
Train_HitMiss_sampled <- Train_HitMiss_sampled %>% group_by(Hit) %>% dplyr::mutate(id = row_number())

#long format
Train_HitMiss_long <- melt(setDT(Train_HitMiss_sampled), id.vars = c("id","Hit"), variable.name = "Vars")

Alpha <- ifelse(Train_HitMiss_long$Vars %in% RELIEF_selected, "1", "0")
levels(Train_HitMiss_long$Hit) <- list("Non Grooming"="0","Grooming"="1")

## heatmap
ggplot(Train_HitMiss_long, aes(id, Vars)) + geom_tile(aes(fill = value),colour = "white", na.rm = T) +
  scale_fill_viridis() +  
  facet_wrap(~Hit) +
  #guides(fill=guide_legend(title="Total Incidents")) +
  theme_bw() + theme_minimal() + 
  labs(title = "Heatmap of movement variables during a Grooming interaction",
       x = "interaction N", y = "") +
  theme(strip.text.x = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title.position = "plot")

## heatmap2
ggplot(Train_HitMiss_long, aes(id, Vars)) + geom_tile(aes(fill = value,alpha=Alpha),colour = "white", na.rm = T) +
  scale_fill_viridis() +  
  facet_wrap(~Hit) +
  #guides(fill=guide_legend(title="Total Incidents")) +
  theme_bw() + theme_minimal() + 
  labs(title = "Heatmap of movement variables during a Grooming interaction",
       x = "interaction N", y = "") +
  theme(strip.text.x = element_text(size = 18),
        axis.text.y = element_text(colour = a),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title.position = "plot")

+
  scale_x_discrete(labels=c("0" = "Non Grooming", "1" = "Grooming"))

}



####################################################################################
################ SCRAPS ############################################################


# ## COUNTS
# Xpos <- barplot( Count ~ period , Counts_AUTO_MEAN, beside=T, xlab="", ylab=" ", ylim=c(0,30)
#                  ,main="Auto classified")
# 
# ##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
# segments(x0 = Xpos[1,], 
#          x1 = Xpos[1,], 
#          y0 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="pre"], 
#          y1 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="pre"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="pre"],
#          lwd=2)
# 
# segments(x0 = Xpos[2,], 
#          x1 = Xpos[2,], 
#          y0 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"], 
#          y1 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="post"],
#          lwd=2)
# # 
# # text(x = ((Xpos[1,]+Xpos[2,])/2),
# #      y = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="post"]+15,
# #      stars.pval(posthoc_FREQ_summary$p.value))
# mtext("comparison of detected grooming", line=0, side=3, outer=T, cex=1.5)
# 




# inferred$T_start_UNIX <- as.POSIXct(inferred$T_start_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
# inferred$T_stop_UNIX <- as.POSIXct(inferred$T_stop_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
# 
# #split REPS
# 
# inferred_R3SP <- inferred[which(inferred$REPLICATE=="R3SP"),]
# inferred_R9SP <- inferred[which(inferred$REPLICATE=="R9SP"),]
# 
# #bins of an hour
# inf_R3SP_bin <- table(cut(inferred_R3SP$T_start_UNIX, breaks="hour"))
# inf_R9SP_bin <- table(cut(inferred_R9SP$T_start_UNIX, breaks="hour"))
# 
# 
# barplot(inf_R9SP_bin)

# 
# end <- inferred_R3SP[which(inferred_R3SP$T_stop_UNIX==max(inferred_R3SP$T_stop_UNIX)),"T_stop_sec"]
# 
# myFrame <- as.data.frame(table(myTable))
# 
# ggplot(inf_R3SP_bin, aes(, Y)) + geom_point() + geom_vline(xintercept = as.Date("2020-07-01"))
# 
# 
# 
# #
# barplot(inf_R3SP_bin) + abline(v=end)
# 
# +  abline(v =  ymd_hms(max(inferred_R3SP$T_start_UNIX)-4*3600, tz="GMT"))
# 
# 
# + abline(v = as.POSIXct(strptime("2021-03-17 07:56:42 GMT", format="%Y-%m-%d %H:%M:%OS")))
# 
# + abline(v =  max(inferred_R3SP$T_start_UNIX)-24*3600)
# 
# 
# exp_R3SP_time <-  as.POSIXct( "2021-03-15 12:11:00 GMT"  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" ) #manual
# #exp_R3SP_time <-max(inferred_R3SP$T_start_UNIX)-24*3600 #more formally correct
# exp_R9SP_time <-  as.POSIXct( "2021-04-26 11:26:00 GMT"  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" ) #manual
# #exp_R9SP_time <-max(inferred_R9SP$T_start_UNIX)-24*3600 #more formally correct

# inferred_R3SP$PERIOD_new <- NA
# inferred_R3SP$PERIOD_new <- ifelse(inferred_R3SP$T_start_UNIX < exp_R3SP_time, "pre", "post")
# #n occurrences
# 

# 
# ##############
# # Load manual annotations
# 
# annotations <- read.csv(paste(DATADIR,"/annotations_TRAINING_DATASET.csv",sep = ""), sep = ",")
# #transform zulu time in GMT
# annotations$T_start_UNIX <- as.POSIXct(annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
# annotations$T_stop_UNIX  <- as.POSIXct(annotations$T_stop,  format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
# #assign time in sec to avoid issues on time management and matching
# 
# #SELECT ONLY the exposed nurses
# annotations_R3 <-  annotations[which(annotations$treatment_rep=="R3SP" & annotations$Receiver %in% c(5,17)),]
# annotations_R9 <-  annotations[which(annotations$treatment_rep=="R9SP" & annotations$Receiver %in% c(23,29,32)),]
# annotations <- rbind(annotations_R3,annotations_R9)
# 
# ## count the number of observations of each behaviour - WARNING; some behavs not observed e.g. before the treatment (period), so will need to account for that (next step)
# Counts_by_Behaviour_CLEAN    <- aggregate(Actor ~ Behaviour + period + treatment_rep, FUN=length, na.action=na.pass, annotations); colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
# ## calculate mean durations for each behaviour
# Durations_by_Behaviour_CLEAN <- aggregate(duration ~ Behaviour + period + treatment_rep, FUN=mean, na.rm=T, na.action=na.pass, annotations)
# ## merge counts & durations carefully
# Counts_by_Behaviour_CLEAN    <-  plyr::join(x=Counts_by_Behaviour_CLEAN, y=Durations_by_Behaviour_CLEAN, type = "full", match = "all")
# 
# ## create a data frame with all combinations of the conditioning variables
# all_combos <- expand.grid ( Behaviour=unique(annotations$Behaviour), period=unique(annotations$period), treatment_rep=unique(annotations$treatment_rep))
# 
# ## add the missing cases
# Counts_by_Behaviour_AllCombos <- plyr::join (x = Counts_by_Behaviour_CLEAN , y=all_combos, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            
# 
# ## Focus only on a few important behaviours
# Counts_by_Behaviour_AllCombos$Behaviour <- as.character(Counts_by_Behaviour_AllCombos$Behaviour)  ## naughty R
# Counts_by_Behaviour_AllCombos <- Counts_by_Behaviour_AllCombos[which(Counts_by_Behaviour_AllCombos$Behaviour %in% c("G")),]
# 
# ## replace the NAs with 0 counts            
# Counts_by_Behaviour_AllCombos$Count[which(is.na(Counts_by_Behaviour_AllCombos$Count))] <- 0
# ## finally, get the mean & S.E. for each behav before/after  for barplots
# Counts_by_Behaviour_MEAN  <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)
# Counts_by_Behaviour_SE    <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)
# 
# 
# 
# ## show the mean counts for each behav | stage
# pdf(file=paste(DATADIR,"Grooming_Auto_Man__pre-post.pdf", sep = ""), width=5, height=8)
# par(mfrow=c(1,2), family="serif", mai=c(0.4,0.5,0.5,0.1), mgp=c(1.3,0.3,0), tcl=-0.2,oma=c(0,0,2,0))
# 
# ## COUNTS
# Counts_by_Behaviour_MEAN$period <- factor(Counts_by_Behaviour_MEAN$period , levels = c("pre","post"))
# 
# Xpos <- barplot( Count ~ period , Counts_by_Behaviour_MEAN, beside=T, xlab="", ylab="Behaviour count" , ylim=c(0,30)
#                  ,main="Manual annotation")
# ##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
# segments(x0 = Xpos[1,], 
#          x1 = Xpos[1,], 
#          y0 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="pre"], 
#          y1 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="pre"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="pre"],
#          lwd=2)
# 
# segments(x0 = Xpos[2,], 
#          x1 = Xpos[2,], 
#          y0 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"], 
#          y1 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="post"],
#          lwd=2)
# # 
# # text(x = ((Xpos[1,]+Xpos[2,])/2),
# #      y = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="post"]+15,
# #      stars.pval(posthoc_FREQ_summary$p.value))
# 
# 
# 
# ######### 
