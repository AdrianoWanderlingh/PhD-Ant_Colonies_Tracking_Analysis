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


WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data"

#########################################################################################
######THIS PART OF THE SCRIPT HAS BEEN DEACTIVATED AS THE PRELIMINARY OUTPUT HAS BEEN SAVED ON A SPECIFIC THIRD FILE
####################################################################################

# 
# ## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED!
# window_shift <- 60*20 #approx N of minutes that where given at the end as leeway, minutes can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition
# #window_shift_UNIX <- as.POSIXct(window_shift,origin = "1970-01-01",tz = "GMT")
# 
# 
# #### FUNCTIONS
# #list files recursive up to a certain level (level defined by "n" parameter)
# list.dirs.depth.n <- function(p, n) {
#   res <- list.dirs(p, recursive = FALSE)
#   if (n > 1) {
#     add <- list.dirs.depth.n(res, n-1)
#     c(res, add)
#   } else {
#     res
#   }
# }
# 
# ## DIRECTORIES
# 
# 
# REP1to7_SmallOnly <- read.table(paste(WORKDIR,"/Data/inferred_groomings_REP1to7_SmallOnly.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
# REP8to14_SmallOnly <- read.table(paste(WORKDIR,"/Data/inferred_groomings_REP8to14_SmallOnly.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
# 
# #List all text files in the working directory
# filenames <- list.files(DATADIR,pattern = '\\whole_experiment.txt$')
# filesdir <- paste(DATADIR,filenames, sep="/")
# #Read every text file with header, skipping the 1st row. 
# #Keep only the 5th column after reading the data. 
# result <- lapply(filesdir, function(x) read.table(x,header=T,stringsAsFactors = F, sep=","))
# inferred_LargeCols <- do.call(rbind, result)
# 
# inferred <- rbind(REP1to7_SmallOnly,REP8to14_SmallOnly,inferred_LargeCols)
# 
# 
# inferred$T_start_UNIX <- as.POSIXct(inferred$T_start_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
# inferred$T_stop_UNIX <- as.POSIXct(inferred$T_stop_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
# 
# 
# 
# #base file info
# #TREATMENT
# inferred$TREATMENT <- substr(inferred$REPLICATE,(nchar(inferred$REPLICATE)+1)-2,nchar(inferred$REPLICATE))
# #get rep-treat to match with exp filenames
# inferred$REP_treat <- sub(".*\\/", "", inferred$REPLICATE)
# 
# #
# #inferred$return_time <- NA
# inferred$end_time <- NA
# #inferred$time_since_start <- NA
# 
# ### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!
# 
# #list subdirectories in parent folder EXPERIMENT_DATA
# files_list <- list.dirs.depth.n("/media/cf19810/DISK3/ADRIANO/EXPERIMENT_DATA", n = 1)
# #select REP folders
# files_list <- files_list[grep("REP",files_list)]
# 
# 
# Reps_N_exposed <- data.frame(REP_treat= unique(inferred$REP_treat), N_ants = NA)
# 
# # replicate folder
# for (REP.n in 1:length(files_list)) {
#   # REP.n <- 1    #temp
#   REP.folder      <- files_list[REP.n]
#   REP.files       <- list.files(REP.folder, pattern = "CapDef3.myr")
#   REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
# 
#   #replicate file
#   for (REP.FILES in REP.filefolder) {
# 
#     #get substring in variable until R9BS_
#     REP_treat <- sub("\\_.*", "", basename(REP.FILES))
#     #treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))
# 
#     # REP.FILES <-  REP.filefolder[1]   #temp
#     print(REP.FILES) ##}}
#     #open experiment
#     exp <- fmExperimentOpen(REP.FILES)
#     # exp.Ants <- exp$ants
#     exp_end <- fmQueryGetDataInformations(exp)$end
# 
#     ########## GET EXPOSED ANTS # AW 17June2022
#     e.Ants <- exp$ants
#     Exposed_list <- vector()
#     if (REP_treat %in% Reps_N_exposed$REP_treat) {
#     for (ant in e.Ants){
#       if (TRUE %in% ant$getValues("Exposed")[,"values"]) {
#         exposed <-ant$ID
#         Exposed_list <- c(Exposed_list, exposed) }
#     }
#     explist<- data.frame(Exposed_list,stringsAsFactors = F)
#     #Collapse
#     explist2 <- data.frame(val=paste0(explist$Exposed_list,collapse = ', '),stringsAsFactors = F)
# 
# 
#       Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP_treat) ,"N_ants"] <- explist2
#     }
# 
#     for (ROW in 1:nrow(inferred)) {
#       if (inferred[ROW,"REP_treat"]==REP_treat) {
#         #add end time minus
#         inferred[ROW,"end_time"] <-  as.POSIXct(exp_end ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
#         # inferred[ROW,"return_time"] <-  as.POSIXct( exp_end  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
#         #assign a new common time for analyses
#         #inferred[ROW,"time_since_start"] <-
#       }
#     }
#     rm(list=(c("exp"))) #remove experiment
#   }}
# 
# #subtract end window
# inferred$end_time <- inferred$end_time - window_shift
# ## force unix_time
# inferred$end_time <- as.POSIXct(inferred$end_time,  origin="1970-01-01", tz="GMT" )
# 
# inferred$return_time <- inferred$end_time - 24*3600
# 
# 
# ## add a common time to all rows
# # Negative before exposure, positive after
# inferred$time_stop_since_treat <- as.numeric(difftime(inferred$T_stop_UNIX,inferred$return_time, units = "secs"))
# 
# 
# inferred$PERIOD <- NA
# # Assign Pre and Post labels
# for (ROW in 1:nrow(inferred)) {
# if(inferred[ROW,"time_stop_since_treat"] > 0){ inferred[ROW,"PERIOD"] <- "post"
# }else if ( inferred[ROW,"time_stop_since_treat"] < (-3*3600) & inferred[ROW,"time_stop_since_treat"] > (-27*3600)) {  inferred[ROW,"PERIOD"] <- "pre"  } else{  inferred[ROW,"PERIOD"] <- "EXPOSURE_GAP"}
# 
# }
# 
# 
# ### Given that this operation requires the HD to be plugged in, save the output
# write.table(inferred,file=paste(WORKDIR,"/Data/inferred_groomings_ALL_withCommonStart.txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")
# write.table(Reps_N_exposed,file=paste(WORKDIR,"/Data/N_ants_exposed_xREP.txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")

#######################################################
#start from here


time.break <- c("hour","10min")
tokeep <- NA

for (TIME in time.break) {
  #clean all
  rm(list=setdiff(ls(),c("WORKDIR","DATADIR","TIME","time.break",tokeep)))
  #load
  inferred <- read.table(paste(WORKDIR,"/Data/inferred_groomings_ALL_withCommonStart.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
  Reps_N_exposed <- read.table(paste(WORKDIR,"/Data/N_ants_exposed_xREP.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
  
  #N of occurences!
  table(str_sub( unique(inferred$REP_treat),-2,-1))
  #add count column
  inferred$Count <- 1

if (TIME=="hour") {
  
  inferred$timespan <- inferred$time_stop_since_treat/3600
  
}else{ inferred$timespan <- inferred$time_stop_since_treat/600
}

  
  inferred$timespan <- round(inferred$timespan,0)
  




#N of ants per nest
#inferred_Npairs <- inferred %>% group_by(REP_treat,hour) %>% mutate(N_ants= n_distinct(pair))

# #perform data binning on SECONDS variable
# N_bins <- 51
# inferred_1h <- inferred %>% dplyr::mutate(hour = ntile(time_stop_since_treat, n=N_bins))
# #REMOVE EXP_GAP data
inferred <- inferred[which(inferred$PERIOD!="EXPOSURE_GAP"),]

#inferred_1h_xAnt <- FIX!!!!!!!!!!!!!! 
 

## calculate MEAN durations for timespan bin
inferred_dur_bin_summary    <- aggregate(duration ~ PERIOD + TREATMENT + timespan + REP_treat, FUN=mean, na.action=na.omit, inferred) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
## calculate SUM durations for timespan bin
##IT IS POSSIBLE THAT BEH DETECTION OUTPUT IS FRAGMENTED, SO THE TOTAL SUM COULD BE A MORE RELIABLE MEASURE
inferred_SUMdur_bin_summary    <- aggregate(duration ~ PERIOD + TREATMENT + timespan + REP_treat, FUN=sum, na.action=na.omit, inferred) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
names(inferred_SUMdur_bin_summary)[names(inferred_SUMdur_bin_summary) == 'duration'] <- 'SUM.duration'
#calculate mean by group
inferred_bin_summary    <- aggregate(Count ~ PERIOD + TREATMENT + timespan + REP_treat, FUN=length, na.action=na.omit, inferred) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
## merge counts & durations carefully
inferred_bin_summary    <-  plyr::join(x=inferred_bin_summary, y=inferred_dur_bin_summary, type = "full", match = "all")
inferred_bin_summary    <-  plyr::join(x=inferred_bin_summary, y=inferred_SUMdur_bin_summary, type = "full", match = "all")


##################
# need to normalise by n of ants and n of reps!!!!!!!!!!
###CHECK IF ANTS IN EXP ARE IN INFERRED LIST
Reps_N_exposed$N_received <- NA
for (REP.TREAT in unique(Reps_N_exposed$REP_treat) ) {
  AntList <- as.numeric(unlist(strsplit(Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP.TREAT),"N_ants"],",")))
  all_ants <- c(inferred[which(inferred$REP_treat==REP.TREAT),"ant1"], inferred[which(inferred$REP_treat==REP.TREAT),"ant2"])
  # How many exposed ants received grooming
  Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP.TREAT),"N_received"] <- length(intersect(unique(all_ants),AntList))
}


#divide by N of ants
inferred_bin_summary$Count_byAnt <- NA
for (REP.TREAT in unique(inferred_bin_summary$REP_treat)) {
  inferred_bin_summary[which(inferred_bin_summary$REP_treat==REP.TREAT),"Count_byAnt"] <- inferred_bin_summary[which(inferred_bin_summary$REP_treat==REP.TREAT),"Count"]/ Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP.TREAT),"N_received"] 
  }


#inferred_bin_summary$Count_per_ant <- inferred_bin_summary$Count



## create a data frame with all combinations of the conditioning variables
all_combos1 <- expand.grid (TREATMENT=unique(inferred_bin_summary$TREATMENT) , timespan= unique(inferred_bin_summary$timespan)) #), PERIOD=unique(inferred_bin$PERIOD)

## add the missing cases
Counts_by_Behaviour_AllCombos1 <- plyr::join (x = inferred_bin_summary , y=all_combos1, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )

## RENAME VAR
#names(Counts_by_Behaviour_AllCombos1)[names(Counts_by_Behaviour_AllCombos1) == 'PERIOD_new'] <- 'period'


## replace the NAs with 0 counts
Counts_by_Behaviour_AllCombos1$Count[which(is.na(Counts_by_Behaviour_AllCombos1$Count))] <- 0
## finally, get the mean & S.E. for each behav before/after  for barplots
infer_bin_Count_MEAN  <- aggregate(cbind(Count_byAnt,duration)  ~ PERIOD + timespan + TREATMENT,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)
infer_bin_Count_SE    <- aggregate(cbind(Count_byAnt,duration) ~ PERIOD + timespan + TREATMENT,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)

#get also the SUM of Durations
infer_bin_Dur_SUM  <- aggregate(SUM.duration  ~ PERIOD + timespan + TREATMENT,                 FUN=sum,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)
#normalise SUM.duration by the N of reps present per group

#assign REP.weight: mean/N.reps
REP.WEIGHT <- as.data.frame(mean(table(str_sub( unique(inferred$REP_treat),-2,-1)))/table(str_sub( unique(inferred$REP_treat),-2,-1)))

## add weights to dataframe
##add var to SUM.duration as SUM.dur.weighted
infer_bin_Dur_SUM$SUM.dur.weighted <- REP.WEIGHT$Freq[match(infer_bin_Dur_SUM$TREATMENT, REP.WEIGHT$Var1)] * infer_bin_Dur_SUM$SUM.duration


if (TIME=="hour") {
  #add break
  GAP <- expand.grid(PERIOD= "pre", timespan=c(-2,-1), TREATMENT=c("SP","SS","BP","BS"),Count_byAnt=NA, duration=NA)
  GAP2 <- expand.grid(PERIOD= "pre", timespan=c(-2,-1), TREATMENT=c("SP","SS","BP","BS"), SUM.duration=NA, SUM.dur.weighted=NA)
  
  infer_1h_Count_MEAN <- rbind(infer_bin_Count_MEAN,GAP)
  infer_1h_Count_SE <- rbind(infer_bin_Count_SE,GAP)
  infer_1h_Dur_SUM <-  rbind(infer_bin_Dur_SUM,GAP2)
  
  tokeep <- c("tokeep","infer_1h_Count_MEAN","infer_1h_Count_SE","infer_1h_Dur_SUM")
  
}else{
  #add break
  #redo proper breaks, it should be a repeat between -12 and -6 I guess
  #likely not needed anyway as it will be a barplot
  GAP <- expand.grid(PERIOD= "pre", timespan=c(-2*6,-1*6), TREATMENT=c("SP","SS","BP","BS"),Count_byAnt=NA, duration=NA)
  GAP2 <- expand.grid(PERIOD= "pre", timespan=c(-2*6,-1*6), TREATMENT=c("SP","SS","BP","BS"), SUM.duration=NA, SUM.dur.weighted=NA)
  
  infer_10min_Count_MEAN <- rbind(infer_bin_Count_MEAN,GAP)
  infer_10min_Count_SE <- rbind(infer_bin_Count_SE,GAP)
  infer_10min_Dur_SUM <-  rbind(infer_bin_Dur_SUM,GAP2)
  
  tokeep <- c("tokeep","infer_10min_Count_MEAN","infer_10min_Count_SE","infer_10min_Dur_SUM")
}

} #TIME
############### END LOOP HERE, SAVE THE TWO OBJ AND PLOT!


############# PLOTS ################
##### LINE PLOTS FOR 1H BINS ##

#style
STYLE <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw())

#FREQUENCY
ggplot(infer_1h_Count_MEAN) + 
  aes(x = timespan, y = Count_byAnt, group = TREATMENT,color = TREATMENT) + 
  geom_line(size=1)+
  STYLE +
  geom_ribbon(aes(y = infer_1h_Count_MEAN$Count_byAnt, ymin = infer_1h_Count_MEAN$Count_byAnt - infer_1h_Count_SE$Count_byAnt, ymax = infer_1h_Count_MEAN$Count_byAnt + infer_1h_Count_SE$Count_byAnt, fill = TREATMENT), alpha = .2, colour=NA) +
  # theme(legend.key = element_blank()) +
  geom_vline(xintercept = 0,color = "red")+
  labs(title = "Frequency of Grooming in Sham and Pathogen treated colonies",
       x = "time from treatment", y = "Freq by Hour by Ant") 
 
#MEAN DURATION
ggplot(infer_1h_Count_MEAN) + 
  aes(x = timespan, y = duration, group = TREATMENT,color = TREATMENT) + 
  geom_line(size=1)+
  STYLE +
  geom_ribbon(aes(y = infer_1h_Count_MEAN$duration, ymin = infer_1h_Count_MEAN$duration - infer_1h_Count_SE$duration, ymax = infer_1h_Count_MEAN$duration + infer_1h_Count_SE$duration, fill = TREATMENT), alpha = .2, colour=NA) +
  # theme(legend.key = element_blank()) +
  geom_vline(xintercept = 0,color = "red")+
  labs(title = "Duration of Grooming in Sham and Pathogen treated colonies",
       x = "time from treatment", y = "Mean duration by Hour")

#SUM DURATION
ggplot(infer_1h_Dur_SUM) + 
  aes(x = timespan, y = SUM.dur.weighted, group = TREATMENT,color = TREATMENT) + 
  geom_line(size=1)+
  STYLE +
  # theme(legend.key = element_blank()) +
  geom_vline(xintercept = 0,color = "red")+
  labs(title = "Duration of Grooming in Sham and Pathogen treated colonies",
       subtitle = "adjusted by N of reps",
       x = "time from treatment", y = "Total duration (sec) by Hour")


##### BAR PLOTS FOR 5 mins BINS ##

#FREQUENCY
ggplot(infer_10min_Count_MEAN) + 
  aes(x = timespan, y = Count_byAnt, group = TREATMENT,color = TREATMENT,fill=TREATMENT,color = TREATMENT) + 
  geom_col(size=0.5)+
  STYLE +
  # theme(legend.key = element_blank()) +
  geom_vline(xintercept = 0,color = "red")+
  labs(title = "Frequency of Grooming in Sham and Pathogen treated colonies",
       x = "time from treatment", y = "Freq by 10 min bin by Ant") 

#MEAN DURATION
ggplot(infer_10min_Count_MEAN) + 
  aes(x = timespan, y = duration, group = TREATMENT,color = TREATMENT) + 
  geom_col(size=0.5)+
  STYLE +
  # theme(legend.key = element_blank()) +
  geom_vline(xintercept = 0,color = "red")+
  labs(title = "Duration of Grooming in Sham and Pathogen treated colonies",
       x = "time from treatment", y = "Mean duration by 10 min bin")

#SUM DURATION
ggplot(infer_10min_Dur_SUM) + 
  aes(x = timespan, y = SUM.dur.weighted, group = TREATMENT,fill=TREATMENT ,color = TREATMENT ) + 
  geom_col(size=0.5)+
  STYLE +
  # theme(legend.key = element_blank()) +
  geom_vline(xintercept = 0,color = "red")+
  labs(title = "Duration of Grooming in Sham and Pathogen treated colonies",
       subtitle = "adjusted by N of reps",
       x = "time from treatment", y = "Total duration (sec) by 10 min bin")






### ALSO: IF ANT IS DEAD, EXCLUDE HER (SHE MAY BIAS THE COUNT!!!!!!!!!)



# ###################################################
# # TimeDiff <- difftime(exp_end, To, units = "hours")
# # if(TimeDiff < 24){ PERIOD <- "POST"
# # }else if ( TimeDiff >= 27 & TimeDiff < 51) { PERIOD <- "PRE"  } else{ PERIOD <- "EXPOSURE_GAP"}
# 
# 
# #calculate mean by group_by
# #sum pre post freq
# inferred_by_Rep_Per    <- aggregate(Count ~ PERIOD_new + REPLICATE, FUN=length, na.action=na.omit, inferred_cut) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
# 
# 
# ## create a data frame with all combinations of the conditioning variables
# all_combos1 <- expand.grid ( PERIOD_new=unique(annotations$period), REPLICATE=unique(annotations$treatment_rep))
# 
# ## add the missing cases
# Counts_by_Behaviour_AllCombos1 <- plyr::join (x = inferred_by_Rep_Per , y=all_combos1, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            
# 
# ## RENAME VAR
# names(Counts_by_Behaviour_AllCombos1)[names(Counts_by_Behaviour_AllCombos1) == 'PERIOD_new'] <- 'period'
# 
# 
# ## replace the NAs with 0 counts            
# Counts_by_Behaviour_AllCombos1$Count[which(is.na(Counts_by_Behaviour_AllCombos1$Count))] <- 0
# ## finally, get the mean & S.E. for each behav before/after  for barplots
# Counts_AUTO_MEAN  <- aggregate(Count ~ period,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)
# Counts_AUTO_SE    <- aggregate(Count ~ period,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)
# 
# 
# 
# # mean of the two reps 
# #Counts_by_period  <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)
# 
# Counts_AUTO_MEAN$period  = factor(Counts_AUTO_MEAN$period, levels=c("pre", "post"))
# 
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
# mtext("comparison of detected grooming", line=0, side=3, outer=TRUE, cex=1.5)
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
# Counts_by_Behaviour_CLEAN    <- aggregate(Actor ~ Behaviour + period + treatment_rep, FUN=length, na.action=na.omit, annotations); colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
# ## calculate mean durations for each behaviour
# Durations_by_Behaviour_CLEAN <- aggregate(duration ~ Behaviour + period + treatment_rep, FUN=mean, na.action=na.omit, annotations)
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







# 
# 
# T_start_R3SP <- min(annotations[which(annotations$treatment_rep=="R3SP"),"T_start_UNIX"])
# T_stop_R3SP <- max(annotations[which(annotations$treatment_rep=="R3SP"),"T_start_UNIX"])
# T_start_R9SP <- min(annotations[which(annotations$treatment_rep=="R9SP"),"T_start_UNIX"])
# T_stop_R9SP <- max(annotations[which(annotations$treatment_rep=="R9SP"),"T_start_UNIX"])
# 
# #split by Rep
# inferred_R3SP <- inferred[which(inferred$REPLICATE=="R3SP"),]
# inferred_R9SP <- inferred[which(inferred$REPLICATE=="R9SP"),]
# 
# #cut according to boundaries
# inferred_R3SP <- inferred_R3SP [ which(inferred_R3SP$T_start_UNIX >= T_start_R3SP & inferred_R3SP$T_stop_UNIX <= T_stop_R3SP),]
# inferred_R9SP <- inferred_R9SP [ which(inferred_R9SP$T_start_UNIX >= T_start_R9SP & inferred_R9SP$T_stop_UNIX <= T_stop_R9SP),]
# 
# #recombine
# inferred_cut <- rbind(inferred_R3SP,inferred_R9SP)
# 
# 
# #PERIOD
# inferred_cut$PERIOD_new <- NA
# for (ROW in 1:nrow(inferred_cut)) {
# if(inferred_cut[ROW,"REPLICATE"] == "R3SP"){  inferred_cut[ROW,"PERIOD_new"] <- ifelse(inferred_cut[ROW,"T_start_UNIX"] < exp_R3SP_time, "pre", "post")
#                                                  }else if ( inferred_cut[ROW,"REPLICATE"] == "R9SP") { inferred_cut[ROW,"PERIOD_new"] <- ifelse(inferred_cut[ROW,"T_start_UNIX"] < exp_R9SP_time, "pre", "post") } else{ print("ERROR")}
#   
# }
# 


#dev.off()


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
ggplot(Grooming_int_long, aes(frame, Vars)) + geom_tile(aes(fill = value),colour = "white", na.rm = TRUE) +
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



#balanced_sample = NULL
# 
# 
# for (C in unique(Train_HitMiss$Hit)) {
#   tmp_df = Train_HitMiss%>%filter(Hit==C)
#   tmp<-ovun.sample(Click ~ ., data = tmp_df, method = "under", p = 0.5, seed = 5)$data
#   balanced_sample<-rbind(balanced_sample, tmp)
# }


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
ggplot(Train_HitMiss_long, aes(id, Vars)) + geom_tile(aes(fill = value),colour = "white", na.rm = TRUE) +
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
ggplot(Train_HitMiss_long, aes(id, Vars)) + geom_tile(aes(fill = value,alpha=Alpha),colour = "white", na.rm = TRUE) +
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

