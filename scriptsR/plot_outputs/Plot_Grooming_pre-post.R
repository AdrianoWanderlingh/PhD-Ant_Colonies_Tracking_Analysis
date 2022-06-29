library(FortMyrmidon)
library("ggplot2")
library(lubridate)
library(plotrix)
library(scales)
library(car)
library(lme4)
library(Hmisc)

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

inferred <- read.table(paste(WORKDIR,"/Data/inferred_groomings_ALL_withCommonStart.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
Reps_N_exposed <- read.table(paste(WORKDIR,"/Data/N_ants_exposed_xREP.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

inferred$Count <- 1


inferred$hour <- inferred$time_stop_since_treat/3600
inferred$hour <- round(inferred$hour,0)


#N of ants per nest
#inferred_Npairs <- inferred %>% group_by(REP_treat,hour) %>% mutate(N_ants= n_distinct(pair))

# #perform data binning on SECONDS variable
# N_bins <- 51
# inferred_1h <- inferred %>% dplyr::mutate(hour = ntile(time_stop_since_treat, n=N_bins))
# #REMOVE EXP_GAP data
inferred <- inferred[which(inferred$PERIOD!="EXPOSURE_GAP"),]

#inferred_1h_xAnt <- FIX!!!!!!!!!!!!!! 
 
#calculate mean by group
inferred_1h_summary    <- aggregate(Count ~ PERIOD + TREATMENT + hour + REP_treat, FUN=length, na.action=na.omit, inferred) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"

# need to normalise by n of ants!
###CHECK IF ANTS IN EXP ARE IN INFERRED LIST
Reps_N_exposed$N_received <- NA
for (REP.TREAT in unique(Reps_N_exposed$REP_treat) ) {
  AntList <- as.numeric(unlist(strsplit(Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP.TREAT),"N_ants"],",")))
  all_ants <- c(inferred[which(inferred$REP_treat==REP.TREAT),"ant1"], inferred[which(inferred$REP_treat==REP.TREAT),"ant2"])
  # How many exposed ants received grooming
  Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP.TREAT),"N_received"] <- length(intersect(unique(all_ants),AntList))
}


#divide by N of ants
inferred_1h_summary$Count_byAnt <- NA
for (REP.TREAT in unique(inferred_1h_summary$REP_treat)) {
  inferred_1h_summary[which(inferred_1h_summary$REP_treat==REP.TREAT),"Count_byAnt"] <- inferred_1h_summary[which(inferred_1h_summary$REP_treat==REP.TREAT),"Count"]/ Reps_N_exposed[which(Reps_N_exposed$REP_treat==REP.TREAT),"N_received"] 
  }


## create a data frame with all combinations of the conditioning variables
all_combos1 <- expand.grid (TREATMENT=unique(inferred_1h_summary$TREATMENT) , hour= unique(inferred_1h_summary$hour)) #), PERIOD=unique(inferred_1h$PERIOD)

## add the missing cases
Counts_by_Behaviour_AllCombos1 <- plyr::join (x = inferred_1h_summary , y=all_combos1, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )

## RENAME VAR
#names(Counts_by_Behaviour_AllCombos1)[names(Counts_by_Behaviour_AllCombos1) == 'PERIOD_new'] <- 'period'


## replace the NAs with 0 counts
Counts_by_Behaviour_AllCombos1$Count[which(is.na(Counts_by_Behaviour_AllCombos1$Count))] <- 0
## finally, get the mean & S.E. for each behav before/after  for barplots
infer_1h_Count_MEAN  <- aggregate(Count_byAnt ~ PERIOD + hour + TREATMENT,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)
infer_1h_Count_SE    <- aggregate(Count_byAnt ~ PERIOD + hour + TREATMENT,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)

#add break #not working
GAP <- data.frame(PERIOD= "pre", hour=c(-2,-1), TREATMENT="SP",Count_byAnt=NA)
GAP <- rbind(GAP,data.frame(PERIOD= "pre", hour=c(-2,-1), TREATMENT="SS",Count_byAnt=NA))

infer_1h_Count_MEAN <- rbind(infer_1h_Count_MEAN,GAP)

ggplot(infer_1h_Count_MEAN) + 
  aes(x = hour, y = Count_byAnt, group = TREATMENT,color = TREATMENT) + 
  geom_path()+
  geom_ribbon(aes(y = infer_1h_Count_MEAN$Count_byAnt, ymin = infer_1h_Count_MEAN$Count_byAnt - infer_1h_Count_SE$Count_byAnt, ymax = infer_1h_Count_MEAN$Count_byAnt + infer_1h_Count_SE$Count_byAnt, fill = TREATMENT), alpha = .2, colour=NA) +
  theme_bw() +  
  # theme(legend.key = element_blank()) +
  geom_vline(xintercept = 0,color = "red")+
  labs(title = "Frequency of Grooming in Sham and Pathogen treated colonies",
       x = "time from treatment", y = "Freq by Hour by Ant") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
 

#P1 + scale_colour_viridis_d()


### ALSO: IF ANT IS DEAD, EXCLUDE HER (SHE MAY BIAS THE COUNT!!!!!!!!!)



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

Grooming_int <- Grooming_int %>% select(-contains(c(".x",".y",".angle",".zone",".time_sec","tim-interval","dt_FRAME","time_interval")))


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
  scale_fill_viridis() +  
  #guides(fill=guide_legend(title="Total Incidents")) +
  theme_bw() + theme_minimal() + 
  labs(title = "Heatmap of movement variables during a Grooming interaction",
       x = "frame", y = "variable") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",plot.title.position = "plot")





####################################
#### PLOTTING ALL HIT MISS DATA FOR A SUBSET

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

