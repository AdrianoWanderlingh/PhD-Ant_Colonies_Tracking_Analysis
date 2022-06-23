library(FortMyrmidon)
library("ggplot2")


## TIME WINDOW SHIFT. WARNING: THIS IS AN APPROXIMATION. IN THE FUTURE, THE TIME OF EXP ANTS RETURN PER TRACKING SYSTEM SHOULD BE USED! 
window_shift <- 60*20 #approx N of minutes that where given at the end as leeway, minutes can be skipped because of the "end of exp disruption" and because this causes an offset in the PERIOD transition
#window_shift_UNIX <- as.POSIXct(window_shift,origin = "1970-01-01",tz = "GMT")


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

## DIRECTORIES
WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data"


REP1to7_SmallOnly <- read.table(paste(WORKDIR,"/Data/inferred_groomings_REP1to7_SmallOnly.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
REP8to14_SmallOnly <- read.table(paste(WORKDIR,"/Data/inferred_groomings_REP8to14_SmallOnly.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

inferred <- rbind(REP1to7_SmallOnly,REP8to14_SmallOnly)

inferred$T_start_UNIX <- as.POSIXct(inferred$T_start_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )
inferred$T_stop_UNIX <- as.POSIXct(inferred$T_stop_UNIX,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" )



#base file info
#TREATMENT
inferred$TREATMENT <- substr(inferred$REPLICATE,(nchar(inferred$REPLICATE)+1)-2,nchar(inferred$REPLICATE))
#get rep-treat to match with exp filenames
inferred$REP_treat <- sub(".*\\/", "", inferred$REPLICATE)

#
#inferred$return_time <- NA
inferred$end_time <- NA
#inferred$time_since_start <- NA

### GET EXP END TIME FOR EACH REP AND ASSIGN PRE-POST TIME!

#list subdirectories in parent folder EXPERIMENT_DATA
files_list <- list.dirs.depth.n("/media/cf19810/DISK3/ADRIANO/EXPERIMENT_DATA", n = 1)
#select REP folders
files_list <- files_list[grep("REP",files_list)]

# replicate folder
for (REP.n in 1:length(files_list)) {
  # REP.n <- 1    #temp
  REP.folder      <- files_list[REP.n]
  REP.files       <- list.files(REP.folder, pattern = "CapDef3.myr")
  REP.filefolder  <- paste(REP.folder,REP.files,sep="/")
  
  #replicate file
  for (REP.FILES in REP.filefolder) {
    
    #get substring in variable until R9BS_
    REP_treat <- sub("\\_.*", "", basename(REP.FILES))
    #treat_name = substr(REP_treat_name,(nchar(REP_treat_name)+1)-2,nchar(REP_treat_name))

    # REP.FILES <-  REP.filefolder[1]   #temp
    print(REP.FILES) ##}}
    #open experiment
    exp <- fmExperimentOpen(REP.FILES)
    # exp.Ants <- exp$ants
    exp_end <- fmQueryGetDataInformations(exp)$end
    
    
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


## add a common time to all rows 
# Negative before exposure, positive after
inferred$time_stop_since_treat <- as.numeric(difftime(inferred$T_stop_UNIX,inferred$return_time, units = "secs"))


inferred$PERIOD <- NA
# Assign Pre and Post labels
for (ROW in 1:nrow(inferred)) {
if(inferred[ROW,"time_stop_since_treat"] > 0){ inferred[ROW,"PERIOD"] <- "post"
}else if ( inferred[ROW,"time_stop_since_treat"] < (-3*3600) & inferred[ROW,"time_stop_since_treat"] > (-27*3600)) {  inferred[ROW,"PERIOD"] <- "pre"  } else{  inferred[ROW,"PERIOD"] <- "EXPOSURE_GAP"}

}


###### start plottings

ggplot(inferred,aes(x=time_stop_since_treat,y=REP_treat,col=PERIOD))+geom_point()

table(inferred$PERIOD,inferred$TREATMENT)
##### THE VALUES PRESENT IN EXPOSURE GAP ARE CAUSED BY MIS-AGLINMENT OF TIME? GIVE THAT YOU OPEN FILES, 
# GET EXPOSURE TIME FROM WITH_METADATA AS THE LOWEST NUMBER OF AN EXPOSED ANT OF THE COLONY!!!
xxx <- as.data.frame(table(inferred$PERIOD,inferred$TREATMENT))

## FAI PLOT DECENTE USANDO LA BASE DI TOM!!!!!!!!!!!
ggplot(xxx,aes(y=Freq,col=Var1)) + geom_bar() + facet_wrap(xxx$Var2)


### CALCUALATE MEAN FREQUENCY BY ANT!!!! COLONIES MAY HAVE LARGER OR SMALLER N OF EXPOSED ANTS
### ALSO: IF ANT IS DEAD, EXCLUDE HER (SHE MAY BIAS THE COUNT!!!!!!!!!)
### AGGREGATE DATA PER PATHOGEN AND SHAM

























######################################################
inferred_cut$Count <- 1

###################################################
# TimeDiff <- difftime(exp_end, To, units = "hours")
# if(TimeDiff < 24){ PERIOD <- "POST"
# }else if ( TimeDiff >= 27 & TimeDiff < 51) { PERIOD <- "PRE"  } else{ PERIOD <- "EXPOSURE_GAP"}


#calculate mean by group_by
#sum pre post freq
inferred_by_Rep_Per    <- aggregate(Count ~ PERIOD_new + REPLICATE, FUN=length, na.action=na.omit, inferred_cut) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"


## create a data frame with all combinations of the conditioning variables
all_combos1 <- expand.grid ( PERIOD_new=unique(annotations$period), REPLICATE=unique(annotations$treatment_rep))

## add the missing cases
Counts_by_Behaviour_AllCombos1 <- plyr::join (x = inferred_by_Rep_Per , y=all_combos1, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            

## RENAME VAR
names(Counts_by_Behaviour_AllCombos1)[names(Counts_by_Behaviour_AllCombos1) == 'PERIOD_new'] <- 'period'


## replace the NAs with 0 counts            
Counts_by_Behaviour_AllCombos1$Count[which(is.na(Counts_by_Behaviour_AllCombos1$Count))] <- 0
## finally, get the mean & S.E. for each behav before/after  for barplots
Counts_AUTO_MEAN  <- aggregate(Count ~ period,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)
Counts_AUTO_SE    <- aggregate(Count ~ period,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos1)



# mean of the two reps 
#Counts_by_period  <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)

Counts_AUTO_MEAN$period  = factor(Counts_AUTO_MEAN$period, levels=c("pre", "post"))

## COUNTS
Xpos <- barplot( Count ~ period , Counts_AUTO_MEAN, beside=T, xlab="", ylab=" ", ylim=c(0,30)
                 ,main="Auto classified")

##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
segments(x0 = Xpos[1,], 
         x1 = Xpos[1,], 
         y0 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="pre"], 
         y1 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="pre"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="pre"],
         lwd=2)

segments(x0 = Xpos[2,], 
         x1 = Xpos[2,], 
         y0 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"], 
         y1 = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="post"],
         lwd=2)
# 
# text(x = ((Xpos[1,]+Xpos[2,])/2),
#      y = Counts_AUTO_MEAN$Count [Counts_AUTO_MEAN$period=="post"] + Counts_AUTO_SE$Count [Counts_AUTO_SE$period=="post"]+15,
#      stars.pval(posthoc_FREQ_summary$p.value))
mtext("comparison of detected grooming", line=0, side=3, outer=TRUE, cex=1.5)




# 
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

exp_R3SP_time <-  as.POSIXct( "2021-03-15 12:11:00 GMT"  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" ) #manual
#exp_R3SP_time <-max(inferred_R3SP$T_start_UNIX)-24*3600 #more formally correct
exp_R9SP_time <-  as.POSIXct( "2021-04-26 11:26:00 GMT"  ,format = "%Y-%m-%d %H:%M:%OS",  origin="1970-01-01", tz="GMT" ) #manual
#exp_R9SP_time <-max(inferred_R9SP$T_start_UNIX)-24*3600 #more formally correct

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
# #Get time boundaries (30  min window)



#########################################################################################################
#################### TO DOS: 

# ASSIGN 3H WINDOWS, PRE-POST LABEL BY SKIPPING THE 3h OF EXPOSURE TIME ##########

# PLOT PRE-POST AND SHAM PATHOGEN NEXT TO EACH OTHER 





# 
# 
# 
# for (HOUR in seq(from=0, to=48, by=3)){  ## increments by 3 hours for 48 hours
#   
#   #   ## increment the window start & end times by 1 hour
#   #HOUR <- 0 #TEMP!!! REACTIVAT HOUR LOOP
#   From  <- fmQueryGetDataInformations(exp)$end - 51*3600 + (HOUR * TimeWind) - window_shift
#   To    <- fmQueryGetDataInformations(exp)$end - 48*3600 + (HOUR * TimeWind) - window_shift
#   print(paste0("computing hour ",HOUR))
#   print(paste("Time window, from", From, "to", To))
#   start <- fmTimeCreate(offset=From) #end minus 48 hours plus incremental time
#   end   <- fmTimeCreate(offset=To) #end minus 45 hours plus incremental time
#   
#   #base file information
#   #PERIOD
#   TimeDiff <- difftime(exp_end, To, units = "hours")
#   if(TimeDiff < 24){ PERIOD <- "POST"
#   }else if ( TimeDiff >= 27 & TimeDiff < 51) { PERIOD <- "PRE"  } else{ PERIOD <- "EXPOSURE_GAP"}
#   
#   Period_dt<- data.frame(From, To, PERIOD)
#   Period_dataframe <- rbind(Period_dataframe, Period_dt)
#   
#   # }      
#   # 
#   # table(Period_dataframe$PERIOD) # shall be equal!
#   
#   # RUN FOR PRE AND POST (skip the 3h exposure gap)
#   if (!Period_dt$PERIOD=="EXPOSURE_GAP") {
#     
#     #COMPUTE NETWORK
#     print(paste0("Computing networks and properties"))
#     Graph <- compute_G(exp = exp, start = start, end=end)
#     
#     # COMPUTE NETWORK PROPERTIES
#     Network_prop_hour <- NetProperties(graph=Graph)
#     
#     Network_properties <- rbind(Network_properties,Network_prop_hour)
#     
#   } }#REP LOOP











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


dev.off()
