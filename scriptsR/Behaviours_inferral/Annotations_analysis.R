# USER            <- "Adriano"
# if (USER=="Adriano") {WORKDIR <- "/home/cf19810/Dropbox/Ants_behaviour_analysis"}
# if (USER=="Tom")     {WORKDIR <- "/media/tom/MSATA/Dropbox/Ants_behaviour_analysis"}
# DATADIR <- paste(WORKDIR,"Data",sep="/")


####################### MOVE WHERE RELEVANT ###
#######LOAD TABLE WITH ANTS TASK #####
Task_list_R3SP <- read.table("/home/cf19810/Dropbox/Ants_behaviour_analysis/Data/R3SP_13-03-21.0000_TASK_list_1percent.txt")
Task_list_R3SP$treatment_rep <- "R3SP"
Task_ants_num_R3SP <- aggregate(antID ~ AntTask, FUN=NROW, Task_list_R3SP); colnames(Task_ants_num_R3SP) [match("antID",colnames(Task_ants_num_R3SP))] <- "number_ants"

ExposedAntsR3SP <- c(5,17)
ExposedAntsR9SP <- c(23,29,32)

#File specifing the exposed ants
antID <- c(5,17,23,29,32); Exposed <- c("exp","exp","exp","exp","exp"); treatment_rep <- c("R3SP","R3SP","R9SP","R9SP","R9SP"); ExposedAnts <- data.frame(antID,treatment_rep,Exposed)
# Check that exposed ants are nurses (they are as defined during experiments. If some turn out not to be then the Nurse_list Script has to be altered in T-start and T-end)
#ExposedAnts$AntTask  <- Task_list_R3SP$AntTask[which(Task_list_R3SP$antID %in% ExposedAnts$antID, )]

###############################################################################
###### LOAD MANUAL ANNOTATIONS ################################################
###############################################################################
#Behavioural codes explanation
# beh_codes <- read.csv("behavioural_codes.csv", sep= ",")

library(ggplot2)
library(plotrix)
library(gridExtra)
library(gtools) # to convert p.values in stars
library(pbkrtest)

annotations <- read.csv(paste(DATADIR,"/R3SP_R9SP_All_data_dropped_useless_cols.csv",sep = ""), sep = ",")
annotations$Behaviour <- as.character(annotations$Behaviour)
annotations$Actor <- as.character(annotations$Actor)
annotations$Receiver <- as.character(annotations$Receiver)
#call treatment as period
colnames(annotations)[which(names(annotations) == "treatment")] <- "period"


#convert Zulu time to GMT
annotations$T_start_UNIX <- as.POSIXct(annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
annotations$T_stop_UNIX  <- as.POSIXct(annotations$T_stop, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
annotations$duration <- as.numeric(annotations$T_stop_UNIX - annotations$T_start_UNIX)

#see if milliseconds are shown (number of decimals represented by the number after %OS)
format(annotations$T_start_UNIX[3], "%Y-%m-%d %H:%M:%OS6")

#remove duplicates of directed behaviours (Grooming and Aggression) by keeping only the behaviours where the Focal corresponds to the Actor.
#this seems to work very well with Grooming (cuts 50% of the events) and Agrgession (cuts 15 over 31 events) but affects also 4 Trophallaxis events, check why
annotations_drop_G_A <- annotations[which(annotations$Actor==annotations$Focal),]

#remove duplicates of non-directed behaviours as Trophalllaxis based on multiple columns check
annotations_drop_all_dup <- annotations_drop_G_A[!duplicated(annotations_drop_G_A[,c("T_start_UNIX","T_stop_UNIX","Behaviour","treatment_rep")]),]

#Summary counts of behaviour frequency
Counts_by_Behaviour              <- aggregate(treatment_rep ~ Behaviour + period, FUN=length, annotations); colnames(Counts_by_Behaviour) [match("treatment_rep",colnames(Counts_by_Behaviour))] <- "Count"
Counts_by_Behaviour_drop_G_A     <- aggregate(treatment_rep ~ Behaviour + period, FUN=length, annotations_drop_G_A); colnames(Counts_by_Behaviour_drop_G_A) [match("treatment_rep",colnames(Counts_by_Behaviour_drop_G_A))] <- "Count_drop_G_A"
Counts_by_Behaviour_drop_all_dup <- aggregate(treatment_rep ~ Behaviour + period, FUN=length, annotations_drop_all_dup); colnames(Counts_by_Behaviour_drop_all_dup) [match("treatment_rep",colnames(Counts_by_Behaviour_drop_all_dup))] <- "Count_drop_all_dup"
#check how many cases have been removed
Counts_by_Behaviour <- cbind(Counts_by_Behaviour, Count_drop_G_A = Counts_by_Behaviour_drop_G_A$Count_drop_G_A, Count_drop_all_dup = Counts_by_Behaviour_drop_all_dup$Count_drop_all_dup)
Counts_by_Behaviour

#see total final numer of behaviours
Counts_by_Behaviour_tots <- aggregate(Count_drop_all_dup ~ Behaviour, FUN=sum, Counts_by_Behaviour); colnames(Counts_by_Behaviour) [match("period",colnames(Counts_by_Behaviour))] <- "Totals"

#Over-write cleaned dataset - NOTE 'THIS'annotations' IS USED in the trajectory plotting loop below! 
annotations <- annotations_drop_all_dup

#SAVE THE NEW ANNOTATION FILE 
write.csv(annotations,"/home/cf19810/Dropbox/Ants_behaviour_analysis/Data//R3SP_R9SP_All_data_FINAL_script_output.csv")


## count the number of observations of each behaviour - WARNING; some behavs not observed e.g. before the treatment (period), so will need to account for that (next step)
Counts_by_Behaviour_CLEAN    <- aggregate(Actor ~ Behaviour + period + treatment_rep, FUN=length, na.action=na.omit, annotations); colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
## calculate mean durations for each behaviour
Durations_by_Behaviour_CLEAN <- aggregate(duration ~ Behaviour + period + treatment_rep, FUN=mean, na.action=na.omit, annotations)
## merge counts & durations carefully
Counts_by_Behaviour_CLEAN    <-  plyr::join(x=Counts_by_Behaviour_CLEAN, y=Durations_by_Behaviour_CLEAN, type = "full", match = "all")

## create a data frame with all combinations of the conditioning variables
all_combos <- expand.grid ( Behaviour=unique(annotations$Behaviour), period=unique(annotations$period), treatment_rep=unique(annotations$treatment_rep))

## add the missing cases
Counts_by_Behaviour_AllCombos <- plyr::join (x = Counts_by_Behaviour_CLEAN , y=all_combos, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            

## Focus only on a few important behaviours
Counts_by_Behaviour_AllCombos$Behaviour <- as.character(Counts_by_Behaviour_AllCombos$Behaviour)  ## naughty R
Counts_by_Behaviour_AllCombos <- Counts_by_Behaviour_AllCombos[which(Counts_by_Behaviour_AllCombos$Behaviour %in% c("A","G","SG","T","TB")),]

## replace the NAs with 0 counts            
Counts_by_Behaviour_AllCombos$Count[which(is.na(Counts_by_Behaviour_AllCombos$Count))] <- 0
## finally, get the mean & S.E. for each behav before/after  for barplots
Counts_by_Behaviour_MEAN  <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)
Counts_by_Behaviour_SE    <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)

############################################################################
###### STATISTICS ON COUNTS ################################################
############################################################################
library(lme4)
library(blmeco) # check dispersion for glmer
library(emmeans) #post-hoc comparisons
library(e1071) #calc skewness and other stuff
library(lawstat) #for levene test (homogeneity of variance)


#function to test normality of residuals
test_norm <- function(resids){
  print("Testing normality")
  if (length(resids)<=300){
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
  }else{
    print("More than 300 data points so using the skewness and kurtosis
approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}

str(Counts_by_Behaviour_AllCombos)
Counts_by_Behaviour_AllCombos$treatment_rep <- as.factor(Counts_by_Behaviour_AllCombos$treatment_rep)
Counts_by_Behaviour_AllCombos$period<- as.factor(Counts_by_Behaviour_AllCombos$period)
Counts_by_Behaviour_AllCombos$Behaviour <- as.factor(Counts_by_Behaviour_AllCombos$Behaviour)
Counts_by_Behaviour_AllCombos$Count <- as.numeric(Counts_by_Behaviour_AllCombos$Count)
str(Counts_by_Behaviour_AllCombos)

#Counts_by_Behaviour_AllCombos$period = relevel(Counts_by_Behaviour_AllCombos$period, ref="pre")
m1 <- lmerTest::lmer(sqrt(Count) ~ period * Behaviour * (1|treatment_rep), Counts_by_Behaviour_AllCombos)
summary(m1)
test_norm(residuals(m1)) #test residuals' normality. null hypothesis for the Shapiro-Wilk test is that a variable is normally distributed
qqnorm(residuals(m1))

posthoc_FREQ <- emmeans(m1, specs = trt.vs.ctrlk ~ period | Behaviour)
posthoc_FREQ_summary<-summary(posthoc_FREQ$contrasts)

par(mfrow=c(1,2))
plot(m1)
qqnorm(residuals(m1))
qqline(residuals(m1))
hist(residuals(m1)) 

anova(m1)
##or
#emmeans(m5, list(pairwise ~ period | Behaviour), adjust = "tukey")
##or
# contrast(emmeans(m5,~period|Behaviour, type="response"),
#          method="pairwise", adjust="Tukey")  #method = "trt.vs.ctrl"


# m2 <- glmer(Count ~ period * Behaviour + (1|treatment_rep), data = Counts_by_Behaviour_AllCombos, family = "poisson")
# summary(m2)
# test_norm(residuals(m2))
# qqnorm(residuals(m2))

##check poisson distribution overdispersion
##if the ratio of residual deviance to degrees of freedom it's >1 then the data are overdispersed
#dispersion_glmer(m2)

# #overdispersion can be incorporated by including an observation-level random effect (1|obs), to allow observations to be correlated with themselves
# Counts_by_Behaviour_AllCombos$obs=factor(1:nrow(Counts_by_Behaviour_AllCombos))
# m5 <- glmer(Count ~ period + Behaviour + (1|treatment_rep) +(1|obs), data = Counts_by_Behaviour_AllCombos, family = "poisson")
# summary(m5) # is singular fit a problem?
# dispersion_glmer(m5)#model m5, which assumes poisson distribution and accounts for overdispersion by adding an observation-level random effect shows both lower AIC and dispersion ~1
# anova(m2,m5)


######OTHER STUFF################
# #looking into negative binomial
# m6 <- glmer.nb(Count ~ period + Behaviour + (1|treatment_rep), data = Counts_by_Behaviour_AllCombos)
# #m6<- glm.nb(Count ~ period + Behaviour, data = Counts_by_Behaviour_AllCombos)
# m7 <- glmer(Count ~ treatment + Behaviour + (1|treatment_rep), data = Counts_by_Behaviour_AllCombos, family=binomial(link = "logit"))

# #we can use a quasi-family to estimate the dispersion parameter (Bbut without random effect in the model)
# summary(m4) #model with Quasipoisson distribution

#related packages
# #test normality of residuals
# library(e1071); library(lmerTest); library(MuMIn); library(arm); library(lmtest); library(MASS); library(pbkrtest)

############################################################################
###### STATISTICS ON DURATION ##############################################
############################################################################
annotations$treatment_rep <- as.factor(annotations$treatment_rep)
annotations$period<- as.factor(annotations$period)
annotations$Behaviour <- as.character(annotations$Behaviour)
annotations$Focal <- as.character(annotations$Focal)
str(annotations)

## add the missing cases to avoid A to not be estimated
annotations_allCombs <- plyr::join (x = annotations , y=all_combos, type = "right", match = "all")


annotations_allCombs <- annotations_allCombs[which(annotations_allCombs$Behaviour %in% c("A","G","SG","T","TB")),]
annotations_allCombs$Behaviour <- as.factor(annotations_allCombs$Behaviour)

m_dur1 <- lmerTest::lmer(duration ~ period + Behaviour + (1|treatment_rep), annotations_allCombs)
summary(m_dur1)
test_norm(residuals(m_dur1)) #skewed data

#Log transform data
annotations_allCombs[annotations_allCombs$Behaviour=="A"& annotations_allCombs$period=="pre", ]$duration <- 0.0001
annotations_allCombs <- annotations_allCombs[!annotations_allCombs$duration==0, ]
annotations_allCombs$duration_log <- log(annotations_allCombs$duration)

m_dur2 <- lmerTest::lmer(duration_log ~ period * Behaviour + (1|treatment_rep), annotations_allCombs)
summary(m_dur2)
test_norm(residuals(m_dur2)) #good values of kurtosis and skeweness

plot(m_dur2)
par(mfrow=c(1,2))
qqnorm(resid(m_dur2)); qqline(resid(m_dur2)); hist(resid(m_dur2))

## Homogeneity of variance
#if test is ns the variances are homogenous
levene.test(annotations_allCombs$duration_log, annotations_allCombs$period)

anova(m_dur2) #I guess it shows that there are differences between behaviours but not between period

posthoc_DUR <- emmeans(m_dur2, specs = trt.vs.ctrlk ~ period | Behaviour)
posthoc_DUR_summary<-summary(posthoc_DUR$contrast)

###############################################################################
###### PLOTTING  ##############################################################
###############################################################################

## show the mean counts for each behav | stage
pdf(file=paste(DATADIR,"Behaviour_counts_pre-post.pdf", sep = ""), width=5, height=3)
par(mfrow=c(1,2), family="serif", mai=c(0.4,0.5,0.1,0.1), mgp=c(1.3,0.3,0), tcl=-0.2)
## COUNTS
Xpos <- barplot( Count ~ period + Behaviour , Counts_by_Behaviour_MEAN, beside=T, xlab="", ylab="Behaviour count", ylim=c(0,max( Counts_by_Behaviour_MEAN$Count +  Counts_by_Behaviour_SE$Count*1.4 )))
##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
segments(x0 = Xpos[1,], 
         x1 = Xpos[1,], 
         y0 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="pre"], 
         y1 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="pre"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="pre"],
         lwd=2)

segments(x0 = Xpos[2,], 
         x1 = Xpos[2,], 
         y0 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"], 
         y1 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="post"],
         lwd=2)

text(x = ((Xpos[1,]+Xpos[2,])/2),
     y = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$period=="post"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$period=="post"]+15,
     stars.pval(posthoc_FREQ_summary$p.value))

## DURATIONS
Xpos <- barplot( duration ~ period + Behaviour , Counts_by_Behaviour_MEAN, beside=T, xlab="", ylab="Behaviour duration (s)", ylim=c(0,max( Counts_by_Behaviour_MEAN$duration +  Counts_by_Behaviour_SE$duration, na.rm=T)))
##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
segments(x0 = Xpos[1,], 
         x1 = Xpos[1,], 
         y0 = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$period=="pre"], 
         y1 = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$period=="pre"] + Counts_by_Behaviour_SE$duration [Counts_by_Behaviour_SE$period=="pre"],
         lwd=2)
segments(x0 = Xpos[2,], 
         x1 = Xpos[2,], 
         y0 = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$period=="post"], 
         y1 = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$period=="post"] + Counts_by_Behaviour_SE$duration [Counts_by_Behaviour_SE$period=="post"],
         lwd=2)

text(x = ((Xpos[1,]+Xpos[2,])/2),
     y = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$period=="post"] + Counts_by_Behaviour_SE$duration [Counts_by_Behaviour_SE$period=="post"]+3,
     stars.pval(rev(posthoc_DUR_summary$p.value)))

## COUNT BY ACTOR and RECEIVER
## count the number of observations per actor/receiver and behaviour to find the most interacting indivduals

#Change colname of both Actor and Receiver in the dataframes to "ant" for the loop to work
Counts_by_Behaviour_Rec <- aggregate(Focal ~ Behaviour + period + Receiver + treatment_rep, FUN=length, na.action=na.omit, annotations); colnames(Counts_by_Behaviour_Rec) [match("Focal",colnames(Counts_by_Behaviour_Rec))] <- "Count"; colnames(Counts_by_Behaviour_Rec) [match("Receiver",colnames(Counts_by_Behaviour_Rec))] <- "antID"
Counts_by_Behaviour_Act <- aggregate(Focal ~ Behaviour + period + Actor + treatment_rep, FUN=length, na.action=na.omit, annotations); colnames(Counts_by_Behaviour_Act) [match("Focal",colnames(Counts_by_Behaviour_Act))] <- "Count"; colnames(Counts_by_Behaviour_Act) [match("Actor",colnames(Counts_by_Behaviour_Act))] <- "antID"

Durations_by_Behaviour_Rec <- aggregate(duration ~ Behaviour + period + Receiver + treatment_rep, FUN=mean, na.action=na.omit, annotations); colnames(Durations_by_Behaviour_Rec) [match("Focal",colnames(Durations_by_Behaviour_Rec))] <- "Count"; colnames(Durations_by_Behaviour_Rec) [match("Receiver",colnames(Durations_by_Behaviour_Rec))] <- "antID"
Durations_by_Behaviour_Act <- aggregate(duration ~ Behaviour + period + Actor + treatment_rep, FUN=mean, na.action=na.omit, annotations); colnames(Durations_by_Behaviour_Act) [match("Focal",colnames(Durations_by_Behaviour_Act))] <- "Count"; colnames(Durations_by_Behaviour_Act) [match("Actor",colnames(Durations_by_Behaviour_Act))] <- "antID"

#create place-holder column to have name for plotting
Counts_by_Behaviour_Rec$Receiver <- "Receiver"
Counts_by_Behaviour_Act$Actor <- "Actor"
Rec_Act_Counts_list <- list(Counts_by_Behaviour_Rec=Counts_by_Behaviour_Rec,Counts_by_Behaviour_Act=Counts_by_Behaviour_Act)

Durations_by_Behaviour_Rec$Receiver <- "Receiver"
Durations_by_Behaviour_Act$Actor <- "Actor"

for (SUBJECT in Rec_Act_Counts_list) {
  #print(deparse(substitute(Rec_Act_Counts_list)[SUBJECT]))
  #cat(names(Rec_Act_Counts_list)[SUBJECT])
  ANT_role <- names(SUBJECT)[6]
  ## add the missing cases
  SUBJECT_AllCombos <- plyr::join (x = SUBJECT , y=all_combos, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )
  ## Focus only on a few important behaviours
  SUBJECT_AllCombos$Behaviour <- as.character(SUBJECT_AllCombos$Behaviour)  ## naughty R
  SUBJECT_AllCombos <- SUBJECT_AllCombos[which(SUBJECT_AllCombos$Behaviour %in% c("A","G","SG","T","TB")),]
  ## replace the NAs with 0 counts            
  SUBJECT_AllCombos$Count[which(is.na(SUBJECT_AllCombos$Count))] <- 0
  
  SUBJECT_AllCombos_R3SP <- subset(SUBJECT_AllCombos, treatment_rep == "R3SP")
  # one box per behaviour, showing Receivers and divided by colony (to avoid overlaps of antIDs)
  #colony R3SP Ant
  p <- ggplot(SUBJECT_AllCombos_R3SP, aes(x=period, y=Count, fill=period)) + 
    geom_bar(aes(fill = antID),stat='identity',colour="black", fill=NA) +
    geom_text(aes(label=ifelse(Count>6,as.character(antID),''),color=ifelse(antID %in% ExposedAntsR3SP,"red","black")), size = 2, position = position_stack(vjust = 0.9)) +
    scale_colour_manual(values = c("darkgray","black")) +
    facet_grid(~Behaviour, scale="free") +
    theme(text=element_text(family="serif",size=7),legend.position="none") +
    labs(title = paste("Behaviour frequency by",ANT_role, "\n by exposure period (colony:", "R3SP", ")"),
         subtitle = "stacked bar label represents the AntID. \n black label= pathogen treated ant",
         y ="Behaviour count")
  
  SUBJECT_AllCombos_R9SP <- subset(SUBJECT_AllCombos, treatment_rep == "R9SP")
  # one box per behaviour, showing Receivers and divided by colony (to avoid overlaps of antIDs)
  #colony R3SP Ant
  p1 <- ggplot(SUBJECT_AllCombos_R9SP, aes(x=period, y=Count, fill=period)) + 
    geom_bar(aes(fill = antID),stat='identity',colour="black", fill=NA) +    
    geom_text(aes(label=ifelse(Count>6,as.character(antID),''),color=ifelse(antID %in% ExposedAntsR9SP,"red","")), size = 2, position = position_stack(vjust = 0.9)) +
    scale_colour_manual(values = c("darkgray","black")) +
    facet_grid(~Behaviour, scale="free") +
    theme(text=element_text(family="serif",size=7),legend.position="none") +
    labs(title = paste("Behaviour frequency by",ANT_role,"\n by exposure period (colony:", "R9SP", ")"),
         subtitle = "stacked bar label represents the AntID. \n black label= pathogen treated ant",
         y ="Behaviour count")
  
  grid.arrange(p, p1,
               ncol = 2, nrow = 1)
  
}

## close the pdf
dev.off()


####MODIFY FOR COUNTS AND DURATIONS #########


# # this can be eliminated if not modified!!!!!!!!!!!
# #Change colname of both Actor and Receiver in the dataframes to "ant" for the loop to work
# Counts_by_Behaviour_Rec <- aggregate(Focal ~ Behaviour + period + Receiver + treatment_rep, FUN=length, na.action=na.omit, annotations); colnames(Counts_by_Behaviour_Rec) [match("Focal",colnames(Counts_by_Behaviour_Rec))] <- "Count"; colnames(Counts_by_Behaviour_Rec) [match("Receiver",colnames(Counts_by_Behaviour_Rec))] <- "antID"
# Counts_by_Behaviour_Act <- aggregate(Focal ~ Behaviour + period + Actor + treatment_rep, FUN=length, na.action=na.omit, annotations); colnames(Counts_by_Behaviour_Act) [match("Focal",colnames(Counts_by_Behaviour_Act))] <- "Count"; colnames(Counts_by_Behaviour_Act) [match("Actor",colnames(Counts_by_Behaviour_Act))] <- "antID"
# #create place-holder column to have name for plotting
# Counts_by_Behaviour_Rec$Receiver <- "Receiver"
# Counts_by_Behaviour_Act$Actor <- "Actor"

############################


## COUNT BY EXPOSED AND UNEXPOSED PER TASK ##
#rename Actor/Receiver columns to Role before binding
colnames(Counts_by_Behaviour_Rec) [match("Receiver",colnames(Counts_by_Behaviour_Rec))] <- "Role"
colnames(Counts_by_Behaviour_Act) [match("Actor",colnames(Counts_by_Behaviour_Act))] <- "Role"
colnames(Durations_by_Behaviour_Rec) [match("Receiver",colnames(Durations_by_Behaviour_Rec))] <- "Role"
colnames(Durations_by_Behaviour_Act) [match("Actor",colnames(Durations_by_Behaviour_Act))] <- "Role"

Counts_by_Beh_Role <- dplyr::bind_rows(Counts_by_Behaviour_Rec=Counts_by_Behaviour_Rec,Counts_by_Behaviour_Act=Counts_by_Behaviour_Act)
Durations_by_Beh_Role <- dplyr::bind_rows(Durations_by_Behaviour_Rec=Durations_by_Behaviour_Rec,Durations_by_Behaviour_Act=Durations_by_Behaviour_Act)

## merge counts & durations carefully
Counts_by_Beh_Role   <-  plyr::join(x=Counts_by_Beh_Role, y=Durations_by_Beh_Role, type = "full", match = "all")

#ADD THE R9SP TASKS OR YOU LOOSE R9SP INFORMATION!
Counts_by_Beh_Role_Task <- plyr::join (x = Task_list_R3SP[,c("antID","AntTask","treatment_rep")] , y=Counts_by_Beh_Role, type = "full", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            
NROW(Counts_by_Beh_Role_Task)
NROW(Counts_by_Beh_Role)

#assign "exposed" label to nurses that were exposed, based on ExposedAnts file
Counts_by_Beh_Role_Task_Exp <- plyr::join (x = ExposedAnts[,c("antID","treatment_rep","Exposed")] , y=Counts_by_Beh_Role_Task, type = "full", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            
#assign "non-exposed" label to everyone else
Counts_by_Beh_Role_Task_Exp$Exposed[is.na(Counts_by_Beh_Role_Task_Exp$Exposed)] <- "non_exp"
head(Counts_by_Beh_Role_Task_Exp)
#crete new Exposure+task column for plots
Counts_by_Beh_Role_Task_Exp$Exp_task <- paste(Counts_by_Beh_Role_Task_Exp$Exposed,Counts_by_Beh_Role_Task_Exp$AntTask,sep = "_")
unique(Counts_by_Beh_Role_Task_Exp$Exp_task)

## add the missing cases
Counts_by_Beh_Role_Task_Exp_AllComb <- plyr::join (x = Counts_by_Beh_Role_Task_Exp , y=all_combos, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            

## Focus only on a few important behaviours
Counts_by_Beh_Role_Task_Exp_AllComb <- Counts_by_Beh_Role_Task_Exp_AllComb[which(Counts_by_Beh_Role_Task_Exp_AllComb$Behaviour %in% c("A","G","SG","T","TB")),]
## Skip for now the undetermined groups, REMOVE ONCE TASKS ARE DETERMINED FOR R9SP
Counts_by_Beh_Role_Task_Exp_AllComb <- Counts_by_Beh_Role_Task_Exp_AllComb[which(Counts_by_Beh_Role_Task_Exp_AllComb$Exp_task %in% c("exp_nurse","non_exp_forager","non_exp_nurse")),]


## replace the NAs with 0 counts
Counts_by_Beh_Role_Task_Exp_AllComb$Count[which(is.na(Counts_by_Beh_Role_Task_Exp_AllComb$Count))] <- 0
## finally, get the mean & S.E. for each behav before/after  for barplots
Counts_by_Beh_Role_Task_Exp_MEAN  <- aggregate(cbind(Count,duration) ~ Behaviour + period + Role + Exp_task,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Beh_Role_Task_Exp_AllComb)
Counts_by_Beh_Role_Task_Exp_SE    <- aggregate(cbind(Count,duration) ~ Behaviour + period + Role + Exp_task,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Beh_Role_Task_Exp_AllComb)


#CROSSCHECK OPERATIONS WERE RIGHT, 
#THERE SEEMS TO BE A DISCREPANCY BETWEEN Counts_by_Beh_Role_Task_Exp_MEAN AND Counts_by_Behaviour_MEAN
# SAME FOR THE SE MEASURE. Counts_by_Beh_Role_Task_Exp_MEAN AND SE SHOWS
#OF COURSE IT IS NOT MATCHING! IT IS STILL MISSING THE INFORMATION FROM THE R9SP NEST!!!!!!!!!!!!!!


###############################################################################
###### PLOTTING  ##############################################################
###############################################################################

#probably would be convenient to run a loop for actors and receivers and to provide subsets of data per Behaviour
#given that the focus is on specific responses by Exp_task and Role

#test with one behaviour
# Counts_by_Beh_Role_Task_Exp_MEAN <- Counts_by_Beh_Role_Task_Exp_MEAN[which(Counts_by_Beh_Role_Task_Exp_MEAN$Behaviour %in% c("G")),]
# Counts_by_Beh_Role_Task_Exp_SE <- Counts_by_Beh_Role_Task_Exp_SE[which(Counts_by_Beh_Role_Task_Exp_SE$Behaviour %in% c("G")),]

#come fare??????


## show the mean counts for each behav | stage
pdf(file=paste(DATADIR,"Behaviour_counts_pre-post.pdf", sep = ""), width=5, height=3)
par(mfrow=c(1,2), family="serif", mai=c(0.4,0.5,0.1,0.1), mgp=c(1.3,0.3,0), tcl=-0.2)
## COUNTS
Xpos <- barplot( Count ~ period + Exp_task , Counts_by_Beh_Role_Task_Exp_MEAN, beside=T, xlab="", ylab="Behaviour count", ylim=c(0,max( Counts_by_Beh_Role_Task_Exp_MEAN$Count +  Counts_by_Beh_Role_Task_Exp_SE$Count*1.4 )))
##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
segments(x0 = Xpos[1,], 
         x1 = Xpos[1,], 
         y0 = Counts_by_Beh_Role_Task_Exp_MEAN$Count [Counts_by_Beh_Role_Task_Exp_MEAN$period=="pre"], 
         y1 = Counts_by_Beh_Role_Task_Exp_MEAN$Count [Counts_by_Beh_Role_Task_Exp_MEAN$period=="pre"] + Counts_by_Beh_Role_Task_Exp_SE$Count [Counts_by_Beh_Role_Task_Exp_SE$period=="pre"],
         lwd=2)

segments(x0 = Xpos[2,], 
         x1 = Xpos[2,], 
         y0 = Counts_by_Beh_Role_Task_Exp_MEAN$Count [Counts_by_Beh_Role_Task_Exp_MEAN$period=="post"], 
         y1 = Counts_by_Beh_Role_Task_Exp_MEAN$Count [Counts_by_Beh_Role_Task_Exp_MEAN$period=="post"] + Counts_by_Beh_Role_Task_Exp_SE$Count [Counts_by_Beh_Role_Task_Exp_SE$period=="post"],
         lwd=2)

text(x = ((Xpos[1,]+Xpos[2,])/2),
     y = Counts_by_Beh_Role_Task_Exp_MEAN$Count [Counts_by_Beh_Role_Task_Exp_MEAN$period=="post"] + Counts_by_Beh_Role_Task_Exp_SE$Count [Counts_by_Beh_Role_Task_Exp_SE$period=="post"]+15,
     stars.pval(posthoc_FREQ_summary$p.value))

# ## DURATIONS
# Xpos <- barplot( duration ~ period + Behaviour , Counts_by_Beh_Role_Task_Exp_MEAN, beside=T, xlab="", ylab="Behaviour duration (s)", ylim=c(0,max( Counts_by_Beh_Role_Task_Exp_MEAN$duration +  Counts_by_Beh_Role_Task_Exp_SE$duration, na.rm=T)))
# ##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
# segments(x0 = Xpos[1,], 
#          x1 = Xpos[1,], 
#          y0 = Counts_by_Beh_Role_Task_Exp_MEAN$duration [Counts_by_Beh_Role_Task_Exp_MEAN$period=="pre"], 
#          y1 = Counts_by_Beh_Role_Task_Exp_MEAN$duration [Counts_by_Beh_Role_Task_Exp_MEAN$period=="pre"] + Counts_by_Beh_Role_Task_Exp_SE$duration [Counts_by_Beh_Role_Task_Exp_SE$period=="pre"],
#          lwd=2)
# segments(x0 = Xpos[2,], 
#          x1 = Xpos[2,], 
#          y0 = Counts_by_Beh_Role_Task_Exp_MEAN$duration [Counts_by_Beh_Role_Task_Exp_MEAN$period=="post"], 
#          y1 = Counts_by_Beh_Role_Task_Exp_MEAN$duration [Counts_by_Beh_Role_Task_Exp_MEAN$period=="post"] + Counts_by_Beh_Role_Task_Exp_SE$duration [Counts_by_Beh_Role_Task_Exp_SE$period=="post"],
#          lwd=2)
# 
# text(x = ((Xpos[1,]+Xpos[2,])/2),
#      y = Counts_by_Beh_Role_Task_Exp_MEAN$duration [Counts_by_Beh_Role_Task_Exp_MEAN$period=="post"] + Counts_by_Beh_Role_Task_Exp_SE$duration [Counts_by_Beh_Role_Task_Exp_SE$period=="post"]+3,
#      stars.pval(rev(posthoc_DUR_summary$p.value)))
# 







#check for normality visually first
#look at the data histograms
beh_dist_hist <- ggplot(annotations,aes(x=duration))+geom_histogram()+facet_grid(~Behaviour, scales = "free")+theme_bw() +
  ggtitle("Behaviour duration distribution")
beh_dist_hist

##############################################
######CROSS-CHECK ANNOTATIONS ################
##############################################

#The createDataPartition() function is meant to subset a dataset without losing the probability distribution of your target variable.
annotations$RowID <- seq.int(nrow(annotations))
annotations$BEH_AW <- NA
annotations$RowID<- as.character(annotations$RowID)
str(annotations)

library("caret")
my.ids <- createDataPartition(annotations$Behaviour, p = 0.5)
annotations_subset <- annotations[as.numeric(my.ids[[1]]), ]

#create 25% set as first step of Cross Validation agreement, if agreement is high stop here, else continue check on the larger set 
my.ids2 <- createDataPartition(annotations_subset$Behaviour, p = 0.5)
annotations_subset2 <- annotations_subset[as.numeric(my.ids2[[1]]), ]
# 
# #You can check the distribution of your target variable in the population and in your subset.
# par(mfrow = c(1,3))
# barplot(table(annotations$Behaviour), main = "full dataset")
# barplot(table(annotations_subset$Behaviour), main = "subset 50%")
# barplot(table(annotations_subset2$Behaviour), main = "subset 25%")

print("COMPLETED ANALYSING THE MANUALLY ANNOTATED BEHAVIOUR LABELS")


# avoid file overwriting!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# if ()
#   
# write.csv(annotations_subset, "/home/cf19810/Dropbox/Ants_behaviour_analysis/Cross_Validation/annotations_subset_50%_2.csv")
# write.csv(annotations_subset2, "/home/cf19810/Dropbox/Ants_behaviour_analysis/Cross_Validation/annotations_subset_25%_2.csv")
