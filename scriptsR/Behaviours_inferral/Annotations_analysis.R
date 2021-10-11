
###############################################################################
###### LOAD MANUAL ANNOTATIONS ################################################
###############################################################################
#Behavioural codes explanation
# beh_codes <- read.csv("behavioural_codes.csv", sep= ",")



annotations <- read.csv(paste(DATADIR,"/R3SP_R9SP_All_data_dropped_useless_cols.csv",sep = ""), sep = ",")
annotations$Behaviour <- as.character(annotations$Behaviour)
annotations$Actor <- as.character(annotations$Actor)
annotations$Receiver <- as.character(annotations$Receiver)
#call treatment as period
colnames(annotations)[which(names(annotations) == "treatment")] <- "period"


#convert Zulu time to GMT
annotations$T_start <- as.POSIXct(annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
annotations$T_stop  <- as.POSIXct(annotations$T_stop, format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )
annotations$duration <- as.numeric(annotations$T_stop - annotations$T_start)

#see if milliseconds are shown (number of decimals represented by the number after %OS)
format(annotations$T_start[3], "%Y-%m-%d %H:%M:%OS6")

#remove duplicates of directed behaviours (Grooming and Aggression) by keeping only the behaviours where the Focal corresponds to the Actor.
#this seems to work very well with Grooming (cuts 50% of the events) and Agrgession (cuts 15 over 31 events) but affects also 4 Trophallaxis events, check why
annotations_drop_G_A <- annotations[which(annotations$Actor==annotations$Focal),]

#remove duplicates of non-directed behaviours as Trophalllaxis based on multiple columns check
annotations_drop_all_dup <- annotations_drop_G_A[!duplicated(annotations_drop_G_A[,c("T_start","T_stop","Behaviour","treatment_rep")]),]

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

## count the number of observations of each behaviour - WARNING; some behavs not observed e.g. before the treatment (period), so will need to account for that (next step)
Counts_by_Behaviour_CLEAN    <- aggregate(Actor ~ Behaviour + period + treatment_rep, FUN=length, na.action=na.omit, annotations); colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
## calculate mean durations for each behaviour
Durations_by_Behaviour_CLEAN <- aggregate(duration ~ Behaviour + period + treatment_rep, FUN=mean, na.action=na.omit, annotations)
## merge counts & durations carefully
Counts_by_Behaviour_CLEAN    <-  plyr::join(x=Counts_by_Behaviour_CLEAN, y=Durations_by_Behaviour_CLEAN, , type = "full", match = "all")

## create a data frame with all combinations of the conditioning variables
all_combos <- expand.grid ( Behaviour=unique(annotations$Behaviour), period=unique(annotations$period), treatment_rep=unique(annotations$treatment_rep))

## add the missing cases
Counts_by_Behaviour_AllCombos <- plyr::join (x = Counts_by_Behaviour_CLEAN , y=all_combos, type = "right", match = "all")  #, by.x=c("Behaviour","period","treatment_rep"), by.y=c("Behaviour","period","treatment_rep") )            

## Focus only ion a few important behaviours
Counts_by_Behaviour_AllCombos$Behaviour <- as.character(Counts_by_Behaviour_AllCombos$Behaviour)  ## naughty R
Counts_by_Behaviour_AllCombos <- Counts_by_Behaviour_AllCombos[which(Counts_by_Behaviour_AllCombos$Behaviour %in% c("A","G","SG","T","TB")),]


## replace the NAs with 0 counts            
Counts_by_Behaviour_AllCombos$Count[which(is.na(Counts_by_Behaviour_AllCombos$Count))] <- 0
## finally, get the mean & S.E. for each behav before/after  for barplots
Counts_by_Behaviour_MEAN  <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)
Counts_by_Behaviour_SE    <- aggregate(cbind(Count,duration) ~ Behaviour + period,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)

## show the mean counts for each behav | stage
pdf(file=paste(DATADIR,"Behaviour_counts_pre-post.pdf", sep = ""), width=5, height=3)
par(mfrow=c(1,2), family="serif", mai=c(0.4,0.5,0.1,0.1), mgp=c(1.3,0.3,0), tcl=-0.2)
## COUNTS
Xpos <- barplot( Count ~ period + Behaviour , Counts_by_Behaviour_MEAN, beside=T, xlab="", ylab="Behaviour count", ylim=c(0,max( Counts_by_Behaviour_MEAN$Count +  Counts_by_Behaviour_SE$Count)))
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
## close the pdf
dev.off()

## TO DO : statistics on before / after counts & durations

############################################################################
###### STATISTICS ON COUNTS ################################################
############################################################################
library(lme4)
library(blmeco) # check dispersion for glmer
library(emmeans) #post-hoc comparisons
library(e1071) #calc skewness and other stuff

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
m1 <- lmerTest::lmer(Count ~ period + Behaviour + (1|treatment_rep), Counts_by_Behaviour_AllCombos)
summary(m1)
test_norm(residuals(m1)) #test residuals' normality. null hypothesis for the Shapiro-Wilk test is that a variable is normally distributed
m2 <- glmer(Count ~ period + Behaviour + (1|treatment_rep), data = Counts_by_Behaviour_AllCombos, family = "poisson")
summary(m2)
test_norm(residuals(m2))

par(mfrow=c(1,2))
plot(m2)
qqnorm(residuals(m2))
qqline(residuals(m2))
hist(residuals(m2)) 

#check poisson distribution overdispersion
#if the ratio of residual deviance to degrees of freedom it's >1 then the data are overdispersed
dispersion_glmer(m2)

#overdispersion can be incorporated by including an observation-level random effect (1|obs), to allow observations to be correlated with themselves
Counts_by_Behaviour_AllCombos$obs=factor(1:nrow(Counts_by_Behaviour_AllCombos))
m5 <- glmer(Count ~ period + Behaviour + (1|treatment_rep) +(1|obs), data = Counts_by_Behaviour_AllCombos, family = "poisson")
summary(m5) # is singular fit a problem?
dispersion_glmer(m5)#model m5, which assumes poisson distribution and accounts for overdispersion by adding an observation-level random effect shows both lower AIC and dispersion ~1
anova(m2,m5)

#test residuals' normality
test_norm(residuals(m5)) #residuals are non normal! 
par(mfrow=c(1,2))
plot(m5)
#dev.off()
qqnorm(residuals(m5))
qqline(residuals(m5))
hist(residuals(m5)) #distribution of residuals doesn't look very good

#ISSUE: model m2 (poisson) suffers overdispersion but fixing it causes normality assumption not to be met anymore...

#getting post-hoc comparison per behaviour between period (pre/post)
posthoc <- emmeans(m5, specs = trt.vs.ctrlk ~ period | Behaviour) #EVERY CONTRAST HAS SAME ESTIMATES WITH SUPER HIGH SIGNIFICANCE, EVEN WHEN IT SHOULDN'T (E.G. TROPHALLAXIS)
posthoc$contrast
##or
#emmeans(m5, list(pairwise ~ period | Behaviour), adjust = "tukey")
##or
# contrast(emmeans(m5,~period|Behaviour, type="response"),
#          method="pairwise", adjust="Tukey")  #method = "trt.vs.ctrl"



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
annotations$Behaviour <- as.factor(annotations$Behaviour)
annotations$Focal <- as.factor(annotations$Focal)
str(annotations)

#Counts_by_Behaviour_AllCombos$period = relevel(Counts_by_Behaviour_AllCombos$period, ref="pre")
m_dur1 <- lmerTest::lmer(duration ~ period + Behaviour + (1|treatment_rep), annotations)
summary(m_dur1)
test_norm(residuals(m_dur1)) #skewed data

#how to solve this? look into SUMMARY_DATA_MOVEMENT_131219 FOR POINTERS
m_dur2 <- glmer(duration ~ period + Behaviour + (1|treatment_rep), data = annotations, family = "   ")
summary(m_dur2)
test_norm(residuals(m_dur2)) 





print("I'D LOVE TO HAVE COMPLETED ANALYSING THE MANUALLY ANNOTATED BEHAVIOUR LABELS ALREADY!!")

##############################################
######CROSS-CHECK ANNOTATIONS ################
##############################################

#The createDataPartition() function is meant to subset a dataset without losing the probability distribution of your target variable.
annotations$RowID <- seq.int(nrow(annotations))
annotations$BEH_AW <- NA

library("caret")
my.ids <- createDataPartition(annotations$Behaviour, p = 0.5)
annotations_subset <- annotations[as.numeric(my.ids[[1]]), ]

#You can check the distribution of your target variable in the population and in your subset.
par(mfrow = c(1,2))
barplot(table(annotations$Behaviour), main = "full dataset")
barplot(table(annotations_subset$Behaviour), main = "subset")

# avoid file overwriting!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# if ()
#   
#   list.files(path=DATADIR, pattern=REPLICATE, full.names=T); MyrmidonCapsuleFile <- MyrmidonCapsuleFile[grepl("Capsule_defined.myrmidon",MyrmidonCapsuleFile)]  
#write.csv(annotations_subset, "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/Annotation_Cross_Validation/annotations_subset_50%.csv")
