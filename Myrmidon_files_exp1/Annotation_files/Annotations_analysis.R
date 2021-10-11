
###############################################################################
###### LOAD MANUAL ANNOTATIONS ################################################
###############################################################################
#Behavioural codes explanation
# beh_codes <- read.csv("behavioural_codes.csv", sep= ",")



annotations <- read.csv(paste(DATADIR,"/R3SP_R9SP_All_data_dropped_useless_cols.csv",sep = ""), sep = ",")
annotations$Behaviour <- as.character(annotations$Behaviour)
annotations$Actor <- as.character(annotations$Actor)
annotations$Receiver <- as.character(annotations$Receiver)

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
Counts_by_Behaviour              <- aggregate(treatment_rep ~ Behaviour + treatment, FUN=length, annotations); colnames(Counts_by_Behaviour) [match("treatment_rep",colnames(Counts_by_Behaviour))] <- "Count"
Counts_by_Behaviour_drop_G_A     <- aggregate(treatment_rep ~ Behaviour + treatment, FUN=length, annotations_drop_G_A); colnames(Counts_by_Behaviour_drop_G_A) [match("treatment_rep",colnames(Counts_by_Behaviour_drop_G_A))] <- "Count_drop_G_A"
Counts_by_Behaviour_drop_all_dup <- aggregate(treatment_rep ~ Behaviour + treatment, FUN=length, annotations_drop_all_dup); colnames(Counts_by_Behaviour_drop_all_dup) [match("treatment_rep",colnames(Counts_by_Behaviour_drop_all_dup))] <- "Count_drop_all_dup"
#check how many cases have been removed
Counts_by_Behaviour <- cbind(Counts_by_Behaviour, Count_drop_G_A = Counts_by_Behaviour_drop_G_A$Count_drop_G_A, Count_drop_all_dup = Counts_by_Behaviour_drop_all_dup$Count_drop_all_dup)
Counts_by_Behaviour

#see total final numer of behaviours
Counts_by_Behaviour_tots <- aggregate(Count_drop_all_dup ~ Behaviour, FUN=sum, Counts_by_Behaviour); colnames(Counts_by_Behaviour) [match("treatment",colnames(Counts_by_Behaviour))] <- "Totals"



#Over-write cleaned dataset - NOTE 'THIS'annotations' IS USED in the trajectory plotting loop below! 
annotations <- annotations_drop_all_dup

## count the number of observations of each behaviour - WARNING; some behavs not observed e.g. before the treatment, so will need to account for that (next step)
Counts_by_Behaviour_CLEAN    <- aggregate(Actor ~ Behaviour + treatment + treatment_rep, FUN=length, na.action=na.omit, annotations); colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count"
## calculate mean durations for each behaviour
Durations_by_Behaviour_CLEAN <- aggregate(duration ~ Behaviour + treatment + treatment_rep, FUN=mean, na.action=na.omit, annotations)
## merge counts & durations carefully
Counts_by_Behaviour_CLEAN    <-  plyr::join(x=Counts_by_Behaviour_CLEAN, y=Durations_by_Behaviour_CLEAN, , type = "full", match = "all")

## create a data frame with all combinations of the conditioning variables
all_combos <- expand.grid ( Behaviour=unique(annotations$Behaviour), treatment=unique(annotations$treatment), treatment_rep=unique(annotations$treatment_rep))

## add the missing cases
Counts_by_Behaviour_AllCombos <- plyr::join (x = Counts_by_Behaviour_CLEAN , y=all_combos, type = "right", match = "all")  #, by.x=c("Behaviour","treatment","treatment_rep"), by.y=c("Behaviour","treatment","treatment_rep") )            

## Focus only ion a few important behaviours
Counts_by_Behaviour_AllCombos$Behaviour <- as.character(Counts_by_Behaviour_AllCombos$Behaviour)  ## naughty R
Counts_by_Behaviour_AllCombos <- Counts_by_Behaviour_AllCombos[which(Counts_by_Behaviour_AllCombos$Behaviour %in% c("A","G","SG","T","TB")),]


## replace the NAs with 0 counts            
Counts_by_Behaviour_AllCombos$Count[which(is.na(Counts_by_Behaviour_AllCombos$Count))] <- 0
## finally, get the mean & S.E. for each behav before/after  for barplots
Counts_by_Behaviour_MEAN  <- aggregate(cbind(Count,duration) ~ Behaviour + treatment,                 FUN=mean,      na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)
Counts_by_Behaviour_SE    <- aggregate(cbind(Count,duration) ~ Behaviour + treatment,                 FUN=std.error, na.rm=T, na.action=NULL, Counts_by_Behaviour_AllCombos)

## show the mean counts for each behav | stage
pdf(file=paste(DATADIR,"Behaviour_counts_pre-post.pdf", sep = ""), width=5, height=3)
par(mfrow=c(1,2), family="serif", mai=c(0.4,0.5,0.1,0.1), mgp=c(1.3,0.3,0), tcl=-0.2)
## COUNTS
Xpos <- barplot( Count ~ treatment + Behaviour , Counts_by_Behaviour_MEAN, beside=T, xlab="", ylab="Behaviour count", ylim=c(0,max( Counts_by_Behaviour_MEAN$Count +  Counts_by_Behaviour_SE$Count)))
##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
segments(x0 = Xpos[1,], 
         x1 = Xpos[1,], 
         y0 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$treatment=="pre"], 
         y1 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$treatment=="pre"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$treatment=="pre"],
         lwd=2)
segments(x0 = Xpos[2,], 
         x1 = Xpos[2,], 
         y0 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$treatment=="post"], 
         y1 = Counts_by_Behaviour_MEAN$Count [Counts_by_Behaviour_MEAN$treatment=="post"] + Counts_by_Behaviour_SE$Count [Counts_by_Behaviour_SE$treatment=="post"],
         lwd=2)

## DURATIONS
Xpos <- barplot( duration ~ treatment + Behaviour , Counts_by_Behaviour_MEAN, beside=T, xlab="", ylab="Behaviour duration (s)", ylim=c(0,max( Counts_by_Behaviour_MEAN$duration +  Counts_by_Behaviour_SE$duration, na.rm=T)))
##  add SE bars for the left bar in each behaviour - only add the upper SE limit to avoid the possibility of getting negative counts in the error
segments(x0 = Xpos[1,], 
         x1 = Xpos[1,], 
         y0 = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$treatment=="pre"], 
         y1 = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$treatment=="pre"] + Counts_by_Behaviour_SE$duration [Counts_by_Behaviour_SE$treatment=="pre"],
         lwd=2)
segments(x0 = Xpos[2,], 
         x1 = Xpos[2,], 
         y0 = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$treatment=="post"], 
         y1 = Counts_by_Behaviour_MEAN$duration [Counts_by_Behaviour_MEAN$treatment=="post"] + Counts_by_Behaviour_SE$duration [Counts_by_Behaviour_SE$treatment=="post"],
         lwd=2)
## close the pdf
dev.off()

## TO DO : statistics on before / after counts & durations
str(Counts_by_Behaviour_AllCombos)
Counts_by_Behaviour_AllCombos$treatment_rep <- as.factor(Counts_by_Behaviour_AllCombos$treatment_rep)
Counts_by_Behaviour_AllCombos$treatment<- as.factor(Counts_by_Behaviour_AllCombos$treatment)
Counts_by_Behaviour_AllCombos$Behaviour <- as.factor(Counts_by_Behaviour_AllCombos$Behaviour)
Counts_by_Behaviour_AllCombos$Count <- as.numeric(Counts_by_Behaviour_AllCombos$Count)


library(lme4)
m1 <- lmerTest::lmer(Count ~ treatment + Behaviour + (1|treatment_rep), Counts_by_Behaviour_AllCombos)
m2 <- glmer(Count ~ treatment + Behaviour + (1|treatment_rep), data = Counts_by_Behaviour_AllCombos, family = "poisson")
m3 <- glm(Count ~ treatment + Behaviour, data = Counts_by_Behaviour_AllCombos, family = poisson(link = "log"))
m4 <- glm(Count ~ treatment + Behaviour, data = Counts_by_Behaviour_AllCombos, family = quasipoisson(link = "log"))



m5<- glm.nb(Count ~ treatment + Behaviour, data = Counts_by_Behaviour_AllCombos)
#m6 <- glm(Count ~ treatment + Behaviour, data = Counts_by_Behaviour_AllCombos, family=negative.binomial(3.99))

summary(m1)
summary(m2)
summary(m3) #shows that res deviance/df is >>1
summary(m4) #shows that res deviance/df is >>1
summary(m5) 
summary(m6) 


#test poisson distribution overdispersion
library(AER)
#if the ratio of residual deviance to degrees of freedom it's >1 then the data are overdispersed
dispersiontest(m3,trafo=1)

#another measure of overdispersion with Pearson’s Chi-squared 
dp <- sum(residuals(m3,type ="pearson")^2)/m3$df.residual
dp

#We can check how much the coefficient estimations are affected by overdispersion.
summary(m3, dispersion=dp)
#we can use a quasi-family to estimate the dispersion parameter.
summary(m4) #model with uasipoisson distribution


#with summary(m5) overdispersion is still very high! how to correct for it?
#same problem of this https://stats.stackexchange.com/questions/123999/how-to-account-for-overdispersion-in-a-glm-with-negative-binomial-distribution



#test residuals' normality. null hypothesis for the Shapiro-Wilk test is that a variable is normally distributed
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

test_norm(residuals(m1))
test_norm(residuals(m2))








anova(m1,m2)







#test normality of residuals
library(e1071)

library(lmerTest)
library(MuMIn)
library(arm)
library(lmtest)
library(MASS)
install.packages(c("MASS"))



#including both fixed and random effects. R2c
r.squaredGLMM(m2)
lmtest::bptest(m2)  # Breusch-Pagan test, if significant residuals are not homoscedastic
anova(m2)
display(m2)
par(mfrow=c(1,2))
plot(m1)
#dev.off()
qqnorm(resid(m2))
qqline(resid(m2))
hist(resid(m2))



library(emmeans)
library(pbkrtest)
emmeans(m2, list(pairwise ~ treatment), adjust = "tukey")


print("COMPLETED ANALYSING THE MANUALLY ANNOTATED BEHAVIOUR LABELS!!")

