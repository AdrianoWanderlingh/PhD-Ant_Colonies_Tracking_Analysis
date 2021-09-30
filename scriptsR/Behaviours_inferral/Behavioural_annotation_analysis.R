rm(list=ls())
###########
#
#check the analysis route as outlined in https://docs.google.com/document/d/1jxa8D5hHv12OnhDZBANQGnI4pwNfRtQt/edit#heading=h.s0n4rczfptco
#
##########
setwd("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/Annotation_Cross_Validation/")
getwd()

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(parsedate)
library(forcats) #what's it for?

annotations <- read.csv("R3SP_R9SP_All_data_dropped_useless_cols.csv", sep = ",")
annotations$Behaviour <- as.character(annotations$Behaviour)
annotations$Actor <- as.character(annotations$Actor)
annotations$Receiver <- as.character(annotations$Receiver)

#convert Zulu time to GMT
annotations$T_start <- as.POSIXct(annotations$T_start, format = "%Y-%m-%dT%H:%M:%OSZ") 
annotations$T_stop <- as.POSIXct(annotations$T_stop, format = "%Y-%m-%dT%H:%M:%OSZ")

#see if milliseconds are shown (number of decimals represented by the number after %OS)
format(annotations$T_start[3], "%Y-%m-%d %H:%M:%OS6")

#remove duplicates of directed behaviours (Grooming and Aggression) by keeping only the behaviours where the Focal corresponds to the Actor.
#this seems to work very well with Grooming (cuts 50% of the events) and Agrgession (cuts 15 over 31 events) but affects also 4 Trophallaxis events, check why
annotations_nodup <- annotations[which(annotations$Actor==annotations$Focal),]

Counts_by_Behaviour <- aggregate(treatment_rep ~ Behaviour + treatment, FUN=length, annotations); colnames(Counts_by_Behaviour) [match("treatment_rep",colnames(Counts_by_Behaviour))] <- "Count"
Counts_by_Behaviour_nodups <- aggregate(treatment_rep ~ Behaviour + treatment, FUN=length, annotations_nodup); colnames(Counts_by_Behaviour_nodups) [match("treatment_rep",colnames(Counts_by_Behaviour_nodups))] <- "Count_nodups"

#check how many cases have been removed
Counts_by_Behaviour <- cbind(Counts_by_Behaviour, count_nodups = Counts_by_Behaviour_nodups$Count_nodups)

#data<-basedata[!(basedata$Behaviour=="CB" | basedata$Behaviour=="TQ" | basedata$Behaviour=="GQ" | basedata$Behaviour=="FR" | basedata$Behaviour=="CR"),]

############################################
#######DATA PREPARATION#####################
############################################

#we should prepare the data:
#1. Remove double observations: grooming is both made and received, our data should include only one of the two occurrences (i,e. grooming received). Same for trophallaxis
#3. to ease calculations, times  should be normalised as starting from the observation start time (may be an extra column of the dataset) 
#   to help calcluations if making timed observation blocks (see below)
#2. possibly summarise to make ready for analyisis as needed (read below)


#important to keep in mind when deciding/ performing the tests:
#1. each subject (ant) may have repeated measures if performs the behaviour more than once
#2. the N of measures (both total that by individual) will certainly not match pre and post treatment as the N of individuals performing behaviours is going to be different.
#3. will this interfere with tests?

#3a. likely NO: it is possible to compare time series (duration) of unequal data size, see:
#https://www.researchgate.net/post/how_do_we_measure_the_similarity_between_two_time_series_depending_on_magnitude
#if this is a solution, ignore 3b. Similar considerations should be done not only for Duration but also for frequency

#3b. likely YES: Our factor of interest is more properly the nest, not the indiviual ants. 
#   If the observations are summarised by calculating the mean (total per nest-treatment or by time block per nest-treatment, as 5 mins block, 6 blocks per nest-treatment),
#   of time in seconds there would be perfectly paired repeated measures!
#Sum by time could be, for example: blocks of 5 minutes, where in 0-5min go all observations with start time <5:00
#SO, as A paired t-test (or Wilcoxon test) with sample sizes does not make any sense,
# thus if that is the analysis route taken, it will be needed sum data up but the N of groups should not decrease over a threshold as 
# the statistical power will be affected! see https://www.blopig.com/blog/2013/10/wilcoxon-mann-whitney-test-and-a-small-sample-size/

#It may be the case to use bootstrapping on these data, but the topic should be investigated more.
# For example, not wanting to sum up or make the mean of all individual ants, bootstrapping could help!
# Reading from "Practical 3 Programming in R(1).pdf" (folder "UNDERSTANDING DATA Course"):
# "The main point here is that these are skewed data that only take positive values. They might, for example, be the absolute asymmetry in
# the tail streamers of some bird species with elongated display plumage. With such data you:
# . Can't do a t-test on the mean because they data are so skewed.
# . Can't do a one-sample Wilcoxon test on the median, because it assumes a symmetrical
# distribution.
# . Can't do a permutation test for the median being zero, because no permutation of these data 
#   can have a median equal to or less than zero (because no values are less than zero)."

############################################
#######CALCULATE DURATION###################
############################################

annotations$duration  <- as.numeric(annotations$T_stop - annotations$T_start, units.difftime= "seconds")

############################################
#######SOME PLOTTING########################
############################################

#Regarding plots, the consensus is to use boxplots when residuals cannot be normalised and non-parametric tests
#are used, and mean +- standard error plots (either barplots or interval plots) when parametric tests are used.

#plot by removing some of the extra behaviours 

annotations<-annotations[!(annotations$Behaviour=="CB" | annotations$Behaviour=="TQ" | annotations$Behaviour=="GQ" | annotations$Behaviour=="FR" | annotations$Behaviour=="CR"),]
unique(annotations$Behaviour)
#annotationz<-annotations[which(annotations$Behaviour== "SG"),]

#MEAN BEHAVIOUR DURATION plot
mean_dur_plot<- ggplot(data = annotations , aes(x=Behaviour, y=duration, fill=fct_rev(treatment))) + 
  geom_point(aes(colour=treatment)) + #how the hell can I show the datapoints by treatment??
  scale_fill_manual(values=c("skyblue2","tomato1")) +
  theme(legend.position = c(0.85, 0.85),legend.title = element_blank()) +
  ggtitle("Mean behaviour duration \n before and after pathogen exposure")

#boxplot
mean_dur_boxplot <- mean_dur_plot + geom_boxplot()

#violinplot
mean_dur_plot + geom_violin() 

#data_agg <- aggregate(frequency~Behaviour+treatment, annotations, sum)

#make summary for the plotting#
table(annotations$Behaviour)
Counts_by_Behaviour <- aggregate(treatment_rep ~ Behaviour + treatment, FUN=length, annotations); colnames(Counts_by_Behaviour) [match("treatment_rep",colnames(Counts_by_Behaviour))] <- "Count"

#remove duplicates of the encounter behavours
#if (annotations$Behaviour=="G") {annotations_G_nodups<- annotations$Behaviour[!duplicated(annotations$Behaviour),]}

# for (BEH in c("G","T"))
# {
#   annot_BEH <- annotations[which(annotations$Behaviour==BEH),]
#   ## remove doubled allo-grooming interactions
#   if (BEH=="G") {annot_BEH <- annot_BEH[!duplicated(annot_BEH),]}  ## leave NOT to catch possible un-matched rows
#   if (BEH=="T") {print("FIND A WAY TO REMOVE DUPLICATES, MAY THAT BE USING THE FOCAL INDIVIDUAL?")}
# }


#BEHAVIOUR FREQUENCY
#
#to this graph should be added the Standard deviation by nest! 
freq_plot <- ggplot(data = Counts_by_Behaviour , aes(x=Behaviour, y=Count, fill=treatment)) + 
  geom_bar(
    aes(fill = fct_rev(treatment)), stat = "identity", position = position_dodge()
  ) + 
  scale_fill_manual(values=c("skyblue2","tomato1")) +
  theme(legend.position = c(0.85, 0.85),legend.title = element_blank()) +
  ggtitle("Behaviour frequency \n before and after pathogen exposure")

grid.arrange(mean_dur_boxplot, freq_plot,
             ncol = 2, nrow = 1)


############################################
#######DATA SUMMARY#########################
############################################

###summary with different tools:

###FREQUENCY summary with tapply
with(annotations, tapply(occurrence, list(Behaviour,treatment), sum))

###DURATION summary with aggregate (nice but doesn't save correctly, needs correction)
duration_summary <- with(annotations, aggregate(duration ~ Behaviour + treatment, FUN =  function(x) c( SD = sd(x), MN= mean(x), MAX= max(x) )))
duration_summary #the function of aggregate works but the object doesn't save correctly!
#check for outliers, may be not relevant to select them if the test is roust for outliers
with(annotations, aggregate(duration ~ Behaviour + treatment, FUN =  function(x) c( out = boxplot.stats(x)$out )))

###FREQUENCY summary with plyr
library(plyr)

freq_summary <- ddply(annotations, c("Behaviour","treatment"), summarise, grp.tot=sum(occurrence))
freq_summary

###DURATION summary with plyr, used for the following histograms
#all values
dur_summary <- ddply(annotations, c("Behaviour","treatment"), summarise, grp.mean=mean(duration))
dur_summary

########################################################################################
#DO NOT EXCLUDE VALUES! LONGER EVENTS (IE FOR GROOMING) SHOWED TO BE MORE COMPLEX BEHAVIOURS WHEN PLOTTED!
#############################################################################################

#values excluding those 6 times larger than the mean (value chosen to exclude the larger ones and
#visualise graphs more easily. outliers should not be excluded, especially in non-parametric tests).
#If dropping them in parametric tests helps, it may be done if reasonably executed (truncate them instead of elimination?)
data_without_outliers <- annotations[!annotations$duration>6*mean(annotations$duration),]
dur_summary_nolarge <- ddply(data_without_outliers, c("Behaviour","treatment"), summarise, grp.mean=mean(duration))
dur_summary_nolarge

############################################
#######CHECKING NORMALITY###################
############################################

#check for normality visually first
#look at the data histograms
ggplot(annotations,aes(x=duration))+geom_histogram()+facet_grid(~Behaviour)+theme_bw()

###plot divided by behaviour and showing treatment with the mean per group 
#with large values
ggplot(annotations,aes(x=duration, color=treatment))+
  geom_histogram(fill="white",alpha=0.5, position="identity")+
  facet_grid(~Behaviour)+
  geom_vline(data=dur_summary, aes(xintercept=grp.mean, color=treatment),
             linetype="dashed")+
  theme_bw() + theme(legend.position="top")


#large values removed
ggplot(data_without_outliers,aes(x=duration, color=treatment))+
  geom_histogram(fill="white",alpha=0.5, position="identity")+
  facet_grid(~Behaviour)+
  geom_vline(data=dur_summary_nolarge, aes(xintercept=grp.mean, color=treatment),
             linetype="dashed")+
  theme_bw() + theme(legend.position="top")


###check for normality with glm (see Constance Script/summary_data_movement_131219.R from line 133 on)
# by seeing if the residuals have a normal distribution


############################################
#######TESTING##############################
############################################

# If we sum up observations by ant (I think it is the most sensible choice):
#per each behaviour, we will need to compare two paired groups,
#If data is parametric (normal residuals), use	Paired t test; if data is non-parametric, use	Wilcoxon signed rank test on paired samples




# WILCOXON TEST ASSUMPTIONS:
#1. dependent variable should be measured at the ordinal or continuous level: true for both occurrence and time_sec
#2. Your independent variable should consist of two categorical, "related groups" or "matched pairs": yes (pre-post)
#3. The distribution of the differences between the two related groups needs to be symmetrical in shape
#to check the symmetry assumption of the wilcoxon test we should probably follow https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/#signed-rank-test-on-paired-samples

