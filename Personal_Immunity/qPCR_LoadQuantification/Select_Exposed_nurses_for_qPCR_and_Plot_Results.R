
#CLEAN UP
gc()
mallinfo::malloc.trim(0L)


##### LIBRARIES
library(dplyr)
library(broman)
library(ggplot2)
library(stringr)
library(plotrix)
library(mallinfo)
## for stats
library(performance)
library(remotes)
library(remotes)
#install_version("MuMIn", "1.46.0")
library(MuMIn)
library(multcomp)
library(multcompView)
library(blmeco)
library(sjPlot)
library(lsmeans)
library(lme4)
library(blmeco) # check dispersion for glmer
library(emmeans) #post-hoc comparisons
library(e1071) #calc skewness and other stuff
library(lawstat) #for levene test (homogeneity of variance)
library(lmPerm) # for permutations
library(lmerTest)
library(car)
library(nlme) #lme
#library(rcompanion)
#install_version("rcompanion", "2.3.0")


#functions
RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}

#
STYLE <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw(),
              scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

STYLE_CONT <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw()
)

#
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


LoadOriData <- FALSE
WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data"
IMMUNITY_DATA <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Pathogen_Quantification_Data"

#############################################################################################
#############################################################################################

if (LoadOriData) {

  # Select Random Exposed nurses to do qpcr quantif
  ## 2 per colony
  
metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021_2022-10-12.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

metadata <- metadata[which(metadata$Exposed==TRUE),]
metadata <- metadata[which(metadata$IsAlive==TRUE),]

# #ADDED for second selection
# # exlcude already selected ants
# already_selected <- read.table(paste(WORKDIR,"/Personal_Immunity/Pathogen Quantification Data/220809-Adriano-MetIS2-colony-checkup_Analysis_with_Identities.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
# already_selected_sub <- already_selected[,2:4]
# already_selected_sub$Analysed_date <- "08/08/22"
# metadata <- left_join(metadata, already_selected_sub, by = c("REP_treat","antID"))                 # Apply left_join dplyr function
# metadata <- subset(metadata, is.na(metadata$Analysed_date))

# keep 2 ants for the big colonies
big <- c("BP","BS")
small <- c("SP","SS")

metadata_big <-
  metadata[which(metadata$treatment %in% big),] %>% 
  group_by(REP_treat) %>% 
  filter(row_number()==c(1,2))

# keep only 1 ant for the small colonies
metadata_small <-
  metadata[which(metadata$treatment %in% small),] %>% 
  group_by(REP_treat) %>% 
  filter(row_number()==1)

metadata_Selected <- rbind(metadata_big,metadata_small)

# #ADDED for second selection
# metadata_Selected <- metadata_small


metadata_Selected <- metadata_Selected[,c("REP_treat","antID","tagIDdecimal","Exposed","N_exposed")]

### ADD THE MISSING REPS, R9BS and R9SS

extra_missing <- data.frame(REP_treat=c("R9SS","R9BS","R9BS") , antID= c(9,3,34), tagIDdecimal= c(222,16,86), Exposed= c(TRUE,TRUE,TRUE), N_exposed= c(2,14,14) )
metadata_Selected <- rbind(metadata_Selected,extra_missing)

# convert decimal to hex

metadata_Selected$tagIDHEX <- convert2hex(metadata_Selected$tagIDdecimal)
metadata_Selected$tagIDHEX <- paste0("0x",metadata_Selected$tagIDHEX)

metadata_Selected <- metadata_Selected[,c(1,6,2,3,4,5)]


#save it! 

write.table(metadata_Selected,file="/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Select_Exposed_nurses_for_qPCR_8-08-22.txt",append=F,col.names=T,row.names=F,quote=T,sep=",")

# #ADDED for second selection
# write.table(metadata_Selected,file="/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Pathogen Quantification Data/Select_Exposed_nurses_for_qPCR_17-10-22.txt",append=F,col.names=T,row.names=F,quote=T,sep=",")

########################################################################
##### PLOT OUTPUTS #####################################################



Florant_output <- read.csv(paste0(IMMUNITY_DATA,"/220809-Adriano-MetIS2-colony-checkup_Analysis.csv"),header=T,stringsAsFactors = F, sep=",")
colnames(Florant_output)[which(colnames(Florant_output)=="Reduced.quantification..ng.µL..negative.if.below.detection.threshold.0.001.")] <-   "Red_quantif_ng.µL"

info_ants <- read.table(paste0(IMMUNITY_DATA,"/Select_Exposed_nurses_for_qPCR_8-08-22_ANNOTATED.txt"),header=T,stringsAsFactors = F, sep=",")


info_ants$Well <- gsub("^(.{1})(.*)$",         # Apply gsub
                       "\\10\\2",
                       info_ants$NEW_position)
# Print new string

### merge
Meta_all_combs <- list(info_ants,Florant_output)
Meta_all_combs <- Reduce(function(x, y) merge(x, y, all=TRUE), Meta_all_combs)


#Meta_all_combs[which(Meta_all_combs$Red_quantif_ng.µL == "#VALUE!"),] <- NA
No_CT_value <- 0.000001
Meta_all_combs[which(Meta_all_combs$Red_quantif_ng.µL == "#VALUE!"),"NOTE"] <- paste(Meta_all_combs[which(Meta_all_combs$Red_quantif_ng.µL == "#VALUE!"),"NOTE"] ,"No CT is either no DNA or issue in extraction",sep=".")
Meta_all_combs[which(Meta_all_combs$Red_quantif_ng.µL == "#VALUE!"),"Red_quantif_ng.µL"] <- No_CT_value #No CT (either no DNA or issue in extraction)


Meta_all_combs <- Meta_all_combs[which(Meta_all_combs$Red_quantif_ng.µL != "No sample"),]


Meta_all_combs$Treatment <- RIGHT(Meta_all_combs$REP_treat,2)

#save output!
write.table(Meta_all_combs,file="/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Pathogen_Quantification_Data/220809-Adriano-MetIS2-colony-checkup_Analysis_with_Identities.txt",append=F,col.names=T,row.names=F,quote=T,sep=",")


# Rename by name
Meta_all_combs$Treatment <- as.factor(Meta_all_combs$Treatment)
levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="BS"] <- "Big Sham"
levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="BP"] <- "Big Pathogen"
levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="SS"] <- "Small Sham"
levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="SP"] <- "Small Pathogen"

Meta_all_combs$Red_quantif_ng.µL <- as.numeric(Meta_all_combs$Red_quantif_ng.µL)

No_CT_REPs <- toString(Meta_all_combs[which(Meta_all_combs$Red_quantif_ng.µL == No_CT_value),"REP_treat"])

ggplot(Meta_all_combs,
       aes(x = Treatment, y = Red_quantif_ng.µL,group = Treatment,color = Treatment, label = REP_treat)) +
  #geom_jitter(position = position_jitter(seed = 1)) +
  #geom_text(position = position_jitter(seed = 5),fontface = "bold") +
  geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(Red_quantif_ng.µL == No_CT_value, 0.5, 1)))+
  STYLE +
  theme(legend.position = "none") +
labs(title = "Pathogen Quantification Adriano",
subtitle = "2 ants per large colony, 1 per small colony",
 y = " quantification ng/µL",
caption = paste("Threshold cycle (Ct) missing for" , No_CT_REPs) ) #+
#facet_wrap(~ PERIOD) #, labeller = as_labeller(time_of_day,text.add)
}



#check if there are overlapping rows between two datasets using an indicator
# info_test$included_FIRST <- TRUE
# metadata_Selected$included_SECOND <- TRUE
# res <- merge(info_test, metadata_Selected, all=TRUE)




########################################################################################################
### FULL PATHOGEN EXPOSED DATASET! #####################################################################
########################################################################################################

# read csv
##############################
# The file contains the data for the first run (09/08/22) with 2 ants per colony + half of the pathogen colonies (up to the 3/11/22)
DNA_Results <- read.csv(paste(IMMUNITY_DATA,"/02-11-22_Adriano-DNA_Results_Analysis_with_Identities_Path_plus_sampleShams.csv",sep=""),header=T,stringsAsFactors = F, sep=",")
#fix missing labels
DNA_Results$Sample_position <- sub(".*\\-", "", DNA_Results$Code)
DNA_Results$Sample_Plate <- sub("\\-.*", "", DNA_Results$Code)

DNA_Results$Treatment <- RIGHT(DNA_Results$Colony,2)

# Select only pathogen treated colonies
#DNA_Results <- DNA_Results[ which(DNA_Results$Treatment == "Big Pathogen" | DNA_Results$Treatment == "Small Pathogen") , ]
#filter colonies with only handful of individuals
DNA_Results <- DNA_Results %>% 
  group_by(Colony) %>%
  filter(sum(NROW(Colony))>10)

#see which cols have been processed
table(DNA_Results$Colony)


# Rename by name
DNA_Results$Treatment <- as.factor(DNA_Results$Treatment)
levels(DNA_Results$Treatment)[levels(DNA_Results$Treatment)=="BS"] <- "Big Sham"
levels(DNA_Results$Treatment)[levels(DNA_Results$Treatment)=="BP"] <- "Big Pathogen"
levels(DNA_Results$Treatment)[levels(DNA_Results$Treatment)=="SS"] <- "Small Sham"
levels(DNA_Results$Treatment)[levels(DNA_Results$Treatment)=="SP"] <- "Small Pathogen"

DNA_Results$MbruDNA <- as.numeric(DNA_Results$MbruDNA)




### CHECK THE WELLS POSITIONS TO ADD THE ANT IDENTITY
plate_positions_list <- list.files(path= file.path(IMMUNITY_DATA, "QPCR_SAMPLING_TAGSCANNER"),pattern="PLAQUE",full.names = T)

plate_positions <- do.call("rbind", lapply(plate_positions_list, function(x) {
  dat <- read.csv(x, header=TRUE)
  dat$fileName <- tools::file_path_sans_ext(basename(x))
  dat
}))

colnames(plate_positions)[which(colnames(plate_positions)=="Comment")] <- "Sample_position"
#extract col and plaque info
plate_positions$Sample_Plate <- sub("\\_.*", "", plate_positions$fileName)
plate_positions$Sample_Plate <- gsub('[PLAQUE]','',plate_positions$Sample_Plate)

plate_positions$Colony <- sub(".*\\-", "", plate_positions$fileName)


common_col_names <- intersect(names(DNA_Results), names(plate_positions))
DNA_Results_annotated <- merge(DNA_Results, plate_positions, by=common_col_names, all.x=TRUE)

# ## info on missing data
# data_count1 <- aggregate( TagID ~ Colony,                        # Count NA by group
#                           DNA_Results_annotated,
#                          function(x) { sum(is.na(x)) },
#                          na.action = NULL)

metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021_2022-10-12.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
#rename cols to match the DNA results
colnames(metadata)[which(colnames(metadata)=="REP_treat")] <- "Colony"
colnames(metadata)[which(colnames(metadata)=="antID")] <- "AntID"


#Merge info file of plate positions with the DNA results
common_col_names1 <- intersect(names(DNA_Results_annotated), names(metadata))
DNA_Results_annotated <- merge(DNA_Results_annotated, metadata, by=common_col_names1, all.x=TRUE)
#ASSIGN QUEEN LABEL TO QUEEN ISTEAD OF NURSE (SHOULD BE FIXED IN METADATA!!!)
DNA_Results_annotated[which(DNA_Results_annotated$IsQueen==TRUE),"AntTask"] <- "queen"

# DNA_Results <- left_join(DNA_Results, metadata, by = c("REP_treat","antID"))                 # Apply left_join dplyr function

###### filter objects below detection threshold of the qPCR machine.
## threshold calculated as the point of inflection of the values curve (scatterplot of the ordered Cts), getting the values of second derivative. 
# by point 990 (new data only) there is the start of the threshold, corresponding to Ct of 36. Above 36 the detection is non reliable.
DNA_Results_annotated$Ct1 <- as.numeric(DNA_Results_annotated$Ct1)
DNA_Results_annotated$Ct2 <- as.numeric(DNA_Results_annotated$Ct2)

DNA_Results_annotated$Ct_mean <- rowMeans(DNA_Results_annotated[,c("Ct1","Ct2")],na.rm=T)

# plot Ct_mean VS MbDNA (same equation for each plate)
DNA_Results_annotated$qPCR_Plate <- as.factor(DNA_Results_annotated$qPCR_Plate)
ggplot(DNA_Results_annotated,
       aes(x = Ct_mean, y = log(MbruDNA), colour=qPCR_Plate)) + #,group = Exposed,color = Exposed
  #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  geom_point(aes(fill = qPCR_Plate),position = position_jitterdodge(),size=1,alpha=0.4)#+ #,alpha= 0.8,stroke=0
  
plot(DNA_Results_annotated$Ct_mean,DNA_Results_annotated$MbruDNA)

##### NOT USEFUL IF WE RE-DO THE STANDARD CURVE!
# #trying to replicate what we did with florent... No Ct_mean but CT1 or 2?
# plot(order(sort(DNA_Results_annotated$Ct_mean)),sort(DNA_Results_annotated$Ct_mean))
# f_of_x <- splinefun(order(sort(DNA_Results_annotated$Ct_mean)),sort(DNA_Results_annotated$Ct_mean))
# plot(f_of_x(sort(DNA_Results_annotated$Ct_mean), deriv = 2),xlim=c(0,1500),ylim=c(-0.2,0.2),cex=0.1,type="b")
# 
# 
# #-----------------
# x = order(sort(DNA_Results_annotated$Ct_mean))
# y = sort(DNA_Results_annotated$Ct_mean)
# plot(x,y,type="l")
# lo <- loess(y~x)
# xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
# out = predict(lo,xl)
# lines(xl, out, col='red', lwd=2)
# #-----------------



#filter values above detection threshold
DNA_Results_annotated[which(DNA_Results_annotated$Ct_mean >36),"MbruDNA"] <- 0

# #check objects with no values
# No_CT_REPs <- toString(DNA_Results_annotated[which(DNA_Results_annotated$MbruDNA == 0),"Colony"])
# N_ants_NoCt <- length(DNA_Results_annotated[which(DNA_Results_annotated$MbruDNA == 0),"Colony"])

#remove  dead ants
DNA_Results_annotated <- DNA_Results_annotated[which(!is.na(DNA_Results_annotated$AntTask)),]

#extract info from plot for making mean of different colour
p <- ggplot(DNA_Results_annotated,
            aes(x = Exposed, y = MbruDNA, fill=Treatment)) + geom_boxplot()
plotdat <- ggplot_build(p)$data[[1]]



#####################################################################################
##################                 STATS             ################################
#####################################################################################

# maybe low Rsquared not relevant according to this: https://stats.stackexchange.com/questions/509933/is-there-such-a-thing-as-a-too-low-r-squared-when-running-multiple-linear-regres

# inferred_ByAnt$TREATMENT <- as.factor(inferred_ByAnt$TREATMENT)
# inferred_ByAnt$PERIOD <- as.factor(inferred_ByAnt$PERIOD)
# inferred_ByAnt$REP_treat<- as.factor(inferred_ByAnt$REP_treat)
# inferred_ByAnt$Rec_Name<- as.factor(inferred_ByAnt$Rec_Name)

# Relevel Exposed
DNA_Results_annotated$Exposed <- as.factor(DNA_Results_annotated$Exposed)
levels(DNA_Results_annotated$Exposed)[levels(DNA_Results_annotated$Exposed)=="TRUE"] <- "treated"
levels(DNA_Results_annotated$Exposed)[levels(DNA_Results_annotated$Exposed)=="FALSE"] <- "untreated"

# Create new categories!
DNA_Results_annotated$Ant_status <- paste(DNA_Results_annotated$Exposed, DNA_Results_annotated$AntTask)

## TEST DIFFERENCE BETWEEN SIZES INSIDE OF EACH ANTTASK*TREATMENT CONDITION
#create constant to make zeros a small number
constant <- min ( DNA_Results_annotated$MbruDNA[  which(DNA_Results_annotated$MbruDNA!=0)  ]  ,na.rm=T  )/10

m1 <- lmer(log10(MbruDNA+constant) ~ Treatment * Ant_status + (1|Colony), data = DNA_Results_annotated) # the "/" is for the nesting #  + (1|time_of_day) 
test_norm(residuals(m1))
summary(m1)
Anova(m1)
r.squaredGLMM(m1)
tab_model(m1)


### post-hocs
posthoc_Ant_status <- summary(glht( m1, linfct=mcp(Ant_status="Tukey")),test=adjusted("BH"))

### NOTE: No need for post-hocs if there is no significant interaction
# POST-HOCs within Ant_status
posthoc_Treatment <- emmeans(m1, specs = pairwise ~  Treatment | Ant_status) #When the control group is the last group in emmeans we can use trt.vs.ctrlk to get the correct set of comparisons # f1|f2 translates to “compare levels of f1 within each level of f2” 
posthoc_Treatment_summary<-summary(posthoc_Treatment$contrasts)
print(posthoc_Treatment_summary)

# # POST-HOCs for AntTask by treatment
# posthoc_Ant_status <- emmeans(m1, specs = pairwise ~  Ant_status | Treatment) #When the control group is the last group in emmeans we can use trt.vs.ctrlk to get the correct set of comparisons # f1|f2 translates to “compare levels of f1 within each level of f2” 
# posthoc_Ant_status_summary<-summary(posthoc_Ant_status$contrasts)
# print(posthoc_Ant_status_summary)

# add letters to each mean
# model_means_cld <- cld(object = posthoc_Ant_status,
#                        adjust = "Tukey",
#                        Letters = letters,
#                        alpha = 0.05)
model_means_cld <- cld(posthoc_Ant_status)
model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters))); model_means_cld$Ant_status <- row.names(model_means_cld)

# show output
model_means_cld


### TEST DIFFERENCE BETWEEN SIZES ONLY FOR TREATED NURSES
m3 <- lmer(log10(MbruDNA+constant) ~ Treatment + (1|Colony), data = DNA_Results_annotated[which(DNA_Results_annotated$Ant_status=="treated nurse"),]) # the "/" is for the nesting #  + (1|time_of_day) 


### TEST DIFFERENCES BETWEEN EXPOSED AND UNEXPOSED AS COUNTS ("positive)
DNA_Results_annotated$positive<- 1
DNA_Results_annotated[which(DNA_Results_annotated$MbruDNA==0),"positive"] <- 0
m2 <- glmer(positive ~ Treatment * Ant_status + (1|Colony), data = DNA_Results_annotated,family=binomial) # the "/" is for the nesting #  + (1|time_of_day) 
Anova(m2)

compareqqnorm(m1)
#par(mfrow=c(1,2))
plot(m1)
qqnorm(residuals(m1))
qqline(residuals(m1))
hist(residuals(m1))
#anova(m1)


##########################################################################
#### PROPORTION OF CLOSE TO 0 VALUES 




###########################################################################
#### PERMUTATIONS


# Difference in variances instead of means (is there a mixed effect model for that?) To couple with a permutation test to subsample large colonies

number_permutations <- 1000

treated_large <- DNA_Results_annotated[which(DNA_Results_annotated$Treatment=="Big Pathogen"&DNA_Results_annotated$Ant_status=="treated nurse"),]
treated_small <- DNA_Results_annotated[which(DNA_Results_annotated$Treatment=="Small Pathogen"&DNA_Results_annotated$Ant_status=="treated nurse"),]

large_colonies_list <- unique(treated_large$Colony)

number_treated_small <- aggregate (  AntID ~ Treatment + Colony, FUN=length, data=treated_small)

overall_resampled_table <- NULL

pb = txtProgressBar(min = 0, max = number_permutations, initial = 0)  # set progress bar 


for (permutation in 1:number_permutations){
  
  setTxtProgressBar(pb,permutation) 
  
  #####first randomly assign a small size to the large colonies
  sampling_sizes         <- sample ( number_treated_small$AntID,  size = length(number_treated_small$AntID), replace = F  )
  names(sampling_sizes)  <- large_colonies_list
  
  resampled_table <- NULL
  for (colony in large_colonies_list ){
    treated_large_subset <- treated_large[which(treated_large$Colony==colony),]
    sampled_indices      <- sample (  1:nrow(treated_large_subset),size =  sampling_sizes[colony], replace=F  )
    resampled_table      <- rbind(resampled_table,treated_large_subset[sampled_indices,])
  }
  overall_resampled_table <- rbind(overall_resampled_table,data.frame(permutation=permutation,resampled_table))
  
  
  if (permutation == number_permutations) { close(pb) } 
}

standard_deviations_permuted <- aggregate ( log10(MbruDNA + constant) ~ Treatment + Colony + permutation,FUN=sd,data=overall_resampled_table     )
names(standard_deviations_permuted)[which(names(standard_deviations_permuted)=="log10(MbruDNA + constant)")] <- "sd"
mean_standard_deviations_permuted <- aggregate ( sd ~ permutation, FUN=mean, data=standard_deviations_permuted)
hist(mean_standard_deviations_permuted$sd)

observed_standard_deviations <- aggregate ( log10(MbruDNA + constant) ~ Treatment + Colony ,FUN=sd,data=treated_small     )
names(observed_standard_deviations)[which(names(observed_standard_deviations)=="log10(MbruDNA + constant)")] <- "sd"
mean_observed_sd_small <- mean(observed_standard_deviations$sd,na.rm=T)
abline(v=mean_observed_sd_small,col="red")

prop_lower <- length (   which(mean_standard_deviations_permuted$sd< mean_observed_sd_small ))/length(mean_standard_deviations_permuted$sd)
prop_higher<- length (   which(mean_standard_deviations_permuted$sd> mean_observed_sd_small ))/length(mean_standard_deviations_permuted$sd)
prop_more_extreme <- min(prop_lower,prop_higher)
pval <- 2*prop_more_extreme

CI_95 <- quantile(  mean_standard_deviations_permuted$sd, probs=c(0.025,0.975)    )
CI_95
mean_observed_sd_small

##########################################################################
##################### PLOTTING ###########################################
##########################################################################

### HIST OF ANTS THAT RECEIVED A LOAD BY TASK
ggplot(DNA_Results_annotated, aes(x = log(MbruDNA),fill=Exposed)) +
  geom_histogram(position = "identity", bins = 30,alpha=0.7) + facet_wrap(~Treatment)  +
  xlab("LOG MbruDNA")  #+
 # scale_fill_discrete(labels=c('untreated', 'treated'))


# DNA_Results_annotated_MEAN    <- aggregate(MbruDNA ~ Treatment + Exposed , FUN=mean, na.rm=T, na.action=na.pass, DNA_Results_annotated) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count_byAnt"
# DNA_Results_annotated_SE    <- aggregate(MbruDNA ~ Treatment + Exposed , FUN=std.error, na.rm=T, na.action=na.pass, DNA_Results_annotated) ; colnames(DNA_Results_annotated_SE) [match("MbruDNA",colnames(DNA_Results_annotated_SE))] <- "SD_MbruDNA"
# DNA_Results_for_barplot    <-  plyr::join(x=DNA_Results_annotated_MEAN, y=DNA_Results_annotated_SE, type = "full", match = "all")
# 
# ######### barplot
# ggplot(DNA_Results_for_barplot, aes(x=Exposed, y=MbruDNA, fill= Treatment))+
#   #geom_bar(stat='identity',position = position_dodge())+
#   geom_errorbar( aes(x=Exposed,ymin=MbruDNA-SD_MbruDNA, ymax=MbruDNA+SD_MbruDNA),position=position_dodge2(width=0.9, preserve = "single"))+
#   geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
#   STYLE+
#   labs(#title = "Pathogen Quantification Adriano",
#     #subtitle = "All ants",
#     y = "Mean  M. brunneum quantification ng/µL",
#     x="") +
#   scale_x_discrete(labels = c("untreated", "treated")) #+
#   #facet_wrap(~ Exposed)
#   # labs(#title= "Grooming Location",
#   #      subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate")), y = "Mean Freq by ant")+



####### BY Ant_status #####################################################

## BOXPLOT | LOAD VS Ant_status BY TREATMENT
p <- ggplot(DNA_Results_annotated,
       aes(x = Ant_status, y = log(MbruDNA+constant))) + #,group = Exposed,color = Exposed
  #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  geom_boxplot(aes(colour=Treatment),lwd=0.8) + #lwd=0.8
  geom_point(aes(fill = Treatment,colour=Treatment),position = position_jitterdodge(),size=1,alpha=0.4)+ #,alpha= 0.8,stroke=0
  #geom_point() +
  STYLE +
  #theme(legend.position = "none") +
  labs(#title = "Pathogen Quantification Adriano",
    #subtitle = "All ants",
    y = "LOG  M. brunneum quantification ng/µL per ant",
    #x="",
    #caption = paste("Threshold cycle (Ct) missing for", N_ants_NoCt , "ants in cols", No_CT_REPs) 
    )

p + geom_text(data = model_means_cld, aes(x = Ant_status, y = 3,group=Ant_status, label = V1,cex=2,fontface="bold"),show.legend = FALSE)#+ ,position = position_jitterdodge()
  #facet_grid(~ AntTask) +




####### BY LOAD VS DISTANCE TO QUEEN #####################################################

#PLOT and STATS of NETWORK distance to queen (in network: aggregated_distance_to_queen) vs pathogen load to see if big and small cols have different slopes... (see graph in notebook). 
#get average distance by colony removing isqueen =T. 
#Network distance will depend on network size, so we will need to do randomisations to test for significance. plot actual data and re-sampled data (as Yuko ulrich did with 95% conficence interval and mean)


NetworkProp_individual <- read.table(paste(WORKDIR,"/Data/NetworkProp_individual.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
# remove ants with NA exposur status
NetworkProp_individual <- NetworkProp_individual[which(!is.na(NetworkProp_individual$Exposed)),]
#rename vars to match DNA_results
NetworkProp_individual <- NetworkProp_individual %>% 
  rename(
    treatment_condition = treatment,
    treatment = size_treat,
    Colony = REP_treat,
    AntID=antID
  )
# Relevel Exposed
NetworkProp_individual$Exposed <- as.factor(NetworkProp_individual$Exposed)
levels(NetworkProp_individual$Exposed)[levels(NetworkProp_individual$Exposed)=="TRUE"] <- "treated"
levels(NetworkProp_individual$Exposed)[levels(NetworkProp_individual$Exposed)=="FALSE"] <- "untreated"
#ASSIGN QUEEN LABEL TO QUEEN ISTEAD OF NURSE (SHOULD BE FIXED IN METADATA!!!)
NetworkProp_individual[which(NetworkProp_individual$IsQueen==TRUE),"AntTask"] <- "queen"

#keep only colonies that have been subjected to qPCR
NetworkProp_individual <- NetworkProp_individual[NetworkProp_individual$Colony %in% unique(DNA_Results_annotated$Colony), ]

#calculate means by colony for plotting
#for aggregate dist
NetworkProp_ColMean <- NetworkProp_individual %>% group_by(Colony, AntTask, Exposed, treatment) %>% summarize( mean_distance_to_queen = mean(aggregated_distance_to_queen))
# Remove duplicate rows
NetworkProp_ColMean <- NetworkProp_ColMean[!duplicated(NetworkProp_ColMean), ]
#for pathogen load
DNA_Results_ColMean <- DNA_Results_annotated %>% group_by(Colony, AntTask, Exposed, treatment,Treatment, Ant_status) %>% summarize( mean_MbruDNA = mean(MbruDNA))
DNA_Results_ColMean <- DNA_Results_ColMean[!duplicated(DNA_Results_ColMean), ]

#add info from network to DNA_Results_annotated
#Merge network data with the DNA results COLONY MEANS
common_col_names2 <- intersect(names(DNA_Results_ColMean), names(NetworkProp_ColMean))
DNA_Results_ColMeans <- merge(DNA_Results_ColMean, NetworkProp_ColMean, by=common_col_names2, all.x=TRUE)


#calculate means by individual for stats (now in 3h blocks)
#for aggregate dist
NetworkProp_IndMean <- NetworkProp_individual %>% group_by(Colony, AntID, AntTask, Exposed, treatment) %>% summarize( mean_distance_to_queen = mean(aggregated_distance_to_queen))
# Remove duplicate rows
NetworkProp_IndMean <- NetworkProp_IndMean[!duplicated(NetworkProp_IndMean), ]

#Merge network data with the DNA results INDIVIDUAL MEANS
common_col_names3 <- intersect(names(DNA_Results_annotated), names(NetworkProp_IndMean))
DNA_Results_IndMean <- dplyr::left_join(DNA_Results_annotated, NetworkProp_IndMean, by=common_col_names3)
DNA_Results_IndMean <- DNA_Results_IndMean[!duplicated(DNA_Results_IndMean), ]


# scatter Plot the data
ggplot(DNA_Results_ColMeans, aes(x = log10(mean_MbruDNA+constant), y = mean_distance_to_queen, color = Treatment)) + #, group = Treatment
  geom_point(size = 2) +
  STYLE_CONT +
  geom_smooth(data=subset(DNA_Results_ColMeans, Treatment=="Big Pathogen"), method = "lm", formula = y ~ x ) + #+ I(x^2)
  geom_smooth(data=subset(DNA_Results_ColMeans, Treatment=="Small Pathogen"), method = "lm", formula = y ~ x) +
  facet_grid(rows= vars(Ant_status))
  
  
  



#### test correlationbetween LOAD and DISTANCE_to_queen
# remove queens as distance is 0
m.distance <- lmer(log10(MbruDNA+constant) ~ Treatment * mean_distance_to_queen * Ant_status + (1|Colony), data = DNA_Results_IndMean[which(DNA_Results_IndMean$Ant_status!="untreated queen"),]) # the "/" is for the nesting #  + (1|time_of_day) 
#anova(m.distance)

test_norm(residuals(m.distance))
summary(m.distance)
Anova(m.distance) # difference between Type II Wald chisquare tests and Type III Analysis of Variance Table with Satterthwaite's method
r.squaredGLMM(m.distance)
tab_model(m.distance)


### post-hocs
posthoc_Qdist_Ant_status <- summary(glht( m.distance, linfct=mcp(Ant_status="Tukey")),test=adjusted("BH"))
print(cld(posthoc_Qdist_Ant_status))

posthoc_Qdist_treat<- summary(glht(m.distance, linfct=mcp(Treatment="Tukey")), test=adjusted("BH"))
print(cld(posthoc_Qdist_treat))


##### post-hocs for interactions
## following this: https://stats.stackexchange.com/questions/5250/multiple-comparisons-on-a-mixed-effects-model
DNA_Results_IndMean$Treat_Status <- interaction(DNA_Results_IndMean$Treatment, DNA_Results_IndMean$Ant_status)
m.distance_Treat_Status <- lme(log10(MbruDNA+constant) ~ Treat_Status * mean_distance_to_queen, random = 1|Colony, data=DNA_Results_IndMean) 
comp.timegeno <- glht(model, linfct=mcp(Treat_Status="Tukey")) 

#####

# 
# ### NOTE: No need for post-hocs if there is no significant interaction
# # POST-HOCs within Ant_status
# posthoc_Treatment <- emmeans(m.distance, specs = pairwise ~  Treatment | Ant_status) #When the control group is the last group in emmeans we can use trt.vs.ctrlk to get the correct set of comparisons # f1|f2 translates to “compare levels of f1 within each level of f2” 
# posthoc_Treatment_summary<-summary(posthoc_Treatment$contrasts)
# print(posthoc_Treatment_summary)

# # POST-HOCs for AntTask by treatment
# posthoc_Qdist_Ant_status <- emmeans(m.distance, specs = pairwise ~  Ant_status | Treatment) #When the control group is the last group in emmeans we can use trt.vs.ctrlk to get the correct set of comparisons # f1|f2 translates to “compare levels of f1 within each level of f2” 
# posthoc_Qdist_Ant_status_summary<-summary(posthoc_Qdist_Ant_status$contrasts)
# print(posthoc_Qdist_Ant_status_summary)

# add letters to each mean
# model_means_cld <- cld(object = posthoc_Qdist_Ant_status,
#                        adjust = "Tukey",
#                        Letters = letters,
#                        alpha = 0.05)
model_means_cld <- cld(posthoc_Qdist_Ant_status)
model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters))); model_means_cld$Ant_status <- row.names(model_means_cld)

# show output
model_means_cld

# Obtain slopes
#m.interaction$coefficients # only working with lm
# m.lst <- lstrends(m.interaction, "Treatment", var="mean_distance_to_queen")

# Compare slopes
# pairs(m.lst)

# #### test the difference between residual variances
# # you can compare Pearson correlations with a Fisher transform, also called a "Fisher's r-to-z", as follows.  
# library(psych)
# library(data.table)
# iris <- as.data.table(iris)
# # Calculate Pearson's R
# m.correlations <- iris[, cor(Sepal.Length, Petal.Width), by = Species]
# m.correlations
# # Compare R values with Fisher's R to Z
# paired.r(m.correlations[Species=="setosa", V1], m.correlations[Species=="versicolor", V1], 
#          n = iris[Species %in% c("setosa", "versicolor"), .N])
# paired.r(m.correlations[Species=="setosa", V1], m.correlations[Species=="virginica", V1], 
#          n = iris[Species %in% c("setosa", "virginica"), .N])
# paired.r(m.correlations[Species=="virginica", V1], m.correlations[Species=="versicolor", V1], 
#          n = iris[Species %in% c("virginica", "versicolor"), .N]) 








#-------------------------------------------------------------------------

#### test the difference between slopes?
m.interaction <- lm(Sepal.Length ~ Petal.Width*Species, data = iris)
anova(m.interaction)
  
# Obtain slopes
m.interaction$coefficients
m.lst <- lstrends(m.interaction, "Species", var="Petal.Width")

# Compare slopes
pairs(m.lst)
  
#### test the difference between residual variances
# you can compare Pearson correlations with a Fisher transform, also called a "Fisher's r-to-z", as follows.  
  
library(psych)
library(data.table)
iris <- as.data.table(iris)
# Calculate Pearson's R
m.correlations <- iris[, cor(Sepal.Length, Petal.Width), by = Species]
m.correlations
# Compare R values with Fisher's R to Z
paired.r(m.correlations[Species=="setosa", V1], m.correlations[Species=="versicolor", V1], 
         n = iris[Species %in% c("setosa", "versicolor"), .N])
paired.r(m.correlations[Species=="setosa", V1], m.correlations[Species=="virginica", V1], 
         n = iris[Species %in% c("setosa", "virginica"), .N])
paired.r(m.correlations[Species=="virginica", V1], m.correlations[Species=="versicolor", V1], 
         n = iris[Species %in% c("virginica", "versicolor"), .N]) 
  
  

# "We used linear regression to compare the relationship of Sepal Length to Petal Width for each Species. 
# We did not find a significant interaction in the relationships of Sepal Length to Petal Width  ( ANOVA FOR P VAL + m.lst FOR ESTIMATES) for
# I. Setosa (B = 0.9), I. Versicolor (B = 1.4), nor I. Virginica (B = 0.6); F (2, 144) = 1.6, p = 0.19. 
# A Fisher's r-to-z comparison indicated that the Pearson correlation for I. Setosa (r = 0.28) was significantly
# lower (p = 0.02) than I. Versicolor (r = 0.55). Similarly, the correlation for I. Virginica (r = 0.28) was 
# significantly weaker (p = 0.02) than the one observed for I. Versicolor."
  
  








### TO DO'S:
# DO I HAVE MODEL FOR VARIANCE?
# DOES THE PERMUTATION TEST GIVE A SIG DIFFERENCE? ANNOTATE SCRIPT nd save results!


#compare distance to queen VS load SLOPES




# ### keep ant task
# DNA_Results_annotated_MEAN_T    <- aggregate(MbruDNA ~ Treatment + Exposed +AntTask, FUN=mean, na.rm=T, na.action=na.pass, DNA_Results_annotated) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count_byAnt"
# DNA_Results_annotated_SE_T    <- aggregate(MbruDNA ~ Treatment + Exposed +AntTask, FUN=std.error, na.rm=T, na.action=na.pass, DNA_Results_annotated) ; colnames(DNA_Results_annotated_SE_T) [match("MbruDNA",colnames(DNA_Results_annotated_SE_T))] <- "SD_MbruDNA"
# DNA_Results_for_barplot_T    <-  plyr::join(x=DNA_Results_annotated_MEAN_T, y=DNA_Results_annotated_SE_T, type = "full", match = "all")
# 
# ######### barplot
# ggplot(DNA_Results_for_barplot_T, aes(x=Exposed, y=MbruDNA, fill= Treatment))+
#   #geom_bar(stat='identity',position = position_dodge())+
#   geom_errorbar( aes(x=Exposed,ymin=MbruDNA-SD_MbruDNA, ymax=MbruDNA+SD_MbruDNA),position=position_dodge2(width=0.9, preserve = "single"))+
#   geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
#   STYLE+
#   labs(#title = "Pathogen Quantification Adriano",
#     #subtitle = "All ants",
#     y = "Mean  M. brunneum quantification ng/µL",
#     x="") +
#   scale_x_discrete(labels = c("untreated", "treated")) +
#   facet_grid(~ AntTask)#+
# #facet_wrap(~ Exposed)
# # labs(#title= "Grooming Location",
# #      subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate")), y = "Mean Freq by ant")+

