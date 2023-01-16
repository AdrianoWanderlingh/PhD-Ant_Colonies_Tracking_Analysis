
##### CLEAN UP
gc()
mallinfo::malloc.trim(0L)


##### LIBRARIES
library(dplyr)
library(broman)
library(ggplot2)
library(stringr)
library(plotrix)
library(mallinfo)
library(scales) # for colours

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
library(afex)
library(DHARMa) # residual diagnostics (here used for binomial)
#library(rcompanion)
#install_version("rcompanion", "2.3.0")

#to avoid conflicts between lmer and anova results, set the contrasts globally: 
afex::set_sum_contrasts()

#it is necessary to set the contrasts option in R. Because the multi-way ANOVA model is over-parameterised, it is necessary to choose a contrasts setting that sums to zero, otherwise the ANOVA analysis will give incorrect results with respect to the expected hypothesis. (The default contrasts type does not satisfy this requirement.)
options(contrasts = c("contr.sum","contr.poly"))




##### PLOT STYLES

#Create a custom color scale FOR COLONIES + treatments
FullPal <- scales::viridis_pal(option = "D")(20)
myColorsSmall  <- tail(FullPal,5)
myColorsLarge  <- head(FullPal,5) 
Cols_treatments <- c("#440154FF","#FDE725FF") #"#31688EFF"
myColors      <- c(myColorsLarge,myColorsSmall, Cols_treatments)
names(myColors) <- c("R3BP","R5BP","R7BP","R8BP","R12BP","R1SP", "R2SP", "R5SP", "R7SP","R11SP","Big Pathogen","Small Pathogen")
colScale <- scale_colour_manual(name = "Colony",values = myColors,drop=TRUE)
fillScale <- scale_fill_manual(name = "Colony",values = myColors,drop=TRUE)


# ggplot PLOT STYLE
STYLE <- list(colScale, fillScale,
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw(),
              scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

STYLE_CONT <- list(colScale, fillScale,
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw()
)


##### FUNCTIONS

#right function
RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}

#function to test normality of residuals
test_norm <- function(resids){
  print("Testing normality")
  if (length(resids)<=300){
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
    print("below 0.05, the data significantly deviate from a normal distribution")
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

#function to report a model output
output_lmer <- function(model){
  print("------------RESIDUALS NORMALITY------------")
  test_norm(residuals(model))
  print("------------SUMMARY------------")
  print(summary(model))
  print("------------ANOVA------------")
  print(Anova(model))
  print("------------RSQUARED------------")
  print(r.squaredGLMM(model))
  tab_model(model)
}

## extras that could be added to function above
# plot(m2)
# qqnorm(residuals(m1))
# qqline(residuals(m1))
# hist(residuals(m1))
#anova(m1)

# function to perform posthocs
posthoc_list <- list()
compute_posthocs <- function(model){
  #check taht there are significant outputs
  if (length(row.names(Anova(model)[Anova(model)$'Pr(>Chisq)' < 0.05,]))==0){
    print("there are no significant vars.")}else{
  for (SIG.VAR in row.names(Anova(model)[Anova(model)$'Pr(>Chisq)' < 0.05,])) {
    if (grepl(":", SIG.VAR)) {
      warning(paste0(SIG.VAR, "is an interaction, currently this situation is not handled by the function."))
    }else{
    #check if the variable is not numeric . to do so, we need to access the dataframe from the model
    if (!is.numeric(get(gsub("\\[.*","",as.character(model@call)[3]))[,SIG.VAR])) {
    print(paste0("Performing posthocs for the significant var: ",SIG.VAR))
    arg <- list("Tukey")
    names(arg) <- SIG.VAR
    # glht (mcp) : General linear hypotheses Testing (glht) and multiple comparisons for parametric models (mcp)
    cmp <- do.call(mcp, arg)
    posthoc_SIG.VAR <- summary(glht(model, linfct = cmp),test=adjusted("BH"))
    # Set up a compact letter display of all pair-wise comparisons
    model_means_cld <- cld(posthoc_SIG.VAR)
    #create dataframe usable with ggplot geom_text
    model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters)))
    # add column name
    model_means_cld$newcol <- NA
    colnames(model_means_cld)[which(names(model_means_cld) == "newcol")] <- SIG.VAR
    model_means_cld[,SIG.VAR] <- row.names(model_means_cld)
    rownames(model_means_cld) <- NULL
    # add to list
    posthoc_list <- c(posthoc_list, list(model_means_cld))
    names(posthoc_list)[length(posthoc_list)] <- paste(deparse(substitute(model)),SIG.VAR,sep="_")
    print(paste(deparse(substitute(model)),SIG.VAR,sep="_"))
    print(model_means_cld)
  } # SIG.VAR LOOP
    } #check if is an interaction
  } #check if numeric
  print("call posthoc_list to get posthocs")
  } # if significant vars exist
  return(posthoc_list)
}


#function for POOLED standard deviation (to calculate the SE of the difference between means)
# from here: https://rpubs.com/brouwern/SEdiff2means
## Note the formulas squares SD to get variance
var.pooled <- function(N1,N2,SD1,SD2){
  (N1*SD1^2 + N2*SD2^2)/(N1+N2)
}
# Standard error of difference
## Note that this uses sample size, NOT degrees of freedom (N)
SE.diff <- function(var.pool, n1,n2){
  sqrt(var.pool*(1/n1 + 1/n2))
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
  ##### OUTPUTS #####################################################
  
  Florent_output <- read.csv(paste0(IMMUNITY_DATA,"/221219-Adriano-DNA-220809-newplate_Analysis.csv"),header=T,stringsAsFactors = F, sep=",")
  colnames(Florent_output)[which(colnames(Florent_output)=="CT.19.12")] <-   "Ct1"
  colnames(Florent_output)[which(colnames(Florent_output)=="Mbru.DNA.19.12")] <-   "MbruDNA"
  Florent_output$CT.09.08 <- NULL
  Florent_output$Mbru.DNA.09.08 <- NULL
  
  info_ants <- read.table(paste0(IMMUNITY_DATA,"/Select_Exposed_nurses_for_qPCR_8-08-22_ANNOTATED.txt"),header=T,stringsAsFactors = F, sep=",")
  
  #format well name with 0 after letter
  info_ants$Well <- gsub("^(.{1})(.*)$",         # Apply gsub
                         "\\10\\2",
                         info_ants$NEW_position)
  ### merge
  Meta_all_combs <- list(info_ants,Florent_output)
  Meta_all_combs <- Reduce(function(x, y) merge(x, y, all=TRUE), Meta_all_combs)
  
  
  #Meta_all_combs[which(Meta_all_combs$MbruDNA == "#DIV/0!"),] <- NA
  No_CT_value <- 0.000001
  Meta_all_combs[which(Meta_all_combs$MbruDNA == "#DIV/0!"),"NOTE"] <- paste(Meta_all_combs[which(Meta_all_combs$MbruDNA == "#DIV/0!"),"NOTE"] ,"No CT is either no DNA or issue in extraction",sep=".")
  Meta_all_combs[which(Meta_all_combs$MbruDNA == "#DIV/0!"),"MbruDNA"] <- No_CT_value #No CT (either no DNA or issue in extraction)
  
  
  Meta_all_combs <- Meta_all_combs[which(Meta_all_combs$MbruDNA != "No sample"),]
  
  
  Meta_all_combs$Treatment <- RIGHT(Meta_all_combs$REP_treat,2)
  
  #save output!
  write.table(Meta_all_combs,file="/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Pathogen_Quantification_Data/220809-Adriano-MetIS2-colony-checkup_Analysis_with_Identities.txt",append=F,col.names=T,row.names=F,quote=T,sep=",")
  
  
  # Rename by name
  Meta_all_combs$Treatment <- as.factor(Meta_all_combs$Treatment)
  levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="BS"] <- "Big Sham"
  levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="BP"] <- "Big Pathogen"
  levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="SS"] <- "Small Sham"
  levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="SP"] <- "Small Pathogen"
  
  Meta_all_combs$MbruDNA <- as.numeric(Meta_all_combs$MbruDNA)
  
  No_CT_REPs <- toString(Meta_all_combs[which(Meta_all_combs$MbruDNA == No_CT_value),"REP_treat"])
  
  ggplot(Meta_all_combs,
         aes(x = Treatment, y = MbruDNA,group = Treatment,color = Treatment, label = REP_treat)) +
    #geom_jitter(position = position_jitter(seed = 1)) +
    #geom_text(position = position_jitter(seed = 5),fontface = "bold") +
    geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == No_CT_value, 0.5, 1)))+
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




#######################################################################################################
### DATA PREP: FULL PATHOGEN EXPOSED DATASET ##########################################################
#######################################################################################################

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

plot(DNA_Results_annotated$Ct_mean,DNA_Results_annotated$MbruDNA)

ggplot(DNA_Results_annotated,
       aes(x = Ct_mean, y = log(MbruDNA), colour=qPCR_Plate)) + #,group = Exposed,color = Exposed
  #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  geom_point(aes(fill = qPCR_Plate),position = position_jitterdodge(),size=1,alpha=0.4)#+ #,alpha= 0.8,stroke=0
  
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

# #extract info from plot for making mean of different colour
# p <- ggplot(DNA_Results_annotated,
#             aes(x = Exposed, y = MbruDNA, fill=Treatment)) + geom_boxplot()
# plotdat <- ggplot_build(p)$data[[1]]

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


#####################################################################################
##################              STATS & PLOTS        ################################
#####################################################################################

### HIST OF ANTS THAT RECEIVED A LOAD BY TASK
ggplot(DNA_Results_annotated, aes(x = log(MbruDNA),fill=Exposed)) +
  geom_histogram(position = "identity", bins = 30,alpha=0.7) + facet_wrap(~Treatment)  +
  xlab("LOG MbruDNA")  #+
# scale_fill_discrete(labels=c('untreated', 'treated'))


# maybe low Rsquared not relevant according to this: 
# https://stats.stackexchange.com/questions/509933/is-there-such-a-thing-as-a-too-low-r-squared-when-running-multiple-linear-regres

#create constant to make zeros a small number
constant <- min ( DNA_Results_annotated$MbruDNA[  which(DNA_Results_annotated$MbruDNA!=0)  ]  ,na.rm=T  )/10

#####################################################################################
##### TEST DIFFERENCE OF MbruDNA.LOAD BETWEEN SIZES FOR EACH Ant_status (ANTTASK*TREATMENT CONDITION)

m1 <- lmer(log10(MbruDNA+constant) ~ Treatment * Ant_status + (1|Colony), data = DNA_Results_annotated) # the "/" is for the nesting #  + (1|time_of_day) 
output_lmer(m1)
# for reporting: report sig. results from the Anova: " N=1040,  P<0.0001 (LMM, Ant_status)" 

#plot the interaction terms
#interaction.plot(DNA_Results_annotated$Treatment, DNA_Results_annotated$Ant_status, DNA_Results_annotated$MbruDNA, legend = TRUE, xlab = "Treatment", ylab = "MbruDNA", pch = 19)
# produce post-hocs for significant 
##IMPROVE: HOW TO SAVE IN posthoc_list WITHOUT EXPLICITLY CALLING IT? 
posthoc_list <- compute_posthocs(m1)

### NOTE: No need for post-hocs if there is no significant interaction

# ### POST-HOCs
# # within ant-status
# posthoc_Ant_status <- summary(glht( m1, linfct=mcp(Ant_status="Tukey")),test=adjusted("BH"))
# 
# model_means_cld <- cld(posthoc_Ant_status)
# model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters))); model_means_cld$Ant_status <- row.names(model_means_cld)
# 
# # show output
# model_means_cld

############
#### BOXPLOT | MbruDNA VS Ant_status BY TREATMENT
p <- ggplot(DNA_Results_annotated,
            aes(x = Ant_status, y = log(MbruDNA+constant))) + #,group = Exposed,color = Exposed
  #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  geom_boxplot(aes(colour=Treatment),lwd=0.8) + #lwd=0.8
  geom_point(aes(fill = Treatment,colour=Colony),position = position_jitterdodge(),size=1,alpha=0.4)+ #,alpha= 0.8,stroke=0
  #geom_point() +
  STYLE +
  #theme(legend.position = "none") +
  labs(#title = "Pathogen Quantification Adriano",
    #subtitle = "All ants",
    y = "LOG  M. brunneum quantification ng/µL per ant",
    #x="",
    #caption = paste("Threshold cycle (Ct) missing for", N_ants_NoCt , "ants in cols", No_CT_REPs) 
  )+ 
  geom_text(data = posthoc_list$m1_Ant_status, aes(x = Ant_status, y = 3,group=Ant_status, label = V1,cex=2,fontface="bold"),show.legend = FALSE)#+ ,position = position_jitterdodge()


#####################################################################################
##### TEST DIFFERENCE OF MbruDNA.LOAD BETWEEN SIZES ONLY FOR TREATED NURSES (ANTTASK*TREATMENT CONDITION)
m3 <- lmer(log10(MbruDNA+constant) ~ Treatment + (1|Colony), data = DNA_Results_annotated[which(DNA_Results_annotated$Ant_status=="treated nurse"),]) # the "/" is for the nesting #  + (1|time_of_day) 
output_lmer(m3)
posthoc_list <- compute_posthocs(m3)


#####################################################################################
#### PERMUTATION TEST to check for difference in VARIANCE of MbruDNA.LOAD in TREATED NURSES by TREATMENT
### Is the difference in variance depending only on the difference in sample size between Big and Small colonies?
## permutation run in large colonies (~180) subsampled to small cols size (~30)
number_permutations <- 1000

## HERE ONLY FOR TREATED NURSES
for (STATUS in c("treated nurse")) { # "untreated forager","untreated nurse"  ,
  #subset data by colony size
  STATUS_large <- DNA_Results_annotated[which(DNA_Results_annotated$Treatment=="Big Pathogen"&DNA_Results_annotated$Ant_status==STATUS),]
  STATUS_small <- DNA_Results_annotated[which(DNA_Results_annotated$Treatment=="Small Pathogen"&DNA_Results_annotated$Ant_status==STATUS),]
  
  large_colonies_list <- unique(STATUS_large$Colony)
  
  #get the N of treated nurses in the small colonies
  number_STATUS_small <- aggregate (  AntID ~ Treatment + Colony, FUN=length, data=STATUS_small)
  
  #set empty resampled df
  overall_resampled_table <- NULL
  
  # set visual progress bar for permutations
  pb = txtProgressBar(min = 0, max = number_permutations, initial = 0)  # set progress bar
  
  print(paste0("permutations for status: ", STATUS))
  for (permutation in 1:number_permutations){
    setTxtProgressBar(pb,permutation)
    #####first randomly assign a small size to the large colonies
    sampling_sizes         <- sample ( number_STATUS_small$AntID,  size = length(number_STATUS_small$AntID), replace = F  )
    names(sampling_sizes)  <- large_colonies_list
    
    #resample from the large colonies using the randomly assigned size
    resampled_table <- NULL
    for (colony in large_colonies_list ){
      STATUS_large_subset <- STATUS_large[which(STATUS_large$Colony==colony),]
      sampled_indices      <- sample (  1:nrow(STATUS_large_subset),size =  sampling_sizes[colony], replace=F  )
      resampled_table      <- rbind(resampled_table,STATUS_large_subset[sampled_indices,])
    }
    #generate permutation output
    overall_resampled_table <- rbind(overall_resampled_table,data.frame(permutation=permutation,resampled_table))
    # close progress bar
    if (permutation == number_permutations) { close(pb) }
  }
  
  #calculate SD per each permutation
  standard_deviations_permuted <- aggregate ( log10(MbruDNA + constant) ~ Treatment + Colony + permutation,FUN=sd,data=overall_resampled_table     )
  names(standard_deviations_permuted)[which(names(standard_deviations_permuted)=="log10(MbruDNA + constant)")] <- "sd"
  mean_standard_deviations_permuted <- aggregate ( sd ~ permutation, FUN=mean, data=standard_deviations_permuted)
  #produce histogram
  hist(mean_standard_deviations_permuted$sd)
  
  observed_standard_deviations <- aggregate ( log10(MbruDNA + constant) ~ Treatment + Colony ,FUN=sd,data=STATUS_small     )
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
  
  #report results
  print(paste0("N. permutations=",number_permutations," pval=",pval))
  if (unname(CI_95["2.5%"]) < mean_observed_sd_small & mean_observed_sd_small < unname(CI_95["97.5%"])) {
    print("the standard deviation of the SMALL colonies is INSIDE the distribution of permuted mean SD of the BIG colonies, therefore the observed difference in SD is only caused by the difference in sample sizes between Small and Big.")
  }else{
    print("the standard deviation of the SMALL colonies is OUTSIDE the distribution of permuted mean SD of the BIG colonies, therefore the observed difference in SD is REAL.")
  }
}# STATUS LOOP

#####################################################################################
##### TEST DIFFERENCE OF MbruDNA.LOAD BETWEEN SIZES FOR EACH Ant_status (ANTTASK*TREATMENT CONDITION)
##### MbruDNA.LOAD AS BINARY FACTOR (1 = POSITIVE = non-zero load)
DNA_Results_annotated$positive<- 1
DNA_Results_annotated[which(DNA_Results_annotated$MbruDNA==0),"positive"] <- 0
m2 <- glmer(positive ~ Treatment * Ant_status + (1|Colony), data = DNA_Results_annotated,family=binomial) # the "/" is for the nesting #  + (1|time_of_day) 
output_lmer(m2)

# binomial model diagnostics with DHARMa
simulationOutput <- simulateResiduals(fittedModel = m2)
# Main plot function from DHARMa, which gives 
# Left: a qq-plot to detect overall deviations from the expected distribution
# Right: a plot of the residuals against the rank-transformed model predictions
plot(simulationOutput)
# Plotting standardized residuals against the predicted value
plotResiduals(simulationOutput, fitted(m1), xlab = "fitted(m1)")

### POST-HOCs
# within ant-status
posthoc_list <- compute_posthocs(m2)
warning("posthoc_list shows only 1 letter in cld ")

# posthoc_Ant_status2 <- summary(glht(m2, linfct=mcp(Ant_status="Tukey")),test=adjusted("BH"))
# model_means_cld2 <- cld(posthoc_Ant_status2)
# model_means_cld2 <- as.data.frame(t(t(model_means_cld2$mcletters$Letters))); model_means_cld2$Ant_status <- row.names(model_means_cld2)
# # show output
# model_means_cld2

############
#### BARPLOT | PROPORTION OF ZERO LOAD VS Ant_status BY TREATMENT
table(DNA_Results_annotated$positive,DNA_Results_annotated$Treatment)

DNA_Results_ColPropPos <- DNA_Results_annotated %>% group_by(AntTask, Exposed, treatment,Treatment, Ant_status) %>% summarise(count = n(), positive = sum(positive == 1))
DNA_Results_ColPropPos$negative <- DNA_Results_ColPropPos$count - DNA_Results_ColPropPos$positive 
#prop of zero 
DNA_Results_ColPropPos$propZeros <- (DNA_Results_ColPropPos$negative / DNA_Results_ColPropPos$count)*100
#select only relevant cols
DNA_Results_ColPropPos <- DNA_Results_ColPropPos[,c("Treatment","Ant_status","propZeros")]

# #bad way to get table, it works as there are single vals per each condition...
# #table useless? not possible to test percentages...
# DNA_Results_ColPropPosTABLE <- tapply(DNA_Results_ColPropPos$propZeros, list(DNA_Results_ColPropPos$Treatment, DNA_Results_ColPropPos$Ant_status), mean)

ggplot(DNA_Results_ColPropPos, aes(fill=Treatment, y=propZeros/100, x=Ant_status)) + 
  geom_bar(position="dodge", stat="identity") +
  STYLE +
  labs(#title = "Pathogen Quantification Adriano",
    #subtitle = "All ants",
    y = "proportion of close to 0 values",
    #x="",
    #caption = paste("Threshold cycle (Ct) missing for", N_ants_NoCt , "ants in cols", No_CT_REPs) 
  )


#####################################################################################
##### TEST DIFFERENCE OF MbruDNA.LOAD VS DISTANCE TO QUEEN BETWEEN SIZES FOR EACH Ant_status (ANTTASK*TREATMENT CONDITION)

# load network data

#get average distance by colony removing isqueen =T. 
#Network distance will depend on network size, so we will need to do randomisations to test for significance. plot actual data and re-sampled data (as Yuko ulrich did with 95% conficence interval and mean)
NetworkProp_individual <- read.table(paste(WORKDIR,"/Data/NetworkProp_individual.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
# add space Usage info
SpaceUsage <- read.table(paste(WORKDIR,"/Data/Space_Usage.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

#Merge network data with the DNA results INDIVIDUAL MEANS
common_col_names_NetSpace <- intersect(names(NetworkProp_individual), names(SpaceUsage))
NetworkProp_individual <- dplyr::left_join(NetworkProp_individual, SpaceUsage, by=common_col_names_NetSpace)

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
## order PERIOD and time_of _day levels
NetworkProp_individual$period = factor(NetworkProp_individual$period, levels=c("pre","post"))
#NetworkProp_individual$time_of_day = factor(NetworkProp_individual$time_of_day, levels=sort(as.numeric(unique(DNA_Results_IndNet$time_of_day))))
NetworkProp_individual$time_of_day = factor(NetworkProp_individual$time_of_day, levels=as.character(c(12,15,18,21,0,3,6,9)))


#ASSIGN QUEEN LABEL TO QUEEN INSTEAD OF NURSE (SHOULD BE FIXED IN METADATA!!!)
NetworkProp_individual[which(NetworkProp_individual$IsQueen==TRUE),"AntTask"] <- "queen"
#keep only colonies that have been subjected to qPCR
NetworkProp_individual <- NetworkProp_individual[NetworkProp_individual$Colony %in% unique(DNA_Results_annotated$Colony), ]

# #calculate means by colony for plotting
# #for aggregate dist
# NetworkProp_ColMean <- NetworkProp_individual %>% group_by(Colony, AntTask, Exposed, treatment) %>% summarize( mean_distance_to_queen = mean(aggregated_distance_to_queen))
# # Remove duplicate rows
# NetworkProp_ColMean <- NetworkProp_ColMean[!duplicated(NetworkProp_ColMean), ]
# #for pathogen load
# DNA_Results_ColMean <- DNA_Results_annotated %>% group_by(Colony, AntTask, Exposed, treatment,Treatment, Ant_status) %>% summarize( mean_MbruDNA = mean(MbruDNA))
# DNA_Results_ColMean <- DNA_Results_ColMean[!duplicated(DNA_Results_ColMean), ]
# 
# #add info from network to DNA_Results_annotated
# #Merge network data with the DNA results COLONY MEANS
# common_col_names2 <- intersect(names(DNA_Results_ColMean), names(NetworkProp_ColMean))
# DNA_Results_ColMeans <- merge(DNA_Results_ColMean, NetworkProp_ColMean, by=common_col_names2, all.x=TRUE)

##### OVERALL DISTANCE TO QUEEN
#calculate means by individual for stats (now in 3h blocks)
#for aggregate dist
NetworkProp_IndMean <- NetworkProp_individual %>% group_by(period, Colony, AntID, AntTask, Exposed, treatment) %>% summarize( mean_distance_to_queen = mean(aggregated_distance_to_queen),  mean_distance_to_treated = mean(mean_aggregated_distance_to_treated))
# Remove duplicate rows
NetworkProp_IndMean <- NetworkProp_IndMean[!duplicated(NetworkProp_IndMean), ]

#Merge network data with the DNA results INDIVIDUAL MEANS
common_col_names3 <- intersect(names(DNA_Results_annotated), names(NetworkProp_IndMean))
DNA_Results_IndMean <- dplyr::left_join(DNA_Results_annotated, NetworkProp_IndMean, by=common_col_names3)
# DNA_Results_IndMean <- DNA_Results_IndMean[!duplicated(DNA_Results_IndMean), ]

#CREATE BINS for distance to Q and to Treated nurses
for (   col_id in unique(DNA_Results_IndMean$Colony)  ){
  DNA_Results_IndMean[which(DNA_Results_IndMean$Colony==col_id),"MbruDNA_binned_per_colony"] <- as.numeric(  cut(log10(DNA_Results_IndMean[which(DNA_Results_IndMean$Colony==col_id),"MbruDNA"]+constant),breaks=100,ordered_result = T))
  DNA_Results_IndMean[which(DNA_Results_IndMean$Colony==col_id&DNA_Results_IndMean$Ant_status!="untreated queen"),"mean_distance_to_queen_ordered"] <- as.numeric(  cut(DNA_Results_IndMean[which(DNA_Results_IndMean$Colony==col_id&DNA_Results_IndMean$Ant_status!="untreated queen"),"mean_distance_to_queen"],breaks=100,ordered_result = T))
  DNA_Results_IndMean[which(DNA_Results_IndMean$Colony==col_id&DNA_Results_IndMean$Ant_status!="treated nurse"),"mean_distance_to_treated_ordered"] <- as.numeric(  cut(DNA_Results_IndMean[which(DNA_Results_IndMean$Colony==col_id&DNA_Results_IndMean$Ant_status!="treated nurse"),"mean_distance_to_treated"],breaks=100,ordered_result = T))
}


#####################################################################################
##### TEST correlation between MbruDNA.LOAD and DISTANCE_to_queen for PERIOD POST
##### PLOT and STATS of NETWORK distance to queen (in network: aggregated_distance_to_queen) vs pathogen load to see if big and small cols have different slopes... (see graph in notebook). 
# # remove queens as distance is 0 + load only makes sense for post
## data subset NO Q, NO T
DNA_Results_IndMean_NoQ_post <- DNA_Results_IndMean[which(DNA_Results_IndMean$Ant_status!="untreated queen" & DNA_Results_IndMean$period=="post"),]

m.loadVSdistanceQ <- lmer(log10(MbruDNA+constant) ~ Treatment * mean_distance_to_queen * Ant_status + (1|Colony), data = DNA_Results_IndMean_NoQ_post) # the "/" is for the nesting #  + (1|time_of_day) 
output_lmer(m.loadVSdistanceQ)

posthoc_list <- compute_posthocs(m.loadVSdistanceQ)

#####################################################################################
##### TEST correlation between MbruDNA.LOAD and BINNED DISTANCE_to_queen for PERIOD POST

#model
m.loadVSdistanceQ_BIN <- lmer(log10(MbruDNA+constant) ~ Treatment * mean_distance_to_queen_ordered * Ant_status + (1|Colony), data = DNA_Results_IndMean_NoQ_post) # the "/" is for the nesting #  + (1|time_of_day) 
output_lmer(m.loadVSdistanceQ_BIN)

posthoc_list <- compute_posthocs(m.loadVSdistanceQ_BIN)

##### post-hocs for interactions
## following this: https://stats.stackexchange.com/questions/5250/multiple-comparisons-on-a-mixed-effects-model
DNA_Results_IndMean_NoQ_post$Treat_Status <- interaction(DNA_Results_IndMean_NoQ_post$Treatment, DNA_Results_IndMean_NoQ_post$Ant_status)
m.loadVSdistanceQ_Treat_Status <- lmer(log10(MbruDNA+constant) ~ Treat_Status * mean_distance_to_queen_ordered + (1|Colony), data = DNA_Results_IndMean_NoQ_post) 
output_lmer(m.loadVSdistanceQ_Treat_Status)

##perform post-hocs
posthoc_list <- compute_posthocs(m.loadVSdistanceQ_Treat_Status)

# ### post-hocs
# comp.Treat_Status <- glht(m.loadVSdistanceQ_Treat_Status, linfct=mcp(Treat_Status="Tukey"),test=adjusted("BH")) 
# print(cld(comp.Treat_Status))
# # add letters to each mean
# Treat_Status_means_cld <- cld(comp.Treat_Status)
# Treat_Status_means_cld <- as.data.frame(t(t(Treat_Status_means_cld$mcletters$Letters))); Treat_Status_means_cld$Ant_status <- row.names(Treat_Status_means_cld)

#modify to split the interaction terms
posthoc_list$m.loadVSdistanceQ_Treat_Status <- tidyr::separate(data = posthoc_list$m.loadVSdistanceQ_Treat_Status, col = Treat_Status, into = c("Treatment", "Ant_status"), sep = "\\.")

DNA_Results_CLD_means <- DNA_Results_IndMean_NoQ_post %>% group_by(Ant_status, treatment) %>% summarize(mean_distance_to_queen_ordered = mean(mean_distance_to_queen_ordered), MbruDNA = mean(MbruDNA))
DNA_Results_CLD_means <- cbind(posthoc_list$m.loadVSdistanceQ_Treat_Status ,DNA_Results_CLD_means)
# Remove Duplicate Column Names
DNA_Results_CLD_means <- DNA_Results_CLD_means[!duplicated(colnames(DNA_Results_CLD_means))]

## DATA SUBSETS FOR PLOTTING
## data subset NO Q, NO T
DNA_Results_IndMean_NoQ_NoT <- DNA_Results_IndMean[which(DNA_Results_IndMean$Ant_status!="untreated queen"&DNA_Results_IndMean$Ant_status!="treated nurse"),]
## data subset NO Q
DNA_Results_IndMean_NoQ <- DNA_Results_IndMean[which(DNA_Results_IndMean$Ant_status!="untreated queen"),]

for (DISTANCE_TO_Q in c("mean_distance_to_queen","mean_distance_to_queen_ordered" )) {
  warning("DNA_Results_CLD_means now only works for binned data!")
  warning("FIX DNA_Results_CLD_means to include empty lettering for the pre period")
  # scatter Plot of mean_distance_to_queen VS MbruDNA per INDIVIDUAL (regression line per TREATMENT)
 print( ggplot(DNA_Results_IndMean_NoQ_NoT, 
               aes(x = DNA_Results_IndMean_NoQ_NoT[,DISTANCE_TO_Q] , y = log10(MbruDNA+constant))) + #, group = Treatment
    #stat_ellipse(aes(colour=Colony)) +
    geom_point(aes(colour=Colony),size=1, alpha=0.6) +
    geom_smooth(aes(colour=Treatment),data=subset(DNA_Results_IndMean[which(DNA_Results_IndMean$Ant_status!="untreated queen"&DNA_Results_IndMean$Ant_status!="treated nurse"),]), method = "lm", formula = y ~ x ) + #+ I(x^2)
    STYLE_CONT +
    facet_grid(Ant_status ~ period) + 
      labs(subtitle = "for period PRE there is NO LOAD, it is here reported \n to see differences in behavioral patterns",
        caption = "regression line per Treatment",
        x = DISTANCE_TO_Q) +
  geom_text(data = DNA_Results_CLD_means[which(DNA_Results_CLD_means$Ant_status!="treated nurse"),], aes(x = log10(MbruDNA+constant)+0.5, group=Ant_status, label = V1,cex=2,fontface="bold"),show.legend = FALSE)#+ ,position = position_jitterdodge()
 )
  
  # scatter Plot of mean_distance_to_queen VS MbruDNA per INDIVIDUAL (regression line per COLONY)
  print( ggplot(DNA_Results_IndMean_NoQ_NoT, 
         aes(x = DNA_Results_IndMean_NoQ_NoT[,DISTANCE_TO_Q] , y = log10(MbruDNA+constant))) + #, group = Treatment
    geom_point(aes(colour=Colony),size=1,alpha=0.6) +
    geom_smooth(aes(colour=Colony),data=subset(DNA_Results_IndMean[which(DNA_Results_IndMean$Ant_status!="untreated queen"&DNA_Results_IndMean$Ant_status!="treated nurse"),]), method = "lm", formula = y ~ x ) + #+ I(x^2)
    STYLE_CONT +
    facet_grid(Ant_status~ period + Treatment) +
      labs(subtitle = "for period PRE there is NO LOAD, it is here reported \n to see differences in behavioral patterns",
           caption = "regression line per Colony, faceted by Col",
        x = DISTANCE_TO_Q)
  )
  
    ## Treatment:mean_distance_to_queen
    ## there is an obvious difference caused by the different size of the network. permutations should be made BEFORE the mean_distance_to_queen is calculated, in the Network script.
    print( ggplot(data= DNA_Results_IndMean_NoQ,
                  aes(x = Ant_status, y = DNA_Results_IndMean_NoQ[,DISTANCE_TO_Q])) + #,group = Exposed,color = Exposed
             #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
             geom_point(aes(fill = Treatment,colour=Colony),position = position_jitterdodge(),size=1,alpha=0.6)+ #,alpha= 0.8,stroke=0
             geom_boxplot(aes(colour=Treatment),lwd=0.8,alpha = 0.3) + #lwd=0.8
             STYLE +
             facet_grid(. ~ period)+
             labs(#caption = "regression line per Colony",
                  y = DISTANCE_TO_Q)
    )
}



#####################################################################################
##### TEST:
# - correlation between DISTANCE_to_queen and period overall and per time_of_day
# - correlation between time_spent_outside and period overall and per time_of_day
##### 3H BLOCK DIST TO QUEEN
#Merge network data with the DNA results INDIVIDUAL
common_col_names4 <- intersect(names(DNA_Results_annotated), names(NetworkProp_individual))
DNA_Results_IndNet <- dplyr::left_join(DNA_Results_annotated, NetworkProp_individual, by=common_col_names4)
# DNA_Results_IndNet <- DNA_Results_IndNet[!duplicated(DNA_Results_IndNet), ]

#CREATE BINS for distance to Q and to Treated nurses
for (   col_id in unique(DNA_Results_IndNet$Colony)  ){
  DNA_Results_IndNet[which(DNA_Results_IndNet$Colony==col_id),"MbruDNA_binned_per_colony"] <- as.numeric(  cut(log10(DNA_Results_IndNet[which(DNA_Results_IndNet$Colony==col_id),"MbruDNA"]+constant),breaks=100,ordered_result = T))
  DNA_Results_IndNet[which(DNA_Results_IndNet$Colony==col_id&DNA_Results_IndNet$Ant_status!="untreated queen"),"aggregated_distance_to_queen_ordered"] <- as.numeric(  cut(DNA_Results_IndNet[which(DNA_Results_IndNet$Colony==col_id&DNA_Results_IndNet$Ant_status!="untreated queen"),"aggregated_distance_to_queen"],breaks=100,ordered_result = T))
  DNA_Results_IndNet[which(DNA_Results_IndNet$Colony==col_id&DNA_Results_IndNet$Ant_status!="treated nurse"),"mean_distance_to_treated_ordered"] <- as.numeric(  cut(DNA_Results_IndNet[which(DNA_Results_IndNet$Colony==col_id&DNA_Results_IndNet$Ant_status!="treated nurse"),"mean_aggregated_distance_to_treated"],breaks=100,ordered_result = T))
}

#delta pre-post/pre (normalised difference as it is divided by pre period)
### CALCULATING THE RELATIVE PERCENTAGE PERIOD DIFFERENCE FOR:
# - distance_to_queen
# - time spent outside


### TO BE ADDED!!!


#use the functions above to calculate the pooled SD and the N samples!
DNA_Results_IndNet_DELTA <- DNA_Results_IndNet %>%
  group_by(Treatment,time_of_day, Colony, AntID, AntTask, Ant_status, Exposed) %>%
  dplyr::summarise(# delta 
                   Norm_Qdist_postpre_diff = ((aggregated_distance_to_queen[match("post", period)] - aggregated_distance_to_queen[match("pre", period)]) / aggregated_distance_to_queen[match("pre", period)])*100,
                   #Norm_Qdist_ord_postpre_diff = (aggregated_distance_to_queen_ordered[match("post", period)] - aggregated_distance_to_queen_ordered[match("pre", period)]) / aggregated_distance_to_queen_ordered[match("pre", period)]
  )
# convert from class() "grouped_df","tbl_df","tbl"  to  data.frame
DNA_Results_IndNet_DELTA <- as.data.frame(DNA_Results_IndNet_DELTA)

#CREATE BINS for RELATIVE DIFF pre-post/pre distance to Q
for (   col_id in unique(DNA_Results_IndNet_DELTA$Colony)  ){
  DNA_Results_IndNet_DELTA[which(DNA_Results_IndNet_DELTA$Colony==col_id&DNA_Results_IndNet_DELTA$Ant_status!="untreated queen"),"Norm_Qdist_postpre_diff_ordered"] <- as.numeric(  cut(DNA_Results_IndNet_DELTA[which(DNA_Results_IndNet_DELTA$Colony==col_id&DNA_Results_IndNet_DELTA$Ant_status!="untreated queen"),"Norm_Qdist_postpre_diff"],breaks=100,ordered_result = T))
}

# We can show the data 2 ways:
# Norm_Qdist_postpre_diff: is the delta of the non-binned distance, normalisation at the colony level is guaranteed by calculating it as pre-post/pre
# Norm_Qdist_postpre_diff_ordered: the previous measure is binned as performed before. this messes up the zero but the relationships between variables are unchanged, meaning that it is likely unnecessary.


#####################################################################################
##### TEST correlation between Norm_Qdist_postpre_diff and Status
#model
## FIX DISTRIBUTION

hist(DNA_Results_IndNet_DELTA[which(DNA_Results_IndNet_DELTA$Ant_status!="untreated queen"),"Norm_Qdist_postpre_diff"])
# try the following: https://stats.stackexchange.com/questions/338868/fitting-a-distribution-to-skewed-data-with-negative-values

PercDistQ <- lmer(Norm_Qdist_postpre_diff ~ Treatment * Ant_status + (1|Colony), data = DNA_Results_IndNet_DELTA[which(DNA_Results_IndNet_DELTA$Ant_status!="untreated queen"),]) # the "/" is for the nesting #  + (1|time_of_day) 
output_lmer(PercDistQ)

posthoc_list <- compute_posthocs(PercDistQ)



# PREP DATA FOR BARPLOT INSTEAD OF BOXPLOT
###### AGGREGATE VALUES FOR PRE.POST | FOR PLOTS #####
###### AGGREGATE ALL, not time
# Aggregate the data by Colony
DNA_Results_IndNet_summ <- DNA_Results_IndNet %>% 
  group_by(Treatment, Colony, AntTask, Ant_status, Exposed, period) %>% 
  summarise(MEAN_aggregated_distance_to_queen = mean(aggregated_distance_to_queen))

# remove the colony level
#calculate mean and SD by colony
DNA_Results_IndNet_summ <- DNA_Results_IndNet_summ %>% 
  group_by(Treatment, AntTask, Ant_status, Exposed, period) %>% 
  summarise(N_Count_REP = length(MEAN_aggregated_distance_to_queen),
            SD_aggregated_distance_to_queen = sd(MEAN_aggregated_distance_to_queen),
            MEAN_aggregated_distance_to_queen = mean(MEAN_aggregated_distance_to_queen))

## CALCULATING THE DIFFERENCES POST MINUS AFTER REQUIRES SUBTRACTING MEANS, THEREFORE SE HAS TO BE CALCULATED ACCORDINGLY
# https://rpubs.com/brouwern/SEdiff2means

### CALCULATING THE DIFFERENCES POST MINUS AFTER (DELTA)
#use the functions above to calculate the pooled SD and the N samples!
## NOTE: NOT SURE THIS IS THE RIGHT WAY TO CALCULATE THE SE HERE, GIVEN THAT THIS IS A PRE-POST DIFF DIVIDED BY PRE AND NOT JUST A PRE-POST DIFF.  
DNA_Results_IndNet_summ_DELTA <- DNA_Results_IndNet_summ %>%
  group_by(Treatment, Ant_status) %>%
  dplyr::summarise(Mean_Norm_Qdist_postpre_diff = ((MEAN_aggregated_distance_to_queen[match("post", period)] - MEAN_aggregated_distance_to_queen[match("pre", period)]) / MEAN_aggregated_distance_to_queen[match("pre", period)])*100 ,
                   SE_Norm_Qdist_postpre_diff = (sqrt(var.pooled(N_Count_REP[match("post", period)], N_Count_REP[match("pre", period)], SD_aggregated_distance_to_queen[match("post", period)], SD_aggregated_distance_to_queen[match("pre", period)])*(1/N_Count_REP[match("post", period)] + 1/N_Count_REP[match("pre", period)])))*100,
  )
# convert from class() "grouped_df","tbl_df","tbl"  to  data.frame
DNA_Results_IndNet_summ_DELTA <- as.data.frame(DNA_Results_IndNet_summ_DELTA)

#### NOTE: repetition of the above but it is hard to pass an argument like list(c("Treatment", "Ant_Status"),c("time_of_day","Treatment", "Ant_Status")) to group_by()....
# PREP DATA FOR BARPLOT INSTEAD OF BOXPLOT
###### AGGREGATE VALUES FOR PRE.POST | FOR PLOTS #####
###### AGGREGATE TO time_of_day LEVEL 
# Aggregate the data by Colony
DNA_Results_IndNet_Time <- DNA_Results_IndNet %>% 
  group_by(Treatment,time_of_day, Colony, AntTask, Ant_status, Exposed, period) %>% 
  summarise(MEAN_aggregated_distance_to_queen = mean(aggregated_distance_to_queen))

# remove the colony level
#calculate mean and SD by colony
DNA_Results_IndNet_Time <- DNA_Results_IndNet_Time %>% 
  group_by(Treatment,time_of_day, AntTask, Ant_status, Exposed, period) %>% 
  summarise(N_Count_REP = length(MEAN_aggregated_distance_to_queen),
            SD_aggregated_distance_to_queen = sd(MEAN_aggregated_distance_to_queen),
            MEAN_aggregated_distance_to_queen = mean(MEAN_aggregated_distance_to_queen))

## CALCULATING THE DIFFERENCES POST MINUS AFTER REQUIRES SUBTRACTING MEANS, THEREFORE SE HAS TO BE CALCULATED ACCORDINGLY
# https://rpubs.com/brouwern/SEdiff2means

### CALCULATING THE DIFFERENCES POST MINUS AFTER (DELTA)
#use the functions above to calculate the pooled SD and the N samples!
## NOTE: NOT SURE THIS IS THE RIGHT WAY TO CALCULATE THE SE HERE, GIVEN THAT THIS IS A PRE-POST DIFF DIVIDED BY PRE AND NOT JUST A PRE-POST DIFF.  
DNA_Results_IndNet_Time_DELTA <- DNA_Results_IndNet_Time %>%
  group_by(Treatment, time_of_day, Ant_status) %>%
  dplyr::summarise(Mean_Norm_Qdist_postpre_diff = ((MEAN_aggregated_distance_to_queen[match("post", period)] - MEAN_aggregated_distance_to_queen[match("pre", period)]) / MEAN_aggregated_distance_to_queen[match("pre", period)])*100 ,
                   SE_Norm_Qdist_postpre_diff = (sqrt(var.pooled(N_Count_REP[match("post", period)], N_Count_REP[match("pre", period)], SD_aggregated_distance_to_queen[match("post", period)], SD_aggregated_distance_to_queen[match("pre", period)])*(1/N_Count_REP[match("post", period)] + 1/N_Count_REP[match("pre", period)])))*100,
  )
# convert from class() "grouped_df","tbl_df","tbl"  to  data.frame
DNA_Results_IndNet_Time_DELTA <- as.data.frame(DNA_Results_IndNet_Time_DELTA)

## BOXPLOT
## there is an obvious difference caused by the different size of the network. permutations should be made BEFORE the mean_distance_to_queen is calculated, in the Network script.
ggplot(data= DNA_Results_IndNet_DELTA[which(DNA_Results_IndNet_DELTA$Ant_status!="untreated queen"),],
       aes(x = time_of_day, y = Norm_Qdist_postpre_diff)) + #,group = Exposed,color = Exposed
  #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  #geom_point(aes(fill = Treatment,colour=Colony),position = position_jitterdodge(),size=1,alpha=0.6)+ #,alpha= 0.8,stroke=0
  geom_hline(aes(yintercept = 0, colour = "red"))+
  geom_boxplot(aes(colour=Treatment),lwd=0.8,alpha = 0.3) + #lwd=0.8
  STYLE +
  facet_grid(. ~ Ant_status)+
  labs(#caption = "regression line per Colony",
    y = "Norm_Qdist_ord_postpre_diff")

# BARPLOT BY ALL_TIME
ggplot(DNA_Results_IndNet_summ_DELTA[which(DNA_Results_IndNet_summ_DELTA$Ant_status!="untreated queen"),], aes(x=Ant_status, y=Mean_Norm_Qdist_postpre_diff, fill= Treatment))+
  geom_errorbar( aes(x=Ant_status,ymin=Mean_Norm_Qdist_postpre_diff-SE_Norm_Qdist_postpre_diff, ymax=Mean_Norm_Qdist_postpre_diff+SE_Norm_Qdist_postpre_diff),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  STYLE +
  #facet_grid(. ~ Ant_status) +
  labs(caption = "Percentage difference in distance to Q, calculated as ((post-pre)/pre)*100",
       y = "% difference in distance to Q after exposure")



# BARPLOT BY time_of_day
ggplot(DNA_Results_IndNet_Time_DELTA[which(DNA_Results_IndNet_Time_DELTA$Ant_status!="untreated queen"),], aes(x=time_of_day, y=Mean_Norm_Qdist_postpre_diff, fill= Treatment))+
  geom_errorbar( aes(x=time_of_day,ymin=Mean_Norm_Qdist_postpre_diff-SE_Norm_Qdist_postpre_diff, ymax=Mean_Norm_Qdist_postpre_diff+SE_Norm_Qdist_postpre_diff),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  #facet_wrap(~ time_of_day) + #, labeller = as_labeller(time_of_day,text.add)
  STYLE +
  facet_grid(. ~ Ant_status) +
  labs(caption = "Percentage difference in distance to Q, calculated as ((post-pre)/pre)*100",
    y = "% difference in distance to Q after exposure")





#############################################################
#### TEST correlation between DISTANCE_to_treated and DISTANCE_to_queen for PERIOD POST
# # remove queens as distance is 0

##### ADD MODEL. DIFFERNT FROM PREVIOUS AS SHOULD INCLUDE PERIOD IN FACTORS!

# scatter Plot of mean_distance_to_queen VS MbruDNA per INDIVIDUAL (regression line per TREATMENT)
print( ggplot(DNA_Results_IndMean_NoQ_NoT, 
              aes(x = mean_distance_to_queen_ordered , y = mean_distance_to_treated_ordered)) + #, group = Treatment
         #stat_ellipse(aes(colour=Colony)) +
         geom_point(aes(colour=Colony),size=1, alpha=0.4) +
         geom_smooth(aes(colour=Treatment),data=subset(DNA_Results_IndMean[which(DNA_Results_IndMean$Ant_status!="untreated queen"&DNA_Results_IndMean$Ant_status!="treated nurse"),]), method = "lm", formula = y ~ x ) + #+ I(x^2)
         STYLE_CONT +
         facet_grid(Ant_status~ period) #+
         # labs(caption = "regression line per Treatment",
         #      x = DISTANCE_TO_Q)
       #geom_text(data = DNA_Results_CLD_means, aes(x = log10(MbruDNA+constant)+0.5, group=Ant_status, label = V1,cex=2,fontface="bold"),show.legend = FALSE)#+ ,position = position_jitterdodge()
)













# 
# ###########################################################################
# #### PERMUTATION TEST to check for difference in distance to the QUEEN BY ANT_STATUS IN BIG AND SMALL COLONIES
# 
# ## Treatment:mean_distance_to_queen
# ## there is an obvious difference caused by the different size of the network. permutations should be made BEFORE the mean_distance_to_queen is calculated, in the Network script.
# 
# number_permutations <- 1000
# 
# ## HERE ONLY FOR TREATED NURSES, LOOP FOR ALL THREE CATEGORIES!!
# for (STATUS in c("untreated forager","untreated nurse"  ,"treated nurse")) {
#   STATUS_large <- DNA_Results_IndMean[which(DNA_Results_IndMean$Treatment=="Big Pathogen"&DNA_Results_IndMean$Ant_status==STATUS),]
#   STATUS_small <- DNA_Results_IndMean[which(DNA_Results_IndMean$Treatment=="Small Pathogen"&DNA_Results_IndMean$Ant_status==STATUS),]
#   
#   large_colonies_list <- unique(STATUS_large$Colony)
#   
#   number_STATUS_small <- aggregate (  AntID ~ Treatment + Colony, FUN=length, data=STATUS_small)
#   
#   overall_resampled_table <- NULL
#   
#   pb = txtProgressBar(min = 0, max = number_permutations, initial = 0)  # set progress bar 
#   
#   print(paste0("permutations for status: ", STATUS))
#   for (permutation in 1:number_permutations){
#     
#     setTxtProgressBar(pb,permutation) 
#     
#     #####first randomly assign a small size to the large colonies
#     sampling_sizes         <- sample ( number_STATUS_small$AntID,  size = length(number_STATUS_small$AntID), replace = F  )
#     names(sampling_sizes)  <- large_colonies_list
#     
#     resampled_table <- NULL
#     for (colony in large_colonies_list ){
#       STATUS_large_subset <- STATUS_large[which(STATUS_large$Colony==colony),]
#       sampled_indices      <- sample (  1:nrow(STATUS_large_subset),size =  sampling_sizes[colony], replace=F  )
#       resampled_table      <- rbind(resampled_table,STATUS_large_subset[sampled_indices,])
#     }
#     overall_resampled_table <- rbind(overall_resampled_table,data.frame(permutation=permutation,resampled_table))
#     
#     
#     if (permutation == number_permutations) { close(pb) } 
#   }
#   
#   MEANS_permuted <- aggregate ( mean_distance_to_queen ~ Treatment + Colony + permutation,FUN=mean,data=overall_resampled_table     )
#   names(MEANS_permuted)[which(names(MEANS_permuted)=="mean_distance_to_queen")] <- "mean"
#   mean_MEANS_permuted <- aggregate ( mean ~ permutation, FUN=mean, data=MEANS_permuted)
#   hist(mean_MEANS_permuted$mean)
#   
#   observed_MEANS <- aggregate ( mean_distance_to_queen ~ Treatment + Colony ,FUN=mean,data=STATUS_small     )
#   names(observed_MEANS)[which(names(observed_MEANS)=="mean_distance_to_queen")] <- "mean"
#   mean_observed_mean_small <- mean(observed_MEANS$mean,na.rm=T)
#   abline(v=mean_observed_mean_small,col="red")
#   
#   prop_lower <- length (   which(mean_MEANS_permuted$mean< mean_observed_mean_small ))/length(mean_MEANS_permuted$mean)
#   prop_higher<- length (   which(mean_MEANS_permuted$mean> mean_observed_mean_small ))/length(mean_MEANS_permuted$mean)
#   prop_more_extreme <- min(prop_lower,prop_higher)
#   pval <- 2*prop_more_extreme
#   
#   CI_95 <- quantile(  mean_MEANS_permuted$mean, probs=c(0.025,0.975)    )
#   CI_95
#   mean_observed_mean_small
#   
#   #report results
#   print(paste0("permutations for status: ", STATUS))
#   print(paste0("N. permutations=",number_permutations," pval=",pval))
#   if (mean_observed_mean_small > unname( CI_95["2.5%"]) & mean_observed_mean_small < unname(CI_95["97.5%"])) {
#     print("the mean of the SMALL colonies is INSIDE the distribution of permuted MEANS of the BIG colonies, therefore the observed difference in MEANS is only caused by the difference in sample sizes between Small and Big.")
#   }else{
#     print("the mean of the SMALL colonies is OUTSIDE the distribution of permuted MEANS of the BIG colonies, therefore the observed difference in MEANS is REAL.")
#   }
# }# STATUS LOOP





#-------------------------------------------------------------------------

  




### TO DO'S:
# LOAD BY STATUS:
# DO I HAVE MODEL FOR VARIANCE OF LOAD BY STATUS? (no but likely not needed)

# LOAD VS DISTANCE:
# do PERMUTATION TEST, plotting it in a single cool plot if possible, or separate plots as before with Nath (useless as there are no observable differences)
# compare distance to queen VS load SLOPES instead of just means (NO DIFF IN SLOPES...)
# add SIG letters to lineplot X 

# ANNOTATE AND RE-ORDER SCRIPT
# FIX THRESHOLD DEPENDING ON FLORENT SERIAL DILUTIONS









# EXTRAS ##################################################



#-------------------------------- 
### EMMEANS POSTHOCS

# POST-HOCs within Ant_status
# posthoc_Treatment <- emmeans(m1, specs = pairwise ~  Treatment | Ant_status) #When the control group is the last group in emmeans we can use trt.vs.ctrlk to get the correct set of comparisons # f1|f2 translates to “compare levels of f1 within each level of f2” 
# posthoc_Treatment_summary<-summary(posthoc_Treatment$contrasts)
# print(posthoc_Treatment_summary)

# # POST-HOCs for AntTask by treatment
# posthoc_Ant_status <- emmeans(m1, specs = pairwise ~  Ant_status | Treatment) #When the control group is the last group in emmeans we can use trt.vs.ctrlk to get the correct set of comparisons # f1|f2 translates to “compare levels of f1 within each level of f2” 
# posthoc_Ant_status_summary<-summary(posthoc_Ant_status$contrasts)
# print(posthoc_Ant_status_summary)

# add letters to each mean
# model_means_cld <- cld(object = posthoc_Ant_status,
#                        adjust = "Tukey",
#                        Letters = letters,
#                        alpha = 0.05)



#--------------------------------
# 
# ### post-hocs
# posthoc_Qdist_Ant_status <- summary(glht( m.loadVSdistanceQ, linfct=mcp(Ant_status="Tukey")),test=adjusted("BH"))
# print(cld(posthoc_Qdist_Ant_status))
# 
# posthoc_Qdist_treat<- summary(glht(m.loadVSdistanceQ, linfct=mcp(Treatment="Tukey")), test=adjusted("BH"))
# print(cld(posthoc_Qdist_treat))




#--------------------------------


# ### keep ant task
# DNA_Results_annotated_MEAN_T    <- aggregate(MbruDNA ~ Treatment + Exposed +AntTask, FUN=mean, na.rm=T, na.action=na.pass, DNA_Results_annotated) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count_byAnt"
# DNA_Results_annotated_SE_T    <- aggregate(MbruDNA ~ Treatment + Exposed +AntTask, FUN=std.error, na.rm=T, na.action=na.pass, DNA_Results_annotated) ; colnames(DNA_Results_annotated_SE_T) [match("MbruDNA",colnames(DNA_Results_annotated_SE_T))] <- "SD_MbruDNA"
# DNA_Results_for_barplot_T    <-  plyr::join(x=DNA_Results_annotated_MEAN_T, y=DNA_Results_annotated_SE_T, type = "full", match = "all")