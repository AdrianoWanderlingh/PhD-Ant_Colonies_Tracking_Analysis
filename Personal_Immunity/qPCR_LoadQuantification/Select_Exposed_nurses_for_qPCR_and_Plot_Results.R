# Select Random Exposed nurses to do qpcr quantif
## 2 per colony

######################################### Include the missing reps!!!
library(dplyr)
library(broman)
library(ggplot2)
library(stringr)
library(plotrix)

#functions
RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}

STYLE <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw(),
              scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data"
IMMUNITY_DATA <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Pathogen_Quantification_Data"


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
write.table(Meta_all_combs,file="/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Pathogen Quantification Data/220809-Adriano-MetIS2-colony-checkup_Analysis_with_Identities.txt",append=F,col.names=T,row.names=F,quote=T,sep=",")


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
 y = "Reduced quantification ng/µL",
caption = paste("Threshold cycle (Ct) missing for" , No_CT_REPs) ) #+
#facet_wrap(~ PERIOD) #, labeller = as_labeller(time_of_day,text.add)




#check if there are overlapping rows between two datasets using an indicator
# info_test$included_FIRST <- TRUE
# metadata_Selected$included_SECOND <- TRUE
# res <- merge(info_test, metadata_Selected, all=TRUE)




########################################################################################################
### FULL PATHOGEN EXPOSED DATASET! #####################################################################
########################################################################################################


#### MAKE THIS INDIPENDENT FROM THE ABOVE BIT! MAYBE JUST ADDING AN IF STATEMENT OR SOMETHING ON THE TOP PART

# read csv
##############################
#DESCLAIMER:
# currently (3/11/22) the file contains the data for the first run (09/08/22) with 2 ants per colony + half of the pathogen colonies (up to the 3/11/22)
DNA_Results <- read.csv(paste(IMMUNITY_DATA,"/02-11-22_Adriano-DNA_Results_Analysis_with_Identities_Path_plus_sampleShams.csv",sep=""),header=T,stringsAsFactors = F, sep=",")
#fix missing labels
DNA_Results$Sample_position <- sub(".*\\-", "", DNA_Results$Code)
DNA_Results$Sample_Plate <- sub("\\-.*", "", DNA_Results$Code)

DNA_Results$Treatment <- RIGHT(DNA_Results$Colony,2)

# Rename by name
DNA_Results$Treatment <- as.factor(DNA_Results$Treatment)
levels(DNA_Results$Treatment)[levels(DNA_Results$Treatment)=="BS"] <- "Big Sham"
levels(DNA_Results$Treatment)[levels(DNA_Results$Treatment)=="BP"] <- "Big Pathogen"
levels(DNA_Results$Treatment)[levels(DNA_Results$Treatment)=="SS"] <- "Small Sham"
levels(DNA_Results$Treatment)[levels(DNA_Results$Treatment)=="SP"] <- "Small Pathogen"

DNA_Results$MbruDNA <- as.numeric(DNA_Results$MbruDNA)


# Select only pathogen treated colonies
#DNA_Results <- DNA_Results[ which(DNA_Results$Treatment == "Big Pathogen" | DNA_Results$Treatment == "Small Pathogen") , ]
#filter colonies with only handful of individuals
DNA_Results <- DNA_Results %>% 
                group_by(Colony) %>%
                filter(sum(NROW(Colony))>10)

#see which cols have been processed
table(DNA_Results$Colony)

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

# DNA_Results <- left_join(DNA_Results, metadata, by = c("REP_treat","antID"))                 # Apply left_join dplyr function

#check objects with no values
No_CT_REPs <- toString(DNA_Results_annotated[which(DNA_Results_annotated$MbruDNA == 0),"Colony"])
N_ants_NoCt <- length(DNA_Results_annotated[which(DNA_Results_annotated$MbruDNA == 0),"Colony"])

#### LABELS PLOT
ggplot(DNA_Results_annotated,
       aes(x = Treatment, y = MbruDNA,group = Treatment,color = Treatment, label = Colony)) +
  #geom_jitter(position = position_jitter(seed = 1)) +
  #geom_text(position = position_jitter(seed = 5),fontface = "bold") +
  geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  STYLE +
  theme(legend.position = "none") +
  labs(title = "Pathogen Quantification Adriano",
       subtitle = "All ants",
       y = "Reduced quantification ng/µL",
       caption = paste("Threshold cycle (Ct) missing for", N_ants_NoCt , "ants in cols", No_CT_REPs) ) #+
#facet_wrap(~ PERIOD) #, labeller = as_labeller(time_of_day,text.add)

#remove  dead ants
DNA_Results_annotated <- DNA_Results_annotated[which(!is.na(DNA_Results_annotated$Exposed)),]

#extract info from plot for making mean of different colour
p <- ggplot(DNA_Results_annotated,
            aes(x = Exposed, y = MbruDNA, fill=Treatment)) + geom_boxplot()
plotdat <- ggplot_build(p)$data[[1]]

## BOXPLOT
ggplot(DNA_Results_annotated,
       aes(x = Exposed, y = MbruDNA, colour=Treatment)) + #,group = Exposed,color = Exposed
  #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  geom_boxplot(lwd=0.8) + #lwd=0.8
  geom_point(aes(colour = Treatment),position = position_jitterdodge())+ #,alpha= 0.8,stroke=0
  #geom_point() +
  STYLE +
  #theme(legend.position = "none") +
  labs(#title = "Pathogen Quantification Adriano",
       #subtitle = "All ants",
       y = "Reduced M. brunneum quantification ng/µL per ant",
       x="",
       caption = paste("Threshold cycle (Ct) missing for", N_ants_NoCt , "ants in cols", No_CT_REPs) )+
  scale_x_discrete(labels = c("untreated", "treated")) #+
  #geom_segment(data=plotdat, aes(x=xmin, xend=xmax, y=middle, yend=middle),colour="darkgray",size=1)


DNA_Results_annotated_MEAN    <- aggregate(MbruDNA ~ Treatment + Exposed , FUN=mean, na.rm=T, na.action=na.pass, DNA_Results_annotated) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count_byAnt"
DNA_Results_annotated_SE    <- aggregate(MbruDNA ~ Treatment + Exposed , FUN=std.error, na.rm=T, na.action=na.pass, DNA_Results_annotated) ; colnames(DNA_Results_annotated_SE) [match("MbruDNA",colnames(DNA_Results_annotated_SE))] <- "SD_MbruDNA"
DNA_Results_for_barplot    <-  plyr::join(x=DNA_Results_annotated_MEAN, y=DNA_Results_annotated_SE, type = "full", match = "all")

######### barplot
ggplot(DNA_Results_for_barplot, aes(x=Exposed, y=MbruDNA, fill= Treatment))+
  #geom_bar(stat='identity',position = position_dodge())+
  geom_errorbar( aes(x=Exposed,ymin=MbruDNA-SD_MbruDNA, ymax=MbruDNA+SD_MbruDNA),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  STYLE+
  labs(#title = "Pathogen Quantification Adriano",
    #subtitle = "All ants",
    y = "Mean Reduced M. brunneum quantification ng/µL",
    x="") +
  scale_x_discrete(labels = c("untreated", "treated")) #+
  #facet_wrap(~ Exposed)
  # labs(#title= "Grooming Location",
  #      subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate")), y = "Mean Freq by ant")+



####### nurses vs foragers! #####################################################

## BOXPLOT
ggplot(DNA_Results_annotated,
       aes(x = Exposed, y = MbruDNA, colour=Treatment)) + #,group = Exposed,color = Exposed
  #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  geom_boxplot(lwd=0.8) + #lwd=0.8
  geom_point(aes(colour = Treatment),position = position_jitterdodge())+ #,alpha= 0.8,stroke=0
  #geom_point() +
  STYLE +
  #theme(legend.position = "none") +
  labs(#title = "Pathogen Quantification Adriano",
    #subtitle = "All ants",
    y = "Reduced M. brunneum quantification ng/µL per ant",
    x="",
    caption = paste("Threshold cycle (Ct) missing for", N_ants_NoCt , "ants in cols", No_CT_REPs) )+
  scale_x_discrete(labels = c("untreated", "treated")) +
  facet_grid(~ AntTask)
#geom_segment(data=plotdat, aes(x=xmin, xend=xmax, y=middle, yend=middle),colour="darkgray",size=1)

#-------

### keep ant task
DNA_Results_annotated_MEAN_T    <- aggregate(MbruDNA ~ Treatment + Exposed +AntTask, FUN=mean, na.rm=T, na.action=na.pass, DNA_Results_annotated) #; colnames(Counts_by_Behaviour_CLEAN) [match("Actor",colnames(Counts_by_Behaviour_CLEAN))] <- "Count_byAnt"
DNA_Results_annotated_SE_T    <- aggregate(MbruDNA ~ Treatment + Exposed +AntTask, FUN=std.error, na.rm=T, na.action=na.pass, DNA_Results_annotated) ; colnames(DNA_Results_annotated_SE_T) [match("MbruDNA",colnames(DNA_Results_annotated_SE_T))] <- "SD_MbruDNA"
DNA_Results_for_barplot_T    <-  plyr::join(x=DNA_Results_annotated_MEAN_T, y=DNA_Results_annotated_SE_T, type = "full", match = "all")

######### barplot
ggplot(DNA_Results_for_barplot_T, aes(x=Exposed, y=MbruDNA, fill= Treatment))+
  #geom_bar(stat='identity',position = position_dodge())+
  geom_errorbar( aes(x=Exposed,ymin=MbruDNA-SD_MbruDNA, ymax=MbruDNA+SD_MbruDNA),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  STYLE+
  labs(#title = "Pathogen Quantification Adriano",
    #subtitle = "All ants",
    y = "Mean Reduced M. brunneum quantification ng/µL",
    x="") +
  scale_x_discrete(labels = c("untreated", "treated")) +
  facet_grid(~ AntTask)#+
#facet_wrap(~ Exposed)
# labs(#title= "Grooming Location",
#      subtitle= expression(paste("Mean", phantom(.)%+-%phantom(.), "SE by replicate")), y = "Mean Freq by ant")+


DNA_Results_untreated <- DNA_Results_annotated[which(DNA_Results_annotated$Exposed==FALSE),]
DNA_Results_untreated <- DNA_Results_untreated[which(DNA_Results_untreated$IsAlive==TRUE),]

######### BOXPLOT SUBSET FOR EFFECT ON UNTREATED ANTS ONLY
## BOXPLOT
ggplot(DNA_Results_untreated,
       aes(x = Exposed, y = MbruDNA, colour=Treatment)) + #,group = Exposed,color = Exposed
  #geom_text(position = position_jitter(seed = 5),fontface = "bold",aes(alpha = ifelse(MbruDNA == 0, 0.5, 1)))+
  geom_boxplot(lwd=0.8) + #lwd=0.8
  geom_point(aes(colour = Treatment),position = position_jitterdodge())+ #,alpha= 0.8,stroke=0
  #geom_point() +
  STYLE +
  #theme(legend.position = "none") +
  labs(#title = "Pathogen Quantification Adriano",
    #subtitle = "All ants",
    y = "Reduced M. brunneum quantification ng/µL per ant",
    x ="")+
  scale_x_discrete(labels = c("", "")) +
  facet_grid(~ AntTask)

# 
DNA_Results_for_barplot_untreated    <-  DNA_Results_for_barplot_T[which(DNA_Results_for_barplot_T$Exposed==FALSE),]
  
######### barplot
ggplot(DNA_Results_for_barplot_untreated, aes(x=Exposed, y=MbruDNA, fill= Treatment))+
  #geom_bar(stat='identity',position = position_dodge())+
  geom_errorbar( aes(x=Exposed,ymin=MbruDNA-SD_MbruDNA, ymax=MbruDNA+SD_MbruDNA),position=position_dodge2(width=0.9, preserve = "single"))+
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  STYLE+
  labs(#title = "Pathogen Quantification Adriano",
    #subtitle = "All ants",
    y = "Mean Reduced M. brunneum quantification ng/µL",
    x="") +
  scale_x_discrete(labels = c("", "")) +
  facet_grid(~ AntTask)



#replot these data + stats!

print( ggplot(DNA_Results_annotated, aes(x = MbruDNA)) +
         geom_histogram(position = "identity", bins = 30) + facet_wrap(~Treatment)  +
         xlab("MbruDNA") )  +
  xlim(0.1,15)

#Counts_by_Behaviour_AllCombos1$period = relevel(Counts_by_Behaviour_AllCombos1$period, ref="pre")
print(paste("###################### LMER OF",VAR,"######################"),sep=" ")
#  + (1|time_of_day) is non relevant for both trimmed 4h blocks and full time
m1 <- lmer(sqrt(inferred_ByAnt[,VAR]) ~ PERIOD * TREATMENT + (1|REP_treat/Rec_Name), data = DNA_Results_annotated) # the "/" is for the nesting #  + (1|time_of_day) 

