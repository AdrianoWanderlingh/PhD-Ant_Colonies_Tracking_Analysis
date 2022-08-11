# Select Random Exposed nurses to do qpcr quantif
## 2 per colony

######################################### Include the missing reps!!!
library(dplyr)
library(broman)


WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis"
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data"

metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

metadata <- metadata[which(metadata$Exposed==TRUE),]
metadata <- metadata[which(metadata$IsAlive==TRUE),]

# keep 2 ants for the big colonies
big <- c("BP","BS")
small <- c("SP","SS")

metadata_big <-
  metadata[which(metadata$treatment==big),] %>% 
  group_by(REP_treat) %>% 
  filter(row_number()==c(1,2))

# keep only 1 ant for the small colonies
metadata_small <-
  metadata[which(metadata$treatment==small),] %>% 
  group_by(REP_treat) %>% 
  filter(row_number()==1)

metadata_Selected <- rbind(metadata_big,metadata_small)

metadata_Selected <- metadata_Selected[,c("REP_treat","antID","tagIDdecimal","Exposed","N_exposed")]

### ADD THE MISSING REPS, R9BS and R9SS

extra_missing <- data.frame(REP_treat=c("R9SS","R9BS","R9BS") , antID= c(9,3,34), tagIDdecimal= c(222,16,86), Exposed= c(TRUE,TRUE,TRUE), N_exposed= c(2,14,14) )
metadata_Selected <- rbind(metadata_Selected,extra_missing)

# convert decimal to hex

metadata_Selected$tagIDHEX <- convert2hex(metadata_Selected$tagIDdecimal)
metadata_Selected$tagIDHEX <- paste0("0x0",metadata_Selected$tagIDHEX)

metadata_Selected <- metadata_Selected[,c(1,6,2,3,4,5)]


#save it! 

write.table(metadata_Selected,file="/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data/Select_Exposed_nurses_for_qPCR_8-08-22.txt",append=F,col.names=T,row.names=F,quote=T,sep=",")



########################################################################
##### PLOT OUTPUTS #####################################################


RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}

STYLE <- list(scale_colour_viridis_d(), scale_fill_viridis_d(),
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw(),
              scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

Florant_output <- read.csv("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data/Pathogen Quantification Data/220809-Adriano-MetIS2-colony-checkup_Analysis.csv",header=T,stringsAsFactors = F, sep=",")
colnames(Florant_output)[which(colnames(Florant_output)=="Reduced.quantification..ng.µL..negative.if.below.detection.threshold.0.001.")] <-   "Red_quantif_ng.µL"

info_ants <- read.table("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data/Pathogen Quantification Data/Select_Exposed_nurses_for_qPCR_8-08-22_ANNOTATED.txt",header=T,stringsAsFactors = F, sep=",")


info_ants$Well <- gsub("^(.{1})(.*)$",         # Apply gsub
                       "\\10\\2",
                       info_ants$NEW_position)
# Print new string

### merge
Meta_all_combs <- list(info_ants,Florant_output)
Meta_all_combs <- Reduce(function(x, y) merge(x, y, all=TRUE), Meta_all_combs)


Meta_all_combs[which(Meta_all_combs$Red_quantif_ng.µL == "#VALUE!"),] <- NA

Meta_all_combs <- Meta_all_combs[which(Meta_all_combs$Red_quantif_ng.µL != "No sample"),]


Meta_all_combs$Treatment <- RIGHT(Meta_all_combs$REP_treat,2)

#save output!

write.table(Meta_all_combs,file="/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data/Pathogen Quantification Data/220809-Adriano-MetIS2-colony-checkup_Analysis_with_Identities.txt",append=F,col.names=T,row.names=F,quote=T,sep=",")


# Rename by name
Meta_all_combs$Treatment <- as.factor(Meta_all_combs$Treatment)
levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="BS"] <- "Big Sham"
levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="BP"] <- "Big Pathogen"
levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="SS"] <- "Small Sham"
levels(Meta_all_combs$Treatment)[levels(Meta_all_combs$Treatment)=="SP"] <- "Small Pathogen"

Meta_all_combs$Red_quantif_ng.µL <- as.numeric(Meta_all_combs$Red_quantif_ng.µL)

ggplot(Meta_all_combs,
       aes(x = Treatment, y = Red_quantif_ng.µL,group = Treatment,color = Treatment, label = REP_treat)) +
  #geom_jitter(position = position_jitter(seed = 1)) +
  geom_text(position = position_jitter(seed = 5),fontface = "bold") +
  STYLE +
  theme(legend.position = "none") +
labs(title = "Pathogen Quantification Adriano",
subtitle = "2 ants per large colony, 1 per small colony",
 y = "Reduced quantification ng/µL") #+
facet_wrap(~ PERIOD) #, labeller = as_labeller(time_of_day,text.add)


