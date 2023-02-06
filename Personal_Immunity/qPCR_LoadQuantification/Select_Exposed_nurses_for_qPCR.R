library(dplyr)


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
