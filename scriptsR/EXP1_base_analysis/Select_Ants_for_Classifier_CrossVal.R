###########################################################################
###### RANDOMLY SELECT INDIVIDUALS FOR CLASSIFIER CROSS-VALIDATION ########
###########################################################################

## The selected ants should be followed for: 
# 30 mins after exposure (between 2 and 32 mins post-isolation in 15 mins blocks as done with Vasudha)
# only for grooming

DATADIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Data"

Metadata_Exp1 <- read.table(paste(DATADIR,"/Metadata_Exp1_2021.txt",sep=""),header=T,stringsAsFactors = F, sep=",")

#select only exposed nurses
Metadata_Exp1 <- Metadata_Exp1[which(Metadata_Exp1$Exposed==TRUE),]
Metadata_Exp1 <- Metadata_Exp1[which(Metadata_Exp1$IsAlive==TRUE),]

# RANDOMLY 1 ant per colony (only exposed)
# randomly choose only one row in each Replicate
Metadata_Exp1$Chosen <- 0
Metadata_Exp1[-tapply(-seq_along(Metadata_Exp1$REP_treat),Metadata_Exp1$REP_treat, sample, size=1),]$Chosen <- 1
Metadata_Chosen <- Metadata_Exp1[which(Metadata_Exp1$Chosen==1),]

# -  1 col x treatment - SP, BP, SS, BS - (exclude R3SP and R9SP)
Metadata_Chosen <- Metadata_Chosen[!Metadata_Chosen$REP_treat %in% c("R3SP","R9SP"),]
Metadata_Chosen[-tapply(-seq_along(Metadata_Chosen$treatment),Metadata_Chosen$treatment, sample, size=1),]$Chosen <- 2
Metadata_Chosen <- Metadata_Chosen[which(Metadata_Chosen$Chosen==2),]

#Save output
write.table(Metadata_Chosen,file=paste(DATADIR,"/Select_Ants_for_Classifier_CrossVal_",Sys.Date(),".txt",sep=""),append=F,col.names=T,row.names=F,quote=T,sep=",")
