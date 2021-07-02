getwd()
rm (list=ls())

library("FortMyrmidon")
packageVersion("FortMyrmidon")

define_ants <- F ###change back to TRUE when you do a new colony

setwd("/media/cf19810/DISK4/ADRIANO/NURSEID/R10SS_13-05-21.0000_nurseid_folder/")
myrmidon_file <- "R10SS_13-05-21.0000_nurseid.myrmidon"

##########DEFINE ANTS#############
if (define_ants){source("Define_Ant_Identifications.R")}

# opens an experiment in read-only mode
e <- fmExperimentOpenReadOnly(myrmidon_file)

# Statistics about a fort::myrmidon::TagID in the experiment.
tagStats <- fmQueryComputeTagStatistics(e)
library(data.table)
setDT(tagStats, keep.rownames = "tagID")


#export nurse list: hexadecimal tagIDs corresponding to nurses in your data frame. 
list <- data.frame(tagID=paste("- ", as.character(tagStats$tagID,"tagID" ),sep=""))


# write.table output syntax to be checked in the tracking room
write.table(list, file = paste(e$getDataInformations()[["details"]][["tdd.URI"]],"foragers_list.txt",sep=""), append = FALSE,
            row.names = FALSE, col.names = FALSE,quote = FALSE)
