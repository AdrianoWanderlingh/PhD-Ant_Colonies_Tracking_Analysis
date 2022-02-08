setwd("C://Users//cf19810//Downloads")
library(plyr)

#annotation_val <- read.csv("/home/cf19810/Dropbox/Ants_behaviour_analysis/Cross_Validation/annotations_subset_25%_VALIDATED.csv", sep = ",")
annotation_val <- read.csv("annotations_subset_25%_VALIDATED.csv", sep = ",")



annotation_val$IsEqual <- ifelse(annotation_val$Behaviour == annotation_val$BEH_AW, "1", NA)

annotation_val[which(is.na(annotation_val$IsEqual)),]

Count_equal    <- aggregate(IsEqual ~ Behaviour,FUN = length, na.action=na.omit, annotation_val)
Count_tot    <- aggregate(BEH_AW ~ Behaviour,FUN = length, na.action=na.omit, annotation_val); colnames(Count_tot) [match("BEH_AW",colnames(Count_tot))] <- "Tot.cases"

Counts_25percent <- plyr::join (x = Count_equal , y=Count_tot, match = "all")


Counts_25percent$prop_correct <- Counts_25percent$IsEqual/Counts_25percent$Tot.cases
