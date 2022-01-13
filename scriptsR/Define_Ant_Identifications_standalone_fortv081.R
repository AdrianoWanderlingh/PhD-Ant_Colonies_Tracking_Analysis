rm(list=ls())

#"https://formicidae-tracker.github.io/myrmidon/docs/latest/api/index.html"

#install.packages("adehabitatHR")
#`install.packages("igraph")

######load necessary libraries
#library(adehabitatHR) ####for home range and trajectory analyses
library(FortMyrmidon) ####R bindings
#library(igraph)       ####for network analysis

#####################################################
###### 2. CREATING ANTS #############################
#####################################################
##When you create ants directly in fort-studio, antIDs will be created in the order in which you orient the ants
##This means that if two users were to orient the same experiment, but did it in a slightly different orders, they would end up with different antID lists where different antID would correspond to different tagID
##It also means that if for any reasons you wanted to replicate your results / reorient your ants from scratch, and wanted to compare results to your previous analysis, you may end up with a different antID list and this would make the comparison difficult
##One way to avoid that is to create the ants from R using a fixed code
###To do so, you first need to create a myrmidon file and fill in the first page (defning the space, author, tag size and comment)
###Then, BEFORE DOING ANYTHING ELSE, run the following code in RStudio to create the ants in a set order
###At this step you may already want to apply a filter so that tags that had very few detections compard to others wll be discarded as false detections / low-quality ants


###For this step, you will need to open the experiment as read/write as you will want to overwrite your "blank" myrmidon file with a file that contains antIDs
setwd("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP11/")
tracking_data <- fmExperimentOpen("R11BP_22-05-21_base.myrmidon")

###next you will want to extract the tag statistics to know how many times each tag was detected using the following function:
tag_statistics <- fmQueryComputeTagStatistics(tracking_data)

###create ants
for ( i in 1:nrow(tag_statistics)) {  #####loop over each tag
  if ( tag_statistics[i,"count"] >= 0.001*max(tag_statistics[,"count"],na.rm=T) ) { ### optional: here I decide to create an antID only if the tag detection rate was more than 1/1000 * the best tag detection rate. You can choose your own criterion
    a <- tracking_data$createAnt(); ###this actually creates an antID, i.e. associates a decimal antID number to that particular tagID
    identification <- tracking_data$addIdentification(a$ID,tag_statistics[i,"tagDecimalValue"],fmTimeSinceEver(),fmTimeForever())
    print(identification)
  }
}
tracking_data$save('/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP3/R3SP_13-03-21_automatically_oriented.myrmidon')

##Check: print identifications
ants <- tracking_data$ants
for (a in ants) {
  printf("Ant %s is identified by:\n", fmFormatAntID(a$ID))
  for (i in a$identifications){
    printf(" * %s\n", capture.output(i))
  }
}


###Sace to newly named myrmidon file: "previousname_AntsCreated"
tracking_data$save("/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP11/R11BP_22-05-21_AntsCreated.myrmidon")

###Advantages of creating ants in R
######1. Repeatability (as explained above)
######2. Useful when you want to identify and display ants with specific properties halfway through the experiment, and you don't have the time to orient them manually (e.g. Adriano's experiment requiring to identify nurses)
###Disadvantages: 
####Does not allow to have several tags pointing to the same ant (retag). 
####BUT this can be altered manually in fort-studio afterwards by deleting the second ant and then adding the second tagID as identification to the first ant

rm(list=c("tracking_data"))
gc()

