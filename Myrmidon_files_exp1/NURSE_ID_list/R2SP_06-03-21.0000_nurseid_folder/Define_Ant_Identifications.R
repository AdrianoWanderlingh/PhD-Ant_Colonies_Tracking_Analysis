#Open Experiment
tracking_data <- fmExperimentOpen("R2SP_06-03-21.0000_test2.myrmidon") #if it says it's locked, click Session>Restart R and clear objects from workspace including hidden ones
##GENERAL INFO
tracking_data$getDataInformations()

##user_defined_ant_pose.R
ts <- fmQueryComputeTagStatistics(tracking_data)
  
for ( i in 1:nrow(ts)) {
  # optionnaly maybe skip the tag with some criterion
  if ( ts[i,"count"] >= 0.0001*max(ts[,"count"],na.rm=T) ) {
    a <- tracking_data$createAnt();
    tracking_data$addIdentification(tracking_data$cAnts()$summary[nrow(tracking_data$cAnts()$summary),"antID"],ts[i,"tagDecimalValue"],fmTimeInf(),fmTimeInf())
  }
}

tracking_data$save("/media/adriano/DISK3/ADRIANO/NURSEID/done/R2SP_06-03-21.0000_nurseid_folder/R2SP_06-03-21.0000_test2.myrmidon")
rm(list=c("tracking_data"))
gc()