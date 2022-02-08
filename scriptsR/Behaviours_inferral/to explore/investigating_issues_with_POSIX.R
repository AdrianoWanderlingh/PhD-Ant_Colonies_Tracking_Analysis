
## https://stackoverflow.com/questions/7726034/how-r-formats-posixct-with-fractional-seconds
myformat.POSIXct <- function(x, digits=0) 
{
  x2 <- round(unclass(x), digits)
  attributes(x2) <- attributes(x)
  x <- as.POSIXlt(x2)
  x$sec <- round(x$sec, digits)
  format.POSIXlt(x, paste("%Y-%m-%d %H:%M:%OS",digits,sep=""))
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## extract actor, receiver IDs & start & end times from the hand-annotated data
ACT <- annot_BEH$Actor[ROW]
ENC_TIME_start <- annot_BEH$T_start_UNIX[ROW]
ENC_TIME_stop  <- annot_BEH$T_stop_UNIX[ROW]
Act_Name <- paste("ant",ACT,sep="_")
print(paste("Behaviour:",BEH,"number",ROW,"Actor:",Act_Name,"Receiver:",Rec_Name))
## extract the trajectory for ACT
traj_ACT <-  positions$trajectories[[Act_Name]]
## remove 'time' column as it is confusing - it's not a common time
traj_ACT$time <- NULL
# ## Plot trajectories of both actor & receiver, show on the same panel
# plot  (y ~ x, traj_ACT, pch=".", col=rgb(0,0,1,0.3,1), main=Title, xlim=c(Xmin,Xmax),ylim=c(Ymin,Ymax))
# points(y ~ x, traj_REC, pch=".", col=rgb(1,0,0,0.3,1))
## subset the trajectories of both actor & receiver using the start & end times
traj_ACT <- traj_ACT [ which(traj_ACT$UNIX_time >= ENC_TIME_start & traj_ACT$UNIX_time <= ENC_TIME_stop),]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## For some reason, the 4th decimal place of the UNIX times (thousandths of a sec) for two individuals seen in a given frame are not identical; cbind(format(traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS4"), format(traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS4"))
## so eliminate the 4th-nth decimal place ; retains accuract to a thousandth-of-a-second
# traj_ACT$UNIX_time <- format(traj_ACT$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
# traj_REC$UNIX_time <- format(traj_REC$UNIX_time, "%Y-%m-%d %H:%M:%OS6")
# print(table(substr(x =  format( traj_ACT$UNIX_time, "%OS6"), start = 7, stop = 7)))
# print(table(substr(x =  format( traj_REC$UNIX_time, "%OS6"), start = 7, stop = 7)))

NaRows <- rep(NA,nrow(traj_ACT))
XX <- as.data.frame(matrix(nrow=nrow(traj_ACT), ncol=0))

XX$ORIG             <- traj_ACT$UNIX_time
XX$UNIX_CT          <- as.POSIXct(traj_ACT$UNIX_time , tz = "GMT", origin = "1970-01-01 00:00:00")
XX$UNIX_LT          <- as.POSIXlt(traj_ACT$UNIX_time , tz = "GMT", origin = "1970-01-01 00:00:00")
XX$UNIX_CT_myformat <- as.POSIXct( myformat.POSIXct(XX$UNIX_CT,3), tz = "GMT", origin = "1970-01-01 00:00:00")
XX$LUB_FLOOR        <- lubridate::round_date( traj_ACT$UNIX_time, ".001s")

format( XX[1:5,], "%OS3")
# XX <- as.data.frame(apply(XX,2,as.character))

XX$ORIG <- strptime(x= as.character(XX$ORIG), format="%Y-%m-%d %H:%M:%OS3", tz = "GMT")


