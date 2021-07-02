
trajectory_summary
trajectories_unlist
trajectories_unlist$time <- fmTimeCPtrFromAnySEXP(trajectories_unlist$time)

trajectory_absol_time                 <- function(trajectory){ return (positions$trajectory_summary$start + trajectory$time)}
trajectories_unlist$abs_time <- unlist(lapply(match(positions$trajectory_summary$antID_str,trajectories_unlist$.id),FUN=trajectory_absol_time)) ###the match is once again to ensure we will extract the right information for the right ant
warnings()





###First we would define our own function:
trajectory_duration                    <- function(trajectory){ return (max(trajectory$time,na.rm=T)-min(trajectory$time,na.rm=T))}
###Second let's apply that function to all trajectories and fill in the results
positions$trajectory_summary$duration <- unlist(lapply(positions$trajectories[c(match(positions$trajectory_summary$antID_str,names(positions$trajectories)))],FUN=trajectory_duration)) ###is the only almost working one.....


lapply(list.1, function(x) as.character(data.f$description[match(x, data.f$unit)]))

lapply(positions$trajectories, function(x) trajectory_summary$start[match(x, trajectory_summary$antID_str)])
