merge_trajectory_segments <- function(positions,IF_frames){
  new_summary <- NULL
  new_trajectories <- list()
  
  print(paste("Mergeing trajectory segments for ",length(unique (positions$trajectories_summary$antID)),"ants..."))
  
  for (antID in sort(unique (positions$trajectories_summary$antID))){
     
    ant_indices <- which(positions$trajectories_summary$antID==antID)
    if (length(ant_indices)==1){
      new_summary      <- rbind(new_summary,positions$trajectories_summary[ant_indices,])
      new_trajectories <- c(new_trajectories, list( positions$trajectories[[ant_indices]]      ))
    }else{
      ###the new summary object should contain only one line per ant - and that line should correspond to the first trajectory segment in line
      ant_summary       <- positions$trajectories_summary[ant_indices,]      ### extract all the lines corresponding to that ant in positions$trajectories_summary
      first_segment_idx <- which.min(ant_summary$start)                      ### find which one starts first 
      new_summary      <- rbind(new_summary,ant_summary[first_segment_idx,]) ### copy that line into new_summary 
      
      ###create a single trajectory for that ant 
      ant_trajectory <- NULL ###initialise data frame
      
      ##now loop over the trajectory segments, in correct temporal order order
      for (ant_index in ant_indices[order(ant_summary$start)]){
        traj_segment                               <- positions$trajectories[[ant_index]]
        traj_segment$dt_diff                       <- c(NA,diff(traj_segment$time))
        traj_segment$frame_diff                    <- round(traj_segment$dt_diff*FRAME_RATE) 
        traj_segment[1,"frame"]                    <- IF_frames[match.closest(x=ant_summary[which(ant_index==ant_indices[order(ant_summary$start)]),"start"] , table=as.numeric(IF_frames$time)),"frame_num"]
        traj_segment[2:nrow(traj_segment),"frame"] <- traj_segment[1,"frame"] + cumsum(traj_segment[2:nrow(traj_segment),"frame_diff"])
        traj_segment$time                          <- as.numeric(IF_frames[traj_segment$frame,"time" ]) - ant_summary[first_segment_idx,"start"]
        ant_trajectory       <- rbind(ant_trajectory,traj_segment)
        ###clear trajectory segment to free memory
        positions$trajectories[[ant_index]] <- data.frame()
        rm(list=c("traj_segment"))
        gc()
        mallinfo::malloc.trim(0L)
      }
      ant_trajectory <- ant_trajectory[which(!duplicated(ant_trajectory$frame)),]
      ant_trajectory <-ant_trajectory[which(!names(ant_trajectory)%in%c("dt_diff","frame_diff","frame"))]
      new_trajectories <- c(new_trajectories, list(ant_trajectory))
    }
    
    rm(list=c("ant_indices","ant_summary","first_segment_idx","ant_trajectory","ant_index"))
    gc()
    mallinfo::malloc.trim(0L)
    gc()
    mallinfo::malloc.trim(0L)
    
  }
  
  return(list(trajectories_summary=new_summary,trajectories=new_trajectories))
}

extract_trajectories <- function (e, start ,end , maximumGap ,computeZones=TRUE,showProgress = FALSE, IF_frames){
  ###to avoid crashes, perform queries by chunks of 12 hours
  absolute_end   <- end ###store the absolute end time so we know when to end
  
  ###initialise an object that will store all successive trajectory segments
  trajectories_summary <- data.frame()
  trajectories        <- list()
  
  ###initialise argument for while
  extraction_complete <- F
  
  ###while extraction is not complete, keep querying for 12 hours segments
  while (!extraction_complete){
    ###get IF_frames for start and end
    frame_start <-   IF_frames[ match.closest(x = as.numeric(as.POSIXct(capture.output(print(start)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )), table = as.numeric(IF_frames$time)),"frame_num"]
    frame_end   <-   IF_frames[ match.closest(x = as.numeric(as.POSIXct(capture.output(print(end)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" )), table = as.numeric(IF_frames$time)),"frame_num"]
    
    ###test if remaining time is superior to 12 hours
    if (as.numeric(as.POSIXct(capture.output(print(end)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) ) -   as.numeric(as.POSIXct(capture.output(print(start)), format = "%Y-%m-%dT%H:%M:%OSZ",  origin="1970-01-01", tz="GMT" ) ) > 6 * 60 * 60){
      ###if superior, define new end time as equal to start + 12 hours
      frame_end <- IF_frames[ match.closest(x = IF_frames[which(IF_frames$frame_num==frame_start),"time"]+6*60*60,table = as.numeric(IF_frames$time)),"frame_num"]
      end <- fmTimeCreate(offset=IF_frames[which(IF_frames$frame_num==frame_end) ,"time"   ])
    }else{
      ###otherwise, keep end time and declare extraction to be complete
      extraction_complete <- T
    }

    ###extract trajectory
    print(paste("Querying trajectories between",capture.output(print(start)),"and",capture.output(print(end)),"..."))
    positions <- fmQueryComputeAntTrajectories(e,start = start,end = end,maximumGap = maximumGap,computeZones = computeZones,showProgress = showProgress) #set true to obtain the zone of the ant
    
    ###add to overall objects 
    trajectories_summary <- rbind(trajectories_summary,positions$trajectories_summary)
    trajectories         <- c(trajectories,positions$trajectories)
    
    ###update start and end
    if (!extraction_complete){
      start <- fmTimeCreate(offset=IF_frames[which(IF_frames$frame_num==frame_end-10) ,"time"   ])          ### new start time = end time of last query - a few frames to ensure some overlap to avoid losing data (problems with rounding leading to mismatches between timestamps!)
      end   <- absolute_end ### new end time = absolute_end
     }
    
    ###clear memory
    rm(list=c("positions"))
    gc()
    mallinfo::malloc.trim(0L)
  }
  
  ###create new positions object that contains all the trajectory segments
  positions = list(trajectories_summary=trajectories_summary, "trajectories"=trajectories)
  
  ###clear memory 
  rm(list=c("trajectories_summary","trajectories"))
  gc()
  mallinfo::malloc.trim(0L)
  gc()
  mallinfo::malloc.trim(0L)
  
  ###merge trajectory for each ant so we end up with a single trajectory per ant
  positions <- merge_trajectory_segments (positions,IF_frames)
  return(positions)
  
}