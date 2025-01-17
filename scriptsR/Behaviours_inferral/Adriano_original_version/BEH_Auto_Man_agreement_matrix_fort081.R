##########################################################################################
############## AUTO-MAN AGREEMENT MATRIX #################################################
##########################################################################################

#script dependant on BEH_MAIN_behaviours_analysis_fort081.R

#For previous versions of this script, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral

print(paste("AUTO-MAN AGREEMENT MATRIX",REPLICATE, PERIOD))

# # pointer to list of all the possible ids pairs ordered 
ant_id1 <- rep(1:35);  ant_id2 <- rep(1:35) ## THIS IS SHIT!
allpairs <- expand.grid(ant_id1,ant_id2)
allpairs$pair <- apply(allpairs[,c("Var1","Var2")],1,function(x){paste(sort(x),collapse = "_") })

## exclude  duplicated pairs!
allpairs <- allpairs[which(!duplicated(allpairs$pair)),] 

ids_pairs <- subset(allpairs,select = "pair")

# # initialize interaction matrix each rows represent a binary array, one for each ids pairs, with 1s on the interactions and 0s elsewhere
int_mat_manual <- matrix(0L, nrow = dim(ids_pairs)[1], ncol = (dim(IF_frames)[1] ))
rownames(int_mat_manual) <- ids_pairs$pair
colnames(int_mat_manual) <- c(IF_frames$frame_num)

int_mat_auto <- int_mat_manual
# rownames(int_mat_auto) <- ids_pairs$pair
# colnames(int_mat_auto) <- c(min(IF_frames$frame_num)-1, IF_frames$frame_num,max(IF_frames$frame_num)+1)

# # Manual
for (i in 1:nrow(summary_MAN_REP_PER))
{
  PAIR <- summary_MAN_REP_PER$pair[i]
  int_mat_manual[PAIR,] <- int_mat_manual[PAIR,] + c(rep(0,( summary_MAN_REP_PER$int_start_frame[i]-1)),
                                                     rep(1,(summary_MAN_REP_PER$int_end_frame[i]) - summary_MAN_REP_PER$int_start_frame[i] + 1),
                                                     rep(0,(length(IF_frames$frame_num) - summary_MAN_REP_PER$int_end_frame[i])))
}


## ADRIANO TO CHECK WHY THESE ARE NOT EQUAL: IT IS RIGHT TO NOT BE EQUAL BECAUSE TO KNOW THE REAL LENGTH OF THE INTERACTION WE NEED TO ADD THE END FRAME AND SO WE ADD A N OF ELEMENTS EQUAL TO THE N OF INTERACTIONS
sum(summary_MAN_REP_PER$interaction_length_secs*8)
sum(int_mat_manual)

# # Auto
for (i in 1:nrow(summary_AUTO_REP_PER))
  {  
  PAIR <- summary_AUTO_REP_PER$pair[i]
  int_mat_auto[PAIR,] <- int_mat_auto[PAIR,] + c(rep(0,( summary_AUTO_REP_PER$int_start_frame[i]-1)),
                                                     rep(1,(summary_AUTO_REP_PER$int_end_frame[i]) - summary_AUTO_REP_PER$int_start_frame[i] + 1),
                                                     rep(0,(length(IF_frames$frame_num) - summary_AUTO_REP_PER$int_end_frame[i])))
  }


sum(interacts_AUTO_REP_PER$interactions$Duration) #see above why they differ
sum(summary_AUTO_REP_PER$int_end_frame - summary_AUTO_REP_PER$int_start_frame )
sum(int_mat_auto)

int_mat_err <- int_mat_manual - int_mat_auto

#sum(int_mat_err==2) should always be 0


#Calculate the % disagreement per each interaction
summary_AUTO_REP_PER$disagreement <- NA
for (i in 1:nrow(summary_AUTO_REP_PER))
  {
  PAIR <- summary_AUTO_REP_PER$pair[i]
  Col_Indices <- summary_AUTO_REP_PER$int_start_frame[i] :  summary_AUTO_REP_PER$int_end_frame[i]
  Row_Indices <- which(PAIR == rownames(int_mat_err))
  Overlap <- int_mat_err[ Row_Indices, Col_Indices]
  summary_AUTO_REP_PER$disagreement[i] <-   sum(Overlap)/length(Overlap)
  } 

## visualise the agrement:
## An interaction with mean =0 has total agreement (0)
## An interaction with mean <0 is a false positive (-1)
## An interaction with mean >0 is false negative (+1)

## explore
#plot(disagreement ~ Duration, summary_AUTO_REP_PER); abline(h=0, lty=2)

## APPLY Disagreement THRESHOLDS
summary_AUTO_REP_PER$Hit <- NA
summary_AUTO_REP_PER$Hit [which(abs(summary_AUTO_REP_PER$disagreement) <=  DISAGREEMENT_THRESH)] <- 1 ## 
summary_AUTO_REP_PER$Hit [which(abs(summary_AUTO_REP_PER$disagreement) >   DISAGREEMENT_THRESH)] <- 0

print(paste("AUTO-MAN AGREEMENT MATRIX",REPLICATE, PERIOD,"PERFORMED"))
