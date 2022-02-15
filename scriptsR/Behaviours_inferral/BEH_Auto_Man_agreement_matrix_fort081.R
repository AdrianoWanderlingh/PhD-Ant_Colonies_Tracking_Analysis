print(paste("AUTO-MAN AGREEMENT MATRIX",REPLICATE, PERIOD))

library(Matrix)

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
for (i in 1:nrow(interacts_AUTO_REP_PER$interactions))
  {  
  PAIR <- interacts_AUTO_REP_PER$interactions$pair[i]
  int_mat_auto[PAIR,] <- int_mat_auto[PAIR,] + c(rep(0,( interacts_AUTO_REP_PER$interactions$int_start_frame[i]-1)),
                                                     rep(1,(interacts_AUTO_REP_PER$interactions$int_end_frame[i]) - interacts_AUTO_REP_PER$interactions$int_start_frame[i] + 1),
                                                     rep(0,(length(IF_frames$frame_num) - interacts_AUTO_REP_PER$interactions$int_end_frame[i])))
  }


sum(interacts_AUTO_REP_PER$interactions$Duration) #see above why they differ
sum(int_mat_auto)

int_mat_err <- int_mat_manual - int_mat_auto

#sum(int_mat_err==2) should always be 0


#Calculate the % disagreement per each interaction
interacts_AUTO_REP_PER$interactions$agreement <- NA
for (i in 1:nrow(interacts_AUTO_REP_PER$interactions))
  {
  PAIR <- interacts_AUTO_REP_PER$interactions$pair[i]
  Col_Indices <- interacts_AUTO_REP_PER$interactions$int_start_frame[i] :  interacts_AUTO_REP_PER$interactions$int_end_frame[i]
  Row_Indices <- which(PAIR == rownames(int_mat_err))
  Overlap <- int_mat_err[ Row_Indices, Col_Indices]
  interacts_AUTO_REP_PER$interactions$disagreement[i] <-   sum(Overlap)/length(Overlap)
  } 

## visualise the agrement:
## An interaction with mean =0 has total agreement (0)
## An interaction with mean <0 is a false positive (-1)
## An interaction with mean >0 is false negative (+1)


## explore
#plot(disagreement ~ Duration, interacts_AUTO_REP_PER$interactions); abline(h=0, lty=2)

## APPLY THRESHOLDS
THRESH <- 0.5
interacts_AUTO_REP_PER$interactions$Hit <- NA
interacts_AUTO_REP_PER$interactions$Hit [which(abs(interacts_AUTO_REP_PER$interactions$disagreement) <=  THRESH)] <- 1 ## 
interacts_AUTO_REP_PER$interactions$Hit [which(abs(interacts_AUTO_REP_PER$interactions$disagreement) >   THRESH)] <- 0


# int_err_per_frame.append([(int_mat_err==1).sum() / N_frm, (int_mat_err==-1).sum() / N_frm])

# 
# 1. Decide an OVERLAP RATE threshold to consider the interaction as a true positive to keep (50%, 75%?).  DONE
# 2. Once the True Positives are assigned calculate the False Positives (all AUTO - True positives) and False Negatives (all MAN - True positives). 
# 3. calculate rates: To calculate them, assign a value of 1 to each row figuring as True Positive, and 0 to all others. 
# - row by row function, every time that an interaction detected automatically is ALSO hand labelled (Decide time overlap percentage, full or partial?) assign a value 0 to column “false positive”. Every time that an interaction detected automatically is NOT hand labelled assign a value 1 to column “false positive”.
# - Every time an interaction is hand-labelled and is NOT detected automatically assign a value 1 to column “false negative” (missed): otherwise 0.

