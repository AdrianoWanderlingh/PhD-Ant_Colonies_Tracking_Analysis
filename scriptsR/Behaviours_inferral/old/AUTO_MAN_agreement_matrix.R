library(Matrix)

# # pointer to list of all the possible ids pairs ordered 
ant_id1 <- rep(1:35);  ant_id2 <- rep(1:35) ## THIS IS SHIT!
allpairs <- expand.grid(ant_id1,ant_id2)
allpairs$pair <- apply(allpairs[,c("Var1","Var2")],1,function(x){paste(sort(x),collapse = "_") })

## exclude duplicates; 1-2 == 2-1
# allpairs <- allpairs[which(allpairs$Var1<allpairs$Var2),] ## THIS DOES NOT WORK
## exclude  duplicated pairs!
allpairs <- allpairs[which(!duplicated(allpairs$pair)),] ## THIS DOES WORK

ids_pairs <- subset(allpairs,select = "pair")
# ids_pairs = [id1*10**4 + id2 for id1 in range(1,len(ants_auto)) for id2 in range(id1 + 1,len(ants_auto) + 1)]
# ids_pairs = {k: i for i,k in enumerate(ids_pairs)} 


# # initialize interaction matrix each rows represent a binary array, one for each ids pairs, with 1s on the interactions and 0s elsewhere
int_mat_manual <- matrix(0L, nrow = dim(ids_pairs)[1], ncol = (dim(IF_frames)[1] ))
rownames(int_mat_manual) <- ids_pairs$pair
colnames(int_mat_manual) <- c(IF_frames$frame_num)

int_mat_auto <- int_mat_manual#matrix(0L, nrow = dim(ids_pairs)[1], ncol = (dim(IF_frames)[1] + 2))
# rownames(int_mat_auto) <- ids_pairs$pair
# colnames(int_mat_auto) <- c(min(IF_frames$frame_num)-1, IF_frames$frame_num,max(IF_frames$frame_num)+1)
# int_mat_manual = np.zeros((len(ids_pairs), N_frm + 2), dtype=bool)
# int_mat_auto = np.zeros((len(ids_pairs), N_frm + 2), dtype=bool)


# # Manual
for (i in 1:nrow(summary_MAN_REP_PER))
{
  PAIR <- summary_MAN_REP_PER$pair[i]
  int_mat_manual[PAIR,] <- int_mat_manual[PAIR,] + c(rep(0,( summary_MAN_REP_PER$int_start_frame[i]-1)),
                                                     rep(1,(summary_MAN_REP_PER$int_end_frame[i]) - summary_MAN_REP_PER$int_start_frame[i] + 1),
                                                     rep(0,(length(IF_frames$frame_num) - summary_MAN_REP_PER$int_end_frame[i])))
}


## ADRIANO TO CHECK WHY THESE ARE NOT EQUAL: IT IS RIGHT TO NOT BE EQUAL BECAUSE....
sum(summary_MAN_REP_PER$interaction_length_secs*8)
sum(int_mat_manual)


# for i in fm.Query.ComputeAntInteractions(e_manual,start=start,end=end,maximumGap=fm.Duration(max_gap*10**9))[1]:
#     ids = i.IDs[0]*10**4 + i.IDs[1]
#     int_mat_manual[ids_pairs[ids]] += np.concatenate([np.zeros((1,TimeToFrame[fm.Time.ToTimestamp(i.Start)])), 
#                                                np.ones((1,TimeToFrame[fm.Time.ToTimestamp(i.End)] - TimeToFrame[fm.Time.ToTimestamp(i.Start)] + 1)),
#                                                np.zeros((1,N_frm - TimeToFrame[fm.Time.ToTimestamp(i.End)] + 1))], 1)[0].astype(bool)

# # Auto
for (i in 1:nrow(interacts_AUTO_REP_PER$interactions))
  {  
  PAIR <- interacts_AUTO_REP_PER$interactions$pair[i]
  int_mat_auto[PAIR,] <- int_mat_auto[PAIR,] + c(rep(0,( interacts_AUTO_REP_PER$interactions$int_start_frame[i]-1)),
                                                     rep(1,(interacts_AUTO_REP_PER$interactions$int_end_frame[i]) - interacts_AUTO_REP_PER$interactions$int_start_frame[i] + 1),
                                                     rep(0,(length(IF_frames$frame_num) - interacts_AUTO_REP_PER$interactions$int_end_frame[i])))
  }


sum(interacts_AUTO_REP_PER$interactions$Duration) #already in frames
sum(int_mat_auto)

# int_mat_manual_sparse <- as(int_mat_manual, "dgRMatrix")  # compressed sparse row CSR. CSR is faster for retrieving rows
# int_mat_auto_sparse <- as(int_mat_auto, "dgRMatrix")

## ADRIANO- FIND WHY THERE ARE DUPLICATE ROW NAMES IN BOTH int_mat_auto_sparse & int_mat_man_sparse

int_mat_err <- int_mat_manual - int_mat_auto
# int_mat_err = sparse.csr_matrix(int_mat_manual.astype(int) - int_mat_auto.astype(int))

sum(int_mat_err==2) # why sum(int_mat_err==2) -> 2???


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
plot(disagreement ~ Duration, interacts_AUTO_REP_PER$interactions); abline(h=0, lty=2)

## APPLY THRESHOLDS
THRESH <- 0.5
interacts_AUTO_REP_PER$interactions$Hit <- NA
interacts_AUTO_REP_PER$interactions$Hit [which(abs(interacts_AUTO_REP_PER$interactions$disagreement) <=  THRESH)] <- 1 ## 
interacts_AUTO_REP_PER$interactions$Hit [which(abs(interacts_AUTO_REP_PER$interactions$disagreement) >   THRESH)] <- 0

## visualise the cutoffs
par(mfrow=c(3,1))
hist(log(interacts_AUTO_REP_PER$interactions$disagreement),                                                    main="All disagreements (log)", breaks=seq(-1,1,0.05)); abline(v=0, lty=2)
hist(interacts_AUTO_REP_PER$interactions$disagreement[which(interacts_AUTO_REP_PER$interactions$Hit==1)], main="hits",              breaks=seq(-1,1,0.05), col=2, xlim=c(-1,1)); abline(v=0, lty=2)
hist(interacts_AUTO_REP_PER$interactions$disagreement[which(interacts_AUTO_REP_PER$interactions$Hit==0)], main="misses",            breaks=seq(-1,1,0.05), col=4, xlim=c(-1,1)); abline(v=0, lty=2)

# int_err_per_frame.append([(int_mat_err==1).sum() / N_frm, (int_mat_err==-1).sum() / N_frm])


# 
# 1. Decide an OVERLAP RATE threshold to consider the interaction as a true positive to keep (50%, 75%?).  DONE
# 2. Once the True Positives are assigned calculate the False Positives (all AUTO - True positives) and False Negatives (all MAN - True positives). 
# 3. calculate rates: To calculate them, assign a value of 1 to each row figuring as True Positive, and 0 to all others. 
# - row by row function, every time that an interaction detected automatically is ALSO hand labelled (Decide time overlap percentage, full or partial?) assign a value 0 to column “false positive”. Every time that an interaction detected automatically is NOT hand labelled assign a value 1 to column “false positive”.
# - Every time an interaction is hand-labelled and is NOT detected automatically assign a value 1 to column “false negative” (missed): otherwise 0.

