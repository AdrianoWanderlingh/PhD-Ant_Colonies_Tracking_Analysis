
summary_PCA_vars_NONA <- na.omit(as.matrix(summary_AUTO_transf_vars))
image(as.matrix(summary_PCA_vars_NONA))#

summary_PCA_vars_CorMat <- cor(summary_PCA_vars_NONA)
summary_PCA_vars_Dist   <- as.dist((1 - summary_PCA_vars_CorMat)/2) ## correlation matrix, see ?dist


summary_PCA_vars_Hclust <- hclust(summary_PCA_vars_CorMat)
par(mfrow=c(1,2))
image.plot(summary_PCA_vars_CorMat)

plot(summary_PCA_vars_Hclust )

## prc
summary_PCA_vars_PCA <- prcomp(summary_PCA_vars_NONA, scale=T, center=F)
biplot(summary_PCA_vars_PCA, cex=0.00001)


#redundancy analysi hmisc

library(Hmisc)


?redun

redun(~StDev_angle_ACT.orderNorm + StDev_angle_REC.log_x + mean_abs_turnAngle_ACT.orderNorm + mean_abs_turnAngle_REC.yeojohnson + stDev_turnAngle_ACT.orderNorm + stDev_turnAngle_REC.orderNorm + mean_Mov_Orient_delta_angle_ACT.orderNorm + mean_Mov_Orient_delta_angle_REC.orderNorm +
        moved_distance_px_ACT.orderNorm + moved_distance_px_REC.log_x + mean_speed_pxpersec_ACT.orderNorm + mean_speed_pxpersec_REC.orderNorm + mean_accel_pxpersec2_ACT.orderNorm + mean_accel_pxpersec2_REC.orderNorm + mean_jerk_PxPerSec3_ACT.orderNorm + mean_jerk_PxPerSec3_REC.orderNorm + mean_sqrt_err_px_ACT.log_x + mean_sqrt_err_px_REC.orderNorm + chull_area_ACT.orderNorm + chull_area_REC.orderNorm + int_length_secs.log_x + prop_time_undetected_ACT.sqrt_x + prop_time_undetected_REC.sqrt_x + mean_strghtline_dist_px.orderNorm + mean_pair_orient_diff.orderNorm + mean_ang_Velocity.orderNorm,
      
      data = summary_AUTO_transf_vars, r2=.2)



##########
#TEST: the variable, if highly predictable (in this case x5) is always dropped regardless of the size of 
#the quantile turned into NA (QUANT).
#this means that if a variable is filled with NAs on one side of the distribution, this will not alter its droppability


TOTS <- NULL

set.seed(1)
n <- 100

for (QUANT in c(0.1,0.2,0.4,0.6,0.8)) {
  x1 <- runif(n)
  x2 <- runif(n)
  x2a <- runif(n); x2a[25:30] <- NA
  x3 <- x1 + x2 + runif(n)/10
  x4 <- x1 + x2 + x3 + runif(n)/10; #x4[25:50] <- NA
  x5 <- x1 + x2 + x3 + x4 + runif(n)/10
  x5[which(x5 < quantile(x5, QUANT))] <- NA
  # x5 <- factor(sample(c('a','b','c'),n,replace=TRUE)) #categorical variable
  #x6 <- 1*(x5=='a' | x5=='c')
  REDUND <- redun(~x1+x2+x2a++x3+x4+x5, r2=.8)
  #redun(~x1+x2+x3+x4+x5+x6, r2=.8, allcat=TRUE)
  # x5 is no longer redundant but x6 is
  
  TOTS <- rbind(TOTS,"keep"=REDUND$In,"drop"=REDUND$Out, QUANT)
}
