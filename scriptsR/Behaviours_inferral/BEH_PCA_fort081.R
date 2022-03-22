print(paste("PERFORM LDA ",unique(interaction_MANUAL$PERIOD), unique(interaction_MANUAL$REPLICATE)))

###########################################################################################
# FORECAST VERIFICATION
# https://www.swpc.noaa.gov/sites/default/files/images/u30/Forecast%20Verification%20Glossary.pdf
# https://en.wikipedia.org/wiki/Sensitivity_and_specificity

library(FactoMineR)
library(factoextra)
library(missMDA) #PCA with missing values
library(corrplot)
require(MASS)
require(ggplot2)
require(scales)
require(gridExtra)
library(reshape2)
library(ggbeeswarm)
library(GGally) #plot multicollinearity 

#data structure
#Once the THRESH for $disagreement has been established:
# Hit == 1 is a True Positive
# Hit ==0 is a False Positive
# this is done in BEH_Auto_Man_agreement_matrix_fort081.R

cat(
"
########### LDA Assumptions ########### 
\nThe assumptions of discriminant analysis are the same as those for MANOVA. The analysis is quite sensitive to outliers and the size of the smallest group must be larger than the number of predictor variables.
\n\nMultivariate normality: Independent variables are normal for each level of the grouping variable.
\nHomogeneity of variance/covariance (homoscedasticity): Variances among group variables are the same across levels of predictors. Can be tested with Box's M statistic. It has been suggested, however, that linear discriminant analysis be used when covariances are equal, and that quadratic discriminant analysis may be used when covariances are not equal.
\nMulticollinearity: Predictive power can decrease with an increased correlation between predictor variables.
\nIndependence: Participants are assumed to be randomly sampled, and a participant's score on one variable is assumed to be independent of scores on that variable for all other participants.
\n\nIt has been suggested that discriminant analysis is relatively robust to slight violations of these assumptions, and it has also been shown that discriminant analysis may still be reliable when using dichotomous variables (where multivariate normality is often violated).
         
########### RELEVANT STEPS ########### 

1. ##### Analyse all data and split in train and test only for LDA as done here: https://rpubs.com/zhaojhao/84669  #####

2. ######## Check for collinearity ########
\n --However, LDA cannot be used with data which are highly collinear, such as is the case in the present
study. Therefore, in this investigation PCA was performed prior to LDA to remove the effect of collinearity
and reduce the number of variables in the X data matrix.-- \n\nRef:\nMetabolites 2014, 4, 433-452; doi:10.3390/metabo4020433

Also, is collinearity reduced after variables transformation via BestNormalize? NO

remove correlated var following https://stats.stackexchange.com/questions/121131/removing-collinear-variables-for-lda-qda-in-r

3. ######## Take care of missing cases as they badly affect the analysis given the quantity of NAs ######## 
--To evaluate more precisely the effect of missing values imputation on the
accuracy of the classifier we worked only with the relevant variables in each
dataset. This also sped up the imputation process. The relevant features were
selected using the RELIEF, a filter method for feature selection in supervised
classification, see Acuna et al. (2003) for more details.--

--First we considered the four datasets having missing values. Each of them
was passed through a cleaning process where features with more than 30%
of missing values as well as instances with more than 50% of missing values
were eliminated. We have written a program to perform this task that allows
us to change these percentages as we want. This cleaning process is carried
out in order to have minimize the number of imputations needed. After that
is done we apply the four methods to treat missing values and once that
is finished and we have a complete dataset we compute the 10-fold cross-
  validation estimates of the misclassification error--

READ DISCUSSION, SPECIFICALLY ON HEPATITIS DATASET.

The R functions for all the procedures discussed in this paper are available
in www.math.uprm.eduredgar
\n\nRef:\nAcuna, E., Coaquira, F. and Gonzalez, M. (2003). A Comparison of Feature
Selection Procedures for Classifiers Based on Kernel Density Estimation, in
Proceedings of the International Conference on Computer, Communication and
Control Technologies, Orlando, FL:CCCT'03, Vol I, pp. 468-472 ")

#remove ant names/rep/int/etc in pca (keep only vars)
summary_AUTO_vars <- summary_AUTO[, -match(c("REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","disagreement","Hit"), names(summary_AUTO))] 
#CHECK FOR MULTICOLLINEARITY
#ggpairs(summary_AUTO_vars) #CAREFUL! HUGE PLOT, NEEDS VERY LARGE WINDOW TO BE VISUALISED AND TAKES 30 SEC TO BE BUILT

#PROBLEM! TOO MANY NAs
#summary_AUTO_vars_TEST <- summary_AUTO_vars[rowSums(is.na(summary_AUTO_vars)) > 7,]
##prop missing by row (show only top) - potentially droppable rows?
#KEEP ROWS WITH AT LEAST 20% OF NAs
summary_AUTO_vars_NAOmit <- summary_AUTO_vars[(apply(summary_AUTO_vars, 1, function(x) sum(is.na(x)))/ncol(summary_AUTO_vars)) < 0.05,]

results.cor <- cor(summary_AUTO_vars_NAOmit)
corrplot(results.cor, type="upper", tl.cex= .8)

help("cor.mtest")
#Remove vars (or rows) with too many NAs for the cor.mtest to work
res1 <- cor.mtest(results.cor, conf.level = .95, method="pearson", na.action = "na.omit") #kendall


# corrplot(results.cor, type="upper", order="hclust", tl.cex= .8, 
#          p.mat = res1$p, sig.level = 0.05, insig = "blank", method = "number", number.cex = .7)


###############################################################################
###### NORMALISE VARS BEFORE LDA ##############################################
###############################################################################
library(bestNormalize)
###SPEED UP WITH THE PARALLEL PACKAGE?
#empty base
summary_AUTO_transf <- data.frame()[1:nrow(summary_AUTO), ]
summary_AUTO$ACT <- as.factor(summary_AUTO$ACT)
summary_AUTO$REC <- as.factor(summary_AUTO$REC)
summary_AUTO$Hit <- as.factor(summary_AUTO$Hit)
summary_AUTO$disagreement <- as.factor(summary_AUTO$disagreement)

summary_AUTO$int_start_frame <- as.character(summary_AUTO$int_start_frame)
summary_AUTO$int_end_frame <- as.character(summary_AUTO$int_end_frame)

#transform variables if they are not normal
start.time <- Sys.time()
for (variable in names(summary_AUTO)){
  if (is.numeric(summary_AUTO[,variable])) {
    val <- shapiro.test(summary_AUTO[,variable])
    if (unname(val$p.value)<0.05) {
      print(paste("for [",variable, "] the Shapiro-Wilk Test has p < 0.05, the data is not normal. Transform it.",sep=" "))
      #scale variable (as reccomended for data prep for LDA in Metabolites 2014, 4, 433-452; doi:10.3390/metabo4020433)
      scaled_VAR <- scale(summary_AUTO[,variable])
      #find the best transformation
      #This function currently estimates the Yeo-Johnson, the Box Cox  (if the data is positive), the log_10(x+a), the square-root (x+a), the arcsinh and the ordered quantile normalization
      BNobject <- bestNormalize(scaled_VAR,standardize = FALSE) #with standardize = FALSE also no_tranformation is attempted 
      summary_AUTO_transf$var <- BNobject$x.t; names(summary_AUTO_transf)[names(summary_AUTO_transf) == 'var'] <- paste(variable,class(BNobject$chosen_transform)[1],sep=".")
    }else{print(paste("for [",variable,"] the Shapiro-Wilk Test has p > 0.05, the data is normal. Keep it.",sep=" "))
          summary_AUTO_transf[variable]<- scale(summary_AUTO[,variable])
    }}else{print(paste("non numeric attribute. Pasting [",variable,"] in the new dataset",sep = " "))
           summary_AUTO_transf[variable]<- summary_AUTO[,variable]}
}
end.time <- Sys.time()
(time.taken <- end.time - start.time)

#### CORR PLOT AFTER TRANSFROMATION
summary_AUTO_transf_vars <- summary_AUTO_transf[, -match(c("REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","disagreement","Hit"), names(summary_AUTO_transf))] 
summary_AUTO_transf_vars_NAOmit <- summary_AUTO_transf_vars[(apply(summary_AUTO_transf_vars, 1, function(x) sum(is.na(x)))/ncol(summary_AUTO_transf_vars)) < 0.05,]

results.cor_transf <- cor(summary_AUTO_transf_vars_NAOmit)

pdf(file=paste(DATADIR,"CorrPlots","AUTO","Gap",MAX_INTERACTION_GAP, sep = "_"), width=20, height=5)
par(mfrow=c(1,2),mar = c(4.5, 3.8, 1, 1.1))
corrplot(results.cor, type="upper", tl.cex= .4,title = "\nresults.cor",)
corrplot(results.cor_transf, type="upper", tl.cex= .4,title = "\nresults.cor_transf",)
dev.off()
#it seems that the transformation rises the collinearity instead of reducing it... SOLVE



#summary_AUTO_REP_PER_dget <- dget("/home/cf19810/Documents/Ants_behaviour_analysis/Data/summary_AUTO_REP_PER_16feb22.txt") # load file created with dput 
summary_AUTO_transf$Hit <- as.factor(summary_AUTO_transf$Hit)

#check N of missing values per variable
#TRIM <- c("mean_jerk_PxPerSec3_ACT","mean_jerk_PxPerSec3_REC","mean_accel_PxPerSec2_ACT","mean_accel_PxPerSec2_REC","mean_abs_turnAngle_ACT","mean_abs_turnAngle_REC")
#remove ant names/rep/int/etc in pca (keep only vars)
summary_PCA_vars <- summary_AUTO_transf[, -match(c("REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","disagreement","Hit"), names(summary_AUTO_transf))] 
sapply(summary_PCA_vars, function(x) sum(is.na(x)))
sapply(summary_PCA_vars, function(x) sum(is.infinite(x))) 
summary_PCA_vars_hit <- cbind(summary_PCA_vars,Hit=summary_AUTO_transf$Hit)

#summary_PCA_vars_trim <- summary_AUTO_transf[, -match(c(TRIM,"REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","agreement","disagreement","Hit"), names(summary_AUTO_transf))] 
## eventually trimunwanted vars
#sapply(summary_PCA_vars_trim, function(x) sum(is.na(x))) 


########################################################
######### MISSING DATA #################################
########################################################



####################################
###########plotting ################
####################################
PLOTTING_TRANSF_VAR <- FALSE
if (PLOTTING_TRANSF_VAR) {
  print("plotting vars")

#transform to long format
summary_PCA_long <- reshape2::melt(summary_PCA_vars_hit,id.vars=c("Hit")) #explanation on the warning message https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message


###plot divided by variable and Hit for Grooming
par(mfrow=c(2,3), family="serif" , mar = c(0.1, 0.1, 2.2, 0))

summ_vars_plot <- ggplot(summary_PCA_long, aes(value, fill = Hit)) +
  facet_wrap(variable ~ .,scales="free") +
  theme_bw() +
  theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) #,legend.position="bottom",legend.justification='right'
summ_vars_plot + geom_density(alpha = 0.2) +
  labs(title = "Density plot for movement variables by Hit rate")#,
      #subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))


summ_vars_plot_box <- ggplot(summary_PCA_long, aes(y=value,x=Hit, fill = Hit)) +
  facet_wrap(variable ~ .,scales="free") +
  theme_bw() +
  theme(text=element_text(family="serif",size=12), legend.key.size = unit(0.3, 'cm')) #,legend.position="bottom",legend.justification='right'

summ_vars_plot_box + geom_boxplot(alpha = 0.5)

summ_vars_plot_box + geom_violin(alpha = 0.5) #+ geom_beeswarm() #careful with swarm

}

# #############################
# ####### PCA #################
# #############################
# 
# #PCA with missing values according to Dray & Josse (2015) https://doi.org/10.1007/s11258-014-0406-z
# #imputePCA of the R package missMDA
# #Impute the missing values of a dataset with the Principal Components Analysis model. 
# #Can be used as a preliminary step before performing a PCA on an incomplete dataset.
# #missing values are replaced by random values, and then PCA is applied on the completed data set, and missing values are then updated by the fitted values
# ## Imputation for summary_PCA_vars which contains NAs
# res.comp <- imputePCA(summary_PCA_vars,method="Regularized")
# 
# ## 1. A PCA can be performed on the imputed data for summary_PCA_vars
# res.pca1 <- PCA(res.comp$completeObs)
# 
# ## 2. PCA on the base data for summary_PCA_vars_trim
# RES.PCA <- PCA(summary_PCA_vars_trim) #more explained variance without the trimmeed vars
# 
# #for (RES.PCA in c(res.pca1,res.pca2)) {
# 
# # Extract eigenvalues/variances
# get_eig(res.pca2)
# # Visualize eigenvalues/variances
# fviz_screeplot(RES.PCA, addlabels = TRUE, ylim = c(0, 50))
# # Extract the results for variables
# var <- get_pca_var(RES.PCA)
# var
# 
# corrplot(var$cos2, is.corr = FALSE)
# 
# # Coordinates of variables
# head(var$coord)
# # Contribution of variables
# head(var$contrib)
# # Graph of variables
# # Control variable colors using their contributions
# fviz_pca_var(RES.PCA, col.var="contrib",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE # Avoid text overlapping
# )
# # Contributions of variables to PC1
# fviz_contrib(RES.PCA, choice = "var", axes = 1, top = 10)
# # Contributions of variables to PC2
# fviz_contrib(RES.PCA, choice = "var", axes = 2, top = 10)
# # Extract the results for individuals #GIVE ROWNAME TO DO THIS (I.E. HIT OR INT)
# ind <- get_pca_ind(RES.PCA)
# ind
# # Coordinates of individuals
# head(ind$coord)
# 
# 
# 
# # Use habillage to specify groups for coloring
# fviz_pca_ind(RES.PCA, pointsize =summary_AUTO_transf$disagreement,
#              label = "none", # hide individual labels
#              title = "PCA - indivdual interactions \nPointsize by disagreement rate",
#              habillage = as.factor(summary_AUTO_transf$Hit), # color by groups
#              addEllipses = TRUE
# )
# 
# #}
# 
# #LDA: Linear discriminant analysis also called DFA D.function.A.
# #LDA focuses on finding a feature subspace that maximizes the separability between the groups. 

#############################
####### LDA #################
#############################

# #PCA can only be performed if ImputePCA is performed or if rows with NA (many!!!!) are removed
# pca_prcomp <- prcomp(summary_PCA_vars,
#               center = TRUE,
#               scale. = TRUE) 
# 
# prop.RES.PCA = pca_prcomp$sdev^2/sum(pca_prcomp$sdev^2)


##############################
### LDA ASSUMPTIONS ##########
##############################

#check if I have more predictor vars than the size of the smallest group
table(summary_PCA_vars_hit$Hit)
length(summary_PCA_vars)
#it is a bit risky as they are very close...

lda <- lda(summary_PCA_vars_hit$Hit ~ ., 
           summary_PCA_vars)

## proportion of LDs that encapsulate the variation (in case of 1 Ld, the prop.lda = 1)
#prop.lda = lda$svd^2/sum(lda$svd^2)

## get / compute LDA scores from LDA coefficients / loadings

#PREDICTION SHOULD BE PERFORMED IN THE SECOND HALF OF THE DATASET? (NOT HEREBY ANALYSED)
plda <- predict(object = lda,
                newdata = summary_PCA_vars) #predict(lda_TEST)$x

#### ISSUEEE
#LOTS of missing values affecting the analysis
#N of rows including missing cases:
sum(!complete.cases(summary_PCA_vars))
sum(!complete.cases(plda$class))

#prop of NAs over total in % (missing values)
propNAs_total <- (sum(sapply(summary_PCA_vars, function(x) sum(is.na(x)))) /(nrow(summary_PCA_vars) * ncol(summary_PCA_vars))) *100
##prop missing by variable in %
propNAs_byvar <- (sapply(summary_PCA_vars, function(x) sum(is.na(x)))/nrow(summary_PCA_vars)) * 100

##prop missing by row (show only top) - potentially droppable rows?
prop_Na_row <- (apply(summary_PCA_vars, 1, function(x) sum(is.na(x)))/ncol(summary_PCA_vars))*100
prop_Na_row <- data.frame(prop_missing=sort(prop_Na_row, decreasing=TRUE))
propNAs_byrow_25perc <- ((length(prop_Na_row[which(prop_Na_row>25),]))/nrow(summary_PCA_vars))*100 #show all rows with values over 25% missing
hist(prop_Na_row$prop_missing) + abline(v=20)
#summary of Nas distribution
cat("prop of NAs over total in %",propNAs_total, "\n\nProp missing by variable \n",propNAs_byvar,"\n\nprop of rows with more than 25% missing",propNAs_byrow_25perc)



#create a histogram of the discriminant function values
par(mfrow=c(1,1), oma=c(0,0,2,0),mar = c(4.5, 3.8, 1, 1.1))
ldahist(data = plda$x[,1], g=summary_PCA_vars_hit$Hit,nbins = 25) + mtext("LDA hist", line=0, side=3, outer=TRUE)
## ASSIGN THE PREDICTION TO THE dataframe
summary_PCA_vars_hit_PRED <-  cbind("REPLICATE" = summary_AUTO_transf$REPLICATE,
                                    "PERIOD" = summary_AUTO_transf$PERIOD,
                                    "pair" = summary_AUTO_transf$pair,
                                    "int_start_frame "= summary_AUTO_transf$int_start_frame,
                                    "int_end_frame" = summary_AUTO_transf$int_end_frame,
                                    summary_PCA_vars_hit, 
                                    "Predicted_Hit" =  plda$class)

#########################################
### critical success index (CSI) ########
#########################################

# # CSI is a verification measure of categorical forecast performance
# # calculated as CSI= (hits)/(hits + false alarms + misses)
# # The CSI is not affected by the number of non-event forecasts that verify (correct rejections).

#convert int_start and int_end back to frames
summary_PCA_vars_hit_PRED$int_start_frame <- as.numeric(summary_PCA_vars_hit_PRED$int_start_frame)
summary_PCA_vars_hit_PRED$int_end_frame <- as.numeric(summary_PCA_vars_hit_PRED$int_end_frame)

CSI_val <- NULL

#create matrix
for (REPLICATE in c("R3SP","R9SP")) 
{
  for (PERIOD in c("pre","post"))
  {

trimmed_AUTO_pred <- summary_PCA_vars_hit_PRED[which(summary_PCA_vars_hit_PRED$Predicted_Hit == 1 & summary_PCA_vars_hit_PRED$REPLICATE==REPLICATE & summary_PCA_vars_hit_PRED$PERIOD==PERIOD),]

#assign(paste0("trimmed_AUTO_pred","_", REPLICATE,PERIOD), trimmed_AUTO_pred) 

IF_frames_temp <- get(grep(paste0("IF_frames","_", REPLICATE,PERIOD),ls(),value=TRUE))
int_mat_manual_temp <- get(grep(paste0("int_mat_manual","_", REPLICATE,PERIOD),ls(),value=TRUE))

#this section replicates what is done in the matrix calculation
int_mat_auto_trim <- matrix(0L, nrow = dim(ids_pairs)[1], ncol = (dim(IF_frames_temp)[1] ))
rownames(int_mat_auto_trim) <- ids_pairs$pair
colnames(int_mat_auto_trim) <- c(IF_frames_temp$frame_num)


## Trimmed Auto
if (nrow(trimmed_AUTO_pred)>0) {
for (i in 1:nrow(trimmed_AUTO_pred))
{  
  PAIR <- trimmed_AUTO_pred$pair[i]
  int_mat_auto_trim[PAIR,] <- int_mat_auto_trim[PAIR,] + c(rep(0,( trimmed_AUTO_pred$int_start_frame[i]-1)),
                                                 rep(1,(trimmed_AUTO_pred$int_end_frame[i]) - trimmed_AUTO_pred$int_start_frame[i] + 1),
                                                 rep(0,(length(IF_frames_temp$frame_num) - trimmed_AUTO_pred$int_end_frame[i])))
}#matrix
  

# new category of AUTO  (trimmed prediciton) after eliminating those predicted to be 0 from the LDA prediction and calculate CSI again
# recalculate frame by frame diff matrix using LDA output, but calculate things in a different way as:

 # ROW BY ROW COMPARISON OF DFS
TruePositive  <- sum(int_mat_manual_temp + int_mat_auto_trim==2 & int_mat_manual_temp-int_mat_auto_trim ==0)
TrueNegative  <- sum(int_mat_manual_temp + int_mat_auto_trim==0 & int_mat_manual_temp-int_mat_auto_trim ==0) #not exactly right... the calc should be constrained at certain parts? i.e. for pairs present in the Dataset
FalsePositive <- sum(int_mat_manual_temp-int_mat_auto_trim ==-1)
FalseNegative <- sum(int_mat_manual_temp-int_mat_auto_trim ==1)

#calculate frame by frame CSI
CSI <- TruePositive/(TruePositive+FalseNegative+FalsePositive)

CSI_val <- c(CSI_val,CSI)
CSI <- NULL
}else{ # if DF contains any Hit 
  CSI_val <- c(CSI_val,0)
  
}

  }## REPLICATE
    }##PERIOD

perc_CSI <- sum(CSI_val)/length(CSI_val)*100


######INNCLUDE PCA in outer loop (described as outer of MAIN).
# THE OUTPUT SHOULD INCLUDE:
# VALUES FOR ALL VARS, SAMPLE SIZES, INFO ON CAPS USED (MAYBE), FINAL CSI SCORE


# Effect size
# Some suggest the use of eigenvalues as effect size measures, however, this is generally not supported.[9] Instead, the canonical correlation is the preferred measure of effect size. It is similar to the eigenvalue, but is the square root of the ratio of SSbetween and SStotal. It is the correlation between groups and the function.[9] Another popular measure of effect size is the percent of variance[clarification needed] for each function. This is calculated by: (λx/Σλi) X 100 where λx is the eigenvalue for the function and Σλi is the sum of all eigenvalues. This tells us how strong the prediction is for that particular function compared to the others
# 

cat(paste(" from https://en.wikipedia.org/wiki/Linear_discriminant_analysis#LDA_for_two_classes
\nComparison to logistic regression

\nDiscriminant function analysis is very similar to logistic regression, and both can be used to answer the same research questions.[9] Logistic regression does not have as many assumptions and restrictions as discriminant analysis. However, when discriminant analysis’ assumptions are met, it is more powerful than logistic regression.[28] Unlike logistic regression, discriminant analysis can be used with small sample sizes. It has been shown that when sample sizes are equal, and homogeneity of variance/covariance holds, discriminant analysis is more accurate.[7] Despite all these advantages, logistic regression has none-the-less become the common choice, since the assumptions of discriminant analysis are rarely met.[8][7] 
          ",sep = " "))

  
### PLOTTING PCA AND LDA TOGETHER


#dataset = data.frame(Hit = summary_PCA_vars_hit[,"Hit"],
#                     pca = pca_prcomp$x, lda = plda$x)

# # p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = Hit, shape = Hit), size = 2.5) + 
# #   labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
# #        y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))
# LDmeans <- plyr::ddply(dataset, "Hit", summarise, grp.mean=mean(LD1))
# 
# 
# p1 <- ggplot(dataset, aes(LD1, fill = Hit)) + geom_density(alpha = 0.2)+ 
#       geom_vline(data=LDmeans, aes(xintercept=grp.mean, color=Hit), linetype="dashed") +
#       labs(title = "Density plot for predicted values from model function, by Hit rate",
#            x = paste("proportion of discriminability explained LD1 (", percent(prop.lda[1]), ")", sep=""))
# 
# 
# p2 <- ggplot(dataset) + geom_point(aes(pca.PC1, pca.PC2, colour = Hit, shape = Hit), size = 2.5) +
#   labs(x = paste("PC1 (", percent(prop.RES.PCA[1]), ")", sep=""),
#        y = paste("PC2 (", percent(prop.RES.PCA[2]), ")", sep=""))
# 
# grid.arrange(p1, p2)




###EXTRA #############################################
######################################################


### LDA with data partitioning ####
### PArtition plots ... ####
###ETC ETC

# #------------------------
# library(caret)
# data(mdrr)
# mdrrDescr <- mdrrDescr[, -nearZeroVar(mdrrDescr)]
# mdrrDescr <- mdrrDescr[, -findCorrelation(cor(mdrrDescr), .8)]
# set.seed(1)
# inTrain <- createDataPartition(mdrrClass, p = .75, list = FALSE)[,1]
# train <- mdrrDescr[ inTrain, ]
# test  <- mdrrDescr[-inTrain, ]
# trainClass <- mdrrClass[ inTrain]
# testClass  <- mdrrClass[-inTrain]
# 
# set.seed(2)
# ldaProfile <- rfe(train, trainClass,
#                   sizes = c(1:10, 15, 30),
#                   rfeControl = rfeControl(functions = ldaFuncs, method = "cv"))
# 
# 
# postResample(predict(ldaProfile, test), testClass)
# ldaProfile$optVariables
# 
# #-------------------
# #https://stackoverflow.com/questions/68307682/r-lda-linear-discriminant-analysis-how-to-get-compute-lda-scores-from-lda-co
# #or use
# 
# #https://www.andreaperlato.com/mlpost/linear-discriminant-analysis/
# 
# #https://rpubs.com/ifn1411/LDA
# 
# 
# #----------------------
# # lda$svd explains the ratio of between- and within-group variation.
# lda_TEST$svd
# 
# #--------------------
# #If there are many groups, or if you quickly want to find the maximum probability for each sample, this command will help:
#   
# round(apply(lda_TEST_pred$posterior, MARGIN=1, FUN=max), 2)
# 
# #For this data set, most of these maximum probabilities are large (>0.9), indicating that most samples are 
# #confidently assigned to one group. If most of these probabilities are large, the overall classification is 
# #successful.
# 
# #The scaling value of the lda object gives the loadings (also called the slopes, coefficients, or weights) of each variable on each discriminant function.
# round_scaling <- round(lda_TEST$scaling, 2)
# # round_scaling$var  <- row.names(round_scaling)
# names <- rownames(round_scaling)
# rownames(round_scaling) <- NULL
# round_scaling <- as.data.frame(cbind(names,round_scaling))
# round_scaling <- round_scaling[order(round_scaling$LD1),]
# write.table(round_scaling,paste0(DATADIR,"round_scaling.txt"),sep="\t",row.names=FALSE)
# 
# #-------------------------
# 
# #####
# #From the Partition Plot above, we can see classification for eachof observation in the training dataset based 
# #on the Linear Discriminant Analysis Model, and for every combination of two variables.
# 
# X11(width=15, height=15)
# # Adjust the numbers for the option mar, if needed; see ?partimat and ?par 
# # for more help on margins
# par(mfrow=c(2,3))
# library(klaR)
# klaR::partimat(Hit~mean_movement_angle_diff+prop_time_undetected_REC+StDev_angle_ACT, data=summary_PCA_vars_hit, method="qda",
#                main = "\nPartition Plot for Quadratic Discriminant Analysis Model \nred=incorrectly classified",mar=c(5,4,2,2))


# par(mfrow = c(4,2))
# for(i in 3:10)
#   drawparti(mydata[,1], mydata[,2], mydata[,i], method = "rda",
#             gamma=0, lamda=1)


# look at 
# http://rstudio-pubs-static.s3.amazonaws.com/84669_cd15214061d44e1493ffee69c5d55925.html
# https://rpubs.com/ifn1411/LDA
# https://www.andreaperlato.com/mlpost/linear-discriminant-analysis/
# https://towardsdatascience.com/linear-discriminant-analysis-explained-f88be6c1e00b

#WHAT ARE THE ASSUMPTIONS FOR LDA ANYWAY? (HOMOSCHEDASTICITY FOR LDA BUT NOT FOR QDA)

#svm SEE https://drive.google.com/drive/u/0/folders/0B9OJD63YvZ8Jck5rM0h0cTczTFk?resourcekey=0-NG8IXWZQ064vZnrP0tt5OA 
#BUT IN r




#############   extra: ################################
  
  
  #Remove jerk given the high N of missing cases (it causes the proliferation of vars BUT could be an important discriminating factor!)
  #summary_PCA_vars <- summary_PCA_vars[, -( grep("jerk" , colnames(summary_PCA_vars),perl = TRUE) ) ]
  


##ACCURACY OF CLASSIFICATION #NOT TOO RELEVANT
# accuracy  <- xtabs(~summary_PCA_vars_hit$Hit+plda$class)
#sum(accuracy[row(accuracy) == col(accuracy)]) / sum(accuracy)


# names(dimnames(accuracy)) <- c("True_class", "Predicted_class")
# accuracy <- as.data.frame(accuracy)
# 
# accuracy$Cat <- ifelse(accuracy$True_class == 0 & accuracy$Predicted_class == 0 , "TrueNeg", 
#                        ifelse(accuracy$True_class == 0 & accuracy$Predicted_class == 1 , "FalsePos",
#                               ifelse(accuracy$True_class == 1 & accuracy$Predicted_class == 0 , "FalseNegA", "TruePos")))

# TrueNeg   <- accuracy$Freq[which(accuracy$Cat=="TrueNeg")]
# FalsePos  <- accuracy$Freq[which(accuracy$Cat=="FalsePos")]
# FalseNegA <- accuracy$Freq[which(accuracy$Cat=="FalseNegA")]
# TruePos   <- accuracy$Freq[which(accuracy$Cat=="TruePos")]
# #FalseNegB <- N. of real grooming events
#TRUE NEGB HAS TO BE CALCULATED GETTING BACK TO THE MATRIX CUTTING. SEE NATH NOTES (IS THERE A RECORDING?)

