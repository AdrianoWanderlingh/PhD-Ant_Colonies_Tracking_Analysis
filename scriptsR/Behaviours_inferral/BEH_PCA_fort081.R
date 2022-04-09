##########################################################################################
############## BEH CLASSIFICATION ANALYSIS #################################################
##########################################################################################

#script dependant on BEH_MAIN_behaviours_analysis_fort081.R

#For previous versions of this script, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral


print(paste("PERFORM CLASSIFICATION for the Loop_ID",Loop_ID))

# FORECAST VERIFICATION
# https://www.swpc.noaa.gov/sites/default/files/images/u30/Forecast%20Verification%20Glossary.pdf
# https://en.wikipedia.org/wiki/Sensitivity_and_specificity

##########################################################################################
############## LDA ASSUMPTIONS ANALYSIS ##################################################
##########################################################################################
 
# The script starts with the assumptions checking for the LDA/QDA. Data are checked and transformed accordingly.
# the transformed data are later on used to classify using RandomForest algorithms

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

1. ######## Check for collinearity ########
\n --However, LDA cannot be used with data which are highly collinear, such as is the case in the present
study. Therefore, in this investigation PCA was performed prior to LDA to remove the effect of collinearity
and reduce the number of variables in the X data matrix.-- \n\nRef:\nMetabolites 2014, 4, 433-452; doi:10.3390/metabo4020433

Also, is collinearity reduced after variables transformation via BestNormalize? NO

Careful: correlated and collinear is not the same thing!
remove correlated var following https://stats.stackexchange.com/questions/121131/removing-collinear-variables-for-lda-qda-in-r

2. ######## Take care of missing cases as they badly affect the analysis given the quantity of NAs ######## 
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

#data structure
#Once the THRESH for $disagreement has been established:
# Hit == 1 is a True Positive
# Hit ==0 is a False Positive
# this is done in BEH_Auto_Man_agreement_matrix_fort081.R

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
# COMPARE THIS CORRPLOT WITH THE ONE GENERATED AFTER VARIABLES TRANSFORMATION
# corrplot(results.cor, type="upper", tl.cex= .8)

##Remove vars (or rows) with too many NAs for the cor.mtest to work
# res1 <- cor.mtest(results.cor, conf.level = .95, method="pearson", na.action = "na.omit") #kendall


# corrplot(results.cor, type="upper", order="hclust", tl.cex= .8, 
#          p.mat = res1$p, sig.level = 0.05, insig = "blank", method = "number", number.cex = .7)



#######################################################
######### MISSING CASES ###############################
#######################################################

#LOTS of missing values affecting the analysis

#summary_LDA_vars <- summary_AUTO_vars[, -match(c("REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","disagreement","Hit"), names(summary_AUTO_transf))] 
sapply(summary_AUTO_vars, function(x) sum(is.na(x)))
#summary_LDA_vars_hit <- cbind(summary_LDA_vars,Hit=summary_AUTO_transf$Hit)

#N of rows including missing cases:
sum(!complete.cases(summary_AUTO_vars))

#prop of NAs over total in % (missing values)
propNAs_total <- (sum(sapply(summary_AUTO_vars, function(x) sum(is.na(x)))) /(nrow(summary_AUTO_vars) * ncol(summary_AUTO_vars))) *100
##prop missing by variable in %
propNAs_byvar <- round((sapply(summary_AUTO_vars, function(x) sum(is.na(x)))/nrow(summary_AUTO_vars)),2) * 100

##prop missing by row (show only top) - potentially droppable rows?
prop_Na_row <- (apply(summary_AUTO_vars, 1, function(x) sum(is.na(x)))/ncol(summary_AUTO_vars))*100
prop_Na_row <- data.frame(prop_missing=sort(prop_Na_row, decreasing=TRUE))
propNAs_byrow_25perc <- ((length(prop_Na_row[which(prop_Na_row>25),]))/nrow(summary_AUTO_vars))*100 #show all rows with values over 25% missing
#hist(prop_Na_row$prop_missing) #+ abline(v=20, col="blue")
#summary of Nas distribution
cat("prop of NAs over total in %",propNAs_total, "\n\nProp missing by variable \n",propNAs_byvar,"\n\nprop of rows with more than 25% missing",propNAs_byrow_25perc)

### NAs by class
Hit_FULL.yes <- subset(summary_AUTO, Hit == 1)
Hit_FULL.no <- subset(summary_AUTO, Hit == 0)

round((sapply(Hit_FULL.yes, function(x) sum(is.na(x)))/nrow(Hit_FULL.yes)),3)*100
round((sapply(Hit_FULL.no, function(x) sum(is.na(x)))/nrow(Hit_FULL.no)),3)*100

#drop jerk as it has >5% missing cases for the Hit.yes (minority class) and >17.8%  missing cases for the hit.no (majority class). 
# redun and hierarchical clustering have shown anyway that it is a droppable feature (predictable by other dataset variables)

#Also, currently dropping prop_time_undetected_REC/ACT as they have a very difficult distribution (hard to normalise)
TRIM <- c("mean_jerk_PxPerSec3_ACT","mean_jerk_PxPerSec3_REC","prop_time_undetected_REC","prop_time_undetected_ACT")
summary_AUTO.1 <- summary_AUTO[, -match(c(TRIM), names(summary_AUTO))]

##############################################################################
########## DROP NAs ##########################################################
##############################################################################

# Many NAs are present only in the False positives (majority class), therefore dropping them is going to affect only very few cases of the True Postives
summary_AUTO_NAOmit <- summary_AUTO.1[complete.cases(summary_AUTO.1), ]

###############################################################################
###### NORMALISE VARS BEFORE CLASSIFICATION ###################################
###############################################################################
###SPEED UP WITH THE PARALLEL PACKAGE?
#empty base
summary_AUTO_NAOmit_transf <- data.frame()[1:nrow(summary_AUTO_NAOmit), ]
summary_AUTO_NAOmit$ACT <- as.factor(summary_AUTO_NAOmit$ACT)
summary_AUTO_NAOmit$REC <- as.factor(summary_AUTO_NAOmit$REC)
summary_AUTO_NAOmit$Hit <- as.factor(summary_AUTO_NAOmit$Hit)
summary_AUTO_NAOmit$disagreement <- as.factor(summary_AUTO_NAOmit$disagreement)

summary_AUTO_NAOmit$int_start_frame <- as.character(summary_AUTO_NAOmit$int_start_frame)
summary_AUTO_NAOmit$int_end_frame <- as.character(summary_AUTO_NAOmit$int_end_frame)

#transform variables if they are not normal
start.time <- Sys.time()
for (variable in names(summary_AUTO_NAOmit)){
  if (is.numeric(summary_AUTO_NAOmit[,variable])) {
    val <- shapiro.test(summary_AUTO_NAOmit[,variable])
    if (unname(val$p.value)<0.05) {
      print(paste("for [",variable, "] the Shapiro-Wilk Test has p < 0.05, the data is not normal. Transform it.",sep=" "))
      #scale and center variable (as reccomended for data prep for LDA in Metabolites 2014, 4, 433-452; doi:10.3390/metabo4020433)
      scaled_VAR <- scale(summary_AUTO_NAOmit[,variable])
      #find the best transformation
      #This function currently estimates the Yeo-Johnson, the Box Cox  (if the data is positive), the log_10(x+a), the square-root (x+a), the arcsinh and the ordered quantile normalization
      BNobject <- bestNormalize(scaled_VAR,standardize = FALSE) #with standardize = FALSE also no_tranformation is attempted 
      summary_AUTO_NAOmit_transf$var <- BNobject$x.t; names(summary_AUTO_NAOmit_transf)[names(summary_AUTO_NAOmit_transf) == 'var'] <- paste(variable,class(BNobject$chosen_transform)[1],sep=".")
    }else{print(paste("for [",variable,"] the Shapiro-Wilk Test has p > 0.05, the data is normal. Keep it.",sep=" "))
          summary_AUTO_NAOmit_transf[variable]<- scale(summary_AUTO_NAOmit[,variable])
    }}else{print(paste("non numeric attribute. Pasting [",variable,"] in the new dataset",sep = " "))
           summary_AUTO_NAOmit_transf[variable]<- summary_AUTO_NAOmit[,variable]}
}
end.time <- Sys.time()
(time.taken <- end.time - start.time)

#### CORR PLOT AFTER TRANSFORMATION
summary_AUTO_NAOmit_transf_vars <- summary_AUTO_NAOmit_transf[, -match(c("REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","disagreement","Hit"), names(summary_AUTO_NAOmit_transf))] 
#summary_AUTO_transf_vars_NAOmit <- summary_AUTO_transf_vars[(apply(summary_AUTO_transf_vars, 1, function(x) sum(is.na(x)))/ncol(summary_AUTO_transf_vars)) < 0.05,]

results.cor_transf <- cor(summary_AUTO_NAOmit_transf_vars)

pdf(file=paste(DATADIR,"CorrPlots","AUTO","Gap",MAX_INTERACTION_GAP, ".pdf", sep = "_"), width=20, height=5)
par(mfrow=c(1,2),mar = c(4.5, 3.8, 1, 1.1))
corrplot(results.cor, type="upper", tl.cex= .4,title = "\nresults.cor",)
corrplot(results.cor_transf, type="upper", tl.cex= .4,title = "\nresults.cor_transf",)
dev.off()
#it seems that the transformation rises the collinearity instead of reducing it... SOLVE



#summary_AUTO_REP_PER_dget <- dget("/home/cf19810/Documents/Ants_behaviour_analysis/Data/summary_AUTO_REP_PER_16feb22.txt") # load file created with dput 
summary_AUTO_NAOmit_transf$Hit <- as.factor(summary_AUTO_NAOmit_transf$Hit)

#check N of missing values per variable
#remove ant names/rep/int/etc in pca (keep only vars)
summary_LDA_vars <- summary_AUTO_NAOmit_transf[, -match(c("REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","disagreement","Hit"), names(summary_AUTO_NAOmit_transf))] 
sapply(summary_LDA_vars, function(x) sum(is.na(x)))


#strip attributes so that variables are just numeric (now they contain info on center and scaling)
for (var in colnames(summary_LDA_vars)) {
  # attr(summary_LDA_vars[,deparse(as.name(var))], "scaled:center") <- NULL
  # attr(summary_LDA_vars[,deparse(as.name(var))], "scaled:scale") <- NULL
  summary_LDA_vars[,var]    <-     as.numeric(summary_LDA_vars[,var])
}

#Re-add the class Hit
summary_LDA_vars_hit <- cbind(summary_LDA_vars,Hit=summary_AUTO_NAOmit_transf$Hit)


## eventually trim unwanted vars
#sapply(summary_LDA_vars_trim, function(x) sum(is.na(x))) 


########################################################
######### TEST ASSUMPTIONS #############################
########################################################

### The size of the smallest group must be larger than the number of predictor variables ###
#given that likely the N of TP is going to be very close to the N of variable, reduction in the number of used variables would be beneficial.

### Multivariate normality ###
#With the scaling, centering and transformation of the variables, the distribution should result normal.

### Homogeneity of variance/covariance (homoscedasticity) ###
#Ongoing....

### Multicollinearity ###
#Issues in using the PCR to reduce the N of variables as there is a large N of Nas.
#this can maybe fixed by dropping correlated variables, using redun and similar techniques...
# Features could then be selected using a pre-processing algorithm as RELIEF (as suggested in ACUNA, Edgar;et al. A comparison of feature selection procedures for classifiers based on kernel density estimation. 2003). A review on RELIEF is provided by URBANOWICZ, Ryan J., et al. Relief-based feature selection: Introduction and review .2018.

# As reported in Fig 10, For a class imbalance of 0.9 (i.e. 90% class 0, 10% class 1), we observe that
# ReliefF with a large number of neighbors (i.e. 50%) fails to perform,
# ReliefF with 100 NN and SURF demonstrate slight deficits, but all other
# RBAs -as MultiSURF- perform optimally (see URBANOWICZ, Ryan J., et al. Benchmarking relief-based feature selection methods for bioinformatics data mining. Journal of biomedical informatics, 2018)



########################################################
######### RELIEF feature selection #####################
########################################################

#Feature evaluation algorithms available for classification problems (ReliefF)
# The core idea behind Relief algorithms is to estimate the quality of attributes 
# on the basis of how well the attribute can distinguish between instances that are near to each other.

# Relief calculates a feature score for each feature which can then be applied to rank and select top scoring features for feature selection
# given a randomly selected instance Ri, Relief searches for its
# two nearest neighbors: one from the same class, called nearest hit H, and the
# other from the different class, called nearest miss M. It updates the
# quality estimation W [A] for all attributes A depending on their values for Ri,
# M, and H 


###APPLY IT TO TRANSFORMED VARIABLES?

##### try RELIEF ######
estReliefF <- attrEval("Hit", summary_LDA_vars_hit, estimator="ReliefFequalK", # with ReliefFexpRank the output is very similar
                       ReliefIterations=30)

estReliefF <- as.data.frame(estReliefF)
estReliefF$Variable <- rownames(estReliefF); rownames(estReliefF) <- NULL
#missing cases by var
MissByVar <- round((sapply(summary_LDA_vars, function(x) sum(is.na(x)))/nrow(summary_LDA_vars)) * 100, 2)
MissByVar <- as.data.frame(MissByVar)
MissByVar$Variable <- rownames(MissByVar); rownames(MissByVar) <- NULL
#add information of missing cases to the Relief feature score
estReliefF <- plyr::join(x=estReliefF, y=MissByVar, type = "full", match = "all")
#sort by score
estReliefF[order(estReliefF$estReliefF,decreasing = TRUE),]

# DO the same but adding statistical significance 
# a statistical inferential formalism is needed to avoid imposing arbitrary thresholds to select the most important features
# LE, Trang T., et al. Statistical inference Relief (STIR) feature selection. Bioinformatics, 2019, 35.8: 1358-1365. 

#RF.method = "multisurf" # SURF identifies nearest neighbors (both hits and misses) based on a distance threshold from the target instance defined by the average distance between all pairs of instances in the training data.[20] Results suggest improved power to detect 2-way epistatic interactions over ReliefF. 
metric <- "manhattan"# manhattan #using it instead of euclidean as it is the most common choice as sugegsted by Le et al.2019 
#BUT manhattan is suggested for high-dimensional data (N of features > )

#this takes a few seconds to compute

#Multisurf 
neighbor.idx.observed <- find.neighbors(summary_LDA_vars, summary_LDA_vars_hit$Hit, k = 0, method = "multisurf")
results.list <- stir(summary_LDA_vars, neighbor.idx.observed, k = k, metric = metric, method = "multisurf")
t_sorted_multisurf <- results.list$STIR_T
t_sorted_multisurf$Variable <- rownames(t_sorted_multisurf)
round(t_sorted_multisurf$t.pval.adj,5)

#add information of missing cases to the Relief feature score
estReliefStir <- plyr::join(x=estReliefF, y=t_sorted_multisurf, type = "full", match = "all")


#select only variables with p<0.05
RELIEF_selected <- t_sorted_multisurf[which(t_sorted_multisurf$t.pval.adj<0.01),"Variable"]
summary_LDA_vars_RELIEF <- summary_LDA_vars[, match(c(RELIEF_selected), names(summary_LDA_vars))] 
#check vars correlation
results.cor_RELIEF <- cor(summary_LDA_vars_RELIEF) 
#corrplot(results.cor_RELIEF, type="upper", tl.cex= .4,title = "\nresults.cor",)

#add the Hit classes
summary_LDA_vars_RELIEF_hit <- cbind(summary_LDA_vars_RELIEF,Hit=summary_AUTO_NAOmit_transf$Hit)

##### corelearn can also be used to try other models using CoreModel

# #############################
# ####### PCA #################
# #############################
# 
######### PERFORM PCA TO REMOVE CORRELATED VARIABLES, NOT LINKED WITH RELIEF ########
#
# #PCA with missing values according to Dray & Josse (2015) https://doi.org/10.1007/s11258-014-0406-z
# #imputePCA of the R package missMDA
# #Impute the missing values of a dataset with the Principal Components Analysis model. 
# #Can be used as a preliminary step before performing a PCA on an incomplete dataset.
# #missing values are replaced by random values, and then PCA is applied on the completed data set, and missing values are then updated by the fitted values
# ## Imputation for LDA_VARS which contains NAs
# res.comp <- imputePCA(LDA_VARS,method="Regularized")
# 
# ## 1. A PCA can be performed on the imputed data for LDA_VARS
# res.pca1 <- PCA(res.comp$completeObs)
# 
## 2. PCA on the base data for LDA_VARS_trim
#RES.PCA <- PCA(summary_LDA_vars, scale.unit = FALSE) #scaled when normalised

#with relief selected features
#RES.PCA <- PCA(results.cor_RELIEF, scale.unit = FALSE) #scaled when normalised


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
# Use habillage to specify groups for coloring
# fviz_pca_ind(RES.PCA, pointsize =summary_AUTO_NAOmit_transf$disagreement,
#              label = "none", # hide individual labels
#              title = "PCA - indivdual interactions \nPointsize by disagreement rate",
#              habillage = as.factor(summary_AUTO_NAOmit_transf$Hit), # color by groups
#              addEllipses = TRUE
# )
# 
# #}
# 



#######################################################
######### CHOOSE DATA #################################
#######################################################

LDA_DATA <- "RELIEF" # or "ALL_VARS" 

if(LDA_DATA =="ALL_VARS"){LDA_VARS <- summary_LDA_vars; LDA_VARS_HIT <- summary_LDA_vars_hit}
if(LDA_DATA =="RELIEF"){LDA_VARS <- summary_LDA_vars_RELIEF; LDA_VARS_HIT <- summary_LDA_vars_RELIEF_hit}


#######################################################
######### TEST DATA for LDA/QDA #######################
#######################################################

# #LDA: Linear discriminant analysis also called DFA D.function.A.
# #LDA focuses on finding a feature subspace that maximizes the separability between the groups. 


# #PCA can only be performed if ImputePCA is performed or if rows with NA (many!!!!) are removed
# pca_prcomp <- prcomp(LDA_VARS,
#               center = TRUE,
#               scale. = TRUE) 
# 
# prop.RES.PCA = pca_prcomp$sdev^2/sum(pca_prcomp$sdev^2)


#### Assumptions testing #####
### Homogeneity of variance/covariance (homoscedasticity) ###
#Box’s M Test is extremely sensitive to departures from normality; the fundamental test assumption is that your data is multivariate normally distributed. Therefore, if your samples don’t meet the assumption of normality, you shouldn’t use this test.
# mvn(data = LDA_VARS,subset = LDA_VARS_HIT$Hit,mvnTest="mardia") #crashes R

##### Checking Assumption of Equal Variance-Covariance matrices
#Checking the Assumption of Equal Covariance Ellipse
# heplots::covEllipses(LDA_VARS, 
#                      LDA_VARS_HIT$Hit, 
#                      fill = TRUE, 
#                      pooled = FALSE, 
#                      col = c("blue", "red"), 
#                      variables = c(1:14), 
#                      fill.alpha = 0.05)
#using the BoxM test in order to check our assumption of homogeneity of variance-covariance matrices
# for p<0.05 there is a problem of heterogeneity of variance-covariance matrices
boxm <- heplots::boxM(LDA_VARS, LDA_VARS_HIT$Hit)
boxm #NOT met for the summary_LDA_vars_RELIEF
plot(boxm)

###### Checking Multicollinearity
# One way to detect multicollinearity is to check the Choleski decomposition of the correlation matrix 
#- if there's multicollinearity there will be some diagonal elements that are close to zero
Choleski_decomp <- chol(cor(LDA_VARS))
diag(Choleski_decomp) #looks good

### Checking Assumption of Normality
# Hit.yes <- subset(LDA_VARS_HIT, Hit == 1) 
# Hit.no <- subset(LDA_VARS_HIT, Hit == 0)
## plotting QQplots
# par(mfrow = c(4, 4)) 
# for(i in names(LDA_VARS)) { 
#   qqnorm(Hit.yes[[i]]); qqline(Hit.yes[[i]], col = 2) 
# }
# 
# for(i in names(LDA_VARS)) { 
#   qqnorm(Hit.no[[i]]); qqline(Hit.no[[i]], col = 2) 
# }

# As covariances are not equal, quadratic discriminant analysis will be used.
############# ASSUMPTIONS FOR QDA #################

# - Observation of each class is drawn from a normal distribution (same as LDA).
#   may be strong enough, make sure vars are normalised!

# - QDA assumes that each class has its own covariance matrix (different from LDA).
#   this one is met


####################################
###########plotting ################
####################################
PLOTTING_TRANSF_VAR <- FALSE
if (PLOTTING_TRANSF_VAR) {
  print("plotting vars")

#transform to long format
summary_PCA_long <- reshape2::melt(LDA_VARS_HIT,id.vars=c("Hit")) #explanation on the warning message https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message


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
###########################
###### DATA SAMPLING ######
###########################

#TESTING SMOTE AS IT PROVED TO BE THE HISHEST SCORING FOR F1-SCORE FOR A DATASET (ABALONE) OF VERY SIMILAR CARACTHERISTICS (IMBALANCE, N CASES, N PARAMS)
# IDEA DERIVED FROM TAHIR, Muhammad Atif; KITTLER, Josef; YAN, Fei. Inverse random under sampling for class imbalance problem and its application to multi-label classification. Pattern Recognition, 2012, 45.10: 3738-3750.

# prepare for the display of classification performance 
display <- function(prediction, reference, Hitclass) {
  cm <- caret::confusionMatrix(data=prediction, reference=reference,
                               mode = "sens_spec", positive=Hitclass)
  print(cm$table)
  data.frame(Accuracy=cm$overall[1], Sensitivity=cm$byClass[1], Specificity=cm$byClass[2],
             row.names=NULL)
            }

#dataset
table(LDA_VARS_HIT$Hit)

####################################################
###### CLASSIFIERS ON THE BALANCED DATASET #########
####################################################

#Now we test several sampling algorithms to balance the dataset.
train_sbc   <- SBC(LDA_VARS_HIT,   "Hit") #Under-Sampling Based on Clustering (SBC)
train_ros   <- ROS(LDA_VARS_HIT,   "Hit") #random over-sampling algorithm (ROS)
train_rus   <- RUS(LDA_VARS_HIT,   "Hit") #random under-sampling algorithm (RUS)
train_smote <- SMOTE(LDA_VARS_HIT, "Hit") #synthetic minority over-sampling technique (SMOTE)


##############################
### LDA / QDA ################
##############################

### Quadratic Discriminant Analysis

#train on the class-balanced dataset
fit_QDA_sbc <- qda(Hit ~ ., data = train_sbc)
fit_QDA_ros <- qda(Hit ~ ., data=train_ros)
fit_QDA_rus <- qda(Hit ~ ., data=train_rus)
fit_QDA_smote <- qda(Hit ~ ., data=train_smote)
#predict class on the full dataset
pred_class_QDA_sbc <- predict(fit_QDA_sbc, LDA_VARS_HIT, type="response") #predict class on the full train dataset
pred_class_QDA_ros <- predict(fit_QDA_ros, LDA_VARS_HIT, type="response")
pred_class_QDA_rus <- predict(fit_QDA_rus, LDA_VARS_HIT, type="response")
pred_class_QDA_smote <- predict(fit_QDA_smote, LDA_VARS_HIT, type="response")
# #check performance 
# perf_QDA_sbc <- display(pred_class_QDA_sbc$class, LDA_VARS_HIT$Hit, "1") 
# perf_QDA_ros <- display(pred_class_QDA_ros$class, LDA_VARS_HIT$Hit, "1")  
# perf_QDA_rus <- display(pred_class_QDA_rus$class, LDA_VARS_HIT$Hit, "1")
# perf_QDA_smote <- display(pred_class_QDA_smote$class, LDA_VARS_HIT$Hit, "1")
# QDAPred <- rbind(perf_QDA_ros, perf_QDA_rus, perf_QDA_sbc, perf_QDA_smote)
# rownames(QDAPred) <- c("Balanced by ROS", "Balanced by RUS",
#                               "Balanced by SBC", "Balanced by SMOTE")
# QDAPred

##############################
### RANDOM FOREST ############
##############################

#Random forests is an ensemble learning method for classification that operates by constructing a multitude of decision
#trees at training time. For classification tasks, the output of the random forest is the class selected by most trees.
# to plot RF: https://rpubs.com/markloessi/498787

#### hyperparam tuning
# # NOT USED, doesn't seem to bring significant improvements to the classification
# # RF parameters: 
# # mtry: Number of variables randomly sampled as candidates at each split.
# #test with rus here:
# mtry_rus<- tuneRF(train_rus[1:(length(train_rus)-1)],train_rus$Hit, ntreeTry=500) # The trace specifies whether to print the progress of the search
#                #, trace=TRUE, plot=TRUE) #The plot specifies whether to plot the OOB error as function of mtry
# best.m <- mtry_rus[mtry_rus[, 2] == min(mtry_rus[, 2]), 1]
# print(best.m)
# print(best.m[1]) #first solution
# 
# fit_RF_rus <-randomForest(Hit~.,data=train_rus, mtry=best.m[1])
# print(fit_RF_smote)

#### hyperparam tuning: Random Search
# you can use the trainControl function to specify a number of parameters (including sampling parameters) in your model. 
#The object that is outputted from trainControl will be provided as an argument for RF train
control <- trainControl(method="repeatedcv", number=10, repeats=10, search="random") # Repeated K Fold Cross-Validation

#train on the class-balanced dataset
fit_RF_sbc <- randomForest::randomForest(Hit ~ ., data = train_sbc, tuneLength=20, trControl=control)
fit_RF_ros <- randomForest::randomForest(Hit ~ ., data=train_ros, tuneLength=20, trControl=control)
fit_RF_rus <- randomForest::randomForest(Hit ~ ., data=train_rus, tuneLength=20, trControl=control)
fit_RF_smote <- randomForest::randomForest(Hit ~ ., data=train_smote, tuneLength=20, trControl=control)
#tuneLength : It allows system to tune algorithm automatically. It indicates the number of different values to try for each tunning parameter.
#trControl : looks for the random search parameters established by caret::trainControl
#plot(fit_RF_sbc)

#predict class on the full dataset
pred_class_RF_sbc <- predict(fit_RF_sbc, LDA_VARS_HIT, type="response") #predict class on the full train dataset
pred_class_RF_ros <- predict(fit_RF_ros, LDA_VARS_HIT, type="response")
perf_RF_ros <- display(pred_class_RF_ros, LDA_VARS_HIT$Hit, "1")  
pred_class_RF_rus <- predict(fit_RF_rus, LDA_VARS_HIT, type="response")
pred_class_RF_smote <- predict(fit_RF_smote, LDA_VARS_HIT, type="response")
# #check performance 
# perf_RF_sbc <- display(pred_class_RF_sbc, LDA_VARS_HIT$Hit, "1")
# perf_RF_ros <- display(pred_class_RF_ros, LDA_VARS_HIT$Hit, "1")  
# perf_RF_rus <- display(pred_class_RF_rus, LDA_VARS_HIT$Hit, "1")
# perf_RF_smote <- display(pred_class_RF_smote, LDA_VARS_HIT$Hit, "1")
# RandForestPred <- rbind(perf_RF_ros, perf_RF_rus, perf_RF_sbc, perf_RF_smote)
# rownames(RandForestPred) <- c("Balanced by ROS", "Balanced by RUS",
#                       "Balanced by SBC", "Balanced by SMOTE")
# RandForestPred ##MAYBE ROS AND SMOTE ALWAYS GIVE 1 BECAUSE ARE THE OVERSAMPLING TECHNIQUES AND HAVE ALREADY SEEN THE FLL DATASET?
# #SMOTE is likely very poweful with such data structure but its quality is hard to evaluate on the training data itself

#As we can see, all of these four algorithms helped to improve the performance of random forest more or less. On the premise of ensuring a good accuracy of a classifier, the sensitivity, which is highly related to the minority class that we are more concerned, is improved by our resampling strategies.

## Evaluate variable importance
# importance(fit_RF_smote)
# varImpPlot(fit_RF_smote)


############ RULE EXTRACTION FOR RF ###################
# Stable and Interpretable RUle Set for RandomForests,
# Should output an easily interpretable ruleset for classification of a RandomForest algorithm
# compared to the RandomForest output, the SIRUS rule selection is more stable with respect to data perturbation.
# (Benard et al. 2021)

train_sbc$Hit   <- as.numeric(as.character(train_sbc$Hit))
train_ros$Hit   <- as.numeric(as.character(train_ros$Hit))
train_rus$Hit   <- as.numeric(as.character(train_rus$Hit))
train_smote$Hit <- as.numeric(as.character(train_smote$Hit))

# train SIRUS on the class-balanced dataset
#cross_val_p0 <- sirus.cv(train_smote[1:(length(train_smote)-1)], train_smote$Hit, type = "classif") #takes a while  #Estimate the optimal hyperparameter p0 used to select rules in sirus.fit using cross-validation (Benard et al. 2020, 2021).
fit_RFsirus_sbc   <- sirus.fit(train_sbc[1:(length(train_sbc)-1)], train_sbc$Hit, type = "classif")  # ,p0= cross_val_p0$p0.stab)
fit_RFsirus_ros   <- sirus.fit(train_ros[1:(length(train_ros)-1)], train_ros$Hit, type = "classif") 
fit_RFsirus_rus   <- sirus.fit(train_rus[1:(length(train_rus)-1)], train_rus$Hit, type = "classif") 
fit_RFsirus_smote <- sirus.fit(train_smote[1:(length(train_smote)-1)], train_smote$Hit, type = "classif")

#predict class on the full dataset
class_probs_RFsirus_sbc     <- sirus.predict(fit_RFsirus_sbc, LDA_VARS)
class_probs_RFsirus_ros     <- sirus.predict(fit_RFsirus_ros, LDA_VARS)
class_probs_RFsirus_rus     <- sirus.predict(fit_RFsirus_rus, LDA_VARS)
class_probs_RFsirus_smote   <- sirus.predict(fit_RFsirus_smote, LDA_VARS)

#ASSIGN CLASS (mimick the output of predict(fit_RF_sbc, LDA_VARS_HIT, type="response"))
#create empty vector
pred_class_RFsirus_sbc <- vector(mode="numeric", length=length(class_probs_RFsirus_sbc))
pred_class_RFsirus_ros <- vector(mode="numeric", length=length(class_probs_RFsirus_ros))
pred_class_RFsirus_rus <- vector(mode="numeric", length=length(class_probs_RFsirus_rus))
pred_class_RFsirus_smote <- vector(mode="numeric", length=length(class_probs_RFsirus_smote))
#assign 1 when the class probability is bigger than the proportion of the minority class
pred_class_RFsirus_sbc <- replace(pred_class_RFsirus_sbc, which(class_probs_RFsirus_sbc>fit_RFsirus_sbc$mean), 1)
pred_class_RFsirus_ros <- replace(pred_class_RFsirus_ros, which(class_probs_RFsirus_ros>fit_RFsirus_ros$mean), 1)
pred_class_RFsirus_rus <- replace(pred_class_RFsirus_rus, which(class_probs_RFsirus_rus>fit_RFsirus_rus$mean), 1)
pred_class_RFsirus_smote <- replace(pred_class_RFsirus_smote, which(class_probs_RFsirus_smote>fit_RFsirus_smote$mean), 1)
#ERROR? when looking at the CSI and F1 scores, it seems like the prediction score is very low
# this may be due to the pred_class_RF done in the wrong way? wrong mean probability used?

# table(LDA_VARS_HIT$Hit)
# table(pred_class_RFsirus_sbc)

pred_class_RFsirus_sbc <- as.factor(pred_class_RFsirus_sbc)
pred_class_RFsirus_ros <- as.factor(pred_class_RFsirus_ros)
pred_class_RFsirus_rus <- as.factor(pred_class_RFsirus_rus)
pred_class_RFsirus_smote <- as.factor(pred_class_RFsirus_smote)

#perf_class_RFsirus_sbc <- display(pred_class_RFsirus_sbc, LDA_VARS_HIT$Hit, "1")

#print the rule set (saved at the bottom of MAIN)
SirusRules <- list(SIRUS_rules_sbc =sirus.print(fit_RFsirus_sbc, digits = 5),
              SIRUS_rules_ros =sirus.print(fit_RFsirus_ros, digits = 5),
              SIRUS_rules_rus =sirus.print(fit_RFsirus_rus, digits = 5),
              SIRUS_rules_smote =sirus.print(fit_RFsirus_smote, digits = 5))

#return to factor after SIRUS
train_sbc$Hit <- as.factor(train_sbc$Hit)
train_ros$Hit <- as.factor(train_ros$Hit)
train_rus$Hit <- as.factor(train_rus$Hit)
train_smote$Hit <- as.factor(train_smote$Hit)

######################################################
###### CLASSIFIERS ON THE IMBALANCED DATASET #########
######################################################

####################################
#### LOGISTIC REGRESSION ###########
####################################

model1 = glm(LDA_VARS_HIT$Hit ~ . , data=LDA_VARS, family=binomial) # binomial as there are 2 classes
summary(model1)

#Boolean vector of the coefficients can indeed be extracted by:
toselect.x <- summary(model1)$coeff[-1,4] < 0.05
# select sig. variables
relevant.x <- names(toselect.x)[toselect.x == TRUE] 
# formula with only sig variables
sig.formula <- as.formula(paste("LDA_VARS_HIT$Hit ~",paste(relevant.x, collapse= "+"))) #sig.formula <- as.formula(paste("y ~",relevant.x))
sig.terms <- paste(relevant.x, collapse= "+") #sig.formula <- as.formula(paste("y ~",relevant.x))

#Update model with just the significant variables
model2 <- glm(formula=sig.formula,data=LDA_VARS,family=binomial)

#model2 = update(model1, ~ sig.terms)
summary(model2)

###Predict for training data and find training accuracy
pred.prob = predict(model2, type="response")
pred.prob = ifelse(pred.prob > 0.05, 1, 0)
table(pred.prob, LDA_VARS_HIT$Hit) #terrible...

#Rare event logistic regression?

#maybe  under-sampling nonevents for rare events data as in Logistic Regression for Massive Data with Rare Events
# However, when π0 > 0.1, the performances of the under-sampled
# estimators are almost as good as the full-data estimator.
#see WANG, HaiYing. Logistic regression for massive data with rare events. In: International Conference on Machine Learning. PMLR, 2020. p. 9829-9836.

##############################
### LDA / QDA ################
##############################

lda <- qda(LDA_VARS_HIT$Hit ~ ., 
           LDA_VARS)

## proportion of LDs that encapsulate the variation (in case of 1 Ld, the prop.lda = 1)
#prop.lda = lda$svd^2/sum(lda$svd^2)

## get / compute LDA scores from LDA coefficients / loadings

#train and prediction on the full, imbalanced, dataset
plda <- predict(object = lda,
                newdata = LDA_VARS) #predict(lda_TEST)$x

#create a histogram of the discriminant function values
#par(mfrow=c(1,1), oma=c(0,0,2,0),mar = c(4.5, 3.8, 1, 1.1))
#ldahist(data = plda$x[,1], g=LDA_VARS_HIT$Hit,nbins = 80) + mtext("LDA hist", line=0, side=3, outer=TRUE)

##############################
### RANDOM FOREST ############
##############################

# not included: when RF is trained and tested on the same dataset, it shows always 100% accuracy

###############################################
### kernel discriminant analysis ################
################################################
# Sparse KOS requires the construction of a n × n kernel
# matrix K and is therefore computationally prohibitive
# for large n cases.
###TAKES SEVERAL MINUTES TO COMPUTE ON A SMALL DATATEST OF 4 FEATURES AND 179 OBSERVATIONS
#UNUSABLE.....
###SelectParams( LDA_VARS, LDA_VARS_HIT$Hit, Sigma = NULL, Gamma = NULL)
# 
# SelectParams(Data = Data$TrainData,
#              Cat = Data$CatTrain)
# 
# Predict( Data = Data$TrainData,
#          Cat = Data$CatTrain)
# POSSIBILITY OF USING A SIMPLER VERSION? KOS INSTEAD OF sparseKOS? .....

######## ASSIGN THE PREDICTIONS TO THE dataframe #################
LDA_VARS_HIT_PRED <-  cbind("REPLICATE" = summary_AUTO_NAOmit_transf$REPLICATE,
                            "PERIOD" = summary_AUTO_NAOmit_transf$PERIOD,
                            "pair" = summary_AUTO_NAOmit_transf$pair,
                            "int_start_frame "= summary_AUTO_NAOmit_transf$int_start_frame,
                            "int_end_frame" = summary_AUTO_NAOmit_transf$int_end_frame,
                            LDA_VARS_HIT,
                            # CLASSIFICATION on the imbalanced dataset
                            #NOTE: train and predict done on the training dataset
                            # Quadratic Discriminant analysis
                            "QDA_pred_Hit" =  plda$class,
                            # CLASSIFICATION on dataset balanced with a variety of over and under sampling techniques
                            ## QDA
                            "QDA_SBC_pred_Hit" = pred_class_QDA_sbc$class, #Under-Sampling Based on Clustering (SBC)
                            "QDA_ROS_pred_Hit" = pred_class_QDA_ros$class, #random over-sampling algorithm (ROS)
                            "QDA_RUS_pred_Hit" = pred_class_QDA_rus$class, #random under-sampling algorithm (RUS)
                            "QDA_SMOTE_pred_Hit" = pred_class_QDA_smote$class, #synthetic minority over-sampling technique (SMOTE)
                            # pred_class_QDA_sbc$class
                            ## RandomForest 
                            #NOTE: test done on the FULL training dataset, train done on under or over-sampled dataset
                            "RF_SBC_pred_Hit" = pred_class_RF_sbc, #Under-Sampling Based on Clustering (SBC)
                            "RF_ROS_pred_Hit" = pred_class_RF_ros, #random over-sampling algorithm (ROS)
                            "RF_RUS_pred_Hit" = pred_class_RF_rus, #random under-sampling algorithm (RUS)
                            "RF_SMOTE_pred_Hit" = pred_class_RF_smote, #synthetic minority over-sampling technique (SMOTE)
                            ## Sirus RandomForest
                            "RFsirus_SBC_pred_Hit" = pred_class_RFsirus_sbc, #Under-Sampling Based on Clustering (SBC)
                            "RFsirus_ROS_pred_Hit" = pred_class_RFsirus_ros, #random over-sampling algorithm (ROS)
                            "RFsirus_RUS_pred_Hit" = pred_class_RFsirus_rus, #random under-sampling algorithm (RUS)
                            "RFsirus_SMOTE_pred_Hit" = pred_class_RFsirus_smote #synthetic minority over-sampling technique (SMOTE)
)


#################################################
### Algorithm quality scores ####################
#################################################

### critical success index (CSI)
# # CSI is a verification measure of categorical forecast performance
# # calculated as CSI= (hits)/(hits + false alarms + misses)
# # The CSI is not affected by the number of non-event forecasts that verify (correct rejections).

### F1 score
# # F1 score is a performance metric useful when predicting class labels, where the positive class is more important and FNeg and FPos are equally costly. ###

#convert int_start and int_end back to frames
LDA_VARS_HIT_PRED$int_start_frame <- as.numeric(LDA_VARS_HIT_PRED$int_start_frame)
LDA_VARS_HIT_PRED$int_end_frame <- as.numeric(LDA_VARS_HIT_PRED$int_end_frame)

CSI_val <- NULL
CSI_scores_PRED_REP_PER <- NULL
CSI_scores_ALL <- NULL
CSI_REP_PER_name <- NULL

F1_val <- NULL
F1_scores_PRED_REP_PER <- NULL
F1_scores_ALL <- NULL
F1_REP_PER_name <- NULL

# test all the predictions frame by frame
for (PRED_HIT in grep("pred_Hit", names(LDA_VARS_HIT_PRED), value=TRUE)) {

#create matrices for each REP_PER
for (REPLICATE in c("R3SP","R9SP")) 
{
  for (PERIOD in c("pre","post"))
    {

  IF_frames_temp <- get(grep(paste0("IF_frames","_", REPLICATE,PERIOD),ls(),value=TRUE))
  int_mat_manual_temp <- get(grep(paste0("int_mat_manual","_", REPLICATE,PERIOD),ls(),value=TRUE))
  
  #this section replicates what is done in the matrix calculation
  int_mat_auto_trim <- matrix(0L, nrow = dim(ids_pairs)[1], ncol = (dim(IF_frames_temp)[1] ))
  rownames(int_mat_auto_trim) <- ids_pairs$pair
  colnames(int_mat_auto_trim) <- c(IF_frames_temp$frame_num)
  
  #Create subset of data for the matrix construction
  trimmed_AUTO_pred <- LDA_VARS_HIT_PRED[which(LDA_VARS_HIT_PRED[,PRED_HIT] == 1 & LDA_VARS_HIT_PRED$REPLICATE==REPLICATE & LDA_VARS_HIT_PRED$PERIOD==PERIOD),]
  trimmed_AUTO_Hit  <- LDA_VARS_HIT[which(LDA_VARS_HIT$Hit == 1 & LDA_VARS_HIT_PRED$REPLICATE==REPLICATE & LDA_VARS_HIT_PRED$PERIOD==PERIOD),]
  
  ## Trimmed Auto
    if (nrow(trimmed_AUTO_Hit)>0) {
      if (nrow(trimmed_AUTO_pred)>0) {
        print("pred and Hit >0")
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
    CSI <- (TruePositive/(TruePositive+FalseNegative+FalsePositive))
    CSI_val <- c(CSI_val,CSI)
    
    #F1 Score
    F1 <- TruePositive/(TruePositive+0.5*(FalsePositive+FalseNegative))
    F1_val <- c(F1_val,F1)
    
    CSI <- NULL
    F1 <- NULL
      }else{print("pred = 0 and Hit >0, prediction score is 0 (no matrix has been built)")
            CSI_val <- c(CSI_val,0)
            F1_val <- c(F1_val,0)
      }
    }else{ # if trimmed_AUTO_pred = 0 and trimmed_AUTO_Hit= 0
      print("Hit = 0, there is nothing to predict here")
      CSI_val <- c(CSI_val,NA)
      F1_val <- c(F1_val,NA)
      }
  
  #assign REP_PER name
  #CURRENTLY NOT ASSIGNED IN THE RIGHT ORDER!!
  CSI_REP_PER_name <-paste0("CSI_",MAN_int$REPLICATE, "-", MAN_int$PERIOD)
  F1_REP_PER_name <-paste0("F1_",MAN_int$REPLICATE, "-", MAN_int$PERIOD)
    
  }## REPLICATE
    }##PERIOD
#create DF, useful for the summary DF at the end of the loop
CSI_scores_PRED_REP_PER <- data.frame( REP_PER = CSI_REP_PER_name, Freq = CSI_val, PRED_HIT)
F1_scores_PRED_REP_PER <- data.frame( REP_PER = F1_REP_PER_name, Freq = F1_val, PRED_HIT)

#stack 
CSI_scores_ALL <- rbind(CSI_scores_ALL,CSI_scores_PRED_REP_PER)
F1_scores_ALL <- rbind(F1_scores_ALL,F1_scores_PRED_REP_PER)


#t(column_to_rownames(CSI_scores,"REP_PER"))
CSI_val <- NULL
F1_val <- NULL
}


#calculate overall score
# do it in a more clever way, eg. using AGGREGATE
#perc_CSI <- sum(CSI_val,na.rm = T)/length(CSI_val[!is.na(CSI_val)])
CSI <- setnames(aggregate(CSI_scores_ALL$Freq, list(CSI_scores_ALL$PRED_HIT), FUN=mean,na.rm=T), c("Classifier","CSI_score"))
F1 <- setnames(aggregate(F1_scores_ALL$Freq, list(F1_scores_ALL$PRED_HIT), FUN=mean,na.rm=T), c("Classifier","F1_score"))

CSI$CSI_score <- round(CSI$CSI_score,3)
F1$F1_score <- round(F1$F1_score,3)

#transpose to incorporate in the final output and rename the columns with the quality metric name
tCSI <- setNames(data.frame(t(CSI[,-1])), CSI[,1])
names(tCSI) <- gsub(x = names(tCSI), pattern = "pred_Hit", replacement = "CSI")  
tF1 <- setNames(data.frame(t(F1[,-1])), F1[,1])
names(tF1) <- gsub(x = names(tF1), pattern = "pred_Hit", replacement = "F1")  




##########################################################################
##########################################################################
##### EXTRA STUFF FOR PCA AND LDAs #######################################

## Ongoing work:

# Effect size
# Some suggest the use of eigenvalues as effect size measures, however, this is generally not supported.[9] Instead, the canonical correlation is the preferred measure of effect size. It is similar to the eigenvalue, but is the square root of the ratio of SSbetween and SStotal. It is the correlation between groups and the function.[9] Another popular measure of effect size is the percent of variance[clarification needed] for each function. This is calculated by: (λx/Σλi) X 100 where λx is the eigenvalue for the function and Σλi is the sum of all eigenvalues. This tells us how strong the prediction is for that particular function compared to the others
# 

cat(paste(" from https://en.wikipedia.org/wiki/Linear_discriminant_analysis#LDA_for_two_classes
\nComparison to logistic regression

\nDiscriminant function analysis is very similar to logistic regression, and both can be used to answer the same research questions.[9] Logistic regression does not have as many assumptions and restrictions as discriminant analysis. However, when discriminant analysis’ assumptions are met, it is more powerful than logistic regression.[28] Unlike logistic regression, discriminant analysis can be used with small sample sizes. It has been shown that when sample sizes are equal, and homogeneity of variance/covariance holds, discriminant analysis is more accurate.[7] Despite all these advantages, logistic regression has none-the-less become the common choice, since the assumptions of discriminant analysis are rarely met.[8][7] 
          ",sep = " "))

  
### PLOTTING PCA AND LDA TOGETHER


#dataset = data.frame(Hit = LDA_VARS_HIT[,"Hit"],
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
# klaR::partimat(Hit~mean_movement_angle_diff+prop_time_undetected_REC+StDev_angle_ACT, data=LDA_VARS_HIT, method="qda",
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
  #LDA_VARS <- LDA_VARS[, -( grep("jerk" , colnames(LDA_VARS),perl = TRUE) ) ]
  


##ACCURACY OF CLASSIFICATION #NOT TOO RELEVANT
# accuracy  <- xtabs(~LDA_VARS_HIT$Hit+plda$class)
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

