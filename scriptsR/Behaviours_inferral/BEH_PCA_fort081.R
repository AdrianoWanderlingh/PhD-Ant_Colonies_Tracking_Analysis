##########################################################################################
############## BEH DISCRIMINANT ANALYSIS #################################################
##########################################################################################

#script dependant on BEH_MAIN_behaviours_analysis_fort081.R

#For previous versions of this script, explore: 
# https://github.com/AdrianoWanderlingh/PhD-exp1-data-analysis/tree/main/scriptsR/Behaviours_inferral


print(paste("PERFORM LDA ",unique(interaction_MANUAL$PERIOD), unique(interaction_MANUAL$REPLICATE)))

# FORECAST VERIFICATION
# https://www.swpc.noaa.gov/sites/default/files/images/u30/Forecast%20Verification%20Glossary.pdf
# https://en.wikipedia.org/wiki/Sensitivity_and_specificity

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


###############################################################################
###### NORMALISE VARS BEFORE LDA ##############################################
###############################################################################
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
summary_LDA_vars <- summary_AUTO_transf[, -match(c("REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","disagreement","Hit"), names(summary_AUTO_transf))] 
sapply(summary_LDA_vars, function(x) sum(is.na(x)))
sapply(summary_LDA_vars, function(x) sum(is.infinite(x))) 
summary_LDA_vars_hit <- cbind(summary_LDA_vars,Hit=summary_AUTO_transf$Hit)

#summary_LDA_vars_trim <- summary_AUTO_transf[, -match(c(TRIM,"REPLICATE", "PERIOD","INT","ACT","REC","pair","int_start_frame","int_end_frame","agreement","disagreement","Hit"), names(summary_AUTO_transf))] 
## eventually trimunwanted vars
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
metric <- "euclidean"# manhattan #using it instead of euclidean as it is the most common choice as sugegsted by Le et al.2019 
#BUT manhattan is suggested for high-dimensional data (N of features > )

#this takes a few seconds to compute

#Multisurf 
neighbor.idx.observed <- find.neighbors(summary_LDA_vars, summary_LDA_vars_hit$Hit, k = 0, method = "multisurf")
results.list <- stir(summary_LDA_vars, neighbor.idx.observed, k = k, metric = metric, method = "multisurf")
t_sorted_multisurf <- results.list$STIR_T
t_sorted_multisurf$attribute <- rownames(t_sorted_multisurf)
round(t_sorted_multisurf$t.pval.adj,5)
#select only variables with p<0.05
RELIEF_selected <- t_sorted_multisurf[which(t_sorted_multisurf$t.pval.adj<0.05),"attribute"]
summary_LDA_vars_RELIEF <- summary_LDA_vars[, match(c(RELIEF_selected), names(summary_LDA_vars))] 
#add the classes
summary_LDA_vars_RELIEF_hit <- cbind(summary_LDA_vars_RELIEF,Hit=summary_AUTO_transf$Hit)



#TODOs: test it on untransformed variables, also: why Nas??



##### corelearn can also be used to try other models using CoreModel

#######################################################
######### CHOOSE DATA #################################
#######################################################

LDA_DATA <- "RELIEF" # or "ALL_VARS" 

if(LDA_DATA =="ALL_VARS"){LDA_VARS <- summary_LDA_vars; LDA_VARS_HIT <- summary_LDA_vars_hit}
if(LDA_DATA =="RELIEF"){LDA_VARS <- summary_LDA_vars_RELIEF; LDA_VARS_HIT <- summary_LDA_vars_RELIEF_hit}


#######################################################
######### MISSING CASES ###############################
#######################################################

#### ISSUEEE
#LOTS of missing values affecting the analysis
#N of rows including missing cases:
sum(!complete.cases(LDA_VARS))

#prop of NAs over total in % (missing values)
propNAs_total <- (sum(sapply(LDA_VARS, function(x) sum(is.na(x)))) /(nrow(LDA_VARS) * ncol(LDA_VARS))) *100
##prop missing by variable in %
propNAs_byvar <- (sapply(LDA_VARS, function(x) sum(is.na(x)))/nrow(LDA_VARS)) * 100

##prop missing by row (show only top) - potentially droppable rows?
prop_Na_row <- (apply(LDA_VARS, 1, function(x) sum(is.na(x)))/ncol(LDA_VARS))*100
prop_Na_row <- data.frame(prop_missing=sort(prop_Na_row, decreasing=TRUE))
propNAs_byrow_25perc <- ((length(prop_Na_row[which(prop_Na_row>25),]))/nrow(LDA_VARS))*100 #show all rows with values over 25% missing
hist(prop_Na_row$prop_missing) #+ abline(v=20, col="blue")
#summary of Nas distribution
cat("prop of NAs over total in %",propNAs_total, "\n\nProp missing by variable \n",propNAs_byvar,"\n\nprop of rows with more than 25% missing",propNAs_byrow_25perc)

#######################################################
######### TEST DATA ##################################
#######################################################

### Homogeneity of variance/covariance (homoscedasticity) ###
#Box’s M Test is extremely sensitive to departures from normality; the fundamental test assumption is that your data is multivariate normally distributed. Therefore, if your samples don’t meet the assumption of normality, you shouldn’t use this test.
# mvn(data = LDA_VARS,subset = LDA_VARS_HIT$Hit,mvnTest="mardia") #crashes R

##### Checking Assumption of Equal Variance-Covariance matrices
#Checking the Assumption of Equal Covariance Ellipse
heplots::covEllipses(LDA_VARS, 
                     LDA_VARS_HIT$Hit, 
                     fill = TRUE, 
                     pooled = FALSE, 
                     col = c("blue", "red"), 
                     variables = c(1:3), 
                     fill.alpha = 0.05)
#using the BoxM test in order to check our assumption of homogeneity of variance-covariance matrices
# for p<0.05 there is a problem of heterogeneity of variance-covariance matrices
boxm <- heplots::boxM(LDA_VARS, LDA_VARS_HIT$Hit)
boxm #NOT met for the summary_LDA_vars_RELIEF
plot(boxm)

### Checking Assumption of Normality
Hit.yes <- subset(LDA_VARS_HIT, Hit == 1) 
Hit.no <- subset(LDA_VARS_HIT, Hit == 0)

par(mfrow = c(2, 2)) 
for(i in names(LDA_VARS)) { 
  qqnorm(Hit.yes[[i]]); qqline(Hit.yes[[i]], col = 2) 
}

for(i in names(LDA_VARS)) { 
  qqnorm(Hit.no[[i]]); qqline(Hit.no[[i]], col = 2) 
}
##"prop_time_undetected_REC.sqrt_x is very not normal!


# As covariances are not equal, quadratic discriminant analysis will be used.
############# ASSUMPTIONS FOR QDA #################

# - Observation of each class is drawn from a normal distribution (same as LDA).
#   Just showed to be not valid for prop_time_undetected_REC.sqrt_x!!!

# - QDA assumes that each class has its own covariance matrix (different from LDA).
#   this one is met

####################################
#### LOGISTIC REGRESSION ###########
####################################
paste(names(LDA_VARS),sep = "+")

model1 = glm(LDA_VARS_HIT$Hit ~ . , data=LDA_VARS, family=binomial) # binomial as there are 2 classes
summary(model1)

model2 = update(model1, ~ .-mean_Mov_Orient_delta_angle_REC.orderNorm)
summary(model2)

###Predict for training data and find training accuracy
pred.prob = predict(model2, type="response")
pred.prob = ifelse(pred.prob > 0.5, 1, 0)
table(pred.prob, LDA_VARS_HIT$Hit) #terrible...

#Rare event logistic regression?


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

# #############################
# ####### PCA #################
# #############################
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
# ## 2. PCA on the base data for LDA_VARS_trim
# RES.PCA <- PCA(LDA_VARS_trim) #more explained variance without the trimmeed vars
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
# pca_prcomp <- prcomp(LDA_VARS,
#               center = TRUE,
#               scale. = TRUE) 
# 
# prop.RES.PCA = pca_prcomp$sdev^2/sum(pca_prcomp$sdev^2)


##############################
### LDA / QDA ################
##############################

#check if I have more predictor vars than the size of the smallest group
table(LDA_VARS_HIT$Hit)
length(LDA_VARS)
#it is a bit risky as they are very close...

lda <- lda(LDA_VARS_HIT$Hit ~ ., 
           LDA_VARS)

## proportion of LDs that encapsulate the variation (in case of 1 Ld, the prop.lda = 1)
#prop.lda = lda$svd^2/sum(lda$svd^2)

## get / compute LDA scores from LDA coefficients / loadings

#PREDICTION SHOULD BE PERFORMED IN THE SECOND HALF OF THE DATASET? (NOT HEREBY ANALYSED)
plda <- predict(object = lda,
                newdata = LDA_VARS) #predict(lda_TEST)$x


#create a histogram of the discriminant function values
par(mfrow=c(1,1), oma=c(0,0,2,0),mar = c(4.5, 3.8, 1, 1.1))
ldahist(data = plda$x[,1], g=LDA_VARS_HIT$Hit,nbins = 80) + mtext("LDA hist", line=0, side=3, outer=TRUE)
## ASSIGN THE PREDICTION TO THE dataframe
LDA_VARS_HIT_PRED <-  cbind("REPLICATE" = summary_AUTO_transf$REPLICATE,
                                    "PERIOD" = summary_AUTO_transf$PERIOD,
                                    "pair" = summary_AUTO_transf$pair,
                                    "int_start_frame "= summary_AUTO_transf$int_start_frame,
                                    "int_end_frame" = summary_AUTO_transf$int_end_frame,
                                    LDA_VARS_HIT, 
                                    "Predicted_Hit" =  plda$class)

#########################################
### critical success index (CSI) ########
#########################################

# # CSI is a verification measure of categorical forecast performance
# # calculated as CSI= (hits)/(hits + false alarms + misses)
# # The CSI is not affected by the number of non-event forecasts that verify (correct rejections).

#convert int_start and int_end back to frames
LDA_VARS_HIT_PRED$int_start_frame <- as.numeric(LDA_VARS_HIT_PRED$int_start_frame)
LDA_VARS_HIT_PRED$int_end_frame <- as.numeric(LDA_VARS_HIT_PRED$int_end_frame)

CSI_val <- NULL
CSI_scores <- NULL
CSI_REP_PER_name <- NULL

#create matrix
for (REPLICATE in c("R3SP","R9SP")) 
{
  for (PERIOD in c("pre","post"))
    {
  
  trimmed_AUTO_pred <- LDA_VARS_HIT_PRED[which(LDA_VARS_HIT_PRED$Predicted_Hit == 1 & LDA_VARS_HIT_PRED$REPLICATE==REPLICATE & LDA_VARS_HIT_PRED$PERIOD==PERIOD),]
  
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
    CSI <- (TruePositive/(TruePositive+FalseNegative+FalsePositive))*100
    CSI_val <- c(CSI_val,CSI)
    
    CSI <- NULL
    }else{ # if trimmed_AUTO_pred contains no Hit 
      CSI_val <- c(CSI_val,0)
    }
  #assign REP_PER name
  CSI_REP_PER_name <-paste0("CSI_",MAN_int$REPLICATE, "-", MAN_int$PERIOD)
    }## REPLICATE
    }##PERIOD
#create DF, useful for the summary DF at the end of the loop
CSI_scores <- data.frame( REP_PER = CSI_REP_PER_name, Freq = CSI_val)
#calculate overall score
perc_CSI <- sum(CSI_val)/length(CSI_val)

#t(column_to_rownames(CSI_scores,"REP_PER"))

###########################################################################
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

