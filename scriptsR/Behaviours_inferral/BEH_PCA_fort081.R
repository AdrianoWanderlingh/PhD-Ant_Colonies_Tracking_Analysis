print(paste("PERFORM PCA ",REPLICATE, PERIOD))

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

#summary_AUTO_REP_PER_dget <- dget("/home/cf19810/Documents/Ants_behaviour_analysis/Data/summary_AUTO_REP_PER_16feb22.txt") # load file created with dput 
summary_AUTO_REP_PER_dget <- summary_AUTO_REP_PER
summary_AUTO_REP_PER_dget$Hit <- as.factor(summary_AUTO_REP_PER_dget$Hit)
#data structure
#Once the THRESH for $agreement has been established:
# Hit == 1 is a True Positive
# Hit ==0 is a False Positive
# this is done in BEH_Auto_Man_agreement_matrix_fort081.R


#check N of missing values per variable
TRIM <- c("mean_jerk_PxPerSec3_ant1","mean_jerk_PxPerSec3_ant2","mean_accel_PxPerSec2_ant1","mean_accel_PxPerSec2_ant2","mean_abs_turnAngle_ant1","mean_abs_turnAngle_ant2")
#remove ant names/rep/int/etc in pca (keep only vars)
summary_PCA_vars <- summary_AUTO_REP_PER_dget[, -match(c("REPLICATE", "PERIOD","INT","ant1","ant2","pair","int_start_frame","int_end_frame","mean_prop_time_undetected","agreement","disagreement","Hit"), names(summary_AUTO_REP_PER_dget))] 
sapply(summary_PCA_vars, function(x) sum(is.na(x))) #drop mean_jerk_PxPerSec3 and mean_accel_pxpersec2  (they both don't contribute much)

summary_PCA_vars_trim <- summary_AUTO_REP_PER_dget[, -match(c(TRIM,"REPLICATE", "PERIOD","INT","ant1","ant2","pair","int_start_frame","int_end_frame","mean_prop_time_undetected","agreement","disagreement","Hit"), names(summary_AUTO_REP_PER_dget))] 
summary_PCA_vars_trim_hit <- cbind(summary_PCA_vars_trim,Hit=summary_AUTO_REP_PER_dget$Hit)
sapply(summary_PCA_vars_trim, function(x) sum(is.na(x))) #drop mean_jerk_PxPerSec3 and mean_accel_pxpersec2  (they both don't contribute much)


###########plotting ################

#transform to long format
summary_PCA_long <- melt(summary_PCA_vars_trim_hit,id.vars=c("Hit")) #explanation on the warning message https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message


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

####transform variables??


#############################
####### PCA #################
#############################

#PCA with missing values according to Dray & Josse (2015) https://doi.org/10.1007/s11258-014-0406-z
#imputePCA of the R package missMDA
#Impute the missing values of a dataset with the Principal Components Analysis model. 
#Can be used as a preliminary step before performing a PCA on an incomplete dataset.
#missing values are replaced by random values, and then PCA is applied on the completed data set, and missing values are then updated by the fitted values
## Imputation for summary_PCA_vars which contains NAs
res.comp <- imputePCA(summary_PCA_vars,method="Regularized")

## 1. A PCA can be performed on the imputed data for summary_PCA_vars
res.pca1 <- PCA(res.comp$completeObs)

## 2. PCA on the base data for summary_PCA_vars_trim
RES.PCA <- PCA(summary_PCA_vars_trim) #more explained variance without the trimmeed vars

#for (RES.PCA in c(res.pca1,res.pca2)) {

# Extract eigenvalues/variances
get_eig(res.pca2)
# Visualize eigenvalues/variances
fviz_screeplot(RES.PCA, addlabels = TRUE, ylim = c(0, 50))
# Extract the results for variables
var <- get_pca_var(RES.PCA)
var

corrplot(var$cos2, is.corr = FALSE)

# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Graph of variables
# Control variable colors using their contributions
fviz_pca_var(RES.PCA, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)
# Contributions of variables to PC1
fviz_contrib(RES.PCA, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(RES.PCA, choice = "var", axes = 2, top = 10)
# Extract the results for individuals #GIVE ROWNAME TO DO THIS (I.E. HIT OR INT)
ind <- get_pca_ind(RES.PCA)
ind
# Coordinates of individuals
head(ind$coord)



# Use habillage to specify groups for coloring
fviz_pca_ind(RES.PCA, pointsize =summary_AUTO_REP_PER_dget$disagreement,
             label = "none", # hide individual labels
             title = "PCA - indivdual interactions \nPointsize by disagreement rate",
             habillage = as.factor(summary_AUTO_REP_PER_dget$Hit), # color by groups
             addEllipses = TRUE
)

#}

#LDA: Linear discriminant analysis also called DFA D.function.A.
#LDA focuses on finding a feature subspace that maximizes the separability between the groups. 

pca_prcomp <- prcomp(summary_PCA_vars_trim,
              center = TRUE,
              scale. = TRUE) 

prop.RES.PCA = pca_prcomp$sdev^2/sum(pca_prcomp$sdev^2)


lda <- lda(Hit ~ ., 
           summary_PCA_vars_trim_hit)

prop.lda = lda$svd^2/sum(lda$svd^2)

#get / compute LDA scores from LDA coefficients / loadings

plda <- predict(object = lda,
                newdata = summary_PCA_vars_trim_hit)

dataset = data.frame(Hit = summary_PCA_vars_trim_hit[,"Hit"],
                     pca = pca_prcomp$x, lda = plda$x)
#create a histogram of the discriminant function values
ldahist(data = plda$x[,1], g=summary_PCA_vars_trim_hit$Hit)

#plotting

# p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = Hit, shape = Hit), size = 2.5) + 
#   labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
#        y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))
LDmeans <- plyr::ddply(dataset, "Hit", summarise, grp.mean=mean(LD1))


p1 <- ggplot(dataset, aes(LD1, fill = Hit)) + geom_density(alpha = 0.2)+ 
      geom_vline(data=LDmeans, aes(xintercept=grp.mean, color=Hit), linetype="dashed") +
      labs(title = "Density plot for predicted values from model function, by Hit rate",
           x = paste("proportion of discriminability explained LD1 (", percent(prop.lda[1]), ")", sep=""))


p2 <- ggplot(dataset) + geom_point(aes(pca.PC1, pca.PC2, colour = Hit, shape = Hit), size = 2.5) +
  labs(x = paste("PC1 (", percent(prop.RES.PCA[1]), ")", sep=""),
       y = paste("PC2 (", percent(prop.RES.PCA[2]), ")", sep=""))

grid.arrange(p1, p2)


#------------------------
library(caret)
data(mdrr)
mdrrDescr <- mdrrDescr[, -nearZeroVar(mdrrDescr)]
mdrrDescr <- mdrrDescr[, -findCorrelation(cor(mdrrDescr), .8)]
set.seed(1)
inTrain <- createDataPartition(mdrrClass, p = .75, list = FALSE)[,1]
train <- mdrrDescr[ inTrain, ]
test  <- mdrrDescr[-inTrain, ]
trainClass <- mdrrClass[ inTrain]
testClass  <- mdrrClass[-inTrain]

set.seed(2)
ldaProfile <- rfe(train, trainClass,
                  sizes = c(1:10, 15, 30),
                  rfeControl = rfeControl(functions = ldaFuncs, method = "cv"))


postResample(predict(ldaProfile, test), testClass)
ldaProfile$optVariables

#0-------------------
#https://stackoverflow.com/questions/68307682/r-lda-linear-discriminant-analysis-how-to-get-compute-lda-scores-from-lda-co
#or use
?lda
lda_TEST <- lda(x = summary_PCA_vars_trim_hit[, 1:17], grouping = summary_PCA_vars_trim_hit$Hit)
lda_TEST$scores <- predict(lda_TEST)$x

lda_TEST_pred <- predict(lda_TEST)
accuracy  <- xtabs(~summary_PCA_vars_trim_hit$Hit+lda_TEST_pred$class)
#ACCURACY OF CLASSIFICATION
sum(accuracy[row(accuracy) == col(accuracy)]) / sum(accuracy)



#https://www.andreaperlato.com/mlpost/linear-discriminant-analysis/

#https://rpubs.com/ifn1411/LDA


#----------------------
# lda$svd explains the ratio of between- and within-group variation.
lda_TEST$svd


#--------------------
#If there are many groups, or if you quickly want to find the maximum probability for each sample, this command will help:
  
round(apply(lda_TEST_pred$posterior, MARGIN=1, FUN=max), 2)

#For this data set, most of these maximum probabilities are large (>0.9), indicating that most samples are 
#confidently assigned to one group. If most of these probabilities are large, the overall classification is 
#successful.

#The scaling value of the lda object gives the loadings (also called the slopes, coefficients, or weights) of each variable on each discriminant function.
round_scaling <- round(lda_TEST$scaling, 2)
# round_scaling$var  <- row.names(round_scaling)
names <- rownames(round_scaling)
rownames(round_scaling) <- NULL
round_scaling <- as.data.frame(cbind(names,round_scaling))
round_scaling <- round_scaling[order(round_scaling$LD1),]
write.table(round_scaling,paste0(DATADIR,"round_scaling.txt"),sep="\t",row.names=FALSE)

#-------------------------

#####
#From the Partition Plot above, we can see classification for eachof observation in the training dataset based 
#on the Linear Discriminant Analysis Model, and for every combination of two variables.

X11(width=15, height=15)
# Adjust the numbers for the option mar, if needed; see ?partimat and ?par 
# for more help on margins
par(mfrow=c(2,3))
library(klaR)
klaR::partimat(Hit~mean_movement_angle_diff+prop_time_undetected_ant2+StDev_angle_ant1, data=summary_PCA_vars_trim_hit, method="qda",
               main = "\nPartition Plot for Quadratic Discriminant Analysis Model \nred=incorrectly classified",mar=c(5,4,2,2))


# par(mfrow = c(4,2))
# for(i in 3:10)
#   drawparti(mydata[,1], mydata[,2], mydata[,i], method = "rda",
#             gamma=0, lamda=1)


# look at 
# http://rstudio-pubs-static.s3.amazonaws.com/84669_cd15214061d44e1493ffee69c5d55925.html
# https://rpubs.com/ifn1411/LDA
# https://www.andreaperlato.com/mlpost/linear-discriminant-analysis/
# https://towardsdatascience.com/linear-discriminant-analysis-explained-f88be6c1e00b


#SHOULD I GO FOR LINEAR REGRESSION HAVING ONLY 2 FACTORS VARIABLE?

