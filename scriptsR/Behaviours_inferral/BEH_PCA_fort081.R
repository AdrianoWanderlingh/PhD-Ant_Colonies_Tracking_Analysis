print(paste("PERFORM PCA ",REPLICATE, PERIOD))

library(FactoMineR)
library(factoextra)
library(missMDA) #PCA with missing values

### PCAs ##########
#DO PCA OF summary_AUTO_REP_PER
#remove ant names/rep/int/etc in pca (keep only vars)

#PCA with missing values according to Dray & Josse (2015) https://doi.org/10.1007/s11258-014-0406-z
#imputePCA of the R package missMDA
#Impute the missing values of a dataset with the Principal Components Analysis model. 
#Can be used as a preliminary step before performing a PCA on an incomplete dataset.
#missing values are replaced by random values, and then PCA is applied on the completed data set, and missing values are then updated by the fitted values

summary_AUTO_REP_PER_dget <- dget("/home/cf19810/Documents/Ants_behaviour_analysis/Data/summary_AUTO_REP_PER_16feb22.txt") # load file created with dput 

#data structure
#Once the THRESH for $agreement has been established:
# Hit == 1 is a True Positive
# Hit ==0 is a False Positive
# this is done in BEH_Auto_Man_agreement_matrix_fort081.R

d       <- data.frame(V1 = sample(1:100, 10), V2 = sample(1:100, 10), 
                      V3 = sample(1:100, 10))
result  <- prcomp(d, center = TRUE, scale = TRUE, na.action = na.omit)
result$x                                # $
d$V1[5] <- NA                           # $


summary_PCA_vars <- summary_AUTO_REP_PER_dget[, -match(c("REPLICATE", "PERIOD","INT","ant1","ant2","pair","int_start_frame","int_end_frame","mean_prop_time_undetected","agreement","disagreement","Hit"), names(summary_AUTO_REP_PER_dget))] 

## Imputation
res.comp <- imputePCA(summary_PCA_vars,method="Regularized")

## A PCA can be performed
res.pca <- PCA(res.comp$completeObs)
# Extract eigenvalues/variances
get_eig(res.pca)
# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
# Extract the results for variables
var <- get_pca_var(res.pca)
var

library("corrplot")
corrplot(var$cos2, is.corr = FALSE)

# Coordinates of variables
head(var$coord)
# Contribution of variables
head(var$contrib)
# Graph of variables
# Control variable colors using their contributions
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# Extract the results for individuals #GIVE ROWNAME TO DO THIS (I.E. HIT OR INT)
ind <- get_pca_ind(res.pca)
ind
# Coordinates of individuals
head(ind$coord)



# Use habillage to specify groups for coloring
fviz_pca_ind(res.pca, pointsize =summary_AUTO_REP_PER_dget$interaction_length_secs,
             label = "none", # hide individual labels
             title = "PCA - indivdual interactions \nPointsize by Int Length",
             habillage = as.factor(summary_AUTO_REP_PER_dget$Hit), # color by groups
             addEllipses = TRUE
)

