library(mice)
library(qgcomp)
library(imputeLCMD) # qPCR imputation

# inspired by the package produced in https://ehp.niehs.nih.gov/doi/10.1289/EHP5838 to extend MICE to leftcenslognorm

## from documentation pool {mice}
# The typical sequence of steps to perform a multiple imputation analysis is:
# 1. Impute the missing data by the mice() function, resulting in a multiple imputed data set (class mids);
# 2. Fit the model of interest (scientific model) on each imputed data set by the with() function, resulting an object of class mira;
# 3. Pool the estimates from each model into a single set of estimates and standard errors, resulting in an object of class mipo;
# 4. Optionally, compare pooled estimates from different scientific models by the D1() or D3() functions.
# A common error is to reverse steps 2 and 3, i.e., to pool the multiply-imputed data instead of the estimates. Doing so may severely bias the estimates of scientific interest and yield incorrect statistical intervals and p-values. The pool() function will detect this case.

## you always do the analysis - or extract the statistical information - separately for each imputed dataset. Then you pool (aka combine) the results from each analysis together with specific formulas
# https://www.daviddisabato.com/blog/2021/2/13/analyzing-and-pooling-results-from-multiply-imputed-data

#FLAG
# cut bottom percentile overall or by GROUP
By_Group <- FALSE
warning("if the cutting of the data is performed by group, also the imputation should be performed likewise, with difeerent analysis pipelines from early on!")

# Loading iris dataset 
data("iris")
iris_cut <- iris[,c(1,5)]
iris <- iris[,c(1,5)]
PROB.QUANTILE <- 0.5


if (!By_Group) {
# removing the lower quartile to simulate left censored data.

cutoff <- quantile(iris_cut$Sepal.Length, probs = PROB.QUANTILE)
iris_cut$Sepal.Length[iris_cut$Sepal.Length < cutoff] <- NA
}else{
# Loop through each unique species and remove data points below the cutoff
for(sp in unique(iris_cut$Species)) {
  cutoff <- quantile(iris_cut$Sepal.Length[iris_cut$Species == sp], probs = PROB.QUANTILE)
  iris_cut <- iris_cut %>% mutate(Sepal.Length = ifelse(Sepal.Length < cutoff & Species == sp, NA, Sepal.Length))
}
}


# Calculate the percentage of missing values by species
iris_missing_pct <- iris_cut %>%
  group_by(Species) %>%
  summarise(missing_pct = sum(is.na(Sepal.Length))/n())
# Print the percentage of missing values by species
PercMissing <- round(mean(iris_missing_pct$missing_pct)*100,0)

warning("implement QRILC in this demo!")
# for (GENE in unique(genes_data$gene)) {
#   rel_conc_matrix <- as.matrix(genes_data[which(genes_data$gene==GENE), "rel_conc"])
#   # note that tune.sigma = 1 as it is assumed that the complete data distribution is Gaussian (as shown by the complete gene DEF) .
#   imputed_rel_conc <- impute.QRILC(log(rel_conc_matrix),tune.sigma = 1) # log transform data.. see explanation in paper
#   genes_data[which(genes_data$gene == GENE), "rel_conc_imputed"] <- exp(unlist(imputed_rel_conc[1])) # back-transform with exp
#   sum(is.na(genes_data[is.na(genes_data$rel_conc) & genes_data$gene == GENE, "rel_conc_imputed"])) # should be 0
#   
# }



# Creating the imputation model for MNAR data
imp_model <- mice(
  iris_cut,
  m = 3, #number of imputations
  maxit = 10,
  meth = c(Sepal.Length = "pmm",  # predictive mean matching
           Species = "polyreg"),  # polynomial regression imputation
  seed = 1234
)

# Creating imputed datasets accounting for Left Censoring
#mice.impute.leftcenslognorm
imp_model_LeftCens <- mice(
  iris_cut,
  m = 3, #number of imputations
  method = c("leftcenslognorm", ""),
  lod=c(min(iris_cut$Sepal.Length,na.rm=T), NA), 
  debug=FALSE,
  randseed = 1234)


# Checking the imputed dataset for censored values
as.data.frame(imp_model$imp$Sepal.Length)[as.data.frame(imp_model$imp$Sepal.Length) < cutoff]


warning(" use information in https://www.daviddisabato.com/blog/2021/2/13/analyzing-and-pooling-results-from-multiply-imputed-data to perform stats on MICE objects and pool results ")

#pooled_data <- pool(imp_model)

# Extracting complete cases from original dataset
iris_missing <- iris_cut %>% filter(!is.na(Sepal.Length))
iris_missing$Imputation <- paste0(PercMissing,"% NA")
#adding the original
complete_data <- iris
complete_data$Imputation <- "full"

# Extracting imputed datasets
#Mice
imputed_data <- complete(imp_model, "long")
#Mice leftcenslognorm
imputed_leftCens_data <- complete(imp_model_LeftCens, "long")

# Updating the name of the imputation number column
imputed_data$Imputation <- factor(paste0("mice-",imputed_data$.imp))
imputed_leftCens_data$Imputation <- factor(paste0("miceLC-",imputed_leftCens_data$.imp))

# create a named list with all data frames
data_list <- list(
  iris_missing = iris_missing,
  complete_data = complete_data,
  imputed_data = imputed_data,
  imputed_leftCens_data = imputed_leftCens_data
)

# create a function that adds "name" column to the data frame
add_name_column <- function(df, name) {
  df$name <- name
  return(df)
}

# apply the add_name_column function to each element of the list
named_data_list <- lapply(names(data_list), function(name) {
  add_name_column(data_list[[name]], name)
})

# stack all data frames in the list using bind_rows
full_data <- bind_rows(named_data_list)
# Creating the plot
ggplot(full_data, aes(x = Imputation, y = Sepal.Length)) +
  geom_violin(aes(fill = name, alpha = 0.5), position = "identity") +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, position = position_dodge(width = 0.5)) +
  xlab("Imputation Method") +
  ylab("Sepal Length") +
  scale_fill_discrete(name = "name") +
  scale_color_discrete(name = "name") +
  ggtitle("Comparison of Sepal Length by Imputation method") +
  theme_bw() +
  guides(alpha = FALSE, position = "bottom")
