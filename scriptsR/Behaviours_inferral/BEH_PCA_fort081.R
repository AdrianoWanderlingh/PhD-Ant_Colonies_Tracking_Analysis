print(paste("PERFORM PCA ",REPLICATE, PERIOD))

### PCAs ##########

#USE summary_MANUAL AS A BASE TO START WITH.




#descriptive analysis
#install.packages("DataExplorer")
#install.packages("scales") #may be needed if dataexplorer reports installation issues. Another fix is by checking DataExplorer in the packages tab of Rstudio
library(DataExplorer)
library(corrplot)


str(results)

plot_histogram(results)
plot_density(results)


results.cor = cor(results[,c(5:17)]) #18 for 131219 test

help("cor.mtest")

res1 <- cor.mtest(results.cor, conf.level = .95, method="pearson") #kendall
corrplot(results.cor, type="upper", order="hclust", tl.cex= .8, 
         p.mat = res1$p, sig.level = 0.05, insig = "blank", method = "number", number.cex = .7)
?corrplot

# uncorrelated vars used for 131219:  median_acceleration_mmpersec2, rmse_mm, distance_mm, 
# mean_acceleration_mmpersec2, straightness_index , periodicity_sec

#install.packages("ggpubr")
library("ggpubr")


#results.sub <- results[,c(1:5,8,9,12,16,17)]
results.sub <- results[,c(1:5,8:10,13,16,17)]
str(results.sub)



#--------------------------------


