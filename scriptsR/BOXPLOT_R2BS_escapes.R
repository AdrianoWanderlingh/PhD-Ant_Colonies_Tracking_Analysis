setwd("/home/cf19810/Documents/")
getwd()

library("ggplot2")

data <-read.csv("BOXPLOT_R2BS_escapes.csv",sep = ",")
summary(datav)
head(datav)
data$escapee <- as.factor(data$escapee)
data$nest <- 'R2BS'

ggplot(data = data , aes(x=nest, y=minutes)) + geom_boxplot() + geom_jitter() + ggtitle("time spent outside by escapees. \n tot ants = 36 of 169")
