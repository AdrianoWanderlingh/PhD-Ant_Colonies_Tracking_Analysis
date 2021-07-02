rm(list=ls())

library(adehabitatHR)
library(sp)

width_mm <- 70
# Type column:
# - "1": workers
# - "2": brood
# - "3": queen
# - "4": nest width (70 mm)

setwd("C:/Users/cf19810/Google Drive/Vivi&Adri_Shared/UNI_BRISTOL/EXPERIMENTS/Adriano_raw_data_HOMERANGE_ANALYSIS")

files.list = lapply(list.files(pattern=".csv"), read.delim)

occupied_area_by_ant <- function(data_raw){
  width_pixel <- dist(data_raw[which(data_raw$Type==4),c("X","Y")])
  data_raw[c("X","Y")] <- data_raw[c("X","Y")] * width_mm / width_pixel
  
  ##determine groups ("Type") included in the computation of the occupied area
  test_data <- data_raw[which(data_raw$Type ==1 | data_raw$Type==3),c("Type","X","Y")]
 
  
#test_data <- normalised_file_list[[2]]

#####################make empty coordinate list; will be used for the grid
span <- expand.grid(X=-20:120,Y=-20:120,obs=0) # to be adjusted to be big enough to allow the estimation of home-range. Chosen according to error messages and final graph

###define grid (to do only the first time)
if (!exists("customised_grid")){
  for_grid <- span[,c("X","Y")]###list of all x-y coordinates of the plot (each existing cell (even unvisited, empty cells) defined on 1 line by its x,y coordinates)
  names(for_grid) <- c("x","y")
  customised_grid <- (unique(round(for_grid)))
  coordinates(customised_grid)=c("x","y")
  gridded(customised_grid) <- TRUE
}

###compute home ranges - kernelUD is used to estimate the utilization distribution (UD) of animals
locations_sp <- NULL
locations_sp <- SpatialPoints(test_data[c("X", "Y")])###get the data in the right format for following calculations
#print("locations_sp defined")
if(exists("bandwidth")){rm(list=("bandwidth"))}  
try(bandwidth <- kernelUD(locations_sp, h="href",grid=get("customised_grid"),extent=0)@h$h,silent=T)###get bandwidth (without bounder - because there is a bug when using bounder and href)
#print("bandwidth defined")
if(exists("kud")){rm(list=("kud"))}
try(kud <- kernelUD(locations_sp, h=bandwidth,grid=get("customised_grid"),extent=0),silent=T)###calculate utilization distribution with border definition
#print("kud computed")
###########get home ranges
if(exists("home_range")){rm(list=("home_range"))}
try(home_range <- kernel.area(kud, percent = c(50,75,95),            unin = "m",            unout = "m2", standardize = TRUE),silent=T)######THIS IS THE AREA COVERED IN SQUARED MILLIMITERS for 20 to 95% ranges
print("home range computed")
if(exists("home_range")){
  #centroid <- gCentroid(locations_sp)
  centroid <- data.frame(x=mean(locations_sp$X),y=mean(locations_sp$Y))
  results <- data.frame(home_range_50=home_range["50"], home_range_75=home_range["75"], home_range_95=home_range["95"],centroid_x=centroid$x,centroid_y=centroid$y) #original structure that sums up replicates -> data.frame(replicate=replicate,colony=colony...)
  assign("centroid",centroid)
}

return(results$home_range_95/ sum(complete.cases(test_data)))
#return(results/ sum(complete.cases(test_data))) #to use this, something like this x1 = lapply(x, function(l) l[[1]]) may be needed outside of the function
#obtain different objects in lists, then get overall mean by object, unlisting object

#write.table(get(paste("result_",area,sep="")),file=paste(result_path,"/",area,"/",replicate,"_",colony,"_",treatment,"_",group,"_homerange.dat",sep=""),quote=FALSE,sep="\t",col.names=T,row.names=FALSE,append=T)
#test_data <- data_raw[which(data_raw$Type ==1 | data_raw$Type==3),c("Type","X","Y")] #output including workers and queen for computation of occupied volume
#return(test_data)
}

#How to return multiple objects as home-ranges and occupancy? should be joined in single file?

ant_occupancy_mm2 <- lapply(files.list,FUN=occupied_area_by_ant)

mean_occupied_area_ant95 <- mean(unlist(ant_occupancy_mm2))


small <- 31 #N workers + queen
big <- 181 #N workers + queen
area_small <- mean_occupied_area_ant95 * small
area_big <- mean_occupied_area_ant95 * big

# #plot area occupied by ant
# barplot(unlist(ant_occupancy_mm2))
# abline(h=mean_occupied_area_ant95)
# print(mean_occupied_area_ant95)

#plot nest sizes
require(grid)
plot(c(-10, 60), c(-10,60), type = "n",asp=1, main= paste0("Ant nest size depending on N of ants \n ants density with HR 75% =", round(mean_occupied_area_ant,digits = 2), "mm^2/ant"),
     sub = paste0("Nest area per: ",small-1," workers + queen = ",round(area_small),"mm^2;",
                  "\n",big-1," workers + queen = ",round(area_big),"mm^2;",
                  "\n Science 2018 (fixed size) = ",70*40,"mm^2."),
     xlab="",
     ylab="mm",)
rect( 0, 40, 70, 0, border="gray") #Science 2018 nest size
text(62, 10, labels = "70mm \n x \n 40mm",col = "gray")
rect( 0, sqrt(area_small), sqrt(area_small), 0) #small nest size
text(sqrt(area_small)/2, sqrt(area_small)/2, 
     labels = paste0(round(sqrt(area_small),digits = 2),"mm \n x \n",round(sqrt(area_small),digits = 2),"mm"))
rect( 0, sqrt(area_big), sqrt(area_big), 0) #big nest size
text(sqrt(area_big)/1.5, sqrt(area_big)/1.5, 
     labels = paste0(round(sqrt(area_big),digits = 2),"mm \n x \n",round(sqrt(area_big),digits = 2),"mm"))

warnings()

#squares size comparison, changing n of ants
small25 <- 21 #N workers + queen
big150 <- 121 #N workers + queen
area_small25 <- mean_occupied_area_ant95 * small25
area_big150 <- mean_occupied_area_ant95 * big150


plot(c(-10, 60), c(-10,60), type = "n",asp=1, main= paste0("Ant nest size depending on N of ants \n ants density with HR 95% =", round(mean_occupied_area_ant,digits = 2), "mm^2/ant")
     ,sub = "Size comparisons: \n 180 vs 150 \n 30 vs 25",
     xlab="",
     ylab="mm",)
rect( 0, 40, 70, 0, border="gray") #Science 2018 nest size
text(62, 10, labels = "70mm \n x \n 40mm",col = "gray")

rect( 0, sqrt(area_big), sqrt(area_big), 0, col = "green") #big 180 nest size
rect( 0, sqrt(area_big150), sqrt(area_big150), 0, col = "lightgreen") #big 150 nest size

rect( 0, sqrt(area_small), sqrt(area_small), 0, col = "blue") #small 30 nest size
rect( 0, sqrt(area_small25), sqrt(area_small25), 0, col = "lightblue") #small 25 nest size

#squares size comparison, changing HR size
mean_occupied_area_ant95 <- 19.73358
mean_occupied_area_ant75 <- 10.65043

small <- 31 #N workers + queen
big <- 181 #N workers + queen
area_small95 <- mean_occupied_area_ant95 * small
area_small75 <- mean_occupied_area_ant75 * small

area_big95 <- mean_occupied_area_ant95 * big
area_big75 <- mean_occupied_area_ant75 * big



plot(c(-10, 60), c(-10,60), type = "n",asp=1, main= paste0("Ant nest size depending on N of ants \n Ants density with HR 75% = ", round(mean_occupied_area_ant75,digits = 2), "mm^2/ant (light colour)", "\n and HR 95% = ", round(mean_occupied_area_ant95,digits = 2), "mm^2/ant (full colour)")
     ,sub = paste0("Nest area per: ",small-1," workers + queen = ",round(area_small75),"mm^2;",
                   "\n",big-1," workers + queen = ",round(area_big75),"mm^2;",
                   "\n Science 2018 (fixed size) = ",70*40,"mm^2."),
     xlab="",
     ylab="mm",)
rect( 0, 40, 70, 0, border="gray") #Science 2018 nest size


rect( 0, sqrt(area_big95), sqrt(area_big95), 0, col = "green") #big 180 nest size
rect( 0, sqrt(area_big75), sqrt(area_big75), 0, col = "lightgreen") #big 150 nest size
text(sqrt(area_big75)/1.4, sqrt(area_big75)/1.4, 
     labels = paste0(round(sqrt(area_big75),digits = 2),"mm \n x \n",round(sqrt(area_big75),digits = 2),"mm"))

rect( 0, sqrt(area_small95), sqrt(area_small95), 0, col = "blue") #small 30 nest size
rect( 0, sqrt(area_small75), sqrt(area_small75), 0, col = "lightblue") #small 25 nest size
text(sqrt(area_small75)/2, sqrt(area_small75)/2, 
     labels = paste0(round(sqrt(area_small75),digits = 2),"mm \n x \n",round(sqrt(area_small75),digits = 2),"mm"))


text(65, 10, labels = "70mm \n x \n 40mm",col = "gray")


# #extra plotting for kernelUD and VolumeUD
# ####some plotting
# vud <- getvolumeUD(kud)
# ## Set up graphical parameters
# par(mfrow=c(2,1))
# par(mar=c(0,0,2,0))
# ## The output of kernelUD for the first animal (colony?)
# image(kud)
# points(findmax(kud))
# title("Output of kernelUD")
# ## Convert into a suitable data structure for
# ## the use of contour
# xyz <- as.image.SpatialGridDataFrame(kud)
# contour(xyz, add=TRUE)
# ## and similarly for the output of getvolumeUD
# par(mar=c(0,0,2,0))
# image(vud)
# title("Output of getvolumeUD")
# xyzv <- as.image.SpatialGridDataFrame(vud)
# contour(xyzv, add=TRUE)
# points(centroid$x,centroid$y,col="cyan",pch=16)
# ## store the value of the volume under UD in a vector hr95
# hr95 <- as.data.frame(vud)[,1]
# ## if hr95 is <= 95 then the pixel belongs to the home range
# ## (takes the value 1, 0 otherwise)
# hr95 <- as.numeric(hr95 <= 95)
# ## Converts into a data frame
# hr95 <- data.frame(hr95)
# ## Converts to a SpatialPixelsDataFrame
# coordinates(hr95) <- coordinates(vud)
# gridded(hr95) <- TRUE
# ## display the results
# image(hr95)
# points(centroid$x,centroid$y,col="cyan",pch=16)