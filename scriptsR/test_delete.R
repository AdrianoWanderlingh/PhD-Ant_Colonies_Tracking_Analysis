# Create dummy data
df1<-data.frame(start = c(1,1,1), altro = c(5,6,7))
df1
df3<-data.frame(a = c(10,10,10), time= c(1,2,3))
df4<-data.frame(a = c(12,12,12), time= c(1,2,3))

# Create a list
l1 <- list(data.frame(start = c(1,1,1), altro = c(5,6,7)), df3, df4)
names(l1) <- c("df1","df3", "df4")

lx <- list(l,l1)

# add new column 'b'
# create 'b' values based on column 'a' 
l2<-lapply(lx[which(lx[[2]]=="qualcosa?"),"a"], function(x) 
  cbind(x, b = x$a*4))
#Results in:
  
l2
lx[[2]]

#In your case something like:
  
  filelist<-lapply(filelist, function(x) 
    cbind(x, b = x$SampleID))
  
  
  
  
  
  data2$Cost <-data2$Cost + data1$Cost[match(data2$X, data1$X)]
  
  
  
  
  data2[data2$X %in% data1$X,]$Cost <- data2[data2$X %in% data1$X,]$Cost + data1[data1$X %in% data2$X,]$Cost 
  
  
  
  
  setwd("C:/Users/cf19810/Google Drive/Vivi&Adri_Shared/UNI_BRISTOL/EXPERIMENTS/Adriano_raw_data_HOMERANGE_ANALYSIS")
  
  files.list = lapply(list.files(pattern=".csv"), read.delim)
  
  FUNX <- function(data_raw){
    width_pixel <- dist(data_raw[which(data_raw$Type==4),c("X","Y")])
    data_raw[c("X","Y")] <- data_raw[c("X","Y")] * width_mm / width_pixel
    
    ##determine groups ("Type") included in the computation of the occupied area
    test_data <- data_raw[which(data_raw$Type ==1 | data_raw$Type==3),c("Type","X","Y")]
  

    return(results$home_range_95/ sum(complete.cases(test_data)))
    
    #write.table(get(paste("result_",area,sep="")),file=paste(result_path,"/",area,"/",replicate,"_",colony,"_",treatment,"_",group,"_homerange.dat",sep=""),quote=FALSE,sep="\t",col.names=T,row.names=FALSE,append=T)
    #test_data <- data_raw[which(data_raw$Type ==1 | data_raw$Type==3),c("Type","X","Y")] #output including workers and queen for computation of occupied volume
    #return(test_data)
  }
  
ant_occupancy_mm2 <- lapply(files.list,FUN=FUNX)
