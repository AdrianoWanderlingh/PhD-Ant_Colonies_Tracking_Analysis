##################################################################################################
################################## RT qPCR IMMUNE RELATED GENES ##################################
##################################################################################################

# load the packages
library(tidyverse)
library(platetools) # to plot plate
library(naniar) # show missing points


##################################################################
################## QUALITY CHECK #################################
##################################################################

DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Adriano_RTqPCR_immune_genes"


##### LOAD FILES
# check if the master file has already been created
if (!file.exists(paste(DATADIR, "Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv",sep="/"))){
# open source files
files <- list.files(path = DATADIR, pattern = "*.txt")
#here fread() is needed as the files are very messy (missing closing lines)
data_list <- lapply(files, FUN = function(files) {data.table::fread(paste(DATADIR,files,sep="/"), header = TRUE, sep = "\t")})
names(data_list) <- gsub(".txt", "", files)

# make cols uniform
my_list <- lapply(data_list, function(x) {
  # if(any(grepl("Tm", colnames(x)))) {
  #   parse(grep("Tm", colnames(x), value = TRUE)) <- "Tm_product"
  # } else 
    if(any(grepl("V", colnames(x)))){
    x[, grep("V", colnames(x))] <- NULL
    }
  # if there is no column named "Tm", add it, to make N of cols equal among dfs
  if (!any(grepl("Tm", colnames(x)))) {
    x$Tm <- NA
  }
  #strip off all column names
  colnames(x) <- NULL
  # assign a new name
  colnames(x) <- c("Well", "Well_Type", "Threshold","Ct","Tm_Product")
  return(x)
})

#assign column with the name
for (i in 1:length(my_list)) {
  my_list[[i]]$name <- names(my_list[i])
}

# rbind the list elements in a single dataframe
genes_data <- do.call(rbind, my_list)

genes_data$Well_Type <- NULL

# write the dataframe to a csv file
write.csv(genes_data, paste(DATADIR, "Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv",sep="/"), row.names = FALSE)
}

##### READ STARTING FILE
# read the csv file
genes_data <- read.csv(paste(DATADIR, "Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv",sep="/"))


# combine all the duplicates to have Ct1, Ct2, Tm1, Tm2
class(genes_data$Ct)

# remove all the No_Ct values
genes_data[which(genes_data$Ct == "No Ct"),"Ct"] <- NA
genes_data$Ct <- as.numeric(genes_data$Ct)



# text <- "221111-Adriano-RTqPCR-1-DEF-1"
# pieces <- strsplit(text, "-")[[1]]
# result <- tail(pieces, 3)
# result_string <- paste(result, collapse = "-")

# extract the gene labels
genes_data <- genes_data %>%
  #str_match will cut everything before the 3rd dash % sub will cut everythig after the space
  mutate(gene_rep = sub(" .*", "", stringr::str_match(name, '([^-]+)(?:-[^-]+){3}$')[,1]))

#remove the duplicate plate info
genes_data$gene_rep <- substr(genes_data$gene_rep,1,nchar(genes_data$gene_rep)-2)

head(genes_data)


# combine all the duplicates to have Ct1, Ct2, Tm1, Tm2
summarised_data <- genes_data %>%
  mutate(Ct = as.numeric(Ct)) %>%
  group_by(Well, gene_rep) %>%
  summarise(mean_Ct = mean(Ct,na.rm = TRUE),mean_Tm = mean(Tm_Product,na.rm = TRUE))
# add gene label
summarised_data <- summarised_data %>%
  separate(gene_rep, into = c("Plate", "gene"), sep = "-")

# Using the delta-delta-Ct method
# One common way of analysing qPCR data is to use the “delta-delta-Ct” method. This involves calculating the difference between the Ct of the housekeeping gene and the test gene, then calculating the difference between the treated samples and the control.
# split test genes from housekeeping
test_gene_data <- summarised_data %>%
  filter(gene != "EF1")
# isolate housekeeping and change colnames
ref_gene_data <- summarised_data %>%
  filter(gene == "EF1") %>%
  rename("housekeeping" = "gene","ref_Ct" = "mean_Ct", "ref_Tm" = "mean_Tm")

#recombine data
summarised_data <- left_join(test_gene_data, ref_gene_data, by = c("Well", "Plate"))

#create a new column containing the delta Ct between the housekeeping gene and our gene of interest, and plot the delta Ct for each treatment and replicate.
summarised_data <- mutate(summarised_data, delta_Ct = ref_Ct - mean_Ct)


# plot RAW data
ggplot(summarised_data, aes(x = gene, y = mean_Ct,color = Plate)) +
  geom_boxplot(aes(colour=Plate),lwd=0.8,alpha = 0.3)

# plot delta_Ct
ggplot(summarised_data, aes(x = gene, y = delta_Ct)) +
  geom_boxplot(aes(colour=Plate),lwd=0.8,alpha = 0.3)

# Calculate the mean delta Ct for each treatment.
treatment_summary <- summarised_data %>%
  group_by(gene) %>%
  summarise(mean_delta_Ct = mean(delta_Ct,na.rm = TRUE))




########## TODOS

# calculate the slopes and intercepts per rep (which are the controls?)
# plot stuff










# # https://liz-is.github.io/qpcr-analysis-with-r/aio.html
# # Now we can calculate the delta delta Ct of each replicate compared to the mean of the control sample.
# mean_control <- filter(treatment_summary, RNAi == "Control") %>% pull(mean_delta_Ct)
# 
# summarised_data <- summarised_data %>% 
#   mutate(delta_delta_Ct = mean_control - delta_Ct)
# 
# ggplot(summarised_data, aes(x = RNAi, y = delta_delta_Ct)) +
#   geom_point()
#
# # Calculating relative DNA concentration
# # If you want to calculate the relative DNA concentration, you can use the fact that the amount of cDNA theoretically doubles every cycle.
# combined_data <- combined_data %>%
#   mutate(rel_conc = 2^-delta_delta_Ct)
# # We can now plot the relative concentration.
# ggplot(combined_data, aes(x = RNAi, y = rel_conc)) +
#   geom_point() +
#   scale_y_continuous(labels = scales::percent, limits = c(0, 1.3))


# # should be split by gene and rep 
# plate_id <- rep(c("My Plate"), each = 96)
# 
# platetools::z_grid(data = summarised_data$mean_Ct,
#        well = summarised_data$Well,
#        plate_id = plate_id) +
#   ggtitle("Virus Neutralization Test")





# a = slope of the standard curve
# b = intercept of the standard curve
# genesDNA <- 10^(slope*average(Ct1-Ct2)+intercept))


