##################################################################################################
################################## RT qPCR IMMUNE RELATED GENES ##################################
##################################################################################################

# load the packages
library(tidyverse)
library(platetools) # to plot plate
library(naniar) # show missing points
library(reshape2) # for melt



##################################################################
################## QUALITY CHECK #################################
##################################################################

WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Adriano_RTqPCR_immune_genes"
DATADIR <-  "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data"
IMMUNITY_DATA <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Pathogen_Quantification_Data"

####define to_keep variables to keep clearing memory between runs
to_keep <- c(ls(),c("to_keep"))


##### LOAD FILES
# check if the master file has already been created
if (!file.exists(paste(WORKDIR, "Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv", sep = "/"))) {
  # open source files
  files <- list.files(path = WORKDIR, pattern = "*Data.txt")
  # here fread() is needed as the files are very messy (missing closing lines)
  data_list <- lapply(files, FUN = function(files) {
    data.table::fread(paste(WORKDIR, files, sep = "/"), header = TRUE, sep = "\t")
  })
  names(data_list) <- gsub(".txt", "", files)

  # make cols uniform
  my_list <- lapply(data_list, function(x) {
    # if(any(grepl("Tm", colnames(x)))) {
    #   parse(grep("Tm", colnames(x), value = TRUE)) <- "Tm_product"
    # } else
    if (any(grepl("V", colnames(x)))) {
      x[, grep("V", colnames(x))] <- NULL
    }
    # if there is no column named "Tm", add it, to make N of cols equal among dfs
    if (!any(grepl("Tm", colnames(x)))) {
      x$Tm <- NA
    }
    # strip off all column names
    colnames(x) <- NULL
    # assign a new name
    colnames(x) <- c("Sample_Well", "Well_Type", "Threshold", "Ct", "Tm_Product")
    return(x)
  })

  # assign column with the name
  for (i in 1:length(my_list)) {
    my_list[[i]]$name <- names(my_list[i])
  }

  # rbind the list elements in a single dataframe
  genes_data <- do.call(rbind, my_list)
  genes_data$Well_Type <- NULL
  
  # combine all the duplicates to have Ct1, Ct2, Tm1, Tm2
  
  # remove all the No_Ct values
  genes_data[which(genes_data$Ct == "No Ct"), "Ct"] <- NA
  genes_data$Ct <- as.numeric(genes_data$Ct)
  
  # extract the gene labels
  genes_data <- genes_data %>%
    # str_match will cut everything before the 3rd dash % sub will cut everythig after the space
    mutate(gene_rep = sub(" .*", "", stringr::str_match(name, "([^-]+)(?:-[^-]+){3}$")[, 1]))
  
  # remove the duplicate plate info
  genes_data$gene_rep <- substr(genes_data$gene_rep, 1, nchar(genes_data$gene_rep) - 2)

  # combine all the duplicates to have Ct1, Ct2, Tm1, Tm2
  summarised_data <- genes_data %>%
    mutate(Ct = as.numeric(Ct)) %>%
    group_by(Sample_Well, gene_rep) %>%
    summarise(mean_Ct = mean(Ct, na.rm = TRUE), mean_Tm = mean(Tm_Product, na.rm = TRUE))
  # add gene label
  summarised_data <- summarised_data %>%
    separate(gene_rep, into = c("Sample_Plate", "gene"), sep = "-")
  
  ### load ORIGINAL PLATE POSITIONS reference
  # open source files
  Oriplate <- list.files(path = WORKDIR, pattern = "OriPlates")
  Oriplate_list <- lapply(Oriplate, FUN = function(files) {
    read.csv(paste(WORKDIR, files, sep = "/"), header = FALSE, sep = ",")
  })
  #assign names to list items by removing start of the string and last 4 characters (.csv)
  names(Oriplate_list) <- gsub("Adriano-RTqPCR-OriPlates_references_PLATE", "", Oriplate)
  names(Oriplate_list) <- substr(names(Oriplate_list), 1, nchar(names(Oriplate_list)) - 4)
  
  # plate size:
  n_samples <- 96
  nrow <- 8
  array <- expand.grid(LETTERS[seq_len(nrow)], seq_len(n_samples / nrow))
  array <- array[order(array$Var1),]

  Oriplate_list <- lapply(Oriplate_list, function(x) {
    # make into a column
    x <- as.vector(as.matrix(x))
    # Generating a 96 well plate layout
    samples <- character(n_samples)
    samples[seq_along(x)] <- x
    samples <- matrix(samples, nrow = nrow)
    # colnames(samples) <- seq_len(n_samples / nrow)
    # rownames(samples) <- LETTERS[seq_len(nrow)]
    samples <- reshape2::melt(samples)
    samples$Sample_Well <- paste0(array$Var1, array$Var2)
    colnames(samples)[which(colnames(samples)=="value")] <- "Code"
    x <- samples[,3:4]
    })
  
  # assign column with the name
  for (i in 1:length(Oriplate_list)) {
    Oriplate_list[[i]]$Sample_Plate <- names(Oriplate_list[i])
  }
  
  # rbind the list elements in a single dataframe
  plate_positions_data <- do.call(rbind, Oriplate_list)
  # #Merge info file of plate positions with the DNA results
  common_col_names <- intersect(names(plate_positions_data), names(summarised_data))
  summarised_data <- merge(summarised_data, plate_positions_data, by=common_col_names, all.x=TRUE)
  #fix missing labels
  summarised_data$Ori_Well <- sub(".*\\-", "", summarised_data$Code)
  summarised_data$Ori_Plate <- sub("\\-.*", "", summarised_data$Code)
  
  
  ### CHECK THE WELLS POSITIONS TO ADD THE ANT IDENTITY
  plate_positions_list <- list.files(path= file.path(IMMUNITY_DATA, "QPCR_SAMPLING_TAGSCANNER"),pattern="PLAQUE",full.names = T)
  
  plate_positions <- do.call("rbind", lapply(plate_positions_list, function(x) {
    dat <- read.csv(x, header=TRUE)
    dat$fileName <- tools::file_path_sans_ext(basename(x))
    dat
  }))
  
  colnames(plate_positions)[which(colnames(plate_positions)=="Comment")] <- "Ori_Well"
  #extract col and plaque info
  plate_positions$Ori_Plate <- sub("\\_.*", "", plate_positions$fileName)
  plate_positions$Ori_Plate <- gsub('[PLAQUE]','',plate_positions$Ori_Plate)
  
  plate_positions$Colony <- sub(".*\\-", "", plate_positions$fileName)
  
  common_col_names2 <- intersect(names(summarised_data), names(plate_positions))
  genes_Results_annotated <- merge(summarised_data, plate_positions, by=common_col_names2, all.x=TRUE)
  
  ### clean up of empty wells (NA)
  genes_Results_annotated <- genes_Results_annotated[!is.na(genes_Results_annotated$Ori_Well),]
  # these ants have not been included in this analysis (used for pre-tests), the wells were empty
  missing_ants <- c("1-E1","2-D6","20-D3","20-C4","55-B7","54-H6","59-B1","59-H7","56-B2","56-E2","27-B11","26-E10","46-H11","47-H1","11-H5","11-C7","33-D1","33-A2")
  genes_Results_annotated <- filter(genes_Results_annotated, !(Code %in% missing_ants))
  ### relabel controls (NoT)
  # needed?
  
  ### ADD THE METADATA
  metadata <- read.table(paste(DATADIR,"/Metadata_Exp1_2021_2022-10-12.txt",sep=""),header=T,stringsAsFactors = F, sep=",")
  #rename cols to match the DNA results
  colnames(metadata)[which(colnames(metadata)=="REP_treat")] <- "Colony"
  # MERGE BASED ON TAG ID AS THE PLATE POSTION ANTidS ARE NOT RELIABLE (misalignment issue spotted for Q of R4BP, which had right TagID=228, but wrong antID (153 instead of 152). )
  colnames(metadata)[which(colnames(metadata)=="tagIDdecimal")] <- "TagID"
  
  #Merge info file of plate positions with the DNA results
  common_col_names3 <- intersect(names(genes_Results_annotated), names(metadata))
  genes_Results_annotated <- merge(genes_Results_annotated, metadata, by=common_col_names3, all.x=TRUE)
  #ASSIGN QUEEN LABEL TO QUEEN ISTEAD OF NURSE (SHOULD BE FIXED IN METADATA!!!)
  genes_Results_annotated[which(genes_Results_annotated$IsQueen==TRUE),"AntTask"] <- "queen"
  
  # #remove extra columns
  genes_Results_annotated <- genes_Results_annotated[ , !(names(genes_Results_annotated) %in% c("X.ScanTime","surviv_time","ExpStart","ExpEnd","Comment","TagID","tagIDdecimal","Ori_Plate","Ori_Well","fileName","identifStart","identifEnd"))]
  # #remove  dead ants
  # DNA_Results_annotated <- DNA_Results_annotated[which(!is.na(DNA_Results_annotated$AntTask)),]
  
  
  # write the dataframe to a csv file
  write.csv(genes_Results_annotated, paste(DATADIR, "Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv", sep = "/"), row.names = FALSE)

  
  to_keep <- c(to_keep,"genes_Results_annotated","controls_gene_data")
  #cleaning
  rm(   list   =  ls()[which(!ls()%in%to_keep)]    )
  gc()
  
}

##### READ STARTING FILE
# read the csv file
genes_data <- read.csv(paste(DATADIR, "Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv", sep = "/"))
genes_data$Sample_Plate <- as.factor(genes_data$Sample_Plate)

# save controls aside
controls_gene_data <- genes_data %>%
  filter(is.na(Colony))
# keep housekeeping and test genes data
genes_data <- genes_data %>%
  filter(!is.na(Colony))

# plot RAW data
ggplot(genes_data, aes(x = gene, y = mean_Ct, color = Sample_Plate)) +
  geom_boxplot(aes(colour = Sample_Plate), lwd = 0.8, alpha = 0.3)


# Using the delta-delta-Ct method
# One common way of analysing qPCR data is to use the “delta-delta-Ct” method. This involves calculating the difference between the Ct of the housekeeping gene and the test gene, then calculating the difference between the treated samples and the control.
# split test genes from housekeeping
test_gene_data <- genes_data %>%
  filter(gene != "EF1")
# isolate housekeeping and change colnames
ref_gene_data <- genes_data %>%
  filter(gene == "EF1") %>%
  rename("housekeeping" = "gene", "ref_Ct" = "mean_Ct", "ref_Tm" = "mean_Tm")

# recombine data
genes_data <- left_join(test_gene_data, ref_gene_data, by = c("Sample_Well", "Sample_Plate"))

# create a new column containing the delta Ct between the housekeeping gene and our gene of interest, and plot the delta Ct for each treatment and replicate.
genes_data <- mutate(genes_data, delta_Ct = mean_Ct - ref_Ct)


# ref_Ct

# plot delta_Ct
ggplot(genes_data, aes(x = gene, y = delta_Ct)) +
  geom_boxplot(aes(colour = Sample_Plate), lwd = 0.8, alpha = 0.3)

# Calculate the mean delta Ct for each treatment.
treatment_summary <- genes_data %>%
  group_by(gene) %>%
  summarise(mean_delta_Ct = mean(delta_Ct, na.rm = TRUE))



# Calculate mean and standard deviation
mean_vec <- mean(ref_gene_data$ref_Ct, na.rm = TRUE)
std_vec <- sd(ref_gene_data$ref_Ct, na.rm = TRUE)

# Select values larger than 1 standard deviation from mean
result <- ref_gene_data$ref_Ct[ref_gene_data$ref_Ct > mean_vec + 2*std_vec | ref_gene_data$ref_Ct < mean_vec - 2*std_vec]

# Print result
print(result)




########## TODOS

# plot stuff










# # https://liz-is.github.io/qpcr-analysis-with-r/aio.html
# # Now we can calculate the delta delta Ct of each replicate compared to the mean of the control sample.
# mean_control <- filter(treatment_summary, RNAi == "Control") %>% pull(mean_delta_Ct)
#
# genes_data <- genes_data %>%
#   mutate(delta_delta_Ct = mean_control - delta_Ct)
#
# ggplot(genes_data, aes(x = RNAi, y = delta_delta_Ct)) +
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
# platetools::z_grid(data = genes_data$mean_Ct,
#        well = genes_data$Well,
#        plate_id = plate_id) +
#   ggtitle("Virus Neutralization Test")





# a = slope of the standard curve
# b = intercept of the standard curve
# genesDNA <- 10^(slope*average(Ct1-Ct2)+intercept))
