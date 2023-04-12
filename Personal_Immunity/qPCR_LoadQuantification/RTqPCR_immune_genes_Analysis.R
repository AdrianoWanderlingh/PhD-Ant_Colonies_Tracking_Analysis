##################################################################################################
################################## RT qPCR IMMUNE RELATED GENES ##################################
##################################################################################################

##### CLEAN UP
gc()
mallinfo::malloc.trim(0L)

# load the packages
library(tidyverse)
library(platetools) # to plot plate
library(naniar) # show missing points
library(reshape2) # for melt
library(scales) # for colours
library(ggpubr) # for creating easily publication ready plots
library(ggsci) # for extra colour palettes
library(rstatix) #provides pipe-friendly R functions for easy statistical analyses.
library(gridExtra) #grid of plots
library(cowplot) # 
library(effects)

## for qPCR missing values imputation
# if (!require("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
#   }n
# BiocManager::install("nondetects") # deal with qPCR non-detections (No_Ct)
# BiocManager::install("qpcrNorm")
# BiocManager::install("imputeLCMD")
# library(qpcrNorm)
# library(nondetects)
# library(HTqPCR)
library(imputeLCMD) # qPCR imputation

## for stats
library(performance)
library(remotes)
# install_version("MuMIn", "1.46.0")
library(MuMIn) # multi-model inference
library(multcomp)
library(multcompView)
library(blmeco)
library(sjPlot)
library(lsmeans)
library(lme4)
library(blmeco) # check dispersion for glmer
library(emmeans) # post-hoc comparisons
library(e1071) # calc skewness and other stuff
library(fitdistrplus)
library(lawstat) # for levene test (homogeneity of variance)
library(lmPerm) # for permutations
library(lmerTest)
library(car)
library(nlme) # lme
library(afex)
library(pracma) # calculate shoulders in a curve
#library(gamm4)
library(gamlss) #Generalized Additive Models for Location Scale and Shape


###source function scripts
print("Loading functions and libraries...")
source(paste("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/FUNCTIONS_Analysis_and_Styling.R",sep="/"))
warning("FIX THE PLOTTING COLOURS")

#### MAXIMUM ACCEPTED DIFF IN CT between duplicates WITH 15% PIPETTING ERROR# #Efficiency is between 1.95 and 2.05, differing CT values.
# FUNCTION MODIFIED FROM  # https://rnajournal.cshlp.org/content/23/5/811.full.pdf, BASED ON THE ABOVE
assign_CTDiff_15PipErr <- function(mean_Ct) {
  if (is.na(mean_Ct)) {
    return(0)
  } else if (mean_Ct > 34.5) {return(1.9)
  } else if (mean_Ct > 33.5) {return(1.3)
  } else if (mean_Ct > 32.5) {return(0.9)
  } else if (mean_Ct > 31.5) {return(0.7)  
  } else if (mean_Ct >    0) {return(0.5)
  } else                     {return(0)
  }
}

### fixed PARAMETERS
#Technical_error
#Technical_error <- 3

### TURN plots on/off
EXPLORE_PLOT <- FALSE # exploratory plots
PLOT          <- FALSE # stats plots


# Initialize a data frame to store results
PipelineTesting <- data.frame(gene = character(), P.Status = numeric(), P.Treatment = numeric(),Estim.Status = numeric(),
                              Estim.Treatment = numeric(), T.E. = numeric(), LOD = numeric(), imputation  = character(), Invalids  = character())

WORKDIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Adriano_RTqPCR_immune_genes"
DATADIR <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data"
IMMUNITY_DATA <- "/home/cf19810/Documents/scriptsR/EXP1_base_analysis/Personal_Immunity/Pathogen_Quantification_Data"

#### define to_keep variables to keep clearing memory between runs
to_keep <- c(ls(), c("to_keep"))

##################################################################
#################### CREATE FILE #################################
##################################################################

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

  #assign NA to no_Ct values, then (in a later stage) assign very low value to relative expression
  genes_data[which(genes_data$Ct == "No Ct"), "Ct"] <- NA
  genes_data$Ct <- as.numeric(genes_data$Ct)

  # extract the gene labels
  genes_data <- genes_data %>%
    # str_match will cut everything before the 3rd dash % sub will cut everythig after the space
    mutate(gene_rep = sub(" .*", "", stringr::str_match(name, "([^-]+)(?:-[^-]+){3}$")[, 1]),
           machine_num = sub(" .*", "", stringr::str_match(name, "([^-]+)(?:-[^-]+){1}$")[, 1])
           )

  # remove the duplicate plate info
  genes_data$gene_rep <- substr(genes_data$gene_rep, 1, nchar(genes_data$gene_rep) - 2)

  genes_data$Ct <- as.numeric(genes_data$Ct)
  # add gene label
  genes_data$gene_rep1 <- genes_data$gene_rep
  genes_data <- genes_data %>%
    separate(gene_rep1, into = c("Sample_Plate", "gene"), sep = "-")
   
  ### load ORIGINAL PLATE POSITIONS reference
  # open source files
  Oriplate <- list.files(path = WORKDIR, pattern = "OriPlates")
  Oriplate_list <- lapply(Oriplate, FUN = function(files) {
    read.csv(paste(WORKDIR, files, sep = "/"), header = FALSE, sep = ",")
  })
  # assign names to list items by removing start of the string and last 4 characters (.csv)
  names(Oriplate_list) <- gsub("Adriano-RTqPCR-OriPlates_references_PLATE", "", Oriplate)
  names(Oriplate_list) <- substr(names(Oriplate_list), 1, nchar(names(Oriplate_list)) - 4)

  # plate size:
  n_samples <- 96
  nrow <- 8
  array <- expand.grid(LETTERS[seq_len(nrow)], seq_len(n_samples / nrow))
  array <- array[order(array$Var1), ]

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
    colnames(samples)[which(colnames(samples) == "value")] <- "Code"
    x <- samples[, 3:4]
  })

  # assign column with the name
  for (i in 1:length(Oriplate_list)) {
    Oriplate_list[[i]]$Sample_Plate <- names(Oriplate_list[i])
  }

  # rbind the list elements in a single dataframe
  plate_positions_data <- do.call(rbind, Oriplate_list)
  # #Merge info file of plate positions with the DNA results
  common_col_names <- intersect(names(plate_positions_data), names(genes_data))
  summarised_data <- merge(genes_data, plate_positions_data, by = common_col_names, all.x = TRUE)
  # remove wells that never contained a sample (no code assigned in the battleplan references_PLATE)
  summarised_data <- summarised_data[!is.na(summarised_data$Code),]
  # fix missing labels
  summarised_data$Ori_Well <- sub(".*\\-", "", summarised_data$Code)
  summarised_data$Ori_Plate <- sub("\\-.*", "", summarised_data$Code)

  # ADD EXTRA QUEENS
  ExtraQueens <- read.csv(paste(WORKDIR, "230214_ExtraShamQueens_AW.csv", sep = "/"), header = T, stringsAsFactors = F, sep = ",")
  ExtraQueens$Ct <- as.numeric(ExtraQueens$Ct)
  # uniform names to main dataset
  ExtraQueens$Threshold <- NA
  names(ExtraQueens)[which(names(ExtraQueens)=="Tm")] <- "Tm_Product"
  ExtraQueens$machine_num <- NA
  ExtraQueens$gene_rep <- paste(ExtraQueens$Code,ExtraQueens$gene, sep="-")
  ExtraQueens$name <- paste(ExtraQueens$gene_rep,ExtraQueens$Sample_Plate, sep="-")
  
#   #combine all the duplicates to have meanCts
#   ExtraQueens <- ExtraQueens %>%
#     group_by(Code, gene) %>%
#     summarise(mean_Ct = ifelse(all(is.na(Ct)), NA, mean(Ct, na.rm = TRUE)), #without this, if both NA, NaNs will be created
#               mean_Tm = ifelse(all(is.na(Tm)), NA, mean(Tm, na.rm = TRUE)),
#               Sample_Plate = Sample_Plate,
#               Code = Code,
#               #Sample_Well= NA, # there are 2 distinct sample wells, unlike the rest of the experiment, so they can't be included in the summarise
#               Ori_Well = Ori_Well,
#               Ori_Plate = Ori_Plate)
# #eliminate duplicate rows!
#   ExtraQueens <- ExtraQueens %>% distinct()
  
  #add extra Q
  summarised_data <- rbind(as.data.frame(summarised_data), as.data.frame(ExtraQueens))
  
  ### CHECK THE WELLS POSITIONS TO ADD THE ANT IDENTITY
  plate_positions_list <- list.files(path = file.path(IMMUNITY_DATA, "QPCR_SAMPLING_TAGSCANNER"), pattern = "PLAQUE", full.names = T)

  plate_positions <- do.call("rbind", lapply(plate_positions_list, function(x) {
    dat <- read.csv(x, header = TRUE)
    dat$fileName <- tools::file_path_sans_ext(basename(x))
    dat
  }))

  # REMOVE DUPLICATES!
  plate_positions <- plate_positions[!grepl("DUPLI", plate_positions$Comment), ]
  # REMOVE AntIDs as they can be misleading
  # (misalignment issue spotted for Q of R4BP, which had right TagID=228, but wrong antID (153 instead of 152). )
  # ALSO: some antIDs are missing from the data (i.e. when the myrmidon file was not loaded), while tagIDs are always present
  plate_positions$AntID <- NULL
  
  colnames(plate_positions)[which(colnames(plate_positions) == "Comment")] <- "Ori_Well"
  # extract col and plaque info
  plate_positions$Ori_Plate <- sub("\\_.*", "", plate_positions$fileName)
  plate_positions$Ori_Plate <- gsub("[PLAQUE]", "", plate_positions$Ori_Plate)

  plate_positions$Colony <- sub(".*\\-", "", plate_positions$fileName)

  common_col_names2 <- intersect(names(summarised_data), names(plate_positions))
  genes_Results_annotated <- merge(summarised_data, plate_positions, by = common_col_names2, all.x = TRUE)

  ### clean up of empty wells (NA)
  #genes_Results_annotated <- genes_Results_annotated[!is.na(genes_Results_annotated$Ori_Well), ]
  # these ants have not been included in this analysis (used for pre-tests), the wells were empty
  missing_ants <- c("1-E1", "2-D6", "20-D3", "20-C4", "55-B7", "54-H6", "59-B1", "59-H7", "56-B2", "56-E2", "27-B11", "26-E10", "46-H11", "47-H1", "11-H5", "11-C7", "33-D1", "33-A2")
  genes_Results_annotated <- filter(genes_Results_annotated, !(Code %in% missing_ants))
  
  # general fixes 
  # Plate RT9, well 8H filled with less milliQwater than needed, likely to be discarded for all genes
  genes_Results_annotated <- genes_Results_annotated[which( ! genes_Results_annotated$Code %in% c("27-H8","27-H1","27-H12")),]
  
  ### relabel controls (NoT)
  # needed?

  ### ADD THE METADATA
  metadata <- read.table(paste(DATADIR, "/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")
  # rename cols to match the DNA results
  colnames(metadata)[which(colnames(metadata) == "REP_treat")] <- "Colony"
  # MERGE BASED ON TAG ID AS THE PLATE POSTION ANTidS ARE NOT RELIABLE (misalignment issue spotted for Q of R4BP, which had right TagID=228, but wrong antID (153 instead of 152). )
  # ALSO: some antIDs are missing from the data (i.e. when the myrmidon file was not loaded), while tagIDs are always present
  colnames(metadata)[which(colnames(metadata) == "tagIDdecimal")] <- "TagID"

  # Merge info file of plate positions with the DNA results
  common_col_names3 <- intersect(names(genes_Results_annotated), names(metadata))
  genes_Results_annotated <- merge(genes_Results_annotated, metadata, by = common_col_names3, all.x = TRUE)

  # #remove extra columns
  genes_Results_annotated <- genes_Results_annotated[, !(names(genes_Results_annotated) %in% c("X.ScanTime", "surviv_time", "ExpStart", "ExpEnd", "Comment", "TagID", "tagIDdecimal", "Ori_Plate", "Ori_Well", "fileName", "identifStart", "identifEnd"))]
  # #remove  dead ants
  # genes_data <- genes_data[which(!is.na(genes_data$AntTask)),]

  names(genes_Results_annotated)[which(names(genes_Results_annotated) == "size_treat")] <- "Treatment"
  # Rename by name
  genes_Results_annotated$Treatment <- as.factor(genes_Results_annotated$Treatment)
  levels(genes_Results_annotated$Treatment)[levels(genes_Results_annotated$Treatment) == "BS"] <- "Big Sham"
  levels(genes_Results_annotated$Treatment)[levels(genes_Results_annotated$Treatment) == "SS"] <- "Small Sham"

  # write the dataframe to a csv file
  write.csv(genes_Results_annotated, paste(WORKDIR, "Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv", sep = "/"), row.names = FALSE)

  #to_keep <- c(to_keep, "genes_Results_annotated")
  # cleaning
  rm(list = ls()[which(!ls() %in% to_keep)])
  gc()
}

################################################################################
################################################################################

for (Limit_of_Detection in c(35,36,37)) {

for (Technical_error in c(2,3)) {

for (IMPUTATION in c("QRILC","HM")) {

  for (ALWAYS_DISCARD in c(T)){ #,F ## we want to discard invalid datapoints
    # Limit_of_Detection <- 37
    # Technical_error <- 3
    # IMPUTATION <- "HM"
    
    ##### READ STARTING FILE
    # read the csv file
    genes_data <- read.csv(paste(WORKDIR, "Adriano_RTqPCR_immune_genes_MASTER_REPORT.csv", sep = "/"))
    genes_data$Sample_Plate <- as.factor(genes_data$Sample_Plate)
    #genes_data$mean_Tm <- NULL # problematic as present only in some reps and messes up some functions
    genes_data$gene <- as.factor(genes_data$gene)
    
    ## remove rows which contain NAs for Ct
    #genes_data <- genes_data[!is.na(genes_data$mean_Ct), ]
    
    # Relevel Exposed
    genes_data$Exposed <- as.factor(genes_data$Exposed)
    levels(genes_data$Exposed)[levels(genes_data$Exposed) == "TRUE"] <- "treated"
    levels(genes_data$Exposed)[levels(genes_data$Exposed) == "FALSE"] <- "untreated"
    # Create new Status category
    genes_data$Ant_status <- paste(genes_data$Exposed, genes_data$AntTask)
    genes_data[which(genes_data$Ant_status=="untreated queen" ),"Ant_status"] <- "queen"
    # #REMOVE  MEAN_tm (caused issues with the extra Q)
    # genes_data$mean_Tm <- NULL
    
    # save controls aside
    controls_gene_data <- genes_data %>%
      filter(is.na(Colony))
    # keep housekeeping and test genes data
    genes_data <- genes_data %>%
      filter(!is.na(Colony))
    
    # remove non-categorised ants (DEAD)
    genes_data <- genes_data[which(genes_data$IsAlive == TRUE), ]
    
    # INSTEAD OF SEPARATING THE GROUPS, ADD A GROUP LABEL BASED OFF ANT_STATUS (3 GROUPS) TO RUN THE ANALYSES AND DO THE PLOTS!!
    genes_data$GROUP <- NA
    genes_data[which(genes_data$Ant_status=="queen"),"GROUP"] <- "QUEEN"
    genes_data[which(genes_data$Ant_status=="treated nurse"),"GROUP"] <- "TREATED_W"
    genes_data[which(genes_data$Ant_status %in% c("untreated forager","untreated nurse")),"GROUP"] <- "UNTREATED_W"
    
    
    ###########################################################
    #######       DETERMINE DETECTION THRESHOLD        ########
    ###########################################################
    
    ###### determine detection threshold  (Limit of Detection) of the qPCR machine.
    ## threshold calculated as the shoulder of the values curve (scatterplot of the ordered Cts), which are the peaks of the second derivative.
    
    if (PLOT) {
      #### Find shoulder
      sorted_vector <- sort(genes_data$Ct)
      # Fit a spline curve with a high spar smoothing parameter
      smoothed_curve <- smooth.spline(sorted_vector, spar = 1)
      # Find the second derivative of the smoothed curve
      smoothed_diff <- diff(smoothed_curve$y, differences = 1)
      smoothed_x_diff <- diff(smoothed_curve$x, differences = 1)
      second_derivative <- diff(smoothed_diff) / (head(smoothed_x_diff, -1) * tail(smoothed_x_diff, -1))
      # Compute the absolute value of the second derivative
      abs_second_derivative <- abs(diff(second_derivative))
      sorted_vector <- sorted_vector[-1]
      # Smooth the absolute value curve
      original_second_derivative <- diff(diff(sorted_vector))
      abs_second_derivative <- abs(original_second_derivative)
      # Find the peaks of the smoothed absolute value curve
      peaks <- as.data.frame(pracma::findpeaks(abs_second_derivative))
      peaks <- peaks[which(peaks$V1>=0.2),"V2"]
      # Extract the corresponding x-coordinates from the smoothed curve
      shoulder_points <- data.frame(x = smoothed_curve$x[peaks], y = smoothed_curve$y[peaks])
      # Machine Limit_of_Detection
      Limit_of_Detection <- round(max(shoulder_points$y),0) #conservative threshold, rounds down to 37
      # Plot the curve
      plot(sorted_vector, col = "black", type = "l", main = "Smoothed Curve of raw Ct values",sub=paste("upper Detection Threshold",Limit_of_Detection,sep=" " ),lwd=3)
      points(shoulder_points$x,shoulder_points$y , col = "red", pch = 19)
    }
    ########################################
    
    # # Create a 2 by 2 grid of plots
    # par(mfrow=c(2,2))
    # # Loop through each unique gene and perform the analysis on the subset of data
    # for (GENE in unique(genes_data$gene)) {
    #   # Subset the data by gene
    #   gene_data <- genes_data[genes_data$gene == GENE,]
    #   # Find shoulder
    #   sorted_vector <- sort(gene_data$Ct)
    #   # Fit a spline curve with a high spar smoothing parameter
    #   smoothed_curve <- smooth.spline(sorted_vector, spar = 1)
    #   # Find the second derivative of the smoothed curve
    #   smoothed_diff <- diff(smoothed_curve$y, differences = 1)
    #   smoothed_x_diff <- diff(smoothed_curve$x, differences = 1)
    #   second_derivative <- diff(smoothed_diff) / (head(smoothed_x_diff, -1) * tail(smoothed_x_diff, -1))
    #   # Compute the absolute value of the second derivative
    #   abs_second_derivative <- abs(diff(second_derivative))
    #   sorted_vector <- sorted_vector[-1]
    #   # Smooth the absolute value curve
    #   original_second_derivative <- diff(diff(sorted_vector))
    #   abs_second_derivative <- abs(original_second_derivative)
    #   # Find the peaks of the smoothed absolute value curve
    #   peaks <- as.data.frame(pracma::findpeaks(abs_second_derivative))
    #   peaks <- peaks[which(peaks$V1>=0.2),"V2"]
    #   # Extract the corresponding x-coordinates from the smoothed curve
    #   shoulder_points <- data.frame(x = smoothed_curve$x[peaks], y = smoothed_curve$y[peaks])
    #   # Machine Limit_of_Detection
    #   Limit_of_Detection <- max(shoulder_points$y)
    #   
    #   # Plot the curve
    #   plot(sorted_vector, col = "black", type = "l", main = paste("Smoothed Curve of raw Ct values for gene", GENE), 
    #        sub=paste("upper Detection Threshold",round(Limit_of_Detection,0),sep=" " ),lwd=3)
    #   points(shoulder_points$x,shoulder_points$y , col = "red", pch = 19)
    # }
    
    ###########################################################
    #######   CHECK FOR DISTRIBUTION OF MISSING DATA   ########
    ###########################################################
    
    # INVESTIGATING NON-DETECTS
    # For each gene, we compute the proportion of non-detects and the average Ct value across replicate samples (Fig. 3).
    # There appears to be a strong relationship between the average expression of the genes across replicate samples and the proportion of non-detects.
    
    # Calculate the mean Ct value per group and the proportion of missing values
    mean_data <- genes_data %>%
      group_by(gene) %>%
      summarise(
        propNA = ifelse(all(is.na(Ct)), NA, sum(is.na(Ct)/length(Ct))),
        sum_NA = ifelse(all(is.na(Ct)), NA, sum(is.na(Ct))),
        mean_Ct =  ifelse(all(is.na(Ct)), NA, mean(Ct, na.rm = TRUE)))
    
    # Plot the distribution of NA by mean Ct per group
    ggplot(mean_data, aes(x = mean_Ct, y = propNA, label= gene)) +
      geom_point() +
      geom_smooth() +
      geom_text(nudge_y = 0.1) +
      xlab("Mean Ct per group") +
      ylab("Proportion of missing values in Ct") +
      ggtitle("Distribution of NA in Ct by mean Ct per group")
    # it seems that genes with lower average expression are far more likely to be non-detects.
    # From this we can conclude that the non-detects do not occur completely at random.
    
    #First, the PCR reactions are run for a fixed number of cycles (typically 40), implying that the observed data are censored at the maximum cycle number. This is a type of non-random missingness in which the missing data mechanism depends on the unobserved value. Knowledge of the technology allows us to conclude that the data are at least subject to fixed censoring; however, as we will later show, the qPCR censoring mechanism may actually be a probabilistic function of the unobserved data.
    #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4133581/
    
    mean_data_Col <- genes_data %>%
      group_by(gene,Colony) %>%
      summarise(
        propNA = sum(is.na(Ct)/length(Ct)),
        sum_NA = sum(is.na(Ct)),
        mean_Ct = mean(Ct, na.rm = TRUE),
        Treatment = Treatment
      )
    mean_data_Col <- mean_data_Col %>% distinct()
    
    if (EXPLORE_PLOT) {
      # Plot the distribution of NA by mean Ct per group
      ggplot(mean_data_Col, aes(x = gene, y = propNA,colour=Treatment)) + #, label= Colony
        geom_boxplot() +
        geom_point(position=position_jitterdodge()) +
        # geom_text(nudge_y = 0.1) +
        xlab("Mean Ct per group") +
        ylab("Proportion of missing values in Ct") +
        ggtitle("Distribution of NA in Ct by mean Ct per group")
      
    }
    # ####### CHECK FOR Ct values which are very different from the mean, per each gene:group (unreliable reads) ########################################################
    # 
    # ######## REMOVE CTs which are more than 1 Sd from mean, therefore misreads
    # # Group by gene and calculate mean and sd of Ct values
    # genes_data_grouped <- genes_data %>%
    #   group_by(gene,GROUP) %>%
    #   summarise(gene_Group_mean_ct = ifelse(all(is.na(Ct)), NA, mean(Ct, na.rm = TRUE)),
    #             gene_Group_sd_ct = ifelse(all(is.na(Ct)), NA, sd(Ct, na.rm = TRUE)))
    # 
    # # Join the mean and sd values to the original dataframe
    # genes_data <- left_join(genes_data, genes_data_grouped, by = c("gene","GROUP"))
    # 
    # # Add a new column to the dataframe with a flag to tell if values have to be discarded as too distant from mean
    # # 1. if an induvidual Ct value is smaller or greater than X*standard deviation, assign "discard" flag.
    # # 2. if both values of a specific gene_rep have the flag "discard" flag, assign "keep". (they will be dealt with in the Ct_divergence part)
    # # leave the remaining values unchanged.
    # genes_data <- genes_data %>%
    #   mutate(Excessive_divergence = ifelse(Ct < gene_Group_mean_ct - 2 * gene_Group_sd_ct | Ct > gene_Group_mean_ct + 2 * gene_Group_sd_ct, "discard", "keep"))
    # ##if both values are out of range, keep them - Group by gene_rep again and update Excessive_divergence flags
    # genes_data <- genes_data %>%
    #   group_by(gene,Code) %>%
    #   mutate(Excessive_divergence = ifelse(all(Excessive_divergence == "discard"), "keep", Excessive_divergence)) %>%
    #   ungroup()
    # 
    # #keep reps for which only 1 read was available
    # genes_data[is.na(genes_data$Excessive_divergence),"Excessive_divergence"] <- "keep"
    # #### DISCARD PHASE: assign NA to Ct values which are to be discarded, so that they will not influence the mean value
    # table(genes_data$Excessive_divergence,genes_data$gene, useNA="ifany")
    # warning("here is where files too far away from mean +- 2sd are discarded")
    # genes_data[which(genes_data$Excessive_divergence=="discard"),"Ct"] <- NA
    
    ###########################################################
    #######   STEP 1: REMOVE INVALID HOUSEKEEPING CTs   #######
    ###########################################################
    
    # valid: both Ct<mean+-2*sd, diff.Ct<0.5
    # invalid: everything else
    
    ##### ISOLATE SUSPICIOUS DATAPOINTS (CONTAMINATIONS, BAD EXTRACTIONS, ETC)
    # Calculate mean and standard deviation by gene
    # here we remove NAs to be able to calculate means by gene
    genes_data_byGroup <- genes_data %>%
      group_by(gene) %>%
      summarise(st.dev_mean_Ct = ifelse(all(is.na(Ct)), NA, sd(Ct, na.rm = TRUE)), 
                grand_mean_Ct = ifelse(all(is.na(Ct)), NA, mean(Ct, na.rm = TRUE)))
    
    
    # combine all the duplicates to have meanCts
    # Group by gene and calculate mean and sd of Ct values
    # Note: means will be NA if one of the two values is NA
    genes_data_Dups <- genes_data %>%
      group_by(Code, gene) %>%
      summarise(mean_Ct = ifelse(all(is.na(Ct)), NA, mean(Ct)), #, na.rm = TRUE) #  #without ifelse, if both NA, NaNs will be created
                mean_Tm = ifelse(all(is.na(Tm_Product)), NA, mean(Tm_Product)),
                abs_diff_Ct = ifelse(all(is.na(Ct)), NA, round(abs(diff(Ct)),3)))
    
    # Join the mean and sd values to the original dataframe
    genes_data <- left_join(genes_data, genes_data_Dups, by = c("Code", "gene")) #by duplicate
    genes_data <- left_join(genes_data, genes_data_byGroup, by = "gene") # by gene
    
    if (EXPLORE_PLOT) {#plot dataset abs_diff_Ct
      ggplot(genes_data, aes(mean_Ct, abs_diff_Ct)) +
        geom_point(alpha=0.5) +
        geom_smooth(method = "lm",formula = y ~ x + I(x^2), se = FALSE) +
        xlab("Mean Ct") +
        ylab("Absolute difference in Ct values") +
        facet_wrap(.~Ant_status + Treatment, nrow=2) 
    }
    ##### REMOVE INVALID HOUSEKEEPING GENE'S EXPRESSION
    # EF1alpha is expected to show expression values inside a gaussian range (CHECK WITH FLORENT ACCORDING TO HIS OBSERVATIONS).
    # excluding EF1 values over threshold and as we can't base the expression of the other genes on it
    # ALSO: given that we work with single ants, the closest we are to threshold, the harder it is to quantify small samples (small individuals).
    # valid datapoints: Ct < mean_Ct+-2*st.dev, diff.Ct<0.5, both Ct dups are not NA
    
    # any on the duplicate is NA
    EF1_discards <- NULL
    
    EF1_discards_NA <- genes_data %>%
      filter(gene == "EF1") %>%
      group_by(Code) %>%
      filter(any(is.na(Ct))) %>%
      ungroup()
    
    # abs_diff_Ct is > 0.5
    EF1_discards_05 <-subset(genes_data, gene == "EF1" & abs_diff_Ct > 0.5)
    # mean_Ct-2*st.dev > Ct > mean_Ct+2*st.dev
    #exclude queens from lowerbound as they have higher expression
    keepQ <- subset(genes_data, gene == "EF1" & GROUP == "QUEEN" & (Ct < (grand_mean_Ct - 2*st.dev_mean_Ct) ))
    EF1_discards_dev <-dplyr::anti_join(subset(genes_data, gene == "EF1" & (Ct < (grand_mean_Ct - 2*st.dev_mean_Ct) | Ct > (grand_mean_Ct + 2*st.dev_mean_Ct) )),keepQ )
    #keep all row couples
    EF1_discards_dev <- subset(genes_data, gene == "EF1" & Code %in% EF1_discards_dev$Code)
    
    #add labels for table
    EF1_discards_NA$Category <- "any_is_NA"
    EF1_discards_05$Category <- "diff_Ct>0.5"
    EF1_discards_dev$Category <- "mean+/-sd"
    
    EF1_discards <- rbind(EF1_discards_NA,
                          EF1_discards_05,
                          EF1_discards_dev)
    
    table(EF1_discards$Ant_status,EF1_discards$Category)
    
    if (EXPLORE_PLOT) {
      # create a frequency table of Code and sort it by frequency
      freq_table <- table(EF1_discards$Code)
      sorted_codes <- names(sort(freq_table, decreasing = TRUE))
      
      # explore pattern of missingness to select possible samples to re-run 
      ggplot(data = EF1_discards, aes(x = factor(Code, levels = sorted_codes), fill = Category)) +
        geom_bar() + labs(title = "Code vs Category", x = "Code", y = "Count") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +  facet_wrap(~Ant_status, ncol = 1)
      
      # # explore pattern of missingness for treated nurses by colony size
      # ggplot(data = EF1_discards, aes(x = factor(Code, levels = sorted_codes), fill = Category)) +
      #   geom_bar() + labs(title = "Code vs Category", x = "Code", y = "Count") +
      #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +  facet_wrap(~Ant_status, ncol = 1)
      
      #keep only treated nurses
      df_subset <- subset(EF1_discards, Ant_status == "treated nurse")
      #visualise overlaps of invalidity by gene for the same samples
      ggplot(data = df_subset, aes(x = factor(Code, levels = sorted_codes), fill = Category)) +
        geom_bar() + labs(title = "Code vs Category (only displaying treated nurses)", x = "Code", y = "Count") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +  facet_wrap(~Treatment, ncol = 1)
    }
    #get individual codes
    EF1_discards <- unique(EF1_discards$Code)
    # Remove invalid EF1 datapoints and all the points depending on them (4genes*2dups per code)
    genes_data <- genes_data[which(!genes_data$Code %in% EF1_discards),]
    
    ###########################################################
    #######            STEP 2: CHECK MEAN_Ct            #######
    ###########################################################
    
    #here, several checks are implemented on the Mean_Ct to assign a flag "Category", which is later used to manipulate the dataset.
    
    ##test DF
    # DF <- structure(list(Ct = c(NA, NA, 39.27, NA, NA, 36.99, 35.8,30.74),
    #                      gene = structure(c(3L, 4L, 1L, 3L, 2L, 1L, 2L, 4L),
    #                      .Label = c("DEF","EF1", "HYM", "PO"), class = "factor"), Code = c("11-A7", "11-A7","11-A7", "11-A7", "11-A7", "11-A7", "11-A7", "11-A7"),
    #                      Category = c(NA,NA, NA, NA, NA, NA, NA, NA)), row.names = 1937:1944, class = "data.frame")
    
    ## if mean_Ct is NA (1 or 2 Ct values missing):
    ## 1 NA value:
    ## - other Ct value >(Limit_of_Detection - Technical_error) -> "impute" mean
    ## - other Ct value <(Limit_of_Detection - Technical_error) -> "discard  NA Ct" -> mean_Ct for these points is recalculated as calculated as  mean(Ct, na.rm = TRUE)
    ## 2 NA values:
    ## -> "impute" mean
    
    genes_data$Category <- NA
    genes_data <- genes_data %>%
      group_by(gene, Code) %>%
      mutate(Category = ifelse(sum(is.na(Ct)) == 2, "impute(2NA)", 
                               ifelse(sum(is.na(Ct)) == 1 & sum(!is.na(Ct) & Ct >= (Limit_of_Detection - Technical_error)) == 1, "impute(>LOD-T.E.)", #When a logical vector is used in a numeric context, FALSE is converted to 0 and TRUE is converted to 1. The sum() function then sums up the resulting 0's and 1's to give the total number of TRUE values.
                                      ifelse(sum(is.na(Ct)) == 1 & sum(!is.na(Ct) & Ct < (Limit_of_Detection - Technical_error)) == 1, "discard_NA_Ct", NA))))
    
    ## if mean_Ct does not contain any Na values
    ## if mean_Ct > Detection threshold -> "impute"
    ## if mean_Ct < Detection threshold -> evaluate abs_diff_Ct
    genes_data <- genes_data %>%
      group_by(gene, Code) %>%
      mutate(Category = ifelse(!any(is.na(Ct)), #if there are NAs, reassign Category
                               ifelse(mean_Ct >= Limit_of_Detection, "impute(>DT)", "evaluate_abs_diff_Ct"),
                               Category))
    
    ##Check that all duplicates have the same label
    # Table <- genes_data %>%
    #   group_by(gene, Code) %>%
    #   summarise(same_Category = n_distinct(Category) == 1) %>%
    #   ungroup()
    # 
    # table(Table$same_Category)
    
    ###########################################################
    #######           STEP 3: CHECK Abs_Diff_Ct         #######
    ###########################################################
    
    # https://rnajournal.cshlp.org/content/23/5/811.full.pdf
    # Most qPCR reactions are conducted as replicates. To avoid unnecessary loss of data points (discarding all replicates with >0.5 cycles difference), as discussed
    # above, we used the Poisson distribution to calculate the acceptable Cq range between replicate measurements for different template numbers in the reaction.
    # The acceptable Cq range is defined as the interval in which 95% of the Cq values are expected to be found, given a certain Cq value
    
    # Assign the relevant Poisson Error Threshold "CTDiff_15PipErr" (not a flag replacement but the actual threshold values)
    genes_data$CTDiff_15PipErr <- sapply(genes_data$mean_Ct, assign_CTDiff_15PipErr)
    
    # "evaluate_abs_diff_Ct":
    # if abs_diff_Ct <= "CTDiff_15PipErr", assign "valid"
    # if abs_diff_Ct > 3, assign "discard higher Ct" (Technical Error)
    # if any other abs_diff_Ct, assign "invalid"
    genes_data$Category[genes_data$Category=="evaluate_abs_diff_Ct" & genes_data$abs_diff_Ct <= genes_data$CTDiff_15PipErr] <- "valid"
    if (ALWAYS_DISCARD==F){
      genes_data$Category[genes_data$Category=="evaluate_abs_diff_Ct" & genes_data$abs_diff_Ct > Technical_error] <- "discard_higher_Ct"
    }else{
     #standarrd pipeline 
      genes_data$Category[genes_data$Category=="evaluate_abs_diff_Ct" & genes_data$abs_diff_Ct > Technical_error] <- "invalid"
    
      # #MOD proposed by Florent:
      # # When > Poisson, check if one of the two is > LOD, 
      # # if above → invalid, 
      # # if below → discard higher Ct
      # genes_data$Category[genes_data$Category=="evaluate_abs_diff_Ct" & genes_data$abs_diff_Ct > genes_data$CTDiff_15PipErr] <- "Check_individual_LOD"
      # #Check if one of the 2 is above LOD
      # # if so, Discard the one above (technical error), else invalid (both below, high expression high variation)
      # genes_data <- genes_data %>%
      #   group_by(gene, Code) %>%
      #   mutate(Category = ifelse(Category == "Check_individual_LOD",
      #                            ifelse(sum(Ct >= Limit_of_Detection) == 1, "discard_higher_Ct", "invalid"),
      #                            Category))

    }
    
# table(genes_data$Category,genes_data$gene)
#if any other abs_diff_Ct, assign "invalid"
genes_data$Category[genes_data$Category=="evaluate_abs_diff_Ct" & (genes_data$abs_diff_Ct > genes_data$CTDiff_15PipErr & genes_data$abs_diff_Ct <= Technical_error)] <- "invalid"
    
    # indiviual sample Ants by Category 
    #as.data.frame(table(genes_data$Category,genes_data$gene,genes_data$Ant_status)/2)
    #dput(head(genes_data[,c("Category","gene","Ant_status")]))
    
    genes_sub <- unique(genes_data[which(genes_data$gene!="EF1"),c("gene", "Category","Ant_status","Code")])
    
    if (EXPLORE_PLOT) {
      ## explore all data distribution
      ggplot(genes_sub , aes(x = gene, fill = Category)) +
        geom_bar(position = "stack") +
        geom_text(position = position_stack(vjust = 0.5), aes(label =after_stat(count)), stat='count',size=2.5) +
        facet_wrap(~Ant_status, nrow=4)+
        scale_fill_brewer(palette = "Set2") +
        labs(title = "Frequency of Category by Ant_status by gene", x = "Gene", y = "Frequency") +
        coord_flip() +
        theme(legend.position = "bottom",
              legend.box = "horizontal",
              legend.direction = "horizontal",
              legend.key.size = unit(0.5, "cm"),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12, face = "bold"))+
        scale_y_continuous(limits = c(0, 380), expand = c(0, 0))
      
      ## explore the invalid data distribution
      freq_table <- table(genes_data[which(genes_data$Category=="invalid"),"Code"])
      sorted_codes <- names(sort(freq_table, decreasing = TRUE))
      # create a subset of the data with counts above 2
      df_subset <- subset(genes_data[which(genes_data$Category=="invalid"),], Code %in% names(freq_table[freq_table > 2]))
      #visualise overlaps of invalidity by gene for the same samples
      ggplot(data = df_subset, aes(x = factor(Code, levels = sorted_codes), fill = gene)) +
        geom_bar() + labs(title = "Code vs Category (only displaying Codes with at least 2 invalid genes)", x = "Code", y = "Count") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +  facet_wrap(~Ant_status, ncol = 1)
      
      #keep only treated nurses / invalid
      df_subset <- subset(genes_data, Ant_status == "treated nurse" & Category=="invalid")
      #visualise overlaps of invalidity by gene for the same samples
      ggplot(data = df_subset, aes(x = Code, fill = gene)) +
        geom_bar() + labs(title = "Code vs invalid (only displaying treated nurses)", x = "Code", y = "Count") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +  facet_wrap(~Treatment, ncol = 1)
      
      #keep only treated nurses / impute(>LOD-T.E.)
      df_subset <- subset(genes_data, Ant_status == "treated nurse" & Category=="impute(>LOD-T.E.)")
      #visualise overlaps of invalidity by gene for the same samples
      ggplot(data = df_subset, aes(x = Code, fill = gene)) +
        geom_bar() + labs(title = "Code vs impute(>LOD-T.E.) (only displaying treated nurses)", x = "Code", y = "Count") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +  facet_wrap(~Treatment, ncol = 1)
      
    }
    
    
    ###########################################################
    #######     STEP 4: RE-ASSING VALUES BY CATEGORY    #######
    ###########################################################
    
    ##How to treat the following "Category" labels:
    # "discard_NA_Ct":        Recalculate mean_Ct by duplicate as above but dropping NAs - mean(Ct, na.rm = TRUE)
    # "discard_higher_Ct":    Discard higher Ct value of the pair
    # "valid":                Use Mean_Ct
    # "impute":               Replace by NAs
    # "invalid":              either discard, use mean or re-run sample (810, 405 duplicates. 336 unique samples)
    
    # "discard_NA_Ct":        Recalculate mean_Ct by duplicate as above but dropping NAs - mean(Ct, na.rm = TRUE)
    genes_data <- genes_data %>%
      group_by(gene, Code) %>%
      mutate(mean_Ct = ifelse(Category == "discard_NA_Ct", mean(Ct, na.rm = TRUE), mean_Ct))
    
    # "discard_higher_Ct":    Discard higher Ct value of the pair (Technical Error) - mean_Ct is the min. value
    genes_data <- genes_data %>%
      group_by(gene, Code) %>%
      mutate(mean_Ct = ifelse(Category == "discard_higher_Ct", min(Ct), mean_Ct))
    
    # "impute":               Replace by NAs
    genes_data <- genes_data %>%
      group_by(gene, Code) %>%
      mutate(mean_Ct = ifelse(grepl("impute", Category), NA, mean_Ct))
    
    # "invalid":              either discard, use mean (already calculated) or re-run sample (810, 405 duplicates. 336 unique samples)
    ## DISCARD DATASET
    genes_data_invDiscard <- genes_data[which(genes_data$Category!="invalid"),]
    #MEAN DATASET IS genes_data
    
    ######### LOOP THROUGH THE TWO VARIANTS OF HOW TO TREAT INVALID SAMPLES: ####
    # Combine the data frames into a list
    df_list <- list(Invalids_MEAN=genes_data, Invalids_DISCARD=genes_data_invDiscard)
    
    df_list <- Map(cbind, df_list, name_df = names(df_list))
    
    ### it could be done more cleanly, but for the moment select the DF to use here:
    warning("select the DF to use here. Invalids should be discarded instead of using the means")
    for (INVALIDS in c("Invalids_DISCARD")) { #,"Invalids_MEAN"
      
      CLEAN_DATA <- df_list[[INVALIDS]]
      #CLEAN_DATA <- df_list$Invalids_DISCARD # df_list$Invalids_MEAN
      NAME_DF    <- unique(CLEAN_DATA$name_df)
      # # Loop through the list and perform a summary analysis
      # lapply(names(df_list), function(CLEAN_DATA) {
      #   nrow(df_list[[CLEAN_DATA]])
      #   cat(CLEAN_DATA,"\n")
      #   })
      
      
      ###########################################################
      #######             CALCULATE DELTA-CT              #######
      ###########################################################
      
      # remove extra cols to compress dataframe
      # here is important to drop the Sample_Well as the extra Queens duplicates won't be eliminated if kept.
      CLEAN_DATA <- CLEAN_DATA[,!(names(CLEAN_DATA) %in% c("Threshold","Ct","Tm_Product","name","machine_num","Sample_Well"))] #c("CTDiff_15PipErr","Category")
      CLEAN_DATA <- CLEAN_DATA %>% distinct()
      # # combine all the duplicates to have meanCts
      ## remove columns which are different, then remove duplicated rows
      
      # “delta-Ct”: difference between the Ct of the housekeeping gene and the test gene
      # split test genes from housekeeping
      test_gene_data <- CLEAN_DATA %>%
        filter(gene != "EF1")
      # isolate housekeeping and change colnames
      ref_gene_data <- CLEAN_DATA %>%
        filter(gene == "EF1") %>%
        rename("housekeeping" = "gene", "ref_Ct" = "mean_Ct") %>% #, "ref_Tm" = "mean_Tm"
        #only keep certain columns
        dplyr::select(housekeeping, ref_Ct, Code)
      
      # recombine data
      # common_col_names4 <- intersect(names(ref_gene_data), names(test_gene_data))
      # common_col_names4 <- common_col_names4[ !common_col_names4 == c("mean_Tm","abs_diff_Ct")] #exclude numeric var
      CLEAN_DATA <- left_join(test_gene_data, ref_gene_data, by = "Code")
      
      # create a new column containing the delta Ct between the housekeeping gene and our gene of interest, and plot the delta Ct for each treatment and replicate.
      CLEAN_DATA <- mutate(CLEAN_DATA, delta_Ct = mean_Ct - ref_Ct)
      
      # # plot delta_Ct
      # ggplot(CLEAN_DATA, aes(x = gene, y = delta_Ct, fill = Treatment)) +
      #   geom_boxplot(aes(colour = AntTask), lwd = 0.8, alpha = 0.3)
      
      # # Calculate the mean delta Ct for each treatment.
      # treatment_summary <- CLEAN_DATA %>%
      #   group_by(gene) %>%
      #   summarise(mean_delta_Ct = mean(delta_Ct, na.rm = TRUE))
      
      
      # Calculating relative DNA concentration
      # to calculate the relative DNA concentration, you can use the fact that the amount of cDNA theoretically doubles every cycle.
      # PRIMERS EFFICIENCY
      primers_eff <- data.frame(gene = c("EF1", "HYM", "DEF", "PO"), efficiency = c(1.99, 1.99, 2, 1.94))
      # add primers efficiency data
      CLEAN_DATA <- left_join(CLEAN_DATA, primers_eff, by = c("gene"))
      
      # as EF1 efficiency is 1.99, it is used as baseline for normalisation
      CLEAN_DATA <- CLEAN_DATA %>%
        mutate(rel_conc = (2 * (efficiency / 1.99))^-delta_Ct)
      
      table(CLEAN_DATA$gene,is.na(CLEAN_DATA$rel_conc),useNA = "ifany")
      
      
      
      
      
      #dput(as.data.frame(CLEAN_DATA[which(CLEAN_DATA$Code=="1-F2"),c(5,10,11,24,26,29)]))
      
      
      
      
      
      
      ############
      #### PROPORTION OF ZERO LOAD VS Ant_status BY TREATMENT
      #create flag to show which values are valid
      CLEAN_DATA$Final_Cat <- as.factor(ifelse(is.na(CLEAN_DATA$rel_conc),"to-impute","valid"))
      
      table(CLEAN_DATA$Final_Cat, CLEAN_DATA$Treatment, CLEAN_DATA$Ant_status)
      table(CLEAN_DATA$Final_Cat, CLEAN_DATA$gene, CLEAN_DATA$Ant_status)
      
      ## explore all data distribution
      ggplot(CLEAN_DATA , aes(x = gene, fill = Final_Cat)) +
        geom_bar(position = "stack") +
        geom_text(position = position_stack(vjust = 0.5), aes(label =after_stat(count)), stat='count',size=3) +
        facet_wrap(~Ant_status+Treatment, nrow=4)+
        #scale_fill_brewer(palette = "Set2") +
        labs(title = "Frequency of Category by Ant_status by gene", x = "Gene", y = "Frequency") +
        coord_flip() +
        theme(legend.position = "bottom",
              legend.box = "horizontal",
              legend.direction = "horizontal",
              legend.key.size = unit(0.5, "cm"),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12, face = "bold"))+
        scale_y_continuous(limits = c(0, 430), expand = c(0, 0))
      
      warning("to fix, the selection deletes the ant status, and at a later stage the propZeros")
      # CLEAN_DATA_ColPropPos <- CLEAN_DATA %>%
      #   group_by(Colony, Treatment, gene) %>%
      #   summarise(count = n(), Final_Cat = sum(Final_Cat == 1))
      # CLEAN_DATA_ColPropPos$negative <- CLEAN_DATA_ColPropPos$count - CLEAN_DATA_ColPropPos$Final_Cat
      # # prop of zero
      # CLEAN_DATA_ColPropPos$propZeros <- (CLEAN_DATA_ColPropPos$negative / CLEAN_DATA_ColPropPos$count) * 100
      # #Standard error
      # CLEAN_DATA_ColPropPos <- CLEAN_DATA_ColPropPos %>% group_by(Treatment, gene) %>% summarise(se_propZeros = sqrt(propZeros*(100-propZeros)/count)/100, propZeros=propZeros,Colony=Colony)
      # 
      # # select only relevant cols
      # # CLEAN_DATA_ColPropPos <- CLEAN_DATA_ColPropPos[, c("Treatment", "Ant_status", "propZeros")]
      # 
      # 
      #   group_by(Treatment, gene) %>%
      #   summarise(se_propZeros = sqrt(propZeros * (100 - propZeros) / count) / 100,
      #             propZeros = mean(propZeros))  # add propZeros to the final summary
      # 
      # 
      # 
      # # #bad way to get table, it works as there are single vals per each condition...
      # # #table useless? not possible to test percentages...
      # # CLEAN_DATA_ColPropPosTABLE <- tapply(CLEAN_DATA_ColPropPos$propZeros, list(CLEAN_DATA_ColPropPos$Treatment, CLEAN_DATA_ColPropPos$Ant_status), mean)
      # 
      # ggplot(CLEAN_DATA_ColPropPos, aes(fill = Treatment, y = propZeros, x = Ant_status)) +
      #   geom_bar(position = "dodge", stat = "identity") +
      #   STYLE +
      #   labs( # title = "Pathogen Quantification Adriano",
      #     # subtitle = "All ants",
      #     y = "% of close to 0 values",
      #     # x="",
      #     # caption = paste("Threshold cycle (Ct) missing for", N_ants_NoCt , "ants in cols", No_CT_REPs)
      #   )
      
      
      
      ####################################################################################
      
      #### PERFORM DATA IMPUTATION HERE
      #### WRITE HERE REASONS TO DROP QRILC (MORE THAN 80% OF POINTS MISSING, THEREFORE HALF-MIN IS BETTER)
      
      #create new assign minimum and imputed conc. columns
      CLEAN_DATA$rel_conc_imputed <- CLEAN_DATA$rel_conc
      # CLEAN_DATA$rel_conc_QC_replace_min <- CLEAN_DATA$rel_conc_QC
      
      
      if (IMPUTATION=="HM") {
        # # ASSIGN HALF-MINIMUM
        # #assing value smaller than the minimum to the rel_conc_QC NAs
        print("replacing NA relative concentrations with min/2 (equivalent to mean_Ct/sqrt(2) as each Ct step is a factor 2 doubling of product)")
        for (GENE in unique(CLEAN_DATA$gene)) {
          for (STATUS in unique(CLEAN_DATA$Ant_status)) {
            #assign the min/2 value by gene to missing datapoints
            CLEAN_DATA[is.na(CLEAN_DATA$rel_conc_imputed) & CLEAN_DATA$gene == GENE & CLEAN_DATA$Ant_status == STATUS, "rel_conc_imputed"] <- min(CLEAN_DATA[which(CLEAN_DATA$gene==GENE  & CLEAN_DATA$Ant_status == STATUS),"rel_conc"],na.rm = T)/2
          }
        }
      }
      
      
      if (IMPUTATION=="QRILC") {
        print("replacing NA relative concentrations with QRILC Quantile Regression Imputation of Left Censored Data, proven the most reliable imputation method in Wei, R et al. 2018 ")
        ### ASSIGN IMPUTATED VALS, AT THE STAGE OF THE REL.CONCENTRATION, AS IT FOLLOWS THE PAPER
        ### Wei, R., Wang, J., Su, M. et al. Missing Value Imputation Approach for Mass Spectrometry-based Metabolomics Data. Sci Rep 8, 663 (2018). https://doi.org/10.1038/s41598-017-19120-0
        ### QRILC
        # A missing data imputation method that performs the imputation of
        # left-censored missing data using random draws from a truncated
        # distribution with parameters estimated using quantile
        # regression. Implemented in the `imputeLCMD::impute.QRILC`
        # result <- data %>% log %>% impute.QRILC(., ...) %>% extract2(1) %>% exp
        
        for (GENE in unique(CLEAN_DATA$gene)) {
          for (STATUS in unique(CLEAN_DATA$Ant_status)) {
            rel_conc_matrix <- as.matrix(CLEAN_DATA[which(CLEAN_DATA$gene==GENE &  CLEAN_DATA$Ant_status == STATUS), "rel_conc"])
            # note that tune.sigma = 1 as it is assumed that the complete data distribution is Gaussian (as shown by the complete gene DEF) .
            imputed_rel_conc <- impute.QRILC(log(rel_conc_matrix),tune.sigma = 1) # log transform data.. see explanation in paper
            CLEAN_DATA[which(CLEAN_DATA$gene == GENE  &  CLEAN_DATA$Ant_status == STATUS), "rel_conc_imputed"] <- exp(unlist(imputed_rel_conc[1])) # back-transform with exp
            sum(is.na(CLEAN_DATA[is.na(CLEAN_DATA$rel_conc) & CLEAN_DATA$gene == GENE  &  CLEAN_DATA$Ant_status == STATUS, "rel_conc_imputed"])) # should be 0
            
          }
        }
      }
      
      
      
      if (PLOT) {
      # Create a list to store the plots
      plot_list <- list()
      # Loop through each gene
      for (GENE in unique(CLEAN_DATA$gene)) {
        # Generate the plot and add it to the list
        plot_list[[GENE]] <- ggplot(CLEAN_DATA[which(CLEAN_DATA$gene == GENE),], aes(x = rel_conc_imputed, fill = factor(Final_Cat, levels = c("to-impute", "valid"))  )) + #, fill = factor(Final_Cat, levels = c("to-impute", "valid"))
          geom_histogram(alpha = 0.5, position = "stack") +
          labs(x = "Value", y = "Frequency", fill = "Final_Cat") +
          ggtitle(GENE) +
          scale_x_continuous(trans='log10') +
          theme(legend.position = "none")
      }
      
      # Create a legend plot. Returns a gtable
      plot_list[["legend"]]  <- cowplot::get_legend(plot_list[["HYM"]] + theme(legend.position = "right"))
      
      # Combine the plots using grid.arrange()
      grid.arrange(
        plot_list[[1]], plot_list[[2]],
        plot_list[[3]], plot_list[[4]],
        ncol = 2, nrow = 2,
        widths = c(1, 1),
        heights = c(1, 1),
        bottom = paste("QRILC imputation, stacked histogram. Data=", NAME_DF)
      )
      
      }
      # ## Check that distributions are normal
      # # the imputation assumes that the data distribution is normal, as expected from qPCR data. 
      # # we can test that the imputation does not cause the distributions to diverge, considering that DEF is our complete distribution
      # CLEAN_DATA <- CLEAN_DATA %>% group_by(gene) %>% mutate( rel_conc_imput_LogStand = scale(log(rel_conc_imputed),center=T,scale=T) )
      # 
      # # create the ANOVA model
      # model <- aov(rel_conc_imput_LogStand ~ gene, data = CLEAN_DATA)
      # summary(model)
      # 
      # # check if the ANOVA is significant
      # if (summary(model)[[1]]$"Pr(>F)"[1] < 0.05) {
      #   # if the ANOVA is significant, print the significance value
      #   pvalue <- summary(model)[[1]]$"Pr(>F)"[1]
      #   cat("The ANOVA is significant (p =", pvalue, ")\n")
      # } else {
      #   pvalue <- "ns"
      #   # if the ANOVA is not significant, print "ns"
      #   cat("The ANOVA is not significant\n")
      # }
      # 
      # 
      # # plot standardised 
      # # plot smoothed curve with separate colors for each GENE
      # # note: rel_conc_imputed variable has negative values, this causes a flipping in the standardised distribution
      # print(
      #   ggplot(CLEAN_DATA, aes(x = rel_conc_imput_LogStand, colour = gene)) +
      #     geom_density(alpha = 0.5) +
      #     labs(x = "Log transformed standardised relative concentrations", y = "Density", fill = "gene") +
      #     ggtitle("all genes distribution", subtitle = paste("rel_conc_imputed and standardised.", "Data=",NAME_DF,sep=" ")) +
      #     annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, label = paste("ANOVA, P=", pvalue))
      # )
      # 
      
      ####################################################################################
      if (EXPLORE_PLOT) {
        # CLEAN_DATA[is.na(CLEAN_DATA$rel_conc),"rel_conc"] <- 1e-06
        for (REL_CONC in c("rel_conc","rel_conc_imputed")) { # "rel_conc_replace_min"
          
          # plot the relative concentration.
          print(
            ggplot(CLEAN_DATA, aes(x = Ant_status, y = !!sym(REL_CONC),colour=Treatment)) +
              geom_boxplot(aes(colour = Treatment), lwd = 0.8, alpha = 0.3) +
              colScale_Treatment +
              STYLE +
              scale_y_continuous(labels = scales::percent,trans='log10') + #, limits = c(1e-08,1000)
              ggtitle(REL_CONC, subtitle =paste("Data=",NAME_DF,sep=" ")) +
              facet_grid(. ~ gene)
          )
        }
      }
      # proportion imputed by ant_status
      round(prop.table(table(CLEAN_DATA$Ant_status,CLEAN_DATA$Final_Cat)),2)
      # proportion imputed by gene
      round(prop.table(table(CLEAN_DATA$gene,CLEAN_DATA$Final_Cat)),2)
      
      ###############
      # clean debris
      remove_debris <- c("CLEAN_DATA_byGroup", "outliers", "test_gene_data", "ref_gene_data", "common_col_names4", "outliers_limit", "primers_eff")
      # cleaning
      rm(list = ls()[which(ls() %in% remove_debris)])
      gc()
      
      ################
      if (EXPLORE_PLOT) {
      ggplot(
        data = CLEAN_DATA,
        aes(x = Ant_status, y = rel_conc_imputed, color = Treatment)
      ) + # ,group = Exposed,color = Exposed
        #geom_point(aes(fill = Treatment, colour = Colony), position = position_jitter() , size = 1, alpha = 0.3, show.legend = FALSE) + # # ,alpha= 0.8,stroke=0
        geom_point(position=position_jitterdodge(), size = 1, alpha = 0.3, show.legend = FALSE) +
        #colScale_Colony +
        # new_scale_color() +
        # new_scale_fill() +
        geom_boxplot(aes(colour = Treatment), lwd = 0.8, alpha = 0.5) +
        #geom_violinhalf(aes(fill = Treatment), trim = FALSE, scale = "width", adjust = 1, draw_quantiles = c(0.25, 0.75)) + 
        colScale_Treatment +
        #geom_boxplot(aes( color=alpha("black",0.3)), lwd = 0.8, alpha = 0.3)+
        STYLE +
        scale_y_continuous(labels = scales::percent, trans = "log10")+
        facet_grid(. ~ gene) 
      
      
      
      ggplot(data = CLEAN_DATA, aes(x = Ant_status, y = rel_conc_imputed, fill = Treatment)) +
        geom_violin(aes(fill= Treatment), trim = FALSE,width =1.2) +
        geom_boxplot(aes(fill= Treatment), width = 0.1, alpha = 0.5, position = position_dodge(width = 1.2) )+
        colFill_Treatment +
        STYLE +
        scale_y_continuous(labels = scales::percent, trans = "log10") +
        facet_grid(. ~ gene) 
      
      
      # library(see) # for half violin plots
      
      
      # ggplot(
      #   data = CLEAN_DATA,
      #   aes(x = Ant_status, y = rel_conc_imputed)
      # ) + 
      #   geom_jitter(aes(group = Treatment, colour = Colony), size = 1, alpha = 0.3, width = 0.2, height = 0, show.legend = FALSE) +
      #   colScale_Colony +
      #   STYLE +
      #   scale_y_continuous(labels = scales::percent, trans = "log10") +
      #   facet_grid(. ~ gene) 
      
      
      # for (GROUP in unique(CLEAN_DATA$GROUP)) {
      #   print(
      #     ggplot(
      #       data = CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),],
      #       aes(x = gene, y = rel_conc_imputed)
      #     ) + # ,group = Exposed,color = Exposed
      #       geom_point(aes(fill = Treatment, colour = Colony), position = position_jitterdodge(), size = 1, alpha = 0.3) + # ,alpha= 0.8,stroke=0
      #       geom_boxplot(aes(colour = Treatment), lwd = 0.8, alpha = 0.3) + #
      #       STYLE +
      #       ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
      #       scale_y_continuous(labels = scales::percent, trans = "log10") #+
      #     # facet_grid(. ~ gene)
      #   )
      # }
      
      # ggplot(
      #         data = CLEAN_DATA,
      #         aes(x = Ant_status, y = rel_conc_imputed)
      #       ) + # ,group = Exposed,color = Exposed
      #         geom_point(aes(fill = Treatment, colour = Colony), position = position_jitter(), size = 1, alpha = 0.3) + # ,alpha= 0.8,stroke=0
      #         geom_boxplot( lwd = 0.8, alpha = 0.3) + #
      #         STYLE +
      #         ggtitle(paste(unique(CLEAN_DATA[,"Ant_status"]), collapse = ", ")) +
      #         scale_y_continuous(labels = scales::percent, trans = "log10") +
      #        facet_grid(. ~ gene)
      }
      
      #####################################################################################
      ##################              STATS & PLOTS        ################################
      #####################################################################################
      
      #####################################################################################
      ##### TEST DIFFERENCE OF RELATIVE CONCENTRATION OF IMMUNE GENES TRANSCRIPTS BETWEEN TREATMENTS (AND EVENTUALLY GROUPS FOR untreated workers)
      
      ### How the analysis is conducted
      ## Ant Statuses are not always comparable for the following reasons:
      
      # - Queens
      # Queens' larger amount of tissue makes the Ct values for EF1 much lower (high expression).
      # this is not a problem as each sample housekeeping gene EF1 Ct value is used to normalise the test genes, but:
      # as queens have different metabolic rates (which are known to impact EF1 expression), the Ct values for EF1 will not be comparable,
      # therefore QUEENS HAVE TO BE ANALYSED SEPARATELY
      
      # - Treated nurses
      # the treatment procedure consists in a stressful event that will vary the ant baseline expression,
      # therefore it is not meaningful to compare treated and non-treated ants
      
      # - Untreated nurses and foragers
      # these ants can be compared
      
      CLEAN_DATA$Treatment <- as.factor(CLEAN_DATA$Treatment)
      CLEAN_DATA <- as.data.frame(CLEAN_DATA)
      
      QUEEN_data_OK <- TRUE # queen data may be worth to reprocess as a few get discarded
      warning(paste0("QUEEN_data_OK: ", QUEEN_data_OK))
      
      #create list of significance outputs for Q data
      mod_Q_list <- list()
      #create list of significance outputs for treated data
      mod_T_list <- list()
      ##temporary list for untreated
      mod_UN_list <- list()
      
      
      for (GROUP in unique(CLEAN_DATA$GROUP)) {
        print(GROUP)
        table(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Treatment"],CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"gene"])
        
        for (GENE in unique(CLEAN_DATA$gene)) {
          # select subset for 1 gene
          #GENE_data <- DF[which(DF$gene == GENE), ]
          GENE_data <- CLEAN_DATA[which(CLEAN_DATA$gene==GENE & CLEAN_DATA$GROUP==GROUP),
                                  c("Colony","GROUP","Ant_status","gene","Treatment","rel_conc_imputed")]
          
          # create constant to make values below zeros a small number
          # GENE_cost <- min(GENE_data$rel_conc_imputed, na.rm = T) * 1.1
          # hist(log10(GENE_data$rel_conc_imputed+GENE_cost))
          
          #### TREATED INDIVIDUALS (no status grouping)
          if (unique(GENE_data$GROUP)=="TREATED_W") {
            # First, fit linear models to explain variation in density
            # descdist(GENE_data$rel_conc_imputed)
            
            # # Fit a GAMLSS model with lognormal distribution
            # gamlss_model <- gamlss(rel_conc_imputed ~ Treatment + random(as.factor(Colony)),
            #                        data = GENE_data,
            #                        family =  ZAGA())
            # # Plot residuals vs. fitted values
            # # plot(fitted(gamlss_model), residuals(gamlss_model), xlab = "Fitted values", ylab = "Pearson residuals")
            # # abline(h = 0, lty = 2, col = "red")
            # # QQ plot of residuals
            # Shap <- shapiro.test(residuals(gamlss_model))
            # qqnorm(residuals(gamlss_model),main = paste(GROUP,GENE,"qqnorm","\nshap.test p=", round(Shap$p.value,4),sep=" "))
            # qqline(residuals(gamlss_model))
            # #https://www.gamlss.com/wp-content/uploads/2013/01/gamlss-manual.pdf
            # wp(gamlss_model,xvar=~Treatment)
            # 
            # #extract significance
            # gamlss_sig <- as.data.frame(summary(gamlss_model, coef.only = TRUE))
            # gamlss_sig$p <- round(gamlss_sig$"Pr(>|t|)",3)
            # gamlss_sig$"Pr(>|t|)" <- NULL
            # mod_T_list <- c(mod_T_list, list(gamlss_sig[!grepl("Intercept", rownames(gamlss_sig)), ]))
            # ID_model <- paste(GROUP,GENE,sep="-")
            # names(mod_T_list)[length(mod_T_list)] <- paste(ID_model, sep = "-")
            # mod_T_list[length(mod_T_list)]
            # 
            
            # # Fit the GAMM
            # gamm_model <- gamm4((log10(rel_conc_imputed + GENE_cost)) ~ Treatment, random = ~ (1 | Colony), data = GENE_data)
            # # Extract residuals
            # gamm_residuals <- resid(gamm_model$gam, type = "pearson")
            # # Plot residuals vs. fitted values
            # plot(fitted(gamm_model$gam), gamm_residuals, xlab = "Fitted values", ylab = "Pearson residuals")
            # abline(h = 0, lty = 2, col = "red")
            # # QQ plot of residuals
            # qqnorm(gamm_residuals)
            # qqline(gamm_residuals)
            # # Summary of the model
            # summary(gamm_model$gam)
            # # Summary of the random effects
            # summary(gamm_model$lme)
            # # Check residuals for the GAMM model
            # plot(gamm_model$lme, resid(., type = "pearson"), xlab = "Fitted values", ylab = "Pearson residuals")
            # abline(h = 0, lty = 2, col = "red")
            # qqnorm(resid(gamm_model$lme, type = "pearson"))
            # qqline(resid(gamm_model$lme, type = "pearson"))
            
            
            
            # # Step 3: Generalized linear mixed model (GLMM)
            # # Fit a GLMM with Gamma distribution and log link function
            # glmm_gamma <- glmmTMB(rel_conc_imputed ~ Treatment + (1 | Colony),
            #                       data = GENE_data,
            #                       family = Gamma(link = "log"))
            # qqnorm(resid(glmm_gamma, type = "pearson"))
            # qqline(resid(glmm_gamma, type = "pearson"))
            # 
            # mod1 <- lmer(log10(rel_conc_imputed + GENE_cost) ~ Treatment + (1 | Colony), data = GENE_data) # ,weights = weights
            # print(GENE)
            # output_lmer(mod1)
            # # hist(residuals(mod1))
            # # qqnorm(residuals(mod1))
            # # qqline(residuals(mod1))
            # out.put <- model.sel(mod1)
            # # models sorted from best (top) to worst (bottom). higher weight is better as it shows that the model is more likely to be the best explanation (hypothesis) for variation in rel_conc_imputed
            # sel.table <- nice_print_model_sel(out.put)
            # print(GENE)
            # print(sel.table)
            # # compute posthocs
            # sel_mod <- get(sel.table[which.max(sel.table$weight), "Model"])
            # ID_model <- paste(GROUP,GENE,sep="-") # to assign name to posthoc's list element
            # posthoc_list <- compute_posthocs(sel_mod)
            
            
            # #given the tiny group size, and the non-relevance of the random factor Colony, use a non-parametric test
            # #Wilcoxon rank-sum test (also known as the Mann-Whitney U test)
            # mod_T <- rstatix::wilcox_test(rel_conc_imputed ~ Treatment, data = GENE_data, alternative = "two.sided")
            # #save test information
            # mod_T_list <- c(mod_T_list, list(mod_T))
            # ID_model <- paste(GROUP,GENE,sep="-")
            # names(mod_T_list)[length(mod_T_list)] <- paste(ID_model, sep = "-")
            # mod_T_list[length(mod_T_list)]
            # 
          }
          #### QUEENS (no status grouping, no random)
          else if (unique(GENE_data$GROUP)=="QUEEN") {
            
            if (QUEEN_data_OK) {
              
              #given the tiny group size, and the non-relevance of the random factor Colony, use a non-parametric test
              #Wilcoxon rank-sum test (also known as the Mann-Whitney U test)
              mod_Q <- rstatix::wilcox_test(rel_conc_imputed ~ Treatment, data = GENE_data, alternative = "two.sided")
              #save test information
              mod_Q_list <- c(mod_Q_list, list(mod_Q))
              ID_model <- paste(GROUP,GENE,sep="-")
              names(mod_Q_list)[length(mod_Q_list)] <- paste(ID_model, sep = "-")
              mod_Q_list[length(mod_Q_list)]
              
              #mod1 <- lmer(rel_conc_imputed + GENE_cost ~ Treatment + (1 | Colony), data = GENE_data)
              # mod1 <- glm(log10(rel_conc_imputed + GENE_cost) ~ Treatment, data = GENE_data) # ,weights = weights
              # print(GENE)
              # output_lmer(mod1)
              # # run anova and store the significant variables in the object SIG
              # SIG <- as.data.frame(Anova(mod1)[Anova(mod1)$"Pr(>Chisq)" < 0.05])
              # if (ncol(SIG)==0) {
              #   SIG <- data.frame(LR.Chisq=NA,Df=NA,Pr="n.s.")
              # }else{
              # SIG$Pr <- round(SIG$"Pr(>Chisq)",3)
              # SIG$"Pr(>Chisq)" <- NULL
              # 
              # SIG$LR.Chisq <- SIG$"LR Chisq"
              # SIG$"LR Chisq" <- NULL
              # 
              # }
              
              # #save test information
              #mod_Q_list <- c(mod_Q_list, list(SIG))
              #ID_model <- paste(GROUP,GENE,sep="-")
              #names(mod_Q_list)[length(mod_Q_list)] <- paste(ID_model, sep = "-")
              
            }
          }
          #### UNTREATED INDIVIDUALS
          else{
            # First,fir candidate linear models to explain variation in density
            mod1 <- lmer(log10(rel_conc_imputed) ~ Treatment + Ant_status + (1 | Colony), data = GENE_data)
            #mod2 <- lmer(log10(rel_conc_imputed + GENE_cost) ~ Treatment + Ant_status + (1 | Colony), data = GENE_data) # weights = weights,
            # We can now use the mod.sel to conduct model selection. The default model selection criteria is Akaike’s information criteria (AIC) with small sample bias adjustment, AICc
            # delta AICc, and the model weights
            out.put <- model.sel(mod1) #,mod2
            # models sorted from best (top) to worst (bottom). higher weight is better as it shows that the model is more likely to be the best explanation (hypothesis) for variation in rel_conc_imputed
            sel.table <- nice_print_model_sel(out.put)    
            print(paste(GROUP,GENE,sep=" : "))
            print(sel.table)
            # compute posthocs
            sel_mod <- get(sel.table[which.max(sel.table$weight), "Model"])
            #output_lmer(sel_mod)
            ID_model <- paste(GROUP,GENE,sep="-") # to assign name to posthoc's list element
            #posthoc_list <- compute_posthocs(sel_mod)
            
            # #extract significance
            # mod_Q_list <- c(mod_Q_list, list(mod_Q))
            # ID_model <- paste(GROUP,GENE,sep="-")
            # names(mod_Q_list)[length(mod_Q_list)] <- paste(ID_model, sep = "-")
            # mod_Q_list[length(mod_Q_list)]
            # 
            # 
            # gamlss_sig <- as.data.frame(summary(gamlss_model, coef.only = TRUE))
            # gamlss_sig$p <- round(gamlss_sig$"Pr(>|t|)",3)
            # gamlss_sig$"Pr(>|t|)" <- NULL
            # mod_T_list <- c(mod_T_list, list(gamlss_sig[!grepl("Intercept", rownames(gamlss_sig)), ]))
            # ID_model <- paste(GROUP,GENE,sep="-")
            # names(mod_T_list)[length(mod_T_list)] <- paste(ID_model, sep = "-")
            # mod_T_list[length(mod_T_list)]
            
            ## Reporting  model estimates for fixed effects in Untreated nurses tests
            estimates <- summary(mod1)$coefficients
            fixed_effects <- estimates[grep("Treatment|Ant_status", rownames(estimates)), 1]
            intercept <- estimates[1, 1]
            # Add fixed effects to intercept
            final_estimates <- intercept + fixed_effects
            Estim.Status <- final_estimates[grep("Ant_status",names(final_estimates))]
            Estim.Treatment <-final_estimates[grep("Treatment",names(final_estimates))]
            
            
            ## Reporting Pvalues in Untreated nurses tests
            P.Status <- Anova(mod1)[grep("Ant_status", row.names(Anova(mod1))),"Pr(>Chisq)"]
            P.Treatment <- Anova(mod1)[grep("Treatment", row.names(Anova(mod1))),"Pr(>Chisq)"]
            # Add results to the data frame
            PipelineTesting <- rbind(PipelineTesting, data.frame(gene = GENE,
                                                                 P.Status = P.Status,
                                                                 P.Treatment = P.Treatment,
                                                                 Estim.Status = Estim.Status,
                                                                 Estim.Treatment = Estim.Treatment,
                                                                 T.E. = Technical_error,
                                                                 LOD = Limit_of_Detection,
                                                                 ALWAYS_DISCARD=ALWAYS_DISCARD,
                                                                 imputation  = IMPUTATION,
                                                                 Invalids  = INVALIDS))
            
            rownames(PipelineTesting) <- NULL
            
            # # Create the pairwise comparisons of treatments and ant status:
            # contrasts <- emmeans(sel_mod, pairwise ~ Treatment * Ant_status, adjust = "tukey")
            # lt <- as.data.frame(cld(contrasts, Letters = letters, decreasing = TRUE))
            # # Give the rows meaningful names:
            # mod_UN_list <- c(mod_UN_list,list(lt))
            # names(mod_UN_list)[length(mod_UN_list)] <- paste(ID_model,deparse(substitute(sel_mod)),sep = "-")
            # 
            #   #emmeans is used to calculate and display the predicted marginal means (i.e., average predicted values) 
            #   # for different combinations of the categorical predictor variables Treatment and Ant_status, in the GLM model mod1.
            #   posthocDF = data.frame(emmeans(mod1, ~ Treatment * Ant_status, type="response"))
            #   
            #   if (PLOT) {
            #   print(
            #     ggplot(posthocDF, aes(x=Ant_status, y=response, fill=Treatment)) + 
            #       geom_bar(stat="identity",position ="dodge", width=0.5) +
            #       geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0, size= .5,
            #                     position=position_dodge(.5))+
            #       colFill_Treatment +
            #       colScale_Treatment +
            #       STYLE + #theme_few() +
            #       ggtitle(paste(GENE,"\npredicted marginal means", collapse = ", "),) +
            #       theme(plot.title = element_text(size = 9)) +
            #       ylab("rel. concentration")
            #   )
            # 
            # 
            #   #if the interaction is not significant, the lines in the plot should be approximately parallel, indicating that there is no significant difference in the spread of rel_conc_imputed between Task groups across treatments
            #   ggplot(posthocDF, aes(x = Ant_status, y = response, group = Treatment, color = Treatment)) +
            #     geom_line(size = 1) +
            #     geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Treatment), alpha = 0.2) +
            #     colFill_Treatment +
            #     colScale_Treatment +
            #     STYLE + #theme_few() +
            #     ggtitle(paste(GENE, "\npredicted marginal means", collapse = ", "),) +
            #     theme(plot.title = element_text(size = 9)) +
            #     ylab("rel. concentration") +
            #     scale_x_discrete(expand = c(0, 0))
            #   }
          }
        }
      }
      
      #### PREPARE post hoc LABELS FOR THE PLOTS
      
      # FOR WORKERS  (LMM) 
      # MIGHT BE INCLUDED IN THE POSTHOC FUNCTION
      for (i in seq_along(posthoc_list)) {
        # posthoc_list[[i]]$GROUP <- sub("-.*", "", names(posthoc_list[i]))
        # posthoc_list[[i]]$variable <- sub(".*-", "", names(posthoc_list[i]))
        posthoc_list[[i]][c("GROUP", "gene","variable")] <- as.data.frame(str_split_fixed(names(posthoc_list[i]), '-', 3))
      }
      
      label_status <- bind_rows(posthoc_list[grepl("Ant_status", names(posthoc_list))], .id = "column_label")
      label_treatment <- bind_rows(posthoc_list[grepl("Treatment", names(posthoc_list))], .id = "column_label")
      
      # FOR TREATED (Mann-Whitney U test) -MERGE WITH Q ONE -
      for (i in seq_along(mod_T_list)) {
        mod_T_list[[i]][c("GROUP", "gene")] <- as.data.frame(str_split_fixed(names(mod_T_list[i]), '-', 2))
      }
      mod_T <- do.call(rbind.data.frame, mod_T_list)
      
      # FOR QUEEN (Mann-Whitney U test)
      for (i in seq_along(mod_Q_list)) {
        # posthoc_list[[i]]$GROUP <- sub("-.*", "", names(posthoc_list[i]))
        # posthoc_list[[i]]$variable <- sub(".*-", "", names(posthoc_list[i]))
        mod_Q_list[[i]][c("GROUP", "gene")] <- as.data.frame(str_split_fixed(names(mod_Q_list[i]), '-', 2))
      }
      mod_Q <- do.call(rbind.data.frame, mod_Q_list)
      
      if (PLOT) {
        ### PLOTTING ###
        for (GROUP in unique(CLEAN_DATA$GROUP)) {
          
          if (GROUP == "TREATED_W") {
            ## DO NOT ANALYSE DO NOT PLOT, NOT BASAL EXPRESSION ANYMORE
            # print(
            #   # ggplot(
            #   #   data = CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),],
            #   #   aes(x = gene, y = rel_conc_imputed, color = Treatment)
            #   # ) + 
            #   #   geom_point(position=position_jitterdodge(), size = 1, alpha = 0.3, show.legend = FALSE) +
            #   #   geom_boxplot(aes(colour = Treatment), lwd = 0.8, alpha = 0.2) +
            #   #   colScale_Treatment +
            #   #   STYLE +
            #   #   ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
            #   #   scale_y_continuous(labels = scales::percent, trans = "log10") +
            #   # #geom_text(data = label_treatment[which(label_treatment$GROUP==GROUP),], aes(x = gene, y = 10, group = Treatment, label = V1, fontface = "bold"), position = position_jitterdodge(seed = 2)) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
            #   #   geom_text(data = mod_T[which(mod_T$GROUP==GROUP),], aes(x = gene, y = 1.5, label = p, fontface = "bold") , color = "black") #, position = position_jitterdodge(seed = 2) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
            #   
            #   ggplot(data =  CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),], 
            #          aes(x = gene, y = rel_conc_imputed, color = Treatment)
            #   ) +
            #     geom_violin(aes(fill= Treatment), trim = FALSE,width =1.1, alpha = 0.6) +
            #     geom_boxplot(aes(fill= Treatment),colour ="black", width = 0.1, alpha = 0.6, position = position_dodge(width = 1.1) )+
            #     colFill_Treatment +
            #     colScale_Treatment +
            #     STYLE +
            #     ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
            #     scale_y_continuous(labels = scales::percent, trans = "log10") +
            #     geom_text(data = mod_T[which(mod_T$GROUP==GROUP),], aes(x = gene, y = 100, label = p, fontface = "bold") , color = "black") #, position = position_jitterdodge(seed = 2) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
            # )
          }else if (GROUP == "QUEEN"){
            
            if (QUEEN_data_OK) {
              # separated from the above as the label is assigned differently!
              print(
                # ggplot(
                #   data = CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),],
                #   aes(x = gene, y = rel_conc_imputed, color = Treatment)
                # ) + 
                #   geom_point(position=position_jitterdodge(), size = 1, alpha = 0.3, show.legend = FALSE) +
                #   geom_boxplot(aes(colour = Treatment), lwd = 0.8, alpha = 0.2) +
                #   colScale_Treatment +
                #   STYLE +
                #   ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
                #   scale_y_continuous(labels = scales::percent, trans = "log10") +
                # #geom_text(data = label_treatment[which(label_treatment$GROUP==GROUP),], aes(x = gene, y = 10, group = Treatment, label = V1, fontface = "bold"), position = position_jitterdodge(seed = 2)) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
                #   geom_text(data = mod_Q[which(mod_Q$GROUP==GROUP),], aes(x = gene, y = 1, label = p, fontface = "bold") , color = "black") #, position = position_jitterdodge(seed = 2) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
                
                ggplot(data =  CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),], 
                       aes(x = gene, y = rel_conc_imputed, color = Treatment)
                ) +
                  #geom_violin(aes(fill= Treatment), trim = FALSE,width =1.1, alpha = 0.6) +
                  #geom_boxplot(aes(fill= Treatment),colour ="black", width = 0.1, alpha = 0.6, position = position_dodge(width = 1.1) )+
                  geom_boxplot(aes(fill = Treatment),colour ="black", alpha = 0) +
                  geom_point(position=position_jitterdodge(),aes(fill = "black"), size = 1, alpha = 0.8, show.legend = FALSE) +
                  colFill_Treatment +
                  colScale_Treatment +
                  STYLE +
                  ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
                  scale_y_continuous(labels = scales::percent, trans = "log10") +
                  geom_text(data = mod_Q[which(mod_Q$GROUP==GROUP),], aes(x = gene, y = 1, label = p, fontface = "bold") , color = "black") #, position = position_jitterdodge(seed = 2) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
              )
              
            }
          }else if (GROUP == "UNTREATED_W"){
            # print(
            # ggplot(data =  CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),],
            #        aes(x = Ant_status, y = rel_conc_imputed, color = Treatment)) +
            #   geom_violin(aes(fill= Treatment), trim = FALSE,width =1.1, alpha = 0.6) +
            #   geom_boxplot(aes(fill= Treatment),colour ="black", width = 0.1, alpha = 0.6, position = position_dodge(width = 1.1) )+
            #   colFill_Treatment +
            #   colScale_Treatment +
            #   STYLE +
            #   ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
            #   scale_y_continuous(labels = scales::percent, trans = "log10") +
            #   facet_grid(. ~ gene) 
            #   #   geom_text(data = label_treatment[which(label_treatment$GROUP==GROUP),], aes(x = gene, y = 10, group = Treatment, label = V1, fontface = "bold"), position = position_jitterdodge(seed = 2)) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
            # )  
            
            # PLOT FOR NON-TREATED ANTS
            ## plot by Treatment (size)
            
            p1 <- NULL
            p2 <- NULL
            
            # p1 <- ggplot(
            #   data = CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),],
            #   aes(x = gene, y = rel_conc_imputed, color = Treatment)
            # ) +
            #   geom_violin(aes(fill= Treatment), trim = TRUE,width =1.1, alpha = 0.6) +
            #   geom_boxplot(aes(fill= Treatment),colour ="black", width = 0.1, alpha = 0.6, position = position_dodge(width = 1.1) )+
            #   #geom_point(position=position_jitterdodge(), size = 1, alpha = 0.3, show.legend = FALSE) +
            #   #geom_boxplot(aes(colour = Treatment), lwd = 0.8, alpha = 0.2) +
            #   colFill_Treatment +
            #   colScale_Treatment +
            #   STYLE +
            #   ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
            #   scale_y_continuous(labels = scales::percent, trans = "log10") #+
            # #geom_text(data = label_treatment[which(label_treatment$GROUP==GROUP),], aes(x = gene, y = 10, group = Treatment, label = V1, fontface = "bold"), position = position_jitterdodge(seed = 2)) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
            # 
            
            # PLOT FOR NON-TREATED ANTS
            # ## plot by Treatment (size)
            # p2 <- ggplot(
            #   data = CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),],
            #   aes(x = Ant_status, y = rel_conc_imputed, color = treatment)
            # ) +
            #   #geom_point(aes(colour = Colony),position=position_jitterdodge(), size = 1, alpha = 0.3) +
            #   geom_violin(trim = TRUE,width =1.1, alpha = 0.6) +
            #   geom_boxplot(colour ="black", width = 0.1, alpha = 0.6, position = position_dodge(width = 1.1) )+
            #   colFill_Treatment +
            #   colScale_Treatment +
            #   STYLE +
            #   ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
            #   scale_y_continuous(labels = scales::percent, trans = "log10") +
            #   facet_grid(. ~ gene) #+
            # #geom_text(data = label_status[which(label_status$GROUP==GROUP),], aes(x = Ant_status, y = 10, group = Ant_status, label = V1, fontface = "bold"), position = position_jitterdodge(seed = 2)) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
    
            # PLOT FOR NON-TREATED ANTS
            # ## plot by Treatment (size)
            
            #base model params
            base_model <- PipelineTesting[which(PipelineTesting$T.E.==3 & PipelineTesting$LOD==37 & PipelineTesting$imputation=="HM"& PipelineTesting$Invalids=="Invalids_DISCARD" &PipelineTesting$ALWAYS_DISCARD==T),]
            base_model$P.Status.stars <- sapply(base_model$P.Status, add_star)
            base_model$P.Treatment.stars <- sapply(base_model$P.Treatment, add_star)
            
            
            p1 <- ggplot(
                  data = CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),],
                  aes(x = gene, y = rel_conc_imputed, color = Treatment)
                ) +
                  #geom_violin(aes(fill= Treatment), trim = TRUE,width =1.1, alpha = 0.6) +
                  #geom_boxplot(aes(fill= Treatment),colour ="black", width = 0.1, alpha = 0.6, position = position_dodge(width = 1.1) )+
                  #geom_point(position=position_jitterdodge(), size = 1, alpha = 0.3, show.legend = FALSE) +
                  #geom_boxplot(aes(colour = Treatment), lwd = 0.8, alpha = 0.2) +
              #geom_violin(trim = TRUE,width =1.1, alpha = 0.6) +
              stat_summary(
                fun.data = mean_cl_normal, # The function used to calculate the confidence intervals
                geom = "errorbar",
                width = 0.4,
                size = 1,
                aes(color = Treatment), # Color the errorbar based on treatment
                position = position_dodge(width = 0.4)
                ) +
              stat_summary( # To add a point for the mean
                fun = mean,
                geom = "point",
                size = 1.2, # Customize the size
                aes(color = Treatment, group = Treatment), # Map group aesthetic to treatment
                position = position_dodge(width = 0.4)  
                ) +
                  colFill_Treatment +
                  colScale_Treatment +
                  STYLE +
                  ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) +
                  scale_y_continuous(labels = scales::percent, trans = "log10") +
                #geom_text(data = label_treatment[which(label_treatment$GROUP==GROUP),], aes(x = gene, y = 10, group = Treatment, label = V1, fontface = "bold"), position = position_jitterdodge(seed = 2)) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
            geom_text(data = base_model, aes(x = gene, y = 0.35, label = P.Treatment.stars, fontface = "bold") , color = "black") #, position = position_jitterdodge(seed = 2) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
            
            
            
            # PLOT FOR NON-TREATED ANTS
            # ## plot by Ant_status
            p2 <- ggplot(
              data = CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),],
              aes(x = Ant_status, y = rel_conc_imputed)
            ) +
              #geom_violin(trim = TRUE,width =1.1, alpha = 0.6) +
              stat_summary(
                fun.data = mean_cl_normal, # The function used to calculate the confidence intervals
                geom = "errorbar",
                width = 0.4,
                size = 1
              ) +
              stat_summary( # To add a point for the mean
                fun = mean,
                geom = "point",
                size = 1.2 # Customize the size
              ) +
              colFill_Treatment +
              colScale_Treatment +
              STYLE +
              ggtitle(paste(unique(CLEAN_DATA[which(CLEAN_DATA$GROUP==GROUP),"Ant_status"]), collapse = ", ")) + #,subtitle = "95% C.I."
              scale_y_continuous(labels = scales::percent, trans = "log10") +
              facet_grid(. ~ gene)+
              geom_text(data = base_model, aes(x = 1.5, y = 0.35, label = P.Status.stars, fontface = "bold") , color = "black") #, position = position_jitterdodge(seed = 2) # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
            

            # arrange plots side by side
            grid.arrange(p1, p2, ncol = 2)
            
          }
        }
      } #PLOT?
      
      print(mod_UN_list)
      
      
      #}) ### LOOP BETWEEN THE TWO VARIANTS OF HOW TO TREAT THE INVALIDS (KEEP MEAN OR DISCARD VALUES)
      
      
    } # INVALIDS discard or keep mean loop
    
  }
} # IMPUTATION HM or QRILC
} # Technical_error
} # Limit_of_Detection

warning("CHECK ALL THE STEPS, REPORT % DATAPOINTS PRESENT, UPDATE SLIDES, PRODUCE PLOTS FOR BOTH CONDITIONS - MODIFY PLOT SUBTITLE TO SPECIFY DATASET USED-")

#base model params (WARNING: THIS IS REPEATED IN UNTREATED_w PLOT!)
base_model <- PipelineTesting[which(PipelineTesting$T.E.==3 & PipelineTesting$LOD==37 & PipelineTesting$imputation=="HM"& PipelineTesting$Invalids=="Invalids_DISCARD" &PipelineTesting$ALWAYS_DISCARD==T),]






# Create scatterplot of model estimates VS p-values
sp <- ggplot(data = PipelineTesting, aes_string(x = "Estim.Status", y = "P.Status")) +
  geom_point(aes(color = gene, shape=imputation), alpha = 0.5,size =2,position = position_jitter(width = 0.05, height = 0))  + # CAREFUL: there is a tiny amount of jitter for readability
  geom_point(data =  base_model, aes_string(x = "Estim.Status", y = "P.Status"), shape = 4, size = 3) +
  #geom_rug(aes(color = gene))+
  geom_text(data = base_model, aes_string(x = "Estim.Status", y = "P.Status", label = "gene"), vjust = -2) +
  labs(x = "Model Estimates", y = "p-value") + #, title = "Model Estimates vs. P-values"
  theme_classic()  +
  scale_shape_manual(values = c(15, 17)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")+
  scale_color_jco() +
  theme(aspect.ratio=1) +
  border()  +
  theme(plot.margin = margin(t = -20, r = -20, b = 0, l = 0, unit = "pt"))

  
#theme(legend.position = "bottom")


## Scatter plot colored by groups ("Species")
# sp <- ggscatter(iris, x = "Sepal.Length", y = "Sepal.Width",
#                 color = "Species", palette = "jco",
#                 size = 3, alpha = 0.6)+

# Marginal density plot of x (top panel) and y (right panel)
xplot <- gghistogram(PipelineTesting, "Estim.Status", fill = "gene",color = "transparent",
                   palette = "jco")  +
  theme(plot.margin = margin(t = 0, r = 0, b = -20, l = 0, unit = "pt")) + 
  clean_theme() +
  rremove("legend")
yplot <- gghistogram(PipelineTesting, "P.Status", fill = "gene", color = "transparent",
                   palette = "jco")  +
  theme(plot.margin = margin(t = 0, r = 0, b = -20, l = 0, unit = "pt")) + 
  rotate()+ 
  clean_theme() +
  rremove("legend")
# Cleaning the plots
sp <- sp + rremove("legend")
# Arranging the plot using cowplot
cowplot::plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
                   rel_widths = c(2, 1), rel_heights = c(1, 2)
                   )

### ADD: KEEP LEGEND ASIDE AND READD IT (AS DONE FOR THE QRILC PLOTS)








plot_list <- list()

#   plot_list[[VARIABLE]]

for (VARIABLE in c("P.Status","Estim.Status","P.Treatment","Estim.Treatment")) { 
  
  # plot the relative concentration.
  plot_list[[VARIABLE]] <-  ggplot(PipelineTesting, aes(x = !!sym(VARIABLE), fill = gene)) +
      geom_histogram(alpha = 0.5, bins = 50, position = "identity") +
      scale_fill_manual(values = c("red", "green", "blue")) +
      ggtitle(paste("Values for'",VARIABLE,"' | Pipeline Testing (",nrow(PipelineTesting)/3,"combs)",sep=" "),
              subtitle="arrow for base model") +
      theme_minimal() 
      # Add vertical line if VAR contains a word starting with "P"
      if (length(grep("^P", VARIABLE)) > 0) {
        plot_list[[VARIABLE]] <- plot_list[[VARIABLE]] + geom_vline(xintercept = 0.05, color = "red")
      }
  
  # Add arrows for reference estimates
  for (i in 1:length(base_model)) {
    plot_list[[VARIABLE]] <- plot_list[[VARIABLE]] + annotate(geom = "segment", x = base_model[i,VARIABLE], y = 0, 
                      xend = base_model[i,VARIABLE], yend = 10, 
                      arrow = arrow(length = unit(0.5, "cm")))
    
    plot_list[[VARIABLE]] <- plot_list[[VARIABLE]] + geom_text(x =  base_model[i,VARIABLE], y = 10, 
                       label = base_model[i,"gene"], hjust = -1, angle = 90)
  }
  #print(p)
}

grid.arrange(
  plot_list[[1]], plot_list[[2]],
  plot_list[[3]], plot_list[[4]],
  ncol = 2, nrow = 2,
  widths = c(1, 1),
  heights = c(1, 1)
  #bottom = paste("QRILC imputation, stacked histogram. Data=", NAME_DF)
)



as.data.frame(
PipelineTesting %>%
  group_by(gene) %>%
  summarize(Sig.Status.freq = mean(P.Status < 0.05)*100,Sig.Treat.freq = mean(P.Treatment < 0.05)*100) %>%
  mutate(Sig.Status.freq = paste0(round(Sig.Status.freq, 1), "%"),Sig.Treat.freq = paste0(round(Sig.Treat.freq, 1), "%"))
)


# add the same for the other treatments! and do separate plots
# add column "what"



##### SCRAPS #################################


# # outliers <- genes_data %>%
# #   inner_join(genes_data_byGroup, by = "gene") %>%
# #   filter(mean_Ct < grand_mean_Ct - 2 * st.dev_mean_Ct | mean_Ct > grand_mean_Ct + 2 * st.dev_mean_Ct) %>%
# #   filter(!(AntTask == "queen" & mean_Ct < grand_mean_Ct))
# 
# 
# # # CHECK EF1 with mean_Ct value greater than 2*st.dev mean_Ct value for each group
# # #exclude Queens from outliers (on the lowerbound, they are expected to have a lot of datapoints -low Ct-, but the upperbound should be cleaned)
# # # outliers are not removed or filtered at this stage
# # outliers <- genes_data %>%
# #   inner_join(genes_data_byGroup, by = "gene") %>%
# #   filter(mean_Ct < grand_mean_Ct - 2 * st.dev_mean_Ct | mean_Ct > grand_mean_Ct + 2 * st.dev_mean_Ct) %>%
# #   filter(!(AntTask == "queen" & mean_Ct < grand_mean_Ct))
# # 
# # # plot a scatterplot with lines connecting Ct values between groups
# # ggplot(outliers, aes(x = gene, y = mean_Ct, group = Code)) +
# #   geom_point(aes(color = gene)) +
# #   geom_line(alpha = 0.1) +
# #   ggtitle("Scatterplot of mean_Ct values by group gene with lines connecting values") +
# #   xlab("Group") +
# #   ylab("mean_Ct value") #+
# # # scale_color_discrete(name = "Group", labels = unique(outliers$gene))
# 
# 
# 
# # # check outliers min upper threshold
# outliers_limit <- min(outliers %>%
#                         filter(gene == "EF1") %>%
#                         filter(mean_Ct > grand_mean_Ct + 2 * st.dev_mean_Ct) %>%
#                         pull(mean_Ct))
# #outliers_limit <-  mean(genes_data$mean_Ct,na.rm=T)+2*sd(genes_data$mean_Ct,na.rm=T)
# 
# 
# ## we then discard samples with EF1 mean_ct > 32 as they are probably low quality samples (bad extraction, bad crushing, etc)
# ## filter rows with mean_Ct value greater than 32 for group EF1
# #genes_data <- genes_data[!(genes_data$gene == "EF1") | !(genes_data$mean_Ct >= outliers_limit),]
# genes_data[which(genes_data$gene == "EF1" & genes_data$mean_Ct >= outliers_limit),"mean_Ct"] <- NA
# 
# #####################################
# 
# 
# 
# #base vals
# sum(is.na(genes_data$mean_Ct))
# genes_data$NA_Ct <- ifelse(is.na(genes_data$mean_Ct),NA,"values")
# table(genes_data$NA_Ct,genes_data$gene, useNA = "ifany")
# #Accepted Ct diff (0= non-acceptable, 1= acceptable)
# table(genes_data$CTDiffAcc_15PipErr,genes_data$gene, useNA = "ifany")
# # assign accepted to missing mean_Ct
# genes_data[which(genes_data$abs_diff_Ct==0),"CTDiffAcc_15PipErr"] <- 1
# 
# 
# ##### HOW TO DEAL WITH ABS DIFF = NA? SHOULD THESE BE IMPUTED TOO? OR WE KEEP THEM AS POSITIVE EVEN IF WE DON'T HAVE THEIR COUNTERPART? THESE COULD BE A LOT OF POINTS
# # % of saved datapoints (1) with threshold at 0.5 # see plot too
# genes_data$CTDiffAcc_0.5 <- 0
# genes_data$CTDiffAcc_0.5 <- ifelse((genes_data$abs_diff_Ct<=0.5),1,genes_data$CTDiffAcc_0.5)
# table(genes_data$CTDiffAcc_0.5,useNA = "ifany")/nrow(genes_data)
# 
# # N of saved datapoints (1) with poisson adjusted threshold
# table(genes_data$CTDiffAcc_15PipErr,useNA = "ifany")/nrow(genes_data) ## MORE POINTS RETAINED
# 
# #########Step 3: If mean Ct>critical value (35), recode as nondetect##########
# ##############################################################################
# #If Ct>=35 THEN recode as nondetect (0)
# #If Ct=0 in both duplo's, also recode as nondetect (0).
# genes_data$CtBigger35 <- 1
# genes_data$CtBigger35 <- ifelse((genes_data$mean_Ct>=35 ),0,1) #| genes_data$mean_Ct <1
# table(genes_data$CtBigger35, genes_data$gene, useNA="ifany")
# #if the mean is missing, reassing default value
# genes_data[is.na(genes_data$CtBigger35),"CtBigger35"] <- 1
# #########Step 4: Make final variable##########
# ##############################################
# 
# ###########Ronde et al. 2017 pipeline################
# ## THEY SUGGEST TO USE MULTIPLE IMPUTATION
# # QUALITY CONTROL
# #0=impute (invalid) - to be imputed
# #1=valid
# #2=undetectable - set to minimum
# genes_data$QC <- 0
# genes_data$QC <- ifelse((genes_data$CTDiffAcc_15PipErr==1 & genes_data$CtBigger35==1),1,genes_data$QC) # if the diff between reps is acceptable & lower than 35, it is valid
# genes_data$QC <- ifelse((genes_data$CTDiffAcc_15PipErr==0 & genes_data$CtBigger35==1),0,genes_data$QC) # if the diff between reps is unacceptable and lower than 35, invalid impute
# genes_data$QC <- ifelse((genes_data$CTDiffAcc_15PipErr==1 & genes_data$CtBigger35==0),2,genes_data$QC) # if the diff between reps is acceptable and higher than 35, undetectable, set to min
# # OR: genes_data$QC <- ifelse((genes_data$CTDiffAcc_15PipErr==1 & genes_data$CtBigger35==0),1,genes_data$QC)
# genes_data$QC <- ifelse((genes_data$CTDiffAcc_15PipErr==0 & genes_data$CtBigger35==0),2,genes_data$QC) # if the diff between reps is unacceptable and higher than 35, undetectable, set to min
# #genes_data$QC <- ifelse((is.na(genes_data$mean_Ct)),0,genes_data$QC) # if mean_ct is missing, impute value (not in Ronde et al. 2017, added by me)
# 
# table(genes_data$QC,useNA = "ifany")
# 
# table(genes_data$QC,genes_data$gene, useNA="ifany")
# 






# 
# #### assign min to undetectables
# #0=impute (invalid) - to be imputed
# #1=valid
# #2=undetectable - set to minimum
# for (GENE in unique(genes_data$gene)) {
#   #assign the min/2 value by gene to missing datapoints
#   genes_data[which(genes_data$QC == "2" & genes_data$gene == GENE), "rel_conc_imputed"] <- min(genes_data[which(genes_data$gene==GENE),"rel_conc_replace_min"],na.rm = T)/2
# }
# 



# Create a minimal reproducible sample

min_sample <- GENE_data %>%
  group_by(Ant_status, Treatment) %>%
  slice_sample(n = 1) %>%
  ungroup()
dput(as.data.frame(min_sample[,c("Ant_status", "Treatment", "rel_conc_imputed","Colony")]))
