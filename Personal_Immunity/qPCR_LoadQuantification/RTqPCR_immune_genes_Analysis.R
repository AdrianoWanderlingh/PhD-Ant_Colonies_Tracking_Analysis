##################################################################################################
################################## RT qPCR IMMUNE RELATED GENES ##################################
##################################################################################################

# load the packages
library(tidyverse)
library(platetools) # to plot plate
library(naniar) # show missing points
library(reshape2) # for melt
library(scales) # for colours

## for stats
library(performance)
library(remotes)
#install_version("MuMIn", "1.46.0")
library(MuMIn) # multi-model inference
library(multcomp)
library(multcompView)
library(blmeco)
library(sjPlot)
library(lsmeans)
library(lme4)
library(blmeco) # check dispersion for glmer
library(emmeans) #post-hoc comparisons
library(e1071) #calc skewness and other stuff
library(lawstat) #for levene test (homogeneity of variance)
library(lmPerm) # for permutations
library(lmerTest)
library(car)
library(nlme) #lme
library(afex)

# #Create a custom color scale FOR COLONIES + treatments
# FullPal <- scales::viridis_pal(option = "D")(20)
# myColorsSmall  <- tail(FullPal,5)
# myColorsLarge  <- head(FullPal,5) 
# Cols_treatments <- c("#440154FF","#FDE725FF") #"#31688EFF"
# myColors      <- c(myColorsLarge,myColorsSmall, Cols_treatments)
# names(myColors) <- c("R3BP","R5BP","R7BP","R8BP","R12BP","R1SP", "R2SP", "R5SP", "R7SP","R11SP","Big Pathogen","Small Pathogen")
# colScale <- scale_colour_manual(name = "Colony",values = myColors,drop=TRUE)
# fillScale <- scale_fill_manual(name = "Colony",values = myColors,drop=TRUE)

######## UPDATED VERSION!
#Create a custom color scale FOR COLONIES + treatments
show_col(scales::viridis_pal(option = "D")(4))
# myColorsSmall  <- tail(FullPal,5)
# myColorsLarge  <- head(FullPal,5)
Cols_treatments <- c("#440154FF","#31688EFF","#35B779FF","#FDE725FF") #"#31688EFF"
myColors      <- c(Cols_treatments)
names(myColors) <- c("Big Pathogen","Big Sham","Small Pathogen","Small Sham")

colScale <- scale_colour_manual(name = "Colony",values = myColors,drop=TRUE)
fillScale <- scale_fill_manual(name = "Colony",values = myColors,drop=TRUE)




# ggplot PLOT STYLE
STYLE <- list(colScale, fillScale,
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
              theme_bw(),
              scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

STYLE_CONT <- list(colScale, fillScale,
                   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
                   theme_bw()
)

#function to test normality of residuals
test_norm <- function(resids){
  print("Testing normality")
  if (length(resids)<=300){
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
    print("below 0.05, the data significantly deviate from a normal distribution")
  }else{
    print("More than 300 data points so using the skewness and kurtosis
approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}

#function to report a model output
output_lmer <- function(model){
  print("------------RESIDUALS NORMALITY------------")
  test_norm(residuals(model))
  print("------------SUMMARY------------")
  print(summary(model))
  print("------------ANOVA------------")
  print(Anova(model))
  print("------------RSQUARED------------")
  print(r.squaredGLMM(model))
  tab_model(model)
}


# function to perform posthocs
posthoc_list <- list()
compute_posthocs <- function(model){
  #check taht there are significant outputs
  if (length(row.names(Anova(model)[Anova(model)$'Pr(>Chisq)' < 0.05,]))==0){
    print("there are no significant vars.")}else{
      for (SIG.VAR in row.names(Anova(model)[Anova(model)$'Pr(>Chisq)' < 0.05,])) {
        if (grepl(":", SIG.VAR)) {
          warning(paste0(SIG.VAR, "is an interaction, currently this situation is not handled by the function."))
        }else{
          #check if the variable is not numeric . to do so, we need to access the dataframe from the model
          if (!is.numeric(get(gsub("\\[.*","",as.character(model@call)[3]))[,SIG.VAR])) {
            print(paste0("Performing posthocs for the significant var: ",SIG.VAR))
            arg <- list("Tukey")
            names(arg) <- SIG.VAR
            # glht (mcp) : General linear hypotheses Testing (glht) and multiple comparisons for parametric models (mcp)
            cmp <- do.call(mcp, arg)
            posthoc_SIG.VAR <- summary(glht(model, linfct = cmp),test=adjusted("BH"))
            # Set up a compact letter display of all pair-wise comparisons
            model_means_cld <- cld(posthoc_SIG.VAR)
            #create dataframe usable with ggplot geom_text
            model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters)))
            # add column name
            model_means_cld$newcol <- NA
            colnames(model_means_cld)[which(names(model_means_cld) == "newcol")] <- SIG.VAR
            model_means_cld[,SIG.VAR] <- row.names(model_means_cld)
            rownames(model_means_cld) <- NULL
            # add to list
            posthoc_list <- c(posthoc_list, list(model_means_cld))
            if (exists("ID_model")) {
              names(posthoc_list)[length(posthoc_list)] <- paste(ID_model,deparse(substitute(model)),SIG.VAR,sep="-")
            }else{
              names(posthoc_list)[length(posthoc_list)] <- paste(deparse(substitute(model)),SIG.VAR,sep="-")
            }
            print(paste(deparse(substitute(model)),SIG.VAR,sep="_"))
            print(model_means_cld)
          } # SIG.VAR LOOP
        } #check if is an interaction
      } #check if numeric
      print("call posthoc_list to get posthocs")
    } # if significant vars exist
  return(posthoc_list)
}

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
  # genes_data <- genes_data[which(!is.na(genes_data$AntTask)),]
  
  names(genes_Results_annotated)[which(names(genes_Results_annotated)=="treatment")] <- "Treatment"
  # Rename by name
  genes_Results_annotated$Treatment <- as.factor(genes_Results_annotated$Treatment)
  levels(genes_Results_annotated$Treatment)[levels(genes_Results_annotated$Treatment)=="BS"] <- "Big Sham"
  levels(genes_Results_annotated$Treatment)[levels(genes_Results_annotated$Treatment)=="SS"] <- "Small Sham"
  
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

##remove rows which contain NAs for Ct
genes_data <- genes_data[!is.na(genes_data$mean_Ct), ]

# Relevel Exposed
genes_data$Exposed <- as.factor(genes_data$Exposed)
levels(genes_data$Exposed)[levels(genes_data$Exposed)=="TRUE"] <- "treated"
levels(genes_data$Exposed)[levels(genes_data$Exposed)=="FALSE"] <- "untreated"
# Create new Status category
genes_data$Ant_status <- paste(genes_data$Exposed, genes_data$AntTask)
#remove non-categorised ants
genes_data <- genes_data[which(genes_data$IsAlive==TRUE),]

# save controls aside
controls_gene_data <- genes_data %>%
  filter(is.na(Colony))
# keep housekeeping and test genes data
genes_data <- genes_data %>%
  filter(!is.na(Colony))

# plot RAW data
ggplot(genes_data, aes(x = gene, y = mean_Ct, color = Sample_Plate)) +
  geom_boxplot(aes(colour = Sample_Plate), lwd = 0.8, alpha = 0.3)


##### ISOLATE SUSPICIOUS DATAPOINTS (CONTAMINATIONS, BAD EXTRACTIONS, ETC)
# Calculate standard deviation by group
stdev_by_group <- genes_data %>%
  group_by(gene) %>%
  summarize(st.dev_mean_Ct = sd(mean_Ct,na.rm = T),grand_mean_Ct = mean(mean_Ct,na.rm = T))

# filter rows with mean_Ct value greater than 2*st.dev mean_Ct value for each group
outliers <- genes_data %>%
  inner_join(stdev_by_group, by = "gene") %>%
  filter(mean_Ct < grand_mean_Ct - 2*st.dev_mean_Ct | mean_Ct > grand_mean_Ct + 2*st.dev_mean_Ct)

# plot a scatterplot with lines connecting Ct values between groups
ggplot(outliers, aes(x = gene, y = mean_Ct, group = Code)) +
  geom_point(aes(color = gene)) +
  geom_line(alpha=0.1) +
  ggtitle("Scatterplot of mean_Ct values by group gene with lines connecting values") +
  xlab("Group") +
  ylab("mean_Ct value") #+
  #scale_color_discrete(name = "Group", labels = unique(outliers$gene))

##### DELTA-Ct METHOD
# “delta-Ct”: difference between the Ct of the housekeeping gene and the test gene
# split test genes from housekeeping
test_gene_data <- genes_data %>%
  filter(gene != "EF1")
# isolate housekeeping and change colnames
ref_gene_data <- genes_data %>%
  filter(gene == "EF1") %>%
  rename("housekeeping" = "gene", "ref_Ct" = "mean_Ct", "ref_Tm" = "mean_Tm")

# recombine data
common_col_names4 <- intersect(names(ref_gene_data), names(test_gene_data))
genes_data <- left_join(test_gene_data, ref_gene_data, by = common_col_names4)

# create a new column containing the delta Ct between the housekeeping gene and our gene of interest, and plot the delta Ct for each treatment and replicate.
genes_data <- mutate(genes_data, delta_Ct = mean_Ct - ref_Ct)

# plot delta_Ct
ggplot(genes_data, aes(x = gene, y = delta_Ct,fill=Treatment)) +
  geom_boxplot(aes(colour = AntTask), lwd = 0.8, alpha = 0.3)
# 
# # Calculate the mean delta Ct for each treatment.
# treatment_summary <- genes_data %>%
#   group_by(gene) %>%
#   summarise(mean_delta_Ct = mean(delta_Ct, na.rm = TRUE))


##### REMOVE HOUSEKEEPING GENE'S EXPRESSION OVER THRESHOLD
# EF1alpha is expected () to show expression values ranging between 25+-2 (CHECK WITH FLORENT ACCORDING TO HIS OBSERVATIONS).
# given that we work with single ants, the closest we are to threshold, the harder it is to quantify small samples (small individuals).
# check outliers min upper threshold
outliers_limit <- min(outliers %>%
  filter(gene == "EF1") %>%
  filter(mean_Ct > grand_mean_Ct + 2*st.dev_mean_Ct) %>%
  pull(mean_Ct))

# we then discard samples with EF1 mean_ct > 32 as they are probably low quality samples (bad extraction, bad crushing, etc)
# filter rows with mean_Ct value greater than 32 for group EF1
genes_data <- genes_data %>%
  filter(ref_Ct < outliers_limit)

# Calculating relative DNA concentration
# to calculate the relative DNA concentration, you can use the fact that the amount of cDNA theoretically doubles every cycle.
#PRIMERS EFFICIENCY 
primers_eff <- data.frame(gene = c("EF1","HYM","DEF", "PO"), efficiency = c(1.99 ,1.99 ,2, 1.94))
# add primers efficiency data
genes_data <- left_join(genes_data, primers_eff, by = c("gene"))

# as EF1 efficiency is 1.99, it is used as baseline for normalisation
genes_data <- genes_data %>%
  mutate(rel_conc = (2*(efficiency/1.99))^-delta_Ct)

# # plot the relative concentration.
# ggplot(genes_data, aes(x = gene, y = rel_conc,fill=Treatment)) +
#   geom_boxplot(aes(colour = AntTask), lwd = 0.8, alpha = 0.3) +
#   scale_y_continuous(labels = scales::percent,trans='log10')

ggplot(data= genes_data,
  aes(x = AntTask, y = rel_conc)) + #,group = Exposed,color = Exposed
  geom_point(aes(fill = Treatment,colour=Colony),position = position_jitterdodge(),size=1,alpha=0.3)+ #,alpha= 0.8,stroke=0
  geom_boxplot(aes(colour=Treatment),lwd=0.8,alpha = 0.3) + 
  STYLE_CONT +
  scale_y_continuous(labels = scales::percent,trans='log10') +
  facet_grid(. ~ gene) # 

### How to treat Queens
# Queens' larger amount of tissue makes the Ct values for EF1 much lower (high expression). 
# this is not a problem as each sample housekeeping gene EF1 Ct value is used to normalise the test genes, but: 
# as queens have different metabolic rates (whic are known to impact EF1 expression), the Ct values for EF1 will not be comparable, 
# therefore QUEENS HAVE TO BE ANALYSED SEPARATELY

# Queen data
genes_data_Q <- genes_data %>%
  filter(IsQueen == TRUE)
# Queen data
genes_data_EXP <- genes_data %>%
  filter(Exposed == "treated")
# unexposed data
genes_data <- genes_data %>%
  filter(IsQueen == FALSE & Exposed == "untreated")

#list for ease of plotting
genes_data_list <- list(genes_data_Q=genes_data_Q, genes_data_EXP= genes_data_EXP,genes_data=genes_data)

for (DF in genes_data_list) {
print(
ggplot(data= DF,
       aes(x = gene, y = rel_conc)) + #,group = Exposed,color = Exposed
  geom_point(aes(fill = Treatment,colour=Colony),position = position_jitterdodge(),size=1,alpha=0.3)+ #,alpha= 0.8,stroke=0
  geom_boxplot(aes(colour=Treatment),lwd=0.8,alpha = 0.3) + #
  STYLE +
  ggtitle(paste(unique(DF$Ant_status),collapse=", ")) +
  scale_y_continuous(labels = scales::percent,trans='log10') #+
  #facet_grid(. ~ gene)
)
  }

#####################################################################################
##################              STATS & PLOTS        ################################
#####################################################################################

#####################################################################################
##### TEST DIFFERENCE OF MbruDNA.LOAD BETWEEN SIZES FOR EACH Ant_status (ANTTASK*TREATMENT CONDITION)

for (GENE in unique(genes_data$gene)) {
  #select subset for 1 gene
  GENE_data <- genes_data[which(genes_data$gene==GENE),]
  #create constant to make values below zeros a small number
  GENE_cost <- min ( GENE_data$rel_conc ,na.rm=T ) * 1.1
  
  #hist(log10(GENE_data$rel_conc+GENE_cost))
  
  # #-----------------
  # m.loadVSdistanceQ_BIN <- lmer(log10(MbruDNA+constant) ~ Treatment * mean_distance_to_queen_ordered * Ant_status + (1|Colony), data = DNA_Results_IndMean_NoQ_post) # the "/" is for the nesting #  + (1|time_of_day) 
  # output_lmer(m.loadVSdistanceQ_BIN)
  # 
  # 
  # m.loadVSdistanceQ_BIN <- lmer(log10(MbruDNA+constant) ~  Treatment 
  #                               + mean_distance_to_queen_ordered 
  #                               + Ant_status 
  #                               + Treatment:mean_distance_to_queen_ordered
  #                               + Treatment:Ant_status
  #                               + Ant_status:mean_distance_to_queen_ordered
  #                               + Treatment:mean_distance_to_queen_ordered:Ant_status )
  # #-----------------
  
  #First, fit 4 candidate linear models to explain variation in density
  mod1<-lmer(log10(rel_conc+GENE_cost) ~ Treatment * Ant_status + (1|Colony), data = GENE_data)
  mod2<-lmer(log10(rel_conc+GENE_cost) ~ Treatment + Ant_status + (1|Colony), data = GENE_data)
  
  print(GENE)
  output_lmer(mod1)
  output_lmer(mod2)
  #We can now use the mod.sel to conduct model selection. The default model selection criteria is Akaike’s information criteria (AIC) with small sample bias adjustment, AICc
  out.put<-model.sel(mod1,mod2)
  # delta AICc, and the model weights
  out.put
  
  #models sorted from best (top) to worst (bottom). higher weight is better as it shows that the model is more likely to be the best explanation (hypothesis) for variation in rel_conc 
  # cclean output
  sel.table<-round(as.data.frame(out.put)[-c(1:5)],3)
  # number of parameters (df) should be K
  names(sel.table)[1] = "K"
  sel.table$Model<-rownames(sel.table)
  rownames(sel.table) <- NULL
  # replace Model name with formulas
  for(i in 1:nrow(sel.table)) sel.table$formula[i]<- as.character(formula(get(sel.table$Model[i])))[3]
  print(GENE)
  print(sel.table)
  
  #compute posthocs
  selected_model <- get(sel.table[which.max(sel.table$weight),"Model"])
  ID_model <- GENE # to assign name to posthoc's list element
  posthoc_list <- compute_posthocs(selected_model)
}

#add information to post_hoc list
# MIGHT BE INCLUDED IN THE POSTHOC FUNCTION
for (i in seq_along(posthoc_list)) {
  posthoc_list[[i]]$gene <- sub("-.*", "", names(posthoc_list[i]))
  posthoc_list[[i]]$variable <- sub(".*-", "", names(posthoc_list[i]))
}


label_status <- bind_rows(posthoc_list[grepl("Ant_status", names(posthoc_list))], .id = "column_label")
label_treatment <- bind_rows(posthoc_list[grepl("Treatment", names(posthoc_list))], .id = "column_label")


## plot by Treatment (size)
    ggplot(data= genes_data,
           aes(x = gene, y = rel_conc,colour=Treatment)) + #,group = Exposed,color = Exposed
      geom_point(aes(fill = Treatment,colour=Colony),position = position_jitterdodge(),size=1,alpha=0.3)+ #,alpha= 0.8,stroke=0
      geom_boxplot(lwd=0.8,alpha = 0.3) + #
      STYLE +
      scale_y_continuous(labels = scales::percent,trans='log10') +
      geom_text(data = label_treatment, aes(x = gene , y = 10, group=Treatment,label = V1,fontface = "bold"),position = position_jitterdodge(seed = 2))# jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
    #facet_grid(. ~ gene)
    
    
## plot by Treatment (size)
    ggplot(data= genes_data,
           aes(x = Ant_status, y = rel_conc,colour=Ant_status)) + #,group = Exposed,color = Exposed
      geom_point(aes(colour=Colony),position = position_jitterdodge(),size=1,alpha=0.3)+ #,alpha= 0.8,stroke=0
      geom_boxplot(lwd=0.8,alpha = 0.3) + #
      STYLE +
      scale_y_continuous(labels = scales::percent,trans='log10') +
      geom_text(data = label_status, aes(x = Ant_status , y = 10, group=Ant_status,label = V1,fontface = "bold"),position = position_jitterdodge(seed = 2)) + # jitterdodge in geom_text will need the grouping aesthetic (eg, colour) to be defined inside ggplot() +
    facet_grid(. ~ gene)


    
    ## TODO!!!
    
##### analyze treated nurses and QUEENS!!
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
##### DELTA-DELTA-Ct METHOD
# One common way of analysing qPCR data is to use the “delta-delta-Ct” method. This involves calculating the difference between the Ct of the housekeeping gene and the test gene, then calculating the difference between the treated samples and the control.
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


##### SCRAPS

# 
# # calculate mean of b for group A
# mean_A <- d %>%
#   filter(c == "A") %>%
#   summarize(mean_b = mean(b))
# 
# # filter rows with mean_Ct value greater than 2*st.dev_mean_Ct value for group EF1
# outliers_highlight <- genes_data %>%
#   filter(gene == "EF1") %>%
#   filter(b > mean_A$mean_b)
# 
# # plot a scatterplot with lines connecting b values between groups
# ggplot(genes_data, aes(x = c, y = b, group = id)) +
#   geom_point(aes(color = c)) +
#   geom_line() +
#   geom_point(data = outliers_highlight, shape = 21, size = 5, fill = "red") +
#   ggtitle("Scatterplot of b values by group c with highlighted lines") +
#   xlab("Group") +
#   ylab("b value")
