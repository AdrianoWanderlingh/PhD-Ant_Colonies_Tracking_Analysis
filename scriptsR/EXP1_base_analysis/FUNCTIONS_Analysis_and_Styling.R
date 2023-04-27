########################################################################
##########  ADRIANO EXP1: ANALYSIS AND STYLING FUNCTIONS ###############

require(report)

########## STATS FUNCTIONS ###############

# function to test normality of residuals
test_norm <- function(resids) {
  print("Testing normality")
  if (length(resids) <= 300) {
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
    print("below 0.05, the data significantly deviate from a normal distribution")
  } else {
    print("More than 300 data points so using the skewness and kurtosis
approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}

# function to report a model output
output_lmer <- function(model) {
  print("------------RESIDUALS NORMALITY------------")
  test_norm(residuals(model))
  print("------------SUMMARY------------")
  print(summary(model))
  print("------------ANOVA------------")
  print(Anova(model))
  print("------------RSQUARED------------")
  print(r.squaredGLMM(model))
  print("------------REPORT------------")
  print(report(model))
  #tab_model(model)
}


# function to perform posthocs
posthoc_list <- list()
interactions_to_explore <- list()
compute_posthocs <- function(model) {
  warning("this function has only been tested with lmer()")
  # check taht there are significant outputs
  if (length(row.names(Anova(model)[Anova(model)$"Pr(>Chisq)" < 0.05, ])) == 0) {
    print("there are no significant vars.")
  } else {
    for (SIG.VAR in row.names(Anova(model)[Anova(model)$"Pr(>Chisq)" < 0.05, ])) {
      if (grepl(":", SIG.VAR)) {
        warning(paste0(SIG.VAR, "is an interaction, currently this situation is not handled by the function."))
        interactions_to_explore <- c(interactions_to_explore, list(paste(GENE,GROUP,SIG.VAR,deparse(substitute(model)), sep = "-") ))
        } else {
        # check if the variable is not numeric . to do so, we need to access the dataframe from the model
        if (!is.numeric(get(gsub("\\[.*", "", as.character(model@call)[3]))[, SIG.VAR])) {
          print(paste0("Performing posthocs for the significant var: ", SIG.VAR))
          arg <- list("Tukey")
          names(arg) <- SIG.VAR
          # glht (mcp) : General linear hypotheses Testing (glht) and multiple comparisons for parametric models (mcp)
          cmp <- do.call(mcp, arg)
          posthoc_SIG.VAR <- summary(glht(model, linfct = cmp), test = adjusted("BH"))
          # Set up a compact letter display of all pair-wise comparisons
          model_means_cld <- cld(posthoc_SIG.VAR)
          # create dataframe usable with ggplot geom_text
          model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters)))
          # add column name
          model_means_cld$newcol <- NA
          colnames(model_means_cld)[which(names(model_means_cld) == "newcol")] <- SIG.VAR
          model_means_cld[, SIG.VAR] <- row.names(model_means_cld)
          rownames(model_means_cld) <- NULL
          # add to list
          posthoc_list <- c(posthoc_list, list(model_means_cld))
          if (exists("ID_model")) {
            names(posthoc_list)[length(posthoc_list)] <- paste(ID_model,SIG.VAR,  deparse(substitute(model)),sep = "-")
          } else {
            names(posthoc_list)[length(posthoc_list)] <- paste(SIG.VAR,deparse(substitute(model)), sep = "-")
          }
          print(paste(deparse(substitute(model)), SIG.VAR, sep = "_"))
          print(model_means_cld)
        } # SIG.VAR LOOP
      } # check if is an interaction
    } # check if numeric
    print("call posthoc_list to get posthocs")
  } # if significant vars exist
  return(posthoc_list)
}


# The calculate_weights function takes as input a categorical variable group and an optional argument ratio, which specifies the desired ratio between the weights of the largest and smallest groups
calculate_weights <- function(group) {
  print("function to calculate_weights for imbalanced datasets")
  #make sure that group is a factor
  group <- as.factor(group)
  # Calculate the proportions of each group
  proportions <- table(group) / length(group)
  print(proportions)
  # Calculate the weights for each group
  weights <- numeric(length(group))
  
  for (g in levels(group)) {
    weights[group == g] <- 1 / proportions[g]
  }
  # Normalize the weights so that they sum up to 1
  weights <- weights / sum(unique(weights))
  return(weights)
}

## pretty print the model selection output
nice_print_model_sel <- function(model_output) {
  # clean output
  sel.table <- round(as.data.frame(model_output)[-c(1:5)], 3)
  # number of parameters (df) should be K
  names(sel.table)[1] <- "K"
  sel.table$Model <- rownames(sel.table)
  rownames(sel.table) <- NULL
  # replace Model name with formulas
  for (i in 1:nrow(sel.table)) sel.table$formula[i] <- as.character(formula(get(sel.table$Model[i])))[3]
  return(sel.table)
}


# convert significance levels to stars
add_star <- function(p) {
  if (p<0.001) {
    return('***')
  } else if (p<0.01) {
    return('**')
  } else if (p<0.05) {
    return('*')
  } else {
    return('ns')
  }
}

#Box-cox transformation
Box_Cox <- function(x) {
  library(MASS)
  bc <- boxcox(x ~ 1, plotit = FALSE)
  #computes the log-likelihood for a range of lambda values and returns the lambda that maximizes the log-likelihood
  lambda <- bc$x[which.max(bc$y)]
  (x^lambda - 1) / lambda
}

####################################################
####################################################
####################################################

###########    PLOT SAVING    ###############

library(ggplot2)

  SavePrint_plot <- function(plot_obj, plot_name, dataset_name, save_dir, plot_size = c(4.63, 2.59), dpi = 100, print_plot=F) {
    save_dir_plots <- paste0(save_dir,"/Grooming_plots/")
    
    # Create the directory if it doesn't exist
    if (!dir.exists(save_dir_plots)) {
      dir.create(save_dir_plots, recursive = TRUE)
    }
    
    # Check if the directory is writable
    if (!file.access(save_dir_plots, 2)) {
      # Save plot as png
      ggsave(paste0(save_dir_plots,dataset_name, "_", plot_name,"_", Sys.Date(), ".png"), plot = plot_obj, width = plot_size[1], height = plot_size[2], dpi = dpi)
      # Save plot as pdf
      ggsave(paste0(save_dir_plots,dataset_name, "_", plot_name,"_", Sys.Date(),".pdf"), plot = plot_obj, width = plot_size[1], height = plot_size[2])
      # Add plot to cumulative pdf file
      pdf(paste0(save_dir_plots,"all_plots.pdf"), onefile = TRUE, paper = "a4") #  width = plot_size[1], height = plot_size[2])
      dev.off()
      if (print_plot) {
        print(plot_obj)
      }
      
    } else {
      cat("Error: The directory is not writable.")
    }
  }


########## STYLING FUNCTIONS ###############
library(RColorBrewer)
library(shades)
library(colorspace)
library(plotwidgets)
library(ggplot2)
library(ggnewscale)

# #Create a custom color scale FOR COLONIES + treatments
# FullPal <- scales::viridis_pal(option = "D")(20)
# myColorsSmall  <- tail(FullPal,5)
# myColorsLarge  <- head(FullPal,5)
# Cols_treatments <- c("#440154FF","#FDE725FF") #"#31688EFF"
# myColors      <- c(myColorsLarge,myColorsSmall, Cols_treatments)
# names(myColors) <- c("R3BP","R5BP","R7BP","R8BP","R12BP","R1SP", "R2SP", "R5SP", "R7SP","R11SP","Big Pathogen","Small Pathogen")
# colScale <- scale_colour_manual(name = "Colony",values = myColors,drop=TRUE)
# fillScale <- scale_fill_manual(name = "Colony",values = myColors,drop=TRUE)

#### CREATE CONSISTENT COLORING FOR ALL THE PLOTTING
# GENERATE 4 MAIN COLORS FOR THE 4 TREATMENTS BS,BP,SS,SP + SHADES FOR THE RESPECTIVE COLONIES

# Create a color palette with 4 colors as distant as possible
colors_full <- scales::viridis_pal(option = "D")(100)
# Create a list to store the subsets of colors
color_subsets <- list()
Shades <- list()
# Loop over the 4 colors to get shades of the colour in a +5, -5 range
for (i in c(10, 30, 70, 90)) { # BE CAREFUL: IF CHANGING THIS, ALSO CHANGE Treat_colors
  color_ramp <- colorRampPalette(colors_full[(i-5):(i+5)])
  #color_ramp <- colorRampPalette(colors_full[(i-1):(i+1)])
  color_subsets[[i]] <- color_ramp(12)
}


Treat_colors <- structure(list(Shades = c(colors_full[10],colors_full[30],colors_full[70],colors_full[90]), 
                               Treatment = c("Big Pathogen", "Big Sham", "Small Pathogen", "Small Sham")),
                          row.names = c(NA, 4L), class = "data.frame")

#show_col(c(colors_full[10],colors_full[30],colors_full[70],colors_full[90]))


#clean the output
color_subsets <- color_subsets[lapply(color_subsets,length)>0]

#Darken the colours progressively to increase contrast among the colonies of the same treatment
for (i in 1:4) {
# Define the color gradient
color_shades <- color_subsets[[i]]

# Convert the colors from hexadecimal to HSL (hue, saturation, lightness)

colors_lightGrad <- c()
# Decrease the lightness of each color by an increasing amount
lightness_decrease <- rev(seq(from = 0, to = 0.2, length.out = length(color_shades)))
lightness_increase <- seq(from = 0, to = 0.2, length.out = length(color_shades))

for (j in 1:length(color_shades)) {
  hsl_colors <- col2hsl(color_shades[j])
  hsl_colors[3] <- hsl_colors[3] - lightness_decrease[j]
  hsl_colors[3] <- hsl_colors[3] + lightness_increase[j]
  colors_lightGrad <- c(colors_lightGrad,hsl2col(hsl_colors))
}

Shades[[i]] <- colors_lightGrad
}

# #inspect output
# par(mfrow = c(2, 2))
# for (i in 1:4) {
#   plot(1:12, 1:12,
#        col = Shades[[i]], # color_subsets
#        pch = 19,
#        cex = 5,
#        xaxt = "n",
#        yaxt = "n",
#        xlab = "",
#        ylab = "") 
# }

### ADD THE METADATA
meta.data <- read.table(paste("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")

# create groups to assign the colours
Cols <- list()
# divide each size_treat into a list element with its colonies inside
for (i in 1:length(unique(meta.data$size_treat))) {
  treatment <- unique(meta.data$size_treat)[i]
  treatment_vector <- unique(meta.data$REP_treat[grepl(treatment, meta.data$REP_treat)])
  Cols[[i]] <- treatment_vector
  names(Cols)[i] <- treatment
}
#name list elements according to the favoured pattern (the colour order I have been using since the first plots)
names(Shades)[1] <- "BP"
names(Shades)[2] <- "BS"
names(Shades)[3] <- "SP"
names(Shades)[4] <- "SS"

# dput((list_of_vectors))
# dput((Shades))

# Create an empty dataframe to store the results
colour_palette <- data.frame()

# bind together colours, REP_treats and treatments
# Loop over the list "Cols" to create the dataframe
for (group in names(Cols)) {
  group_cols <- Cols[[group]]
  group_shades <- Shades[[group]]
  # Take a random subset of the shades that matches the length of the cols
  rand_shades <- sample(group_shades, length(group_cols))
  # Create a dataframe with two columns: "Cols" and "Shades"
  group_colour_palette <- data.frame(Cols = group_cols, Shades = rand_shades, Treatment = group)
  # Append the current group dataframe to the overall dataframe
  colour_palette <- rbind(colour_palette, group_colour_palette)
}

# #visualise output
# ggplot(colour_palette, aes(x=Treatment, y=Shades, fill=Shades)) + 
#   geom_tile(colour="white", size=0.5) + 
#   scale_fill_identity() + 
#   theme_void() + 
#   labs(title="Colours by Treatment", x="Treatment", y="Shades") + 
#   theme(axis.text.y=element_text(angle=0, hjust=1)) +
#   facet_wrap(~Treatment, scales = "free_y")

myColors_Colony <- colour_palette$Shades
names(myColors_Colony) <- colour_palette$Cols
myColors_Treatment <- Treat_colors$Shades
names(myColors_Treatment) <- Treat_colors$Treatment

#### DEFINE THE FILL AND COLOR AS STANDALONE

#COLOR
# geom_point(aes(color = Sample)) +
colScale_Colony <- scale_colour_manual(name = "Colony", values = myColors_Colony, drop = TRUE)
#new_scale_color() +
#new_scale_fill() +
# geom_line(aes(color = Sample)) +
colScale_Treatment <- scale_color_manual(name = "Treatment", values = myColors_Treatment) #for lines

####FILL
# geom_point(aes(color = Sample)) +
colFill_Colony <- scale_fill_manual(name = "Colony", values = myColors_Colony, drop = TRUE)
#new_scale_color() +
#new_scale_fill() +
# geom_line(aes(color = Sample)) +
colFill_Treatment <- scale_fill_manual(name = "Treatment", values = myColors_Treatment) #for lines

#### DEFINE REMAINING PLOT STYLE
# ggplot PLOT STYLE
STYLE <- list(
  #colScale, fillScale,
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  theme_bw(),
  scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

STYLE_CONT <- list(
  #colScale, fillScale,
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  theme_bw()
)


remove <- c("color_ramp", "color_shades", "color_subsets", "colors_full", 
            "colors_lightGrad", "colour_palette", "Cols",
            "group", "group_colour_palette", 
            "group_cols", "group_shades", "hsl_colors", "i", "j", "lightness_decrease", 
            "lightness_increase", "meta.data", "myColors_Colony", 
            "rand_shades", 
            "Shades")
# cleaning
rm(list = ls()[which(ls() %in% remove)])
gc()


# #COLOUR SCALES TAKE NAMED VECTORS AS INPUTS
# myColors <- c(Treat_colors$Shades,colour_palette$Shades)
# names(myColors) <- c(Treat_colors$Treatment,colour_palette$Cols)
# colScale <- scale_colour_manual(values = myColors, drop = TRUE) #name = "Colony", 
# fillScale <- scale_fill_manual(values = myColors, drop = TRUE)



