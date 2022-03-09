
###############################################################################
###### NORMALISE VARS BEFORE LDA ##############################################
###############################################################################
library(bestNormalize)


#THE VARS WITH THE SAME NAME (FOR ACT AND REC) SHOULD BE NORMALISED TOGETHER  with the best common transformation per ROLE!
#COULD THAT BE DONE BY RESHAPING AND SPECIFING IN bestNormalize THAT THERE ARE TWO CATEGORIES?
#ONCE DONE, OUTPUT PLOTS FOR summary_MANUAL_transf

#empty base
summary_MANUAL_transf <- data.frame()[1:nrow(summary_MANUAL), ]
summary_MANUAL$ant1 <- as.factor(summary_MANUAL$ant1)
summary_MANUAL$ant2 <- as.factor(summary_MANUAL$ant2)
#transform variables if they are not normal
  for (variable in names(summary_MANUAL)){
    if (is.numeric(summary_MANUAL[,variable])) {
    val <- shapiro.test(summary_MANUAL[,variable])
    if (unname(val$p.value)<0.05) {
      print(paste("for [",variable, "] the Shapiro-Wilk Test has p < 0.05, the data is not normal. Transform it.",sep=" "))
      #find the best transformation
      #This function currently estimates the Yeo-Johnson, the Box Cox  (if the data is positive), the log_10(x+a), the square-root (x+a), the arcsinh and the ordered quantile normalization
      BNobject <- bestNormalize(summary_MANUAL[,variable])
      summary_MANUAL_transf$var <- BNobject$x.t; names(summary_MANUAL_transf)[names(summary_MANUAL_transf) == 'var'] <- paste(variable,class(BNobject$chosen_transform)[1],sep=".")
      
      
    }else{print(paste("for [",variable,"] the Shapiro-Wilk Test has p > 0.05, the data is normal. Keep it.",sep=" "))
      summary_MANUAL_transf[variable]<- summary_MANUAL[,variable]
    }}else{print(paste("non numeric attribute. Pasting [",variable,"] in the new dataset",sep = " "))
      summary_MANUAL_transf[variable]<- summary_MANUAL[,variable]}
}

#PLOT CHOSEN TRANSFORMATION
# MASS::truehist(BNobject$x.t, 
#                main = paste("Best Transformation:", 
#                             class(BNobject$chosen_transform)[1]), nbins = 12)


###############################################################################
###### PARAMETERS PLOTS #######################################################
###############################################################################

#SUMMARY_DATA
#descriptive analysis
str(summary_MANUAL)

#reshape data for plotting. Split by REC and ACT
summary_data_ACT <-summary_MANUAL %>% dplyr::select(contains(c("BEH", "ACT","interaction_length","strghtline","orient_angle_diff","movement_angle_diff"), ignore.case = TRUE))
summary_data_REC <- summary_MANUAL %>% dplyr::select(contains(c("BEH", "REC","interaction_length","strghtline","orient_angle_diff","movement_angle_diff"), ignore.case = TRUE))
#Rename columns to make them match and bind+ melt columns
summ_data_ACT  <- summary_data_ACT %>% rename_with(~str_remove(., c("_ACT|Act_"))); summ_data_ACT$Role <- "ACT"
summ_data_REC  <- summary_data_REC %>% rename_with(~str_remove(., c("_REC|Rec_"))); summ_data_REC$Role <- "REC"
summ_data_bind <- rbind(summ_data_ACT,summ_data_REC)
summ_data_long <- reshape2::melt(summ_data_bind,id.vars=c("BEH","Role","Name")) #explanation on the warning message https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message

#subset data by BEH for plotting
summary_data_G <- summ_data_long[which(summ_data_long$BEH == "G"),]
if (plot_all_BEH) {
summary_data_T <- summ_data_long[which(summ_data_long$BEH == "T"),]
summary_data_FR <- summ_data_long[which(summ_data_long$BEH == "FR"),]
summary_data_CR <- summ_data_long[which(summ_data_long$BEH == "CR"),]
#keep Trophallaxis, Cross Rest and Front REst togheter and fill by BEH when plotting
summary_data_T_FR_CR <- summ_data_long[which(!summ_data_long$BEH == "G"),]
}

#set up plot
pdf(file=paste(DATADIR,"Parameters_plots_post_30minWindow_4March22.pdf", sep = ""), width=6, height=6)
par(mfrow=c(2,3), family="serif" , mar = c(0.1, 0.1, 2.2, 0))

###plot divided by variable and Role for Grooming
vars_plot_G <- ggplot(summary_data_G, aes(value, fill = Role)) +
  facet_wrap(variable ~ .,scales="free") +
  theme_bw() +
  theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) #,legend.position="bottom",legend.justification='right'
vars_plot_G + geom_histogram(colour='black',alpha = 0.2,position="identity")+
  labs(title = "Histogram plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming",
       subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Replicates:",unique(interaction_MANUAL$REPLICATE)))
vars_plot_G + geom_density(alpha = 0.2) +
  labs(title = "Density plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming",
       subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Replicates:",unique(interaction_MANUAL$REPLICATE)))

if (plot_all_BEH) {
  ###plot divided by variable and Role for Trophallaxis
  vars_plot_T <- ggplot(summary_data_T, aes(value)) +
    facet_wrap(variable ~ .,scales="free") +
    theme_bw() +
    theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) 
  vars_plot_T + geom_histogram(colour='black',alpha = 0.2,position="identity") +
    labs(title = "Histogram for movement variables calculated from coordinates \n for interacting ants during Trophallaxis",
         subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))
  vars_plot_T + geom_density(alpha = 0.2) +
    labs(title = "Density plot for movement variables calculated from coordinates \n for interacting ants during Trophallaxis",
         subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))
  
  ###plot divided by variable and Role for Trophallaxis compared to FR and CR
  vars_plot_T_FR_CR <- ggplot(summary_data_T_FR_CR, aes(value, color = forcats::fct_inorder(BEH))) +
    facet_wrap(variable ~ .,scales="free") +
    theme_bw() +
    theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm'))
  #vars_plot_T_FR_CR + geom_histogram(colour='black',alpha = 0.2,position="identity") +
  #  labs(title = "Histogram for movement variables calculated from coordinates \n for interacting ants during Trophallaxis and for non-interacting ants during Front Rest and Cross Rest")
  vars_plot_T_FR_CR + geom_density(alpha = 0.2,aes(linetype=forcats::fct_inorder(BEH))) +
    labs(title = "Density plot for movement variables calculated from coordinates \n for interacting ants during Trophallaxis and for non-interacting ants during Front Rest and Cross Rest",
         subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO))
}

#INTERACTION_DATA


#-------------------------------------------------------------
#use same plots as before but for interactions
#THIS CASE WILL BE JUST CUTTING FOR 1 INTERACTION
interaction_data_19 <- interaction_MANUAL[which(interaction_MANUAL$ROW == 19),]

#TO BE FIXED
#reshape data for plotting. Split by REC and ACT
interaction_data_19$frame <- seq.int(nrow(interaction_data_19)) 
interaction_data_ACT <-interaction_data_19 %>% dplyr::select(contains(c("ROW","BEH", "Act_Name","traj_BOTH.ACT","angle_diff","straightline","frame"), ignore.case = TRUE))
interaction_data_REC <- interaction_data_19 %>% dplyr::select(contains(c("ROW","BEH", "Rec_Name","traj_BOTH.REC","angle_diff","straightline","frame"), ignore.case = TRUE))
#Rename columns to make them match and bind+ melt columns
int_data_ACT  <- interaction_data_ACT %>% rename_with(~str_remove(., c("_ACT|Act_|traj_BOTH.ACT.|traj_BOTH."))); int_data_ACT$Role <- "ACT"
int_data_REC  <- interaction_data_REC %>% rename_with(~str_remove(., c("_REC|Rec_|traj_BOTH.REC.|traj_BOTH."))); int_data_REC$Role <- "REC"
int_data_bind <- rbind(int_data_ACT,int_data_REC)
#int_data_bind$frame <- seq.int(nrow(int_data_bind)) 
int_data_long <- melt(int_data_bind,id.vars=c("BEH","Role","Name","frame")) #explanation on the warning message https://stackoverflow.com/questions/25688897/reshape2-melt-warning-message

#subset data by BEH for plotting and #remove ROW
interaction_data_G <- int_data_long[which(int_data_long$BEH == "G" & !int_data_long$variable == "ROW" ),]
interaction_data_G <- interaction_data_G[which(!is.na(interaction_data_G$value)),]
#cut the receiver as it will be the same for common parameters and add row number as time 
interaction_data_G <- interaction_data_G[which(!interaction_data_G$Role == "REC"),]
interaction_data_G <- interaction_data_G[which(!interaction_data_G$variable == "x"),]
interaction_data_G <- interaction_data_G[which(!interaction_data_G$variable == "y"),]
interaction_data_G <- interaction_data_G[which(!interaction_data_G$variable == "angle"),]
interaction_data_G <- interaction_data_G[which(!is.na(interaction_data_G$value)),]


#set up plot

# pdf(file=paste(DATADIR,"TEST_TEST_2.pdf", sep = ""), width=1.7, height=3.1)
# par(mfrow=c(2,6), family="serif" , mar = c(0.1, 0.1, 2.2, 0))
# #par(mfrow=c(2,3), family="serif" , mar = c(0.1, 0.1, 2.2, 0))
# 
# 
# ggplot(data=interaction_data_G,
#        aes(x=frame, y=value, colour=variable)) + facet_wrap(variable ~ .,scales="free",ncol = 1) + theme(legend.position = "none") +
#   geom_line(size=0.8)
# 
# dev.off()


###plot divided by variable and Role for Grooming
# vars_plot_G <- ggplot(interaction_data_G, aes(frame, fill = BEH)) +
#   facet_wrap(variable ~ .,scales="free") +
#   theme_bw() +
#   theme(text=element_text(family="serif",size=9), legend.key.size = unit(0.3, 'cm')) #,legend.position="bottom",legend.justification='right'
# # vars_plot_G + geom_histogram(colour='black',alpha = 0.2,position="identity",bins = 500)+
# #   labs(title = "Histogram plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming"#,
# #        #subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO)
# #        )
# vars_plot_G + geom_density(alpha = 0.2) +
#   labs(title = "Density plot for movement variables calculated from coordinates \n for Actors and Receivers during Grooming",
#        #subtitle = paste( "Periods:",unique(interaction_MANUAL$PERIOD),". Window:",time_start_ISO,"-",time_stop_ISO)
#        )


#------------------------------------------------------------

###################### TO REVIEW THE SIGNIFICANCE OF THIS ##############################################################
interaction_data_circ <- circular(interaction_MANUAL$traj_BOTH.angle_diff)

##angular_differences plot per interaction
interaction_data_G <- subset(interaction_MANUAL[!is.na(interaction_MANUAL$traj_BOTH.angle_diff),], BEH == "G")
for (row in unique(interaction_data_G$ROW)) {
  single_interaction <-subset(interaction_data_G,ROW == row)
  p <- plot.circular(single_interaction$traj_BOTH.angle_diff, pch = 16, cex = 0.8, stack=TRUE, bins=100, tcl.text	 = 0.5, tol= 0, shrink = 2.5, sep= 0.04, #xlim=c(-1,1), ylim=c(-2.5,2.5), 
                     main = paste("Behaviour: G, \n", "interaction N: ",row, sep=""), sub=NULL) #REPLICATE, ", ", PERIOD, ", \n", "behaviour:", BEH,
  arrows.circular(mean.circular(single_interaction$traj_BOTH.angle_diff))
  arrows.circular(0, col = "red")
}

if (plot_all_BEH) {
interaction_data_T <- subset(interaction_MANUAL[!is.na(interaction_MANUAL$traj_BOTH.angle_diff),], BEH == "T")
  for (row in unique(interaction_data_T$ROW)) {
    single_interaction <-subset(interaction_data_T,ROW == row)
    p <- plot.circular(single_interaction$traj_BOTH.angle_diff, pch = 16, cex = 0.8, stack=TRUE, bins=100, tcl.text	 = 0.5, tol= 0, shrink = 2.5, sep= 0.04, #xlim=c(-1,1), ylim=c(-2.5,2.5), 
                       main = paste("Behaviour: T, \n", "interaction N: ",row, sep=""), sub=NULL)
    arrows.circular(mean.circular(single_interaction$traj_BOTH.angle_diff))
    arrows.circular(0, col = "red")
  }
}

#means
## calculate circular mean angles for each interaction Row
#interaction_data_mean <- aggregate(traj_BOTH.angle_diff ~ ROW + BEH + REPLICATE + PERIOD, FUN=mean, na.action=na.omit, interaction_MANUAL)
interaction_data_circ_mean <- aggregate(circular(traj_BOTH.angle_diff) ~ ROW + BEH + REPLICATE + PERIOD, FUN=mean, na.action=na.omit, interaction_MANUAL)


for (beh in unique(interaction_data_circ_mean$BEH)) {
  single_interaction <-subset(interaction_data_circ_mean,BEH == beh)
  plot.circular(single_interaction$`circular(traj_BOTH.angle_diff)`, pch = 16, cex = 0.8,stack=TRUE, bins=100, tcl.text	 = 0.5, tol= 0, shrink = 2.5, sep= 0.08, xlim=c(-0.4,0.4), ylim=c(-0.4,0.4), 
                main = paste("Mean interaction angle for ", beh, "\n", "Tot N interactions: ",NROW(subset(interaction_data_circ_mean,BEH == beh)), sep=""), sub=NULL)
  #densityline <- density(single_interaction$`circular(traj_BOTH.angle_diff)`, bw=30)
  #lines(densityline, col=2)
  arrows.circular(mean(single_interaction$`circular(traj_BOTH.angle_diff)`))
  arrows.circular(0, col = "red")
}


#close the pdf
dev.off()


#polar coordinates plot alternative
# p <- ggplot(interaction_data_G, aes(x=as.factor(ROW), y=traj_BOTH.angle_diff)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
# geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
# # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
# ylim(-10,300) +
# # Custom the theme: no axis title and no cartesian grid
# theme_minimal() +
# # theme(
# #   axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank()# , plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
# # ) +
# labs(title = paste(REPLICATE, ", ", PERIOD, ", ", BEH,", ", sep="")) +
# coord_polar(start = 0) # This makes the coordinate polar instead of cartesian.
# print(p)





