#######################################################################################################################################

SOFTWARE FOR MANUSCRIPT 'Social network plasticity decreases disease transmission in a eusocial insect'

#######################################################################################################################################


The code.zip file, together with the repository https://github.com/laurentkeller/anttrackingUNIL contains all software necessary 
to understand and assess the conclusions of the manuscript.
All analyses can be replicated using the data provided in repository https://doi.org/10.5281/zenodo.1322669

#######################################################################################################################################

Please follow the following instructions (in order) to run the analyses:

#######################################################################################################################################
PART 1 - AntTracking Tools
1. Install all software included in repository https://github.com/laurentkeller/anttrackingUNIL (follow instructions contained in README.txt file)
2. Ensure that all executables produced in the previous step are contained within a folder entitled 'executables', 
   which will be used in later steps of the analyses
#######################################################################################################################################

#######################################################################################################################################
PART 2 - Download and unzip data
1. Download the data.zip file from Repository https://doi.org/10.5281/zenodo.1322669
2. Unzip the data.zip file into a folder entitled 'Data_Repository'
3. Download the README.txt file from Repository https://doi.org/10.5281/zenodo.1322669 to obtain a detailed description of the data files
#######################################################################################################################################

#######################################################################################################################################
PART 3 - Download and unzip code
1. Download the code.zip file from Repository https://doi.org/10.5281/zenodo.1322676
2. Unzip the code.zip file into a folder entitled 'Code_Repository'
#######################################################################################################################################
    
#######################################################################################################################################
PART 4 - Process data
1. Install R libraries: - Navigate to 'Code_Repository/1_data_post_processing/source'
                        - Open file 'libraries.R'
                        - Install all libraries listed

2. Input folder paths:  - Navigate to 'Code_Repository/1_data_post_processing'
                        - Open each of the files A_Main_experiment.R, B_Age_experiment.R and C_Survival_experiment.R 
                        - Enter your own path for the data (experiment sub-folder within the Data_Repository folder), for the source code (source subfolder within Code_Repository/1_data_post_processing),
                         and for the 'exectutables' folder on lines 6-8, as in the example provided:
data_path  <- "~/Repositories/Data_Repository/main_experiment"
code_path  <- "~/Repositories/Code_Repository/1_data_post_processing/source"
executables_path <- "~/executables"
                        - Save the files

3. Run analyses:        - Navigate to 'Code_Repository/1_data_post_processing'
                        - Run the code files 'A_Main_experiment.R', 'B_Age_experiment.R' and 'C_Survival_experiment.R'
                          Note that each of these three code files are independent and can be run either individually or in parallel.
#######################################################################################################################################

#######################################################################################################################################
PART 5 - Statistical analyses and plots
1. Install R libraries: - Navigate to 'Code_Repository/2_statistics_and_plotting/source'
                        - Open file 'libraries.R'
                        - Install all libraries listed

2. Input folder paths:  - Navigate to 'Code_Repository/2_statistics_and_plotting'
                        - Open the file 'Statistics_and_plots.R'
                        - Enter your own path for the output figure folder, for the data (Data_Repository folder), and for the source code (source subfolder within Code_Repository/2_statistics_and_plotting),
                        as in the example provided:
figurefolder <- "~/figures"
disk_path    <- "~/Repositories/Data_Repository"
source_path  <- "~/Repositories/Code_Repository/2_statistics_and_plotting/source"
                        - Save the file
3. Run statistical analyses and plots: - Navigate to 'Code_Repository/2_statistics_and_plotting'
                                       - Run code file 'Statistics_and_plots.R'
#######################################################################################################################################



