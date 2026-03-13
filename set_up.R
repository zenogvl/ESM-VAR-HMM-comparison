
### Install and load all required packages. 

install.packages(setdiff(c("brms",
                           "BH",
                           "depmixS4",
                           "gtools",
                           "haven", 
                           "imputeTS", 
                           "Metrics",
                           "mgm", 
                           "parallel",
                           "readxl",
                           "rlist",
                           "stats", 
                           "tidyverse",
                           "vars",
                           "wesanderson",
                           "mHMMbayes",
                           "mlVAR"
                           ), rownames(installed.packages())))  
library(tidyverse)


for(f in list.files("Scripts/worker_functions/")){
  source(paste0("Scripts/worker_functions/", f))
}


dataset_names_all <- gsub("Rescaled_|.csv", "", list.files("Data/CleanRescaled/"))
dataset_names_all <- dataset_names_all[!(dataset_names_all %in% "Knouse_2023")]
dataset_names <- gsub("RMSE_|.csv", "", list.files("Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/RMSE/"))
dataset_names <- dataset_names[!(dataset_names %in% c("Rombaoa_2023", "Wendt_2020"))]


variables_list <- readRDS("Data/variables_list.rds")




standard_figure_with <- 5.9
