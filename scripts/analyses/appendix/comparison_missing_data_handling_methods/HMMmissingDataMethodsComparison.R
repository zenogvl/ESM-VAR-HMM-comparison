
install.packages(setdiff(c("tidyverse", "hmmr", "stats", "reshape2", "ggpubr"), rownames(installed.packages())))  
library(tidyverse)
theme_set(theme_classic())

### Paths to data and results 
dayAndMissingPath <- "Output/N1_Results/RescaledResults/HMM_daymissing_split_results/"
pathData <- "Data/CleanRescaled/"
onlyDayPath <- "Output/N1_Results/RescaledResults/HMM_onlyday_split_results/"
kalmanFilterPath <- "Output/N1_Results/RescaledResults/HMM_kalman_filter_results/"
dataNameRemoval <- "Rescaled_|.csv"

### Load in RMSE Data 

variableList <- readRDS("Data/variables_list.rds")

#RMSE
rmseDayMissingSplit <- list()
for(f in list.files(paste0(dayAndMissingPath, "RMSE/"))){
  rmseDayMissingSplit[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0(dayAndMissingPath, "RMSE/", f))
}

rmseOnlyDaySplit <- list()
for(f in list.files(paste0(onlyDayPath, "RMSE/"))){
  rmseOnlyDaySplit[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0(onlyDayPath, "RMSE/", f))
}

rmseKalmanFiler <- list()
for(f in list.files(paste0(kalmanFilterPath, "RMSE/"))){
  rmseKalmanFiler[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0(kalmanFilterPath, "RMSE/", f))
}

datasetNames <- names(rmseKalmanFiler)
nDatasets <- length(datasetNames)

hmmRmseCopmarison <- data.frame(dataset = datasetNames,
                                rateKalmanFilter = rep(0, nDatasets), 
                                rateDayAndMissingSplit = rep(0, nDatasets),
                                rateOnlyDaySplit = rep(0, nDatasets),
                                meanKalmanFilter = rep(0, nDatasets), 
                                meanDayAndMissingSplit = rep(0, nDatasets),
                                meanOnlyDaySplit = rep(0, nDatasets)
                                )


ds <- "Hensel_2022"

lowestRMSE <- data.frame()
for(ds in datasetNames){
  fittedPP <- intersect(intersect(na.omit(rmseDayMissingSplit[[ds]])$ID, 
                                  na.omit(rmseOnlyDaySplit[[ds]])$ID),
                                  na.omit(rmseKalmanFiler[[ds]])$ID)
  
  rmseDayMissingSplitMean <- rmseDayMissingSplit[[ds]] %>%
    select(-X) %>%
    filter(ID %in% fittedPP) %>%
    rowwise() %>%
    mutate(meanDayMissingSplit = mean(c_across(!ID))) %>%
    #ungroup() %>%
    select(ID, meanDayMissingSplit)
  
  rmseOnlyDaySplitMean <- rmseOnlyDaySplit[[ds]] %>%
    select(-X) %>%
    filter(ID %in% fittedPP) %>%
    rowwise() %>%
    mutate(meanOnlyDaySplit = mean(c_across(!ID))) %>%
    #ungroup() %>%
    select(ID, meanOnlyDaySplit)
  
  rmseKalmanFilerMean <- rmseKalmanFiler[[ds]] %>%
    select(-X) %>%
    filter(ID %in% fittedPP) %>%
    rowwise() %>%
    mutate(meanKalmanFilter = mean(c_across(!ID))) %>%
    #ungroup() %>%
    select(ID, meanKalmanFilter)
 
  rmseMeanAll <- rmseDayMissingSplitMean %>%
    inner_join(rmseOnlyDaySplitMean, by = "ID") %>%
    inner_join(rmseKalmanFilerMean, by = "ID") %>%
    mutate(ndistinct = n_distinct(c_across(!ID))) %>% 
    filter(ndistinct == 3) %>%
    select(-ndistinct) %>%
    mutate(lowestRmse = which.min(c_across(!ID))) %>%
    mutate(lowestRmseFactor = factor(lowestRmse, 
                                     levels = c(1,2,3),
                                     labels = c("DayMissingSplit", "OnlyDaySplit", "KalmanFiler")))
  
  lowestRMSE <- rmseMeanAll %>% 
    ungroup() %>%
    count(lowestRmseFactor) %>%
    mutate(dataset = rep(ds, n_distinct(lowestRmseFactor))) %>%
    mutate(lowestRmseRate = n/sum(n)) %>% 
    rbind(lowestRMSE)
  

  
}

lowestRMSE %>% 
  mutate(Studies = gsub("_", " ", dataset)) %>%
  ggplot(aes(fill = lowestRmseFactor, x = Studies, y = lowestRmseRate)) + 
  geom_bar( stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("HMM") + 
  theme(legend.position = "bottom")




