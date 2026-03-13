install.packages(setdiff(c("tidyverse", "hmmr", "stats", "reshape2", "ggpubr"), rownames(installed.packages())))  
library(tidyverse)
theme_set(theme_classic())

### Paths to data and results 


dayAndMissingPath <- "Output/N1_Results/RescaledResults/HMM_daymissing_split_results/"
pathData <- "Data/CleanRescaled/"
onlyDayPath <- "Output/N1_Results/RescaledResults/HMM_onlyday_split_results/"
kalmanFilterPath <- "Output/N1_Results/RescaledResults/HMM_kalman_filter_results/"
dataNameRemoval <- "Rescaled_|.csv"

##########################################################################################################################
#####################################  Separated by Day and Missing ######################################################
##########################################################################################################################


Ntries <- list()
for(f in list.files(paste0(dayAndMissingPath, "Ntries/"))){
  Ntries[[gsub("Ntries_|.csv", "", f)]] <- read.csv(paste0(dayAndMissingPath, "Ntries/", f))
}
StateSelect <- list()
for(f in list.files(paste0(dayAndMissingPath, "Nstates/"))){
  StateSelect[[gsub("Nstates_|.csv", "", f)]] <- read.csv(paste0(dayAndMissingPath, "Nstates/", f))
}

variables_list <- readRDS("Data/variables_list.rds")
NtriesWObs <- list()
succesPercentage <- matrix(NA, 0, 6, 
                           dimnames = list(c(),c("ds", "S2", "S3", "S4", "S5", "nVars"))
)

for(f in list.files(pathData)){
  ds <- gsub(dataNameRemoval,"",f)
  print(ds)
  #Read in the data
  data <- na.omit(read.csv(paste0(pathData, f)))
  #Create a object that shows the number of observations for every pp
  ncount <- data.frame(ID = as.integer(names(table(data$ID))), 
                       n = as.matrix(table(data$ID)),
                       nvars = length(variables_list[[ds]]), 
                       ds = ds)
  #If no pp was successful create a dummy object
  if(is.null(Ntries[[ds]])){
    Ntries[[ds]] <- cbind(ncount$ID, matrix(NA, nrow(ncount), 5))
    colnames(Ntries[[ds]]) <- c("ID", "X", "X2", "X3", "X4", "X5")
    StateSelect[[ds]] <- cbind(ncount$ID, matrix(NA, nrow(ncount), 2))
    colnames(StateSelect[[ds]]) <- c("ID", "X", "nstate")
  }
  
  #Merge this with the Ntries object 
  merged <- Reduce(function(x, y) merge(x, y, all=TRUE), list(ncount, Ntries[[ds]], StateSelect[[ds]]))
  
  NtriesWObs[[ds]] <- merged <- merged %>% 
    dplyr::select(!X) %>%
    rename(S2 = X2, S3 = X3, S4 = X4, S5 = X5, usedState = nstate)
  
  ## Check how often a fit was succesfull
  #Change all NA's to 100 to indicate failure 
  merged <- merged %>% mutate(across(starts_with("S"), ~replace_na(.,100)))
  #Get the percantage that isn't 100. 
  succesPercentage <- rbind(succesPercentage,
        c(ds,(apply(dplyr::select(merged, starts_with("S")) != 100, 2, sum)/nrow(merged)), length(variables_list[[ds]])))
}

#Create workable df's for percentages 
succesPercentage <- as.data.frame(succesPercentage)
succesPercentage_long <- reshape2::melt(succesPercentage, 
                                        id.vars = c("ds", "nVars"), 
                                        variable.name = "Nstate", 
                                        value.name = "percentage") %>% 
  mutate(percentage = as.numeric(percentage))


#Percentage of succesfull fits per state per study 
ggplot(succesPercentage_long, aes(fill = Nstate, x = ds, y = percentage)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) 

#Create workable df's for pp's 
AllNtries <- bind_rows(NtriesWObs)
AllNtries_long <- reshape2::melt(AllNtries, 
                                 id.vars = c("ID", "n", "ds",  "nvars", "usedState"), 
                                 variable.name = "Nstate", 
                                 value.name = "tries") %>%
  dplyr::mutate(tries = replace_na(tries, 100)) %>%
  filter(n > 20)
  
#Plot number of tries against number of obs
ggplot(AllNtries_long, aes(x = n, y = tries)) +
  geom_point(size = 1)  + 
  facet_wrap(~Nstate)


#Create 4 plots for average number of tries against the number of observations 
AllNtries_long %>%
  group_by(n, Nstate) %>%
  summarise(mean = mean(tries)) %>%
  na.omit() %>%
  ggplot(aes(n , mean)) + 
  geom_point() + 
  facet_wrap(~Nstate)


#Create plots with percentage of successful fits against the number of observations. 
p_obsVSsuccesfit <- AllNtries_long %>%
  mutate(fit = as.integer(tries != 100)) %>%
  filter(n < 190) %>%
  group_by(n, Nstate) %>%
  summarise(mean = mean(fit)) %>%
  na.omit() %>%
  ggplot(aes(n , mean)) + 
  geom_point() + 
  geom_smooth(method = "gam", se = FALSE) +
  facet_wrap(~Nstate, nrow = 1)
p_obsVSsuccesfit

#Plots that shows the average number of selected states of every study against the average number of observations and number of variables
p_stateVSnvars <- AllNtries_long %>%
  na.omit() %>% 
  group_by(ds, nvars) %>%
  summarise(AvState = mean(usedState), 
            AvN = mean(n)) %>%
  #filter(nvars < 20) %>%
  ggplot(aes(x = nvars, y = AvState)) +
  geom_point() 

p_stateVSAvN <- AllNtries_long %>%
  na.omit() %>% 
  group_by(ds, nvars) %>%
  summarise(AvState = mean(usedState), 
            AvN = mean(n)) %>%
  #filter(nvars < 20) %>%
  ggplot(aes(x = AvN, y = AvState)) +
  geom_point() +
  ylab("")
  
ggpubr::ggarrange(p_stateVSnvars, p_stateVSAvN)



##########################################################################################################################
#####################################   Only Separated by Day  ###########################################################
##########################################################################################################################

Ntries_sOD <- list()
for(f in list.files(paste0(onlyDayPath, "Ntries/"))){
  Ntries_sOD[[gsub("Ntries_|.csv", "", f)]] <- read.csv(paste0(onlyDayPath, "Ntries/", f))
}
StateSelect_sOD <- list()
for(f in list.files(paste0(onlyDayPath, "Nstates/"))){
  StateSelect_sOD[[gsub("Nstates_|.csv", "", f)]] <- read.csv(paste0(onlyDayPath, "Nstates/", f))
}

variables_list <- readRDS("Data/variables_list.rds")
NtriesWObs_sOD <- list()
succesPercentage_sOD <- matrix(NA, 0, 6, 
                           dimnames = list(c(),c("ds", "S2", "S3", "S4", "S5", "nVars"))
)
for(f in list.files(pathData)){
  ds <- gsub(dataNameRemoval,"",f)
  print(ds)
  #Read in the data
  data <- na.omit(read.csv(paste0(pathData, f)))
  #Create a object that shows the number of observations for every pp
  ncount <- data.frame(ID = as.integer(names(table(data$ID))), 
                       n = as.matrix(table(data$ID)),
                       nvars = length(variables_list[[ds]]), 
                       ds = ds)
  #If no pp was successful create a dummy object
  if(is.null(Ntries_sOD[[ds]])){
    Ntries_sOD[[ds]] <- cbind(ncount$ID, matrix(NA, nrow(ncount), 5))
    colnames(Ntries_sOD[[ds]]) <- c("ID", "X", "X2", "X3", "X4", "X5")
    StateSelect_sOD[[ds]] <- cbind(ncount$ID, matrix(NA, nrow(ncount), 2))
    colnames(StateSelect_sOD[[ds]]) <- c("ID", "X", "nstate")
  }
  
  #Merge this with the Ntries object 
  merged <- Reduce(function(x, y) merge(x, y, all=TRUE), list(ncount, Ntries_sOD[[ds]], StateSelect_sOD[[ds]]))
  
  NtriesWObs_sOD[[ds]] <- merged <- merged %>% 
    dplyr::select(!X) %>%
    rename(S2 = X2, S3 = X3, S4 = X4, S5 = X5, usedState = nstate)
  
  
  
  ## Check how often a fit was succesfull
  #Change all NA's to 100 to indicate failure 
  merged <- merged %>% mutate(across(starts_with("S"), ~replace_na(.,100)))
  #Get the percantage that isn't 100. 
  succesPercentage_sOD <- rbind(succesPercentage_sOD,
                            c(ds,(apply(dplyr::select(merged, starts_with("S")) != 100, 2, sum)/nrow(merged)), length(variables_list[[ds]])))
}

#Create workable df's for percentages 
succesPercentage_sOD <- as.data.frame(succesPercentage_sOD)
succesPercentage_long_sOD <- reshape2::melt(succesPercentage_sOD, 
                                            id.vars = c("ds", "nVars"), 
                                            variable.name = "Nstate", 
                                            value.name = "percentage") %>% 
  mutate(percentage = as.numeric(percentage))

#succesPercentage_long_sOD$percentage <- round(as.numeric(succesPercentage_long_sOD$percentage),10)


#Percentage of succesfull fits per state per study 
ggplot(succesPercentage_long_sOD, aes(fill = Nstate, x = ds, y = percentage)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) 

#Create workable df's for pp's 
AllNtries_sOD <- bind_rows(NtriesWObs_sOD)
AllNtries_long_sOD <- reshape2::melt(AllNtries_sOD, 
                                 id.vars = c("ID", "n", "ds",  "nvars", "usedState"), 
                                 variable.name = "Nstate", 
                                 value.name = "tries") %>%
  dplyr::mutate(tries = replace_na(tries, 100)) %>%
  filter(n > 20)

#Plot number of tries against number of obs
ggplot(AllNtries_long_sOD, aes(x = n, y = tries)) +
  geom_point(size = 1)  + 
  facet_wrap(~Nstate)


#Create 4 plots for average number of tries against the number of observations 
AllNtries_long_sOD %>%
  group_by(n, Nstate) %>%
  summarise(mean = mean(tries)) %>%
  na.omit() %>%
  ggplot(aes(n , mean)) + 
  geom_point() + 
  facet_wrap(~Nstate)


#Create plots with percentage of successful fits against the number of observations. 
p_obsVSsuccesfit_sOD <- AllNtries_long_sOD %>%
  mutate(fit = as.integer(tries != 100)) %>%
  filter(n < 190) %>%
  group_by(n, Nstate) %>%
  summarise(mean = mean(fit)) %>%
  na.omit() %>%
  ggplot(aes(n , mean)) + 
  geom_point() + 
  geom_smooth(method = "gam", se = FALSE) +
  facet_wrap(~Nstate, nrow = 1)
p_obsVSsuccesfit_sOD

#saveRDS(AllNtries_long_sOD, "Scripts/N=1_pipeline/Pres050924/AllNtries_long.rds")


#Plots that shows the average number of selected states of every study against the average number of observations and number of variables
p_stateVSnvars_sOD <- AllNtries_long_sOD %>%
  na.omit() %>% 
  group_by(ds, nvars) %>%
  summarise(AvState = mean(usedState), 
            AvN = mean(n)) %>%
  #filter(nvars < 20) %>%
  ggplot(aes(x = nvars, y = AvState)) +
  geom_point()

p_stateVSAvN_sOD <- AllNtries_long_sOD %>%
  na.omit() %>% 
  group_by(ds, nvars) %>%
  summarise(AvState = mean(usedState), 
            AvN = mean(n)) %>%
  #filter(nvars < 20) %>%
  ggplot(aes(x = AvN, y = AvState)) +
  geom_point() +
  ylab("")

ggpubr::ggarrange(p_stateVSnvars_sOD, p_stateVSAvN_sOD)




##########################################################################################################################
#####################################   Compare Methods   ################################################################
##########################################################################################################################


#Mean Tries 
avTriesDayAndMissing <- AllNtries %>% 
  dplyr::mutate(across(starts_with("S"), ~replace_na(.,100))) %>%
  filter(n > 20) %>%
  select(S2:usedState) %>%
  apply(2, mean, na.rm = TRUE) 
  
avTriesOnlyDay <- AllNtries_sOD %>% 
  dplyr::mutate(across(starts_with("S"), ~replace_na(.,100))) %>%
  filter(n > 20) %>%
  select(S2:usedState) %>%
  apply(2, mean, na.rm = TRUE) 

rbind(DayPlusMissing = avTriesDayAndMissing, OnlyDay = avTriesOnlyDay)

#Percentage fit 

fitRateDayAndMissing <- AllNtries %>% 
  dplyr::mutate(across(starts_with("S"), ~replace_na(.,100))) %>%
  mutate(S2 = as.integer(S2 != 100),
         S3 = as.integer(S3 != 100),
         S4 = as.integer(S4 != 100),
         S5 = as.integer(S5 != 100),) %>%
  filter(n > 20) %>%
  select(S2:usedState) %>%
  apply(2, mean, na.rm = TRUE) 

fitRateOnlyDay <- AllNtries_sOD %>% 
  dplyr::mutate(across(starts_with("S"), ~replace_na(.,100))) %>%
  mutate(S2 = as.integer(S2 != 100),
         S3 = as.integer(S3 != 100),
         S4 = as.integer(S4 != 100),
         S5 = as.integer(S5 != 100),) %>%
  filter(n > 20) %>%
  select(S2:usedState) %>%
  apply(2, mean, na.rm = TRUE) 

round(rbind(DayPlusMissing = fitRateDayAndMissing, OnlyDay = fitRateOnlyDay),3) 


#RMSE



RMSE_sDM <- list()
for(f in list.files(paste0(dayAndMissingPath, "RMSE/"))){
  RMSE_sDM[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0(dayAndMissingPath, "RMSE/", f))
}
RMSE_sOD <- list()
for(f in list.files(paste0(onlyDayPath, "RMSE/"))){
  RMSE_sOD[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0(onlyDayPath, "RMSE/", f))
}


hmmCompRes <- data.frame(
  matrix(NA, length(names(RMSE_sOD)), 7, 
         dimnames = list(NULL,c("ds", "rateOnlyDaybeterFit", "rateDayAndMissingbeterFit", "rateEqualRMSE", "meanRMSEdiff" ,"meanRMSEDayAndMissing", "meanRSMEOnlyDay")))
)
hmmCompRes$ds <- names(RMSE_sOD)


for(ds in names(RMSE_sDM)){
  
  fittedPP <- intersect(na.omit(RMSE_sDM[[ds]])$ID, na.omit(RMSE_sOD[[ds]])$ID)
  #na.omit(RMSE_sDM[[ds]])$ID
  
  #Day & Missing minus Only day:
  sDMminsOD <- RMSE_sDM[[ds]][RMSE_sDM[[ds]]$ID %in% fittedPP,-c(1,2)] - RMSE_sOD[[ds]][RMSE_sOD[[ds]]$ID %in% fittedPP,-c(1,2)]
  
  hmmCompRes[hmmCompRes$ds == ds, "meanRMSEdiff"] <- mean(apply(sDMminsOD, 2, mean))
 
  #For how many pp only-day is higher than day and missing 
  hmmCompRes[hmmCompRes$ds == ds, "rateOnlyDaybeterFit"] <- sum(apply(sDMminsOD, 1, mean) > 0)/nrow(sDMminsOD)
  hmmCompRes[hmmCompRes$ds == ds, "rateDayAndMissingbeterFit"] <- sum(apply(sDMminsOD, 1, mean) < 0)/nrow(sDMminsOD)
  hmmCompRes[hmmCompRes$ds == ds, "rateEqualRMSE"] <- sum(apply(sDMminsOD, 1, mean) == 0)/nrow(sDMminsOD)
  
  hmmCompRes[hmmCompRes$ds == ds, "meanRMSEDayAndMissing"] <- mean(apply(RMSE_sDM[[ds]][RMSE_sDM[[ds]]$ID %in% fittedPP,-c(1,2)], 2, mean))
  hmmCompRes[hmmCompRes$ds == ds, "meanRSMEOnlyDay"] <- mean(apply(RMSE_sOD[[ds]][RMSE_sOD[[ds]]$ID %in% fittedPP,-c(1,2)], 2, mean))

}

hmmCompRes %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  saveRDS("Scripts/N=1_pipeline/Pres050924/hmmCompRes.rds")


p_hmmFitComp <- hmmCompRes  %>% 
  select(1:4) %>%
  reshape2::melt(id.vars = "ds",
                 variable.name = "fit",
                 value.name = "rate") %>%
  ggplot( aes(fill = fit, x = ds, y = rate)) + 
  geom_bar( stat="identity") +
  geom_hline(yintercept = .5, linetype="dashed", color = "red") +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("HMM") + 
  theme(legend.position = "bottom")


p <- ggpubr::ggarrange(p_hmmFitComp, p_varFitComp, nrow = 2)

ggsave("Output/Plots/HmmVarApendixAfitComp.png", p, width = 10, height = 15)

