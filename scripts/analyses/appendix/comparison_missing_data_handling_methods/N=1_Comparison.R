

RMSE_VAR <- list()
for(f in list.files("Output/N1_Results/RescaledResults/VAR_kalman_filter_RMSE/")){
  RMSE_VAR[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0("Output/N1_Results/RescaledResults/VAR_kalman_filter_RMSE/", f))
}

RMSE_HMM <- list()
for(f in list.files("Output/N1_Results/RescaledResults/HMM_kalman_filter_results/RMSE/")){
  RMSE_HMM[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0("Output/N1_Results/RescaledResults/HMM_onlyday_split_results/RMSE/", f))
}




CompRes <- data.frame(
  matrix(NA, length(names(RMSE_HMM)), 7, 
         dimnames = list(NULL,c("ds", "HMMfit", "VARfit", "Equalfit", "meanRMSEdiff" ,"meanRMSEVAR", "meanRSMEHMM")))
)
CompRes$ds <- names(RMSE_HMM)

for(ds in names(RMSE_HMM)){
  
  fittedPP <- intersect(na.omit(RMSE_VAR[[ds]])$ID, na.omit(RMSE_HMM[[ds]])$ID)
  #na.omit(RMSE_sDM[[ds]])$ID
  
  #Day & Missing minus Only day:
  sDMminsOD <- RMSE_VAR[[ds]][RMSE_VAR[[ds]]$ID %in% fittedPP,-c(1,2)] - RMSE_HMM[[ds]][RMSE_HMM[[ds]]$ID %in% fittedPP,-c(1,2)]
  
  CompRes[CompRes$ds == ds, "meanRMSEdiff"] <- mean(apply(sDMminsOD, 2, mean))
  
  #For how many pp only-day is higher than day and missing 
  CompRes[CompRes$ds == ds, "HMMfit"] <- sum(apply(sDMminsOD, 1, mean) > 0)/nrow(sDMminsOD)
  CompRes[CompRes$ds == ds, "VARfit"] <- sum(apply(sDMminsOD, 1, mean) < 0)/nrow(sDMminsOD)
  CompRes[CompRes$ds == ds, "Equalfit"] <- sum(apply(sDMminsOD, 1, mean) == 0)/nrow(sDMminsOD)
  
  CompRes[CompRes$ds == ds, "meanRMSEVAR"] <- mean(apply(RMSE_VAR[[ds]][RMSE_VAR[[ds]]$ID %in% fittedPP,-c(1,2)], 2, mean))
  CompRes[CompRes$ds == ds, "meanRSMEHMM"] <- mean(apply(RMSE_HMM[[ds]][RMSE_HMM[[ds]]$ID %in% fittedPP,-c(1,2)], 2, mean))
  
}


CompRes %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  saveRDS("Scripts/N=1_pipeline/Pres050924/CompRes.rds")


(p <- CompRes %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  select(1:3) %>%
  reshape2::melt(id.vars = "ds",
                 variable.name = "fit",
                 value.name = "rate") %>%
  ggplot( aes(fill = fit, x = ds, y = rate)) + 
  geom_bar( stat="identity") +
  #geom_hline(yintercept = .5, linetype="dashed", color = "red") +
  theme(axis.text.x = element_text(angle=90)) 
)


sum(CompRes$HMMfit > .75)/length(CompRes$HMMfit)

ggsave("Output/Plots/VarVsHmmLowerRmseRate.png", p, width = 10, height = 8)

RMSE_VAR[[1]]




# ds_names <- gsub("RMSE_|.csv", "", list.files("Output/N1_VAR_Results_kalfil/"))
# 
# VARfiles <- list.files("Output/N1_VAR_Results_kalfil/")
# HMMfiles <- list.files("Output/N1_HMM_Results_V2_sep_only_day_split/RMSE/")
# 
# length(list.files("Output/N1_VAR_Results_kalfil/")) == length(list.files("Output/N1_HMM_Results_V2_sep_only_day_split/RMSE/"))
# 
# VARRes <- list()
# HMMRes <- list()
# for(i in 1:){
#   VARRes <-  read.csv(paste0("Output/N1_VAR_Results_kalfil/", ))
#   HMMRes <- list()
#   
# }




variablesList <- readRDS("Data/variables_list.rds")


VARmeanRMSELong <- data.frame()
for(ds in names(RMSE_VAR)){
  x <- data.frame(ID = as.integer(RMSE_VAR[[ds]]$ID), 
                                dataset = rep(ds, nrow(RMSE_VAR[[ds]])), 
                                VAR = apply(RMSE_VAR[[ds]][variablesList[[ds]]] , 1, mean))
  VARmeanRMSELong <- rbind(VARmeanRMSELong, x)
  
}
HMMmeanRMSELong <- data.frame()
for(ds in names(RMSE_HMM)){
  x <- data.frame(ID = RMSE_HMM[[ds]]$ID, 
                                dataset = rep(ds, nrow(RMSE_HMM[[ds]])), 
                                HMM = apply(RMSE_HMM[[ds]][variablesList[[ds]]] , 1, mean))
  HMMmeanRMSELong <- rbind(HMMmeanRMSELong, x)
}




meanRMSE <- full_join(HMMmeanRMSELong, VARmeanRMSELong) %>%
  filter(!is.na(HMM) & !is.na(VAR)) %>% 
  mutate(diff = HMM - VAR) %>%
  reshape2::melt(id.vars = c("ID", "dataset"),
                 variable.name = "method", 
                 value.name = "meanRMSE") %>%
  filter(meanRMSE > 0 & meanRMSE < 1)
  

max(meanRMSE$meanRMSE)

meanRMSE %>% 
  ggplot(aes(x = dataset, y = meanRMSE, col = method)) + 
  geom_point(position = position_dodge(0.5), size = 1, alpha = .3) +
  theme(axis.text.x = element_text(angle=90)) 


meanRMSE

###########################################################################################################
###############################################  OLD  #####################################################
###########################################################################################################



files <- list.files("Output/N1_HMM_Results/")

files <- files[!(files == "RMSE_Bringmann_2013.csv")]
files <- files[!(files == "RMSE_Koval_2015.csv")]
f <- files[1]

VARRes <- list()
HMMRes <- list()

#Load in RMSE results
for(f in files){
  dsName <- gsub("RMSE_", "", gsub(".csv", "", f))
  VARRes[[dsName]]  <- read.csv(paste0("Output/N1_VAR_Results/",f))
  HMMRes[[dsName]] <- read.csv(paste0("Output/N1_HMM_Results/",f))
  
}

datasets <- gsub("RMSE_", "", gsub(".csv", "", files))



#Check percentage of succesfull HMM

succesfulHMMfit <- list()

ds <- datasets[1]

succesfulHMMfit <- vector("numeric", length = length(datasets))
names(succesfulHMMfit) <- datasets
for(ds in datasets){
  temp <- HMMRes[[ds]] %>% 
    is.na() %>%
    apply(1, any) %>%
    table() 
  
  succesfulHMMfit[ds] <- temp[1]/sum(temp) * 100
}

VARResNew <- list()
HMMResNew <- list()

ds

#Remove pp without RMSE
for(ds in datasets){
  ID <- HMMRes[[ds]]$ID[!apply(is.na(HMMRes[[ds]]), 1, any)]
  VARResNew[[ds]] <- VARRes[[ds]][VARRes[[ds]]$ID %in% ID,]
  HMMResNew[[ds]] <- HMMRes[[ds]][HMMRes[[ds]]$ID %in% ID,]
  
  ID <- VARResNew[[ds]]$ID[!apply(is.na(VARResNew[[ds]]), 1, any)]
  VARResNew[[ds]] <- VARResNew[[ds]][VARResNew[[ds]]$ID %in% ID,]
  HMMResNew[[ds]] <- HMMResNew[[ds]][HMMResNew[[ds]]$ID %in% ID,]
}




compar <- matrix(NA, nrow = length(datasets), ncol = 5)
colnames(compar) <- c("Percentage", 
                      "PP VAR higher RMSE" , "PP HMM higher RMSE", "PP equal RMSE",
                      "How much RMSE VAR above HMM")
rownames(compar) <- datasets

for(ds in datasets){
  #Check if higher or lower for every variable separte 
  compar[ds,1] <- sum(VARResNew[[ds]][,-(1:2)] > HMMResNew[[ds]][,-(1:2)])/(nrow(VARResNew[[ds]][,-(1:2)]) * ncol(VARResNew[[ds]][,-(1:2)]))
  
  #Check for how many individuals VAR outperforms HMM
  
  compar[ds,2] <- sum(apply(VARResNew[[ds]][,-(1:2)] > HMMResNew[[ds]][,-(1:2)], 1, sum) > ncol(VARResNew[[ds]][,-(1:2)])/2)
  compar[ds,3] <- sum(apply(VARResNew[[ds]][,-(1:2)] < HMMResNew[[ds]][,-(1:2)], 1, sum) > ncol(VARResNew[[ds]][,-(1:2)])/2)
  compar[ds,4] <- sum(apply(VARResNew[[ds]][,-(1:2)] < HMMResNew[[ds]][,-(1:2)], 1, sum) == ncol(VARResNew[[ds]][,-(1:2)])/2)
  
  #Mean difference between VAR and HMM
  compar[ds,5] <- mean(unlist(VARResNew[[ds]][,-(1:2)] - HMMResNew[[ds]][,-(1:2)]))
  
}

compar <- as.data.frame( na.omit(compar))

compar$Study <- gsub("_", " ", rownames(compar))

ggplot(compar, aes(x = Study, y = Percentage)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = .5, linetype="dashed", color = "red") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90)) 






















