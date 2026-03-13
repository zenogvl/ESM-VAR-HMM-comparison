
install.packages(setdiff(c("tidyverse", "hmmr", "stats", "reshape2", "ggpubr"), rownames(installed.packages())))  
library(tidyverse)
theme_set(theme_classic())





f <- "RMSE_koval_2015.csv"

RMSE_LD$koval_2015

RMSE_LD <- list()
for(f in list.files("Output/N1_Results/RescaledResults/VAR_listwise_deletion_RMSE/")){
  RMSE_LD[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0("Output/N1_Results/RescaledResults/VAR_listwise_deletion_RMSE/", f))
}
RMSE_KF <- list()
for(f in list.files("Output/N1_Results/RescaledResults/VAR_kalman_filter_RMSE/")){
  RMSE_KF[[gsub("RMSE_|.csv", "", f)]] <- read.csv(paste0("Output/N1_Results/RescaledResults/VAR_kalman_filter_RMSE/", f))
}




varCompRes <- data.frame(
  matrix(NA, length(names(RMSE_KF)), 7, 
         dimnames = list(NULL,c("ds", "rateKFbeterFit", "rateNothingbeterFit", "rateEqualRMSE", "meanRMSEdiff" ,"meanRMSEnothing", "meanRSMEkalmanFilter")))
)
varCompRes$ds <- names(RMSE_KF)

for(ds in names(RMSE_KF)){
  
  fittedPP <- intersect(na.omit(RMSE_LD[[ds]])$ID, na.omit(RMSE_KF[[ds]])$ID)
  #na.omit(RMSE_sDM[[ds]])$ID
  
  #Nothing minus kalman filter 
  NminKF <- RMSE_LD[[ds]][RMSE_LD[[ds]]$ID %in% fittedPP,-c(1,2)] - RMSE_KF[[ds]][RMSE_KF[[ds]]$ID %in% fittedPP,-c(1,2)]
  
   #Mean RMSE diff per ds 
  varCompRes[varCompRes$ds == ds, "meanRMSEdiff"] <- mean(apply(NminKF, 2, mean))

  #For how many pp only-day is higher than day and missing 
  varCompRes[varCompRes$ds == ds, "rateKFbeterFit"] <- sum(apply(NminKF, 1, mean) > 0)/nrow(NminKF)
  varCompRes[varCompRes$ds == ds, "rateNothingbeterFit"] <- sum(apply(NminKF, 1, mean) < 0)/nrow(NminKF)
  varCompRes[varCompRes$ds == ds, "rateEqualRMSE"] <- sum(apply(NminKF, 1, mean) == 0)/nrow(NminKF)

  varCompRes[varCompRes$ds == ds, "meanRMSEnothing"] <- mean(apply(RMSE_LD[[ds]][RMSE_LD[[ds]]$ID %in% fittedPP,-c(1,2)], 2, mean))
  varCompRes[varCompRes$ds == ds, "meanRSMEkalmanFilter"] <- mean(apply(RMSE_KF[[ds]][RMSE_KF[[ds]]$ID %in% fittedPP,-c(1,2)], 2, mean))

}



colnames(varCompRes) <- c("ds", "KFbetter", "Notbetter", "equal", "diff", "MeanNot", "MeanKF")
varCompRes %>% mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  saveRDS("Scripts/N=1_pipeline/Pres050924/varCompRes.rds")





p_varFitComp <- varCompRes %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
  filter(!ds == "Knouse_2023") %>%
  select(1:4) %>%
  reshape2::melt(id.vars = "ds",
                 variable.name = "fit",
                 value.name = "rate") %>%
  ggplot( aes(fill = fit, x = ds, y = rate)) + 
  geom_bar( stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("VAR") + 
  theme(legend.position = "bottom")
