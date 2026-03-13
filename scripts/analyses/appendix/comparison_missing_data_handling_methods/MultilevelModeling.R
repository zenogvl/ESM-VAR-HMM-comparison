

install.packages(setdiff(c("tidyverse", "car", "broom", "stats", "gtools", "broom.mixed", "ggeffects", "lme4", "lmerTest", "MLMusingR", "modelsummary"), rownames(installed.packages())))  
library(tidyverse)
theme_set(theme_classic())
options(scipen = 999)

# library(lme4)
# library(lmerTest)
# library(car)
?modelsummary::modelsummary

RmseHmmPath <- "Output/N1_Results/RescaledResults/HMM_onlyday_split_results/RMSE/"
RmseVarPath <- "Output/N1_Results/RescaledResults/VAR_kalman_filter_RMSE/"

# 
# f <- "Rombaoa_2023"
# f <- ds
# f <- gsub("RMSE_|.csv", "", list.files(RmseHmmPath))[3]


variablesList <- readRDS("Data/variables_list.rds")
studiesInfo <- read.csv("studiesInfo2.csv")
studiesInfo <- studiesInfo %>% 
  mutate(scaleDiff = ScaleMax - ScaleMin)

multiLevelRMSEdata <- data.frame()

for(f in gsub("RMSE_|.csv", "", list.files(RmseHmmPath))){
  print(f)
  
  HMM <-  read.csv(paste0(RmseHmmPath, "RMSE_", f, ".csv"))
  HMM$HMMmeanRMSE <- apply(HMM[,variablesList[[f]]] , 1, mean)
  
  VAR <- read.csv(paste0(RmseVarPath, "RMSE_", f, ".csv"))
  VAR$VARmeanRMSE <- apply(VAR[,variablesList[[f]]] , 1, mean)
  
  meanRMSE <- inner_join( na.omit(HMM), na.omit(VAR), by = "ID") %>%
    select(ID, HMMmeanRMSE, VARmeanRMSE)
 
  nPP <- meanRMSE %>% 
    pull(ID) %>%
    n_distinct()
  
  diffRMSE <- read.csv(paste0("Data/CleanRescaled/Rescaled_", f, ".csv")) %>% 
    na.omit() %>%
    count(ID) %>% 
    inner_join(meanRMSE, by = "ID") %>%
    mutate(RMSEdifference = VARmeanRMSE - HMMmeanRMSE,
           dataset = rep(f, nPP),
           nvars = rep(studiesInfo[studiesInfo$datasets == f, "nvars"],nPP),
           population = rep(studiesInfo[studiesInfo$datasets == f, "Population"],nPP), 
           nBeeps = rep(studiesInfo[studiesInfo$datasets == f, "nBeeps"],nPP),
           scale = rep(studiesInfo[studiesInfo$datasets == f, "scaleDiff"],nPP)) 

  multiLevelRMSEdata <- rbind(multiLevelRMSEdata, diffRMSE)

  
}  

multiLevelRMSEdata


#Proportion pp with lower RMSE for HMM. 
sum(multiLevelRMSEdata$RMSEdifference > 0)/nrow(multiLevelRMSEdata)


#Plot to compare the mean RMSE value 
meanRMSEcomparison <- multiLevelRMSEdata %>%
  reshape2::melt(id.vars = c("ID", "n", "RMSEdifference", "dataset", "nvars", "population", "nBeeps", "scale"),
                 variable.name = "Model",
                 value.name = "RMSE") %>%
  mutate(Model = gsub("meanRMSE", "", Model)) %>%
  filter(RMSE < 1) 

orderMeanRmsePerStudy <- meanRMSEcomparison %>%
  group_by(dataset) %>%
  dplyr::summarize(mean = mean(RMSE)) %>%
  arrange(desc(mean))

factor(meanRMSEcomparison$dataset, as.vector(orderMeanRmsePerStudy[,1])$dataset)

(p <- meanRMSEcomparison %>% 
  ggplot(aes(x = factor(dataset, as.vector(orderMeanRmsePerStudy[,1])$dataset), 
             y = RMSE,
             col = Model))  + 
  geom_point(position = position_dodge(0.5), size = 1, alpha = .3) +
  theme(axis.text.x = element_text(angle=90))  +
  xlab("Study") 
)  

ggsave("Output/Plots/meanRmseIndividuals.png", p, width = 10, height = 8)



#Get descriptives of RMSE 
meanRMSEcomparison %>%
  group_by(Model) %>%
  dplyr::summarize(mean = mean(RMSE),
                   sd = sd(RMSE))

t.test(RMSE ~ Model, meanRMSEcomparison)


#Mean RMSE per dataset
multiLevelRMSEdata %>%
  group_by(dataset) %>%
  summarize(mean(RMSEdifference))



######## Multilevel models ########################



#### Plots for ML: 

(p <- multiLevelRMSEdata %>%
  filter(HMMmeanRMSE < 1 & VARmeanRMSE < 1) %>%
  filter(RMSEdifference < 0.3) %>%
  ggplot(aes(x = n, y = RMSEdifference)) +
  geom_point(shape = 20) +
  geom_smooth(aes(color = dataset), method="lm", se=FALSE) + 
  theme(legend.position = "none")
)

ggsave("Output/Plots/RmseDifferanceVsN.png", p, width = 8, height = 6)

  
  
###   Single level model

#Centering variables 
multiLevelRMSEdata <- multiLevelRMSEdata %>% 
  mutate(nCenter = scale(n, center = TRUE, scale = FALSE),
         nVarsCenter = scale(nvars, center = TRUE, scale = FALSE),
         nBeepsCenter = scale(nBeeps, center = TRUE, scale = FALSE),
         scaleCenter = scale(scale, center = TRUE, scale = FALSE),
         populationFactor = factor(population, c("General", "Students", "Clinical"))
  )


mod0 <- lm(RMSEdifference ~1, data = multiLevelRMSEdata)
summary(mod0)


######  Random intercept model
m <- mod2

getVPC <- function(m){
  sigma2_u <- broom.mixed::tidy(m) %>%
    filter(effect == "ran_pars", group == "dataset") %>%
    pull(estimate) %>%
    .^2
  
  sigma2_e <- broom.mixed::tidy(m) %>%
    filter(effect == "ran_pars", group == "Residual") %>%
    pull(estimate) %>%
    .^2
  sigma2_u/(sigma2_u + sigma2_e)
  
}


mod1 <- lmerTest::lmer(RMSEdifference ~ 1 + (1 | dataset), data =  multiLevelRMSEdata)
summary(mod1)

getVPC(mod1)

anova(mod1,mod0)

broom.mixed::tidy(mod1)

#####     Random intercept with variables 




mod2 <- lmerTest::lmer(RMSEdifference ~ 1 + 
                         nCenter + nVarsCenter + populationFactor + nBeepsCenter + scaleCenter +
                         (1 | dataset), data =  multiLevelRMSEdata)


mod2
summary(mod2)
getVPC(mod2)
anova(mod2,mod1)

#Random slope models 

mod3 <-  lmerTest::lmer(RMSEdifference ~ 1 + 
                          nCenter + nVarsCenter + populationFactor + nBeepsCenter + scaleCenter  +
                          (1 + n | dataset), data =  multiLevelRMSEdata)
summary(mod3)
anova(mod3,mod1)


broom.mixed::tidy(mod3) %>% 
  filter(effect == "fixed")

broom.mixed::tidy(mod2) %>% 
  filter(effect == "ran_pars")


#Plots: 
multiLevelRMSEdata %>%
  ggplot(aes(y = RMSEdifference, x = n)) +
  geom_point()



multiLevelRMSEdata %>%
  filter(RMSEdifference > 10)



modelsummary::modelsummary(list(mod1,mod2,mod3), stars = TRUE, gof_omit = 'RMSE|IC|Obs',
                           coef_omit = "SD", output = "latex_tabular")


modelsummary(list("Unconditional" = nullm,
                  "Level-1 variables only" = lev1only,
                  "Full model" = fullm), stars = TRUE, gof_omit = 'RMSE|IC|Obs',
             coef_omit = "SD", title = 'Multilevel Regression Model Results
for Student Engagement.')


