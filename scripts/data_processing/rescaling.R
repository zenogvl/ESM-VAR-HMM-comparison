
cleanFiles <- list.files("Data/Clean/")
variablesList <- readRDS("Data/variables_list.rds")




datasets <- gsub(".csv", "", cleanFiles)

variablesList[datasets[1]]

ds <- "Marian_2023"


N <- length(datasets)
studyInfo <- data.frame(datasets = datasets, 
                        nvars = rep(NA, N),
                        N = rep(NA, N), 
                        ScaleMin = rep(NA, N),
                        ScaleMax = rep(NA, N),
                        nDays= rep(NA, N), 
                        nBeeps= rep(1, N)
                        )

for(ds in datasets){
  data <- read.csv(paste0("Data/Clean/", ds, ".csv"))
  studyInfo$nvars[studyInfo$datasets == ds] <- length(variablesList[[ds]])
  studyInfo$N[studyInfo$datasets == ds] <- length(unique(data$ID))
  studyInfo$ScaleMin[studyInfo$datasets == ds] <- min(data[,variablesList[[ds]]], na.rm = TRUE)
  studyInfo$ScaleMax[studyInfo$datasets == ds] <- max(data[,variablesList[[ds]]], na.rm = TRUE)
  
  if("beep" %in% colnames(data)){
    studyInfo$nBeeps[studyInfo$datasets == ds] <- max(data$beep)
  }
  
  rescaleMin <- studyInfo$ScaleMin[studyInfo$datasets == ds]
  rescaleDenominator <- (studyInfo$ScaleMax[studyInfo$datasets == ds] - studyInfo$ScaleMin[studyInfo$datasets == ds])
  
  data[,variablesList[[ds]]] <- (data[,variablesList[[ds]]] - rescaleMin)/rescaleDenominator
  
  write.csv(data, paste0("Data/CleanRescaled/", "Rescaled_", ds, ".csv"))
}
data$ID




write.csv(studyInfo, "studiesInfo3.csv")



