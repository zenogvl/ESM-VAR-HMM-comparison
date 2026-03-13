

createMetaDataFile <- function(vars, recode = NULL, obsIDvars, scale, author, authorShort = NULL, year, data, etal = TRUE, multipleDS = NULL){
  
  
  if(etal) {
    AuthorData <- paste0("Author Data: ", author, " et al. (", year, ") \n\n", "Citation: \n\n")
  } else {
    AuthorData <- paste0("Author Data: ", author, " (", year, ") \n\n", "Citation: \n\n")
  }
  if(is.null(authorShort)){
    labelData <-  paste0("Label in code: ", author, "_", year, 
                         "\nLabel in text: ", author, " ", year, '\n\n' )
  } else {
    labelData <-  paste0("Label in code: ", authorShort, "_", year, 
                         "\nLabel in text: ", author, " ", year, '\n\n' )
  }
  
  
  format <- paste0("Data format: .csv file \n", dim(data)[1], " Rows, ", dim(data)[2], " Columns \n\nVariables:\n")
  
  topVariablesVector <- c("-ID: participant id number", 
                          "-day: Number describing the day of observation",
                          "-beep: Order of observations on a single day",
                          "-observationID: Unique ID corresponding with the order of observations",
                          "-missing: codes missing observation. 0 for recorded, 1 for missing, 2 for partially") 
  topVars <- paste0(topVariablesVector[c("ID", "day", "beep", "observationID", "missing")  %in% obsIDvars  ], collapse = "\n")
  
  
  if(is.null(recode)){
    emotionVars <- paste0("-", vars, ": Observation on a scale ranging from ", scale, collapse = "\n")
  } else {
    emotionVars <- paste0("-", vars, ": Observation on a scale ranging from ", scale, ". Recoded from: ", recode, collapse = "\n")
  }
  
  
  if(is.null(authorShort)){
    cat(paste0(AuthorData,labelData,format, topVars, "\n",emotionVars), file = paste0("Data/MetaDataCleanData/", "MetaData_", author, "_", year, ".txt"))
  } else {
    cat(paste0(AuthorData,labelData,format, topVars, "\n",emotionVars), file = paste0("Data/MetaDataCleanData/", "MetaData_",authorShort, "_", year, ".txt"))
  }
  
}

