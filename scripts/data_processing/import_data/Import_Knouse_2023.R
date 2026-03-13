

#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Knouse_2023_raw.csv")

#Remove many of the irrelevant columns 
varSelect <- c("ID", "DAY", "SIG","sad", "lonely", "jittery", "nervous", "angry", "hostile", "happy", "excited" )
dataIn <- dataIn[,varSelect]


#Function to get the day and/or beep with rows for missing data
dateToDayBeep <- function(x){
  n_days <- 6
  n_beeps <- 3
  n_observations <- n_days * n_beeps
  
  
  #Get observation ÌD
  x$observationID <- (x$DAY - 1)*n_beeps+x$SIG
  x$missing <- 0
  
  if(nrow(x) < n_observations){
    #Create empty dataframe for missing values 
    missingData <- matrix(NA,nrow = n_observations - nrow(x), ncol = ncol(x), dimnames = list(NULL, colnames(x))) %>% 
      as.data.frame()
    
    #Fill in the relevant information for the missing data
    missingData$ID <- x$ID[1]
    missingData$observationID <- setdiff( 1:n_observations, x$observationID)
    missingData$DAY <- rep(1:n_days, each = n_beeps)[missingData$observationID]
    missingData$SIG <-  rep(1:n_beeps, n_days)[missingData$observationID]
    missingData$missing <- 1
    
    #Combine objects
    out <- rbind(x, missingData)
  } else {
    out <- x
  }
  
  
  #Sort data create output object
  out <- out[order(out$observationID),]
  row.names(out) <- NULL
  
  out
}

    
#Selected variables and new labels 
varSelect <- c("ID", "DAY", "SIG", "observationID", "missing", "sad", "lonely", "jittery", "nervous", "angry", "hostile", "happy", "excited" )
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Sad", "Lonely", "Jittery", "Nervous", "Angry", "Hostile", "Happy", "Excited")
Emotions <- c("Sad", "Lonely", "Jittery", "Nervous", "Angry", "Hostile", "Happy", "Excited")

#Get celan data
data <- dataIn %>% 
  split(dataIn$ID) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
rownames(data) <- 1:nrow(data)

#Write clean data
write.csv(data, "Data/Clean/Knouse_2023.csv")

