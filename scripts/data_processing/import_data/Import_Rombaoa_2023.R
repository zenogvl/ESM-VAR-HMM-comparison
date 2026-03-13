
#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Rombaoa_2023.csv")

#Change ID label
#dataIn$ID <- match(dataIn$XXXX, unique(dataIn$XXXXX), nomatch = NA)


#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  n_days <- 7
  n_beeps <- 4
  n_observations <- n_days * n_beeps
  
  
  #Get Date format 
  x$Date_new <- as.Date(x$StartTime, "%m/%d/%y")
  
  #Get day count 
  x$Day_new <- as.integer(x$Date_new - min(x$Date_new) + 1)

  #Get observation ÌD
  x$observationID <- (x$Day_new - 1)*n_beeps+x$TimePoint
  x$missing <- 0
  
  
  #Remove dubbel beeps 
  if(any(table(x$observationID) > 1)){
    x <- x[-(1:n_observations)[table(x$observationID) > 1],]  
  }
  
   
 
  if(nrow(x) < n_observations){
    #Create empty dataframe for missing values 
    missingData <- matrix(NA,nrow = n_observations - nrow(x), ncol = ncol(x), dimnames = list(NULL, colnames(x))) %>% 
      as.data.frame()
    
    #Fill in the relevant information for the missing data
    missingData$PID <- x$PID[1]
    missingData$observationID <- setdiff( 1:n_observations, x$observationID)
    missingData$Day_new <- rep(1:n_days, each = n_beeps)[missingData$observationID]
    missingData$TimePoint  <-  rep(1:n_beeps, n_days)[missingData$observationID]
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
varSelect <- c("PID", "Day_new", "TimePoint", "observationID", "missing", "Lively", "FullOfEnergy", "Energetic", "Happy", "Pleased", "Helpful", "AtEase", "Calm", "Relaxed", "Sad", "Depressed", "Unhappy", "OnEdge", "Nervous", "Tense", "Hostile", "Resentful", "Angry" )
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Lively", "FullOfEnergy", "Energetic", "Happy", "Pleased", "Helpful", "AtEase", "Calm", "Relaxed", "Sad", "Depressed", "Unhappy", "OnEdge", "Nervous", "Tense", "Hostile", "Resentful", "Angry" )
Emotions <- c("Lively", "FullOfEnergy", "Energetic", "Happy", "Pleased", "Helpful", "AtEase", "Calm", "Relaxed", "Sad", "Depressed", "Unhappy", "OnEdge", "Nervous", "Tense", "Hostile", "Resentful", "Angry" )

#Get clean data
data <- dataIn %>% 
  split(dataIn$PID) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
rownames(data) <- 1:nrow(data)


#Write clean data
write.csv(data, "Data/Clean/Rombaoa_2023.csv")


