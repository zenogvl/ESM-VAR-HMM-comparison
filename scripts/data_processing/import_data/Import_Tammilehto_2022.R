

#Read data 
dataIn <- read.table("Data/FromOriginalPapers/Tammilehto_2022_raw.txt")


#Function to get the day and/or beep
#DTIME describes the time of measurement. There are 14 measurements in a day of which 7 are at night and thus always empty. 
dateToDayBeep <- function(x){
  
  #print(unique(x$KANDI_ID))
  
  #Sort data based on time
  x <- x[order(x$DTIME),]
  
  #### Get day and beep variable ######
  
  #Get the day number by indexing on the known start and end of a day. 
  x$Day <- rep(1:7, each = 14)[x$DTIME]
  
  #Get the beep by indexing on a repeating 1:7 vector. 
  x$beep <- rep(1:7, 14)[x$DTIME]
  
  
  #### Add missing data ####
  n_days <- 7
  n_beeps <- 7
  n_observations <- n_days * n_beeps
  
  
  #Get observation ÌD
  x$observationID <- (x$Day - 1)*n_beeps+x$beep
  x$missing <- 0
  
  if(nrow(x) < n_observations){
    #Create empty dataframe for missing values 
    missingData <- matrix(NA,nrow = n_observations - nrow(x), ncol = ncol(x), dimnames = list(NULL, colnames(x))) %>% 
      as.data.frame()
    
    #Fill in the relevant information for the missing data
    missingData$KANDI_ID <- x$KANDI_ID[1]
    missingData$observationID <- setdiff( 1:n_observations, x$observationID)
    missingData$Day <- rep(1:n_days, each = n_beeps)[missingData$observationID]
    missingData$beep <-  rep(1:n_beeps, n_days)[missingData$observationID]
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
varSelect <- c("KANDI_ID", "Day", "beep", "observationID", "missing", 
               "JOY", "EXCITEME", "PRIDE", "PLEASED", "ANGER", "ANXIETY", "SHAME", "SAD")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", 
                  "Joy", "Exciteme", "Pride", "Pleased", "Anger", "Anxiety", "Shame", "Sad")
Emotions <- c("Joy", "Exciteme", "Pride", "Pleased", "Anger", "Anxiety", "Shame", "Sad")

#Get clean data
data <- dataIn %>% 
  split(dataIn$KANDI_ID) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
rownames(data) <- 1:nrow(data)


#Write clean data
write.csv(data, "Data/Clean/Tammilehto_2022.csv")

