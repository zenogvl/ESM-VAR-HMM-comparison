

#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Wright_2015_raw.csv", sep = ";")


#Get a vector that indicates what row have (partial) missing data
dataIn$partialMissingData <- dataIn %>% 
  select(-(1:2)) %>%
  is.na() %>%
  rowSums() 

#Get a variable for missing data that is 2 when not all observations are missing
dataIn <- dataIn %>% 
  mutate(missing = case_when(partialMissingData == 0  ~ 0,
                             partialMissingData == 10 ~ 1,
                             partialMissingData > 0 & partialMissingData < 10 ~ 2
                             )) %>%
  select(-partialMissingData)


#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  n_days <- 101
  n_observations <- n_days

  
  if(nrow(x) < n_observations){
    #Create empty dataframe for missing values 
    missingData <- matrix(NA,nrow = n_observations - nrow(x), ncol = ncol(x), dimnames = list(NULL, colnames(x))) %>% 
      as.data.frame()
    
    #Fill in the relevant information for the missing data
    missingData$id <- x$id[1]
    missingData$day <- setdiff( 1:n_observations, x$day)
    missingData$missing <- 1
    
    #Combine objects
    out <- rbind(x, missingData)
  } else {
    out <- x
  }
  
  
  #Sort data create output object
  out <- out[order(out$day),]
  row.names(out) <- NULL
  
  out
  
}

#Selected variables and new labels 
varSelect <- c("id",  "day",  "missing", "afraid", "active", "alert", "nervous", "attentive" , "determined", "hostile", "inspired", "ashamed", "upset")
newVarLabels <- c("ID", "day", "missing",  "Afraid", "Active", "Alert", "Nervous", "Attentive" , "Determined", "Hostile", "Inspired", "Ashamed", "Upset")
Emotions <- c("Afraid", "Active", "Alert", "Nervous", "Attentive" , "Determined", "Hostile", "Inspired", "Ashamed", "Upset")


#Get clean data
data <- dataIn %>% 
  split(dataIn$id) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns



#Write clean data
write.csv(data, "Data/Clean/Wright_2015.csv")




