

#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Koval_2019_ds1_raw.csv")


#Change ID label
dataIn$ID <- match(dataIn$person_id, unique(dataIn$person_id), nomatch = NA)

#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  x$Day <- x$daynr + 1
  
  x$beep <- x$Day %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  x
}



#Selected variables and new labels 
varSelect <- c("ID", "Day", "beep", "obs_id", "missed", "angry_m", "afraid_m", "anxious_m", "guilty_m", "embarrassed_m", "ashamed_m","sad_m", "strong_m", "happy_m", "confident_m", "clarity_m")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Angry", "Afraid", "Anxious", "Guilty", "Embarrassed", "Ashamed", "Sad", "Strong", "Happy", "Confident", "Clarity")
Emotions <- c( "Angry", "Afraid", "Anxious", "Guilty", "Embarrassed", "Ashamed", "Sad", "Strong", "Happy", "Confident", "Clarity")


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
write.csv(data, "Data/Clean/Koval_2019_ds1.csv")




