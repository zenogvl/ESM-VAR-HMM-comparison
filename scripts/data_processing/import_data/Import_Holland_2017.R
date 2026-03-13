#Read data
dataIn <- read.csv("Data/FromOriginalPapers/Holland_2017_raw.csv")

#Change ID label
dataIn$ID <- match(dataIn$person_id, unique(dataIn$person_id), nomatch = NA)

#Split into separate df for every individual 
dataIn_split <- split(dataIn, dataIn$person_id)

#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  #Add 1 to day column
  x$day <- x$daynr + 1

  #Get beep 
  x$beep <- x$day %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  x
}

#Create strings to select relevant columns and to rename them  them 
varSelect <- c("ID", "day", "beep", "obs_id", "missed", "angry_m", "anxious_m", "guilty_m", "sad_m", "happy_m", "relaxed_m", "confident_m")
newVarLabels <- c("ID", "day", "beep", "obsNumber", "missed", "Angry", "Anxiety", "Guilty", "Sad", "Happy", "Relaxed", "Confident")

#create day to beep for every pp and rbind them
data <- dataIn %>% 
  split(dataIn$person_id) %>%  #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep 
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Remove row names 
rownames(data) <- NULL

dim(data)
#Print data
write.csv(data, "Data/Clean/Holland_2017.csv")

Emotions <- c("Angry", "Anxiety", "Guilty", "Sad", "Happy", "Relaxed", "Confident")

