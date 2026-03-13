
dataIn <- haven::read_sav("Data/FromOriginalPapers/Grommisch_2020_raw.sav")



#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  #Get beep 
  x$beep <- x$DayNr %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  x
}

#Create strings to select relevant columns and to rename them  them 
varSelect <- c("ppid", "DayNr", "beep", "HAP", "RLX", "CONF", "SAD", "stress", "angry")
newVarLabels <- c("ID", "day", "beep", "Happy", "Relaxed", "Confident", "Sad", "Stressed", "Angry")

#create day to beep for every pp and rbind them
data <- dataIn %>% 
  split(dataIn$ppid) %>%  #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep 
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns



#Write data as csv
write.csv(data, file = "Data/Clean/Grommisch_2020.csv")


Emotions <- c( "Happy", "Relaxed", "Sad", "Stressed", "Angry")