
#Load data
dataIn <- haven::read_sav("Data/FromOriginalPapers/Koval_2015_raw.sav")



#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  #Get beep 
  x$beep <- x$beepday %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  x
}
#Create strings to select relevant columns and to rename them  them 
varSelect <- c("ppid", "beepday", "beep", "sad", "angry", "anx", "depr", "stress", "happy", "relax", "excit")
newVarLabels <- c("ID", "day", "beep", "Sad", "Angry", "Anxiety", "Depressed", "Stressed", "Happy", "Relexed", "Excited")
Emotions <- c("Sad", "Angry", "Anxiety", "Depressed", "Stressed", "Happy", "Relexed", "Excited")

#create day to beep for every pp and rbind them
data <- dataIn %>% 
  split(dataIn$ppid) %>%  #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep 
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Create a vector for missing data
data$missing <- rep(0, nrow(data))
data$missing[apply(is.na(data[,Emotions]),1,any)] <- 1

#Print data 
write.csv(data, "Data/Clean/Koval_2015.csv")



