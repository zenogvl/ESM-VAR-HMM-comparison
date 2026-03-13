
#Read data 
dataIn <- haven::read_sav("Data/FromOriginalPapers/Kalokerinos_2019_raw.sav")

#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  x$Day <- x$beepday - min(x$beepday) + 1
  
  x$beep <- x$Day %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  x
}

#Selected variables and new labels 
varSelect <- c("Participant", "Day", "beep", "beepnum" , "beepdone", "emotion_sad", "emotion_angry", "emotion_disapp", "emotion_ashamed",
               "emotion_anxious", "emotion_stressed", "emotion_proud", "emotion_happy", "emotion_content", "emotion_relief")
newVarLabels <- c("ID", "day", "beep", "obsNumber", "missed", "Sad", "Angry", "Disappointed", "Ashamed",
                  "Anxious", "Stressed", "Proud", "Happy", "Content", "Relieved")
Emotions <- c("Sad", "Angry", "Disappointed", "Ashamed",
              "Anxious", "Stressed", "Proud", "Happy", "Content", "Relieved")

#Get celan data
data <- dataIn %>% 
  mutate(missed = ifelse(beepdone == 1, 0, 1)) %>% #Change the labels for missed surveys from 1 to 0 and 0 to 1. 
  split(dataIn$Participant) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns


#Write clean data
write.csv(data, "Data/Clean/Kalokerinos_2019.csv")



