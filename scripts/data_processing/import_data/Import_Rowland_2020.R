



#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Rowland_2020_raw.csv")


#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  x$observationID <- 1:nrow(x)
  x
}
       

#Selected variables and new labels 
varSelect <- c("subjno", "dayno", "beep", "observationID", "missed", "emo1_m", "emo2_m", "emo3_m", "emo4_m", "emo5_m", "emo6_m", "emo7_m", "emo8_m")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Happy", "Excited", "Relaxed", "Satisfied", "Angry", "Anxious", "Depressed", "Sad")
Emotions <- c("Happy", "Excited", "Relaxed", "Satisfied", "Angry", "Anxious", "Depressed", "Sad")

#Get clean data
data <- dataIn %>% 
  mutate(missed = ifelse(is.na(emo8_m), 1, 0)) %>% #Change the labels for missed surveys from 1 to 0 and 0 to 1. 
  split(dataIn$subjno) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
rownames(data) <- 1:nrow(data)

#Write clean data
write.csv(data, "Data/Clean/Rowland_2020.csv")



