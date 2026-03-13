
#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/vanWoerkom_2022_raw.csv")

#Selected variables and new labels 
varSelect <- c("id", "day.nr", "beep.nr", "notice.nr", "missing.beep", "pa.1", "pa.2", "pa.3", "na.1", "na.2", "na.3", "ph.1", "ph.2")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Cheerful", "Satisfied", "Happy", "Insecure", "Anxious", "Down", "Tired", "Energetic")
Emotions <- c("Cheerful", "Satisfied", "Happy", "Insecure", "Anxious", "Down", "Tired", "Energetic")


#Get clean data
data <- dataIn %>% 
   select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Write clean data
write.csv(data, "Data/Clean/vanWoerkom_2022.csv")



