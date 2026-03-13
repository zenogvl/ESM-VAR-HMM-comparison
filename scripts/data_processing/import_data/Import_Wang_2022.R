
#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Wang_2022_raw.csv")
head(dataIn)


#Selected variables and new labels 
varSelect <- c("ID", "Day", "DailyResponse", "Count", "missed", "Sad", "Anxiety", "Anger")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Sad", "Anxious", "Angry")
Emotions <- c("Sad", "Anxious", "Angry")


#Get clean data
data <- dataIn %>% 
  filter(!(DailyResponse == 7)) %>% #Remove 7th daily observation because it doens't excist 
  filter(!is.na(Count)) %>% #Remove Observations without an obseravtion ID. 
  mutate(missed = ifelse(is.na(Sad), 1, 0)) %>% #Change the labels for missed surveys from 1 to 0 and 0 to 1. 
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns


#Write clean data
write.csv(data, "Data/Clean/Wang_2022.csv")



