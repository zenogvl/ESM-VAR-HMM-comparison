

#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Dejonckheere_2017_raw.csv")

#Create a tamplate for the complete dataset to fill everything to indicate missing data. 
complateDataSetTemplate <- data.frame(person_id = rep(unique(dataIn$person_id),each = 30 ),
                                      obs_id = rep(0:29, length(unique(dataIn$person_id))))

#Selected variables and new labels 
varSelect <- c("ID", "obs_id", "missing", "angry_m", "stressed_m", "anxious_m", "sad_m", "happy_m", "relaxed_m")
newVarLabels <- c("ID", "day", "missing", "Angry", "Stressed", "Anxious", "Sad", "Happy", "Relexed" )
Emotions <- c( "Angry", "Stressed", "Anxious", "Sad", "Happy", "Relexed" )

data <- dataIn %>% 
  full_join(complateDataSetTemplate) %>% 
  mutate(ID = match(person_id, unique(person_id), nomatch = NA)) %>% #Change ID label
  arrange(ID, obs_id) %>% #Sort on ID and day number 
  mutate(missing = case_when(is.na(daynr) ~ 1, !is.na(daynr) ~ 0,)) %>% #Add missing data indicator
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns


#Write clean data
write.csv(data, "Data/Clean/Dejonckheere_2017.csv")


# 
# createMetaDataFile(vars = Emotions, obsIDvars = newVarLabels[1:3], scale = "1-7", 
#                    author = "Dejonckheere", year = "2017", data = data)
