
#Read data 
dataIn <- haven::read_sav("Data/FromOriginalPapers/Stevens_2020_raw.sav")


#Create a tamplate for the complete dataset to fill everything to indicate missing data. 
completeDataSetTemplate <- data.frame(ParticipantID  = rep(unique(dataIn$ParticipantID ),each = 42 ),
                                      NotificationNo  = rep(1:42, length(unique(dataIn$ParticipantID))),
                                      day = rep(rep(1:7, each = 6), length(unique(dataIn$ParticipantID))),
                                      beep = rep(rep(1:6, 7), length(unique(dataIn$ParticipantID)))
                                      )

#Selected variables and new labels 
varSelect <-  c("ID", "day", "beep", "NotificationNo", "Responded_fixed", "Anxious", "Guilty", "Depressed", "Confident", "Happy", "Inspired")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Anxious", "Guilty", "Depressed", "Confident", "Happy", "Inspired")
Emotions <-  c("Anxious", "Guilty", "Depressed", "Confident", "Happy", "Inspired")

data <- dataIn %>% 
  as.data.frame() %>%
  filter(!is.na(NotificationNo)) %>%
  full_join(completeDataSetTemplate) %>% 
  mutate(ID = match(ParticipantID, unique(ParticipantID), nomatch = NA)) %>% #Change ID label
  arrange(ID, NotificationNo) %>%
  filter(!is.na(day)) %>% #Remove some rows with no information 
  slice(-4879) %>% #Remove a buggy duplicate row 
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns


#Write clean data
write.csv(data, "Data/Clean/Stevens_2020.csv")





# 
# createMetaDataFile(vars = Emotions , obsIDvars = newVarLabels, scale = "0-100" , 
#                    author = "Stevens", year = "2020", data = data)
