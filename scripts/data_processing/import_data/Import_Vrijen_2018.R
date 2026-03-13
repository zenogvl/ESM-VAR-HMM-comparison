
#Read .dta data using haven
dataIn <- haven::read_dta("Data/FromOriginalPapers/Vrijen_2018_raw.dta")

#ID with problems with time of day indicators
dataIn[which(rowSums(dataIn[, 8:10])==0),] %>%
  select(ID) %>% 
  unique() %>%
  unlist()

#Remove ID's with time of day problems 
data <- dataIn %>% 
  filter(!row_number() %in% which(rowSums(dataIn[, 8:10])==0)) %>%
  mutate(missed = ifelse(is.na(INTraw), 1, 0)) %>% #create variable for missing values 
  #Select relevant cols
  select(ID, Day, day_part, Time, missed, INTraw:WORraw) 

#Rename columns 
colnames(data) <- c("ID", "day", "beep", "observationID", "missing",
                    "Interested", "Joy", "Sad", "Irritated", "Worried")
Emotions <- c("Interested", "Joy", "Sad", "Irritated", "Worried")

write.csv(data, file = "Data/Clean/Vrijen_2018.csv")

