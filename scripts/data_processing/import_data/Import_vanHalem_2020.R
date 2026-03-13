
#Read data 
dataIn <- readxl::read_xlsx("Data/FromOriginalPapers/vanHalem_2020_raw.xlsx")

#Change ID label
#dataIn$ID <- match(dataIn$XXXX, unique(dataIn$XXXXX), nomatch = NA)

#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  
  #Change to POSIXct format
  x$DateTime <- as.POSIXct(x$Trigger_time)
  
  #Sort by time
  out <- x[order(x$DateTime),]
  
  #Get day counter 
  out$Day <- as.integer(as.Date(out$Trigger_date) - min(as.Date(out$Trigger_date)) + 1)
  
  out$beep <- out$Day %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  out$observationID <- 1:nrow(out)
  
  out
}

#Selected variables and new labels 
varSelect <- c("Participant",  "Day", "beep", "observationID", "missed", 
               "Enthousiast","Ontspannen", "Tevreden", "Prikkelbaar", "Energiek","Kalm","Vrolijk",
               "Geirriteerd","Verveeld","Nerveus","Verdrietig","Boos","Somber", "Onzeker",
               "Angstig","Gelukkig","Bezorgd", "Gestrest")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing",
                  "Enthusiastic", "Relaxed", "Satisfied", "Irritable", "Energetic", "Calm", "Cheerful",
                  "Irritated", "Bored", "Nervous", "Sad", "Angry", "Gloomy",  "Insecure", 
                  "Anxious", "Happy", "Worried", "Stressed")
Emotions <- c("Enthusiastic", "Relaxed", "Satisfied", "Irritable", "Energetic", "Calm", "Cheerful",
              "Irritated", "Bored", "Nervous", "Sad", "Angry", "Gloomy",  "Insecure", 
              "Anxious", "Happy", "Worried", "Stressed")


#Get clean data
data <- dataIn %>% 
  mutate(missed = ifelse(is.na(Missing), 0, 1)) %>% #Change the labels for missed surveys from 1 to 0 and 0 to 1. 
  split(dataIn$Participant) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
#rownames(data) <- 1:nrow(data)

#Write clean data
write.csv(data, "Data/Clean/vanHalem_2020.csv")

