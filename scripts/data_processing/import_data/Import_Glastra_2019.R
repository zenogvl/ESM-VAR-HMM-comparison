

#Read data 
dataIn <- read.table("Data/FromOriginalPapers/Glastra_2019_raw.txt", sep = ",", header = TRUE)

#Change long format data to wide and do some more things 
dataWide <- dataIn %>% 
  select(c(1, 2, 5, 7, 18, 20)) %>% #Only relevant colums 
  filter(!(Prompt.Label == "Instruct") & Responded == 1) %>% #Remove instruction and missed responses
  pivot_wider(names_from = Prompt.Label, values_from = Response) %>%
  mutate(date = as.Date(Notification.Time), #Get date as variable
         time =  format(as.POSIXct(Notification.Time), format = "%H:%M:%S") ) %>% # Get time as variable
  mutate(beep = case_when(time == "11:00:00" ~ 1, #Get beeps based on when notification has been send
                          time == "16:00:00" ~ 2,
                          time == "21:00:00" ~ 3,)) %>% 
  group_by(Participant.ID) %>%
  mutate(day = as.integer(date - min(date) + 1), #Get a day counter for every pp
         totalObservations = sum(as.integer(Responded) )) %>% #Get the number of observations for every pp
  ungroup() %>%
  filter(totalObservations > 14) #Remove pp with less than 1/3 of total observations (Lowers n from 87 to 67)


#Function to add the missing data  
dateToDayBeep <- function(x){
  n_days <- 14
  n_beeps <- 3
  n_observations <- n_days * n_beeps
  
  #Get observation ÌD
  x$observationID <- (x$day - 1)*n_beeps+x$beep
  x$missing <- 0
  
  #If questionaire is partly completed, change missing to 2
  if(any(is.na(x))){
    x$missing[apply(is.na(x), 1, any)] <- 2
  }
  
  if(nrow(x) < n_observations){
    #Create empty dataframe for missing values 
    missingData <- matrix(NA,nrow = n_observations - nrow(x), ncol = ncol(x), dimnames = list(NULL, colnames(x))) %>% 
      as.data.frame()
    
    #Fill in the relevant information for the missing data
    missingData$Participant.ID <- x$Participant.ID[1]
    missingData$observationID <- setdiff( 1:n_observations, x$observationID)
    missingData$day <- rep(1:n_days, each = n_beeps)[missingData$observationID]
    missingData$beep <-  rep(1:n_beeps, n_days)[missingData$observationID]
    missingData$missing <- 1
    
    #Combine objects
    out <- rbind(x, missingData)
  } else {
    out <- x
  }
  
  
  #Sort data create output object
  out <- out[order(out$observationID),]
  row.names(out) <- NULL
  
  out
  
}



#Selected variables and new labels 
varSelect <- c("Participant.ID", "beep", "day", "observationID", "missing", "Energetic", "Cheerful", "Relaxed", "Satisfied", "Nervous", "Irritated", "Bored", "Down")
Emotions  <- c("Energetic", "Cheerful", "Relaxed", "Satisfied", "Nervous", "Irritated", "Bored", "Down")


data <- dataWide %>%
  split(dataWide$Participant.ID) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename(ID = Participant.ID) #Change the name of the ID column


#Write clean data
write.csv(data, "Data/Clean/Glastra_2019.csv")




