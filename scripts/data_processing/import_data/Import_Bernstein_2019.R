

#Read excel file
dataIn <- readxl::read_xlsx("Data/FromOriginalPapers/Bernstein_2019_raw.xlsx")

#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  
  #Remove glitched rows without a time points
  x <- x %>% 
    filter(!is.na(Timepoint))
  
  #Get beep 
  x$beep <- x$Day %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
 

  x
}

varSelect <- c("ID", "Day", "beep", "Timepoint", "Answered", "Cheer", "Fear", "Content", "Sad")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Cheerful", "Anxiety", "Content", "Sad")

#colnames(dataIn)

#x <- split(dataIn, dataIn$ID)[[1]]


data <- dataIn %>% 
  split(dataIn$ID) %>% #Split based on ID
  lapply(dateToDayBeep) %>%
  rlist::list.rbind() %>%
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns


#Write data as csv
write.csv(data, file = "Data/Clean/Bernstein_2019.csv")


#Vector of emotions in study
Emotions <- c("Cheerful", "Anxiety", "Content", "Sad")


#createMetaDataFile(vars = Emotions, obsIDvars = c("ID", "day", "beep", "observationID", "missing"), scale = "0-100", 
