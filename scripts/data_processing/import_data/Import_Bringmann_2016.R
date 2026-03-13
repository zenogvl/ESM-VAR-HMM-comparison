
###### Bringmann et al. 2016 ##############

#Read data
dataIn <- read.csv("Data/FromOriginalPapers/Bringmann_2016/Bringmann_2016_ds1_raw.csv", sep = ";", header = TRUE)
idInfoIn <- read.csv("Data/FromOriginalPapers/Bringmann_2016/id_info.csv", sep = ";")


#Change 9999 to NA. 
dataIn[dataIn == 9999] <- NA

#Create observationID column 
dataIn$observationID <- dataIn$ID %>%
  table() %>%
  as.vector() %>%
  sapply(seq, from = 1) %>%
  unlist()


#Create observationID column 
idInfoIn$observationID <- idInfoIn$PpID %>%
  table() %>%
  as.vector() %>%
  sapply(seq, from = 1) %>%
  unlist()

dateToDayBeep <- function(x){
  
  #Create day month dummy variable 
  monthDayDummy <- as.integer(paste0(x$month, x$day))
  
  #Use match to get the day numbering 
  x$day <- match(monthDayDummy, unique(monthDayDummy), nomatch = NA)
  
  x$beep <- x$day %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  x
}

dayBeepObject <- idInfoIn %>% 
  split(idInfoIn$PpID) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  filter(idInfoIn$PpID %in% unique(dataIn$ID)) %>% #Remove people not in the data file. 
  rename(ID = PpID) %>% 
  select(ID, day, beep, observationID)


data <- dataIn %>% 
  full_join(dayBeepObject) %>% 
  select(ID, observationID, day, beep, KWAAD, DEPRE, DROEV, ANGST, ONTSP, BLIJ) 


data$missing <- rep(0, nrow(data))
data$missing[apply(is.na(data), 1, any)] <- 1

varSelect <- c("ID", "day", "beep", "observationID", "missing", "KWAAD", "DEPRE", "DROEV", "ANGST", "ONTSP", "BLIJ")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Angry", "Depressed", "Sad", "Anxious", "Relaxed", "Happy")

data <- data %>%
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns


#Write data
write.csv(data, file = "Data/Clean/Bringmann_2016.csv")


Emotions <- c("Angry", "Depressed", "Sad", "Anxious", "Relaxed", "Happy")
# 
# createMetaDataFile(vars = Emotions, recode = c("KWAAD", "DEPRE", "DROEV", "ANGST", "ONTSP", "BLIJ"), obsIDvars = newVarLabels[1:5], scale = "1-100", 
#                    author = "Bringmann", year = "2016", data = data)
# 

