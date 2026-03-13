

#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Kroencke_2020_raw.csv")

#Get the emotions from the colnames and by removing part of the string.  
emotions <- gsub("int_", "", colnames(dataIn)[197:210])

#### Emotions are in different columns based on having an interaction or not, merge these columns: ####

#Create empty dataframe to add emotions to
dataComb <- matrix(NA, nrow = nrow(dataIn), ncol = 17, dimnames = list(NULL, c("ID", "Time", "Wave", emotions))) %>%
  as.data.frame()
dataComb$ID <- dataIn$id
dataComb$Time <- dataIn$created_esm
dataComb$Wave <- dataIn$wave

#For every emotion, combine the data 
for(i in emotions){
  dataComb[, i] <- coalesce(dataIn[, paste0("int_",i)], dataIn[, paste0("occup_", i)])
}


#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  
  #### Combine the two waves 
  
  #Create object with numbered days 
  Day <- as.integer(as.Date(x$Time) - min(as.Date(x$Time)) + 1)
  
  if(x$Wave[1] == "Both"){
    #Get the difference in number of days between the rounds 
    DayDistance <- c(Day,max(Day)) - c(0,Day)
    if(max(DayDistance) > 1){
      #Substract the difference between waves from the second wave
      Day[which.max(DayDistance):length(Day)] <- Day[which.max(DayDistance):length(Day)] -  max(DayDistance) + 1
      x$Day <- Day
    } else { # If waves might be just behind each other (to be sure nothing brakes down)
      x$Day <- Day
    }
  } else { # For pp who only did one wave
    x$Day <- Day
  }
  
  x$beep <- x$Day %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  #Get observation ID
  x$observationID <- 1:nrow(x)
  
  x
}


#Selected variables and new labels 
varSelect <- c("ID", "Day", "beep", "observationID", emotions)
Emotions <- str_replace(emotions, "^\\w{1}", toupper)
newVarLabels <- c("ID", "day", "beep", "observationID", Emotions)

#Get celan data
data <- dataComb %>% 
  split(dataComb$ID) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
rownames(data) <- 1:nrow(data)


#Write clean data
write.csv(data, "Data/Clean/Kroencke_2020.csv")
