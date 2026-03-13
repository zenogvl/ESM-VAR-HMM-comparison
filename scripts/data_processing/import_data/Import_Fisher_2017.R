
#Get paths of files without the first one 
path <- "Data/FromOriginalPapers/Fisher_2017_raw/"
files <- paste0(path,list.files(path))[-1]

#Get vector with ID numbers
IDvector <- gsub("_final.RData", "", list.files(path)[-1]) %>%
  gsub(pattern = "P", replacement = "") %>%
  gsub(pattern = "p", replacement = "") %>%
  as.numeric()


dataList <- vector(mode = "list", length = length(files)-1)
#Load pp data, add ID and combine in list of data
for(i in 1:length(files)){
  load(files[i])
  
  #Duplicate columns without colname are removed
  dataList[[i]] <-  cbind(ID = IDvector[i], data[,!is.na(colnames(data))])
}


dateToDayBeep <- function(x){
  
  #Transfrom to time format
  time <- as.POSIXct(x$start, format = "%m/%d/%Y %H:%M")
  
  #Separate date and time 
  x$date <- as.Date(format(time, format = "%Y/%m/%d"))
  x$time <- format(time, format = "%H:%M")
  
  #Create empty df with day and beep vectors
  out <- data.frame(day = vector(mode = "integer", length = nrow(x)),
                    beep = vector(mode = "integer", length = nrow(x)))
  
  #Get count variable of day starting with 1 by substracting the lowest date of all the other dates and transforming it to an integer
  out$day <- as.integer(x$date - min(x$date) + 1)
  
  #Get beep 
  out$beep <- out$day %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  cbind(x,out)
}


#Combine into a single df
data <- NA
data <- lapply(dataList, dateToDayBeep) %>%
  rlist::list.rbind()




varSelect <- c("ID", "day", "beep", "energetic", "enthusiastic", 
               "content", "irritable", "restless", "afraid",
               "angry", "guilty", "worried", "hopeless", 
               "down", "fatigue", "positive", "date", "time")
newVarLabels <- c("ID", "day", "beep", "Energetic", "Enthusiastic",
               "Content", "Irritable", "Restless", "Anxious", 
               "Angry", "Guilty", "Worried", "Hopeless",
               "Depressed", "Fatigue", "Positive", "date", "time")

Emotions <- c("Energetic", "Enthusiastic",
              "Content", "Irritable", "Restless", "Anxious", 
              "Angry", "Guilty", "Worried", "Hopeless",
              "Depressed", "Fatigue", "Positive")



data <- data[,varSelect]
colnames(data) <- newVarLabels

#Write data as csv
write.csv(data, file = "Data/Clean/Fisher_2017.csv")

