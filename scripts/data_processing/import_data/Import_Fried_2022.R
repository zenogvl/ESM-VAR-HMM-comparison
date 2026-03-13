
#Read data
dataIn <- read.csv("Data/FromOriginalPapers/Fried_2022_raw.csv")


#Get a integer ID 
dataIn$ID2 <- match(dataIn$ID, unique(dataIn$ID), nomatch = NA)


#Split data in list based on pp
dataIn_split <- split(dataIn, dataIn$ID)

dateToDayBeep <- function(x){
  #Get count variable of day starting with 1 by substracting the lowest date of all the other dates and transforming it to an integer
  x$day <- as.integer(as.Date(x$Day) - min(as.Date(x$Day)) + 1)
  x
}

#create day to beep for every pp and rbind them
data <- lapply(dataIn_split, dateToDayBeep) %>%
  rlist::list.rbind()

#Select the relevant variables and rename them
varSelect <- c("ID2", "day", "beepvar", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")
newVarLabels  <- c("ID", "day", "beep", "Relaxed", "Irritable", "Worry", "Nervous",  "Future", "Anhedonia", "Tired", "Hungry", "Alone", "Angry")

data <- data[,varSelect]
colnames(data) <- newVarLabels


#Write data as csv
write.csv(data, file = "Data/Clean/Fried_2022.csv")

Emotions <- c("Relaxed", "Irritable", "Worry", "Nervous",  "Future", "Anhedonia", "Tired", "Hungry", "Alone", "Angry")



