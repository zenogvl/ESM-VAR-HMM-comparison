
#Read data
dataIn <- haven::read_sav("Data/FromOriginalPapers/Brans_2013_S2_raw.sav")

dateToDayBeep <- function(x){
  #Create empty df with day and beep vectors
  out <- data.frame(day_new = vector(mode = "integer", length = nrow(x)),
                    beep = vector(mode = "integer", length = nrow(x)))
  
  #Get count variable of day starting with 1 by substracting the lowest date of all the other dates and transforming it to an integer
  out$day_new <- as.integer(x$day - min(x$day) + 1)
  
  #Get beep 
  out$beep <- out$day_new %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  cbind(x,out)
}

#Select the relevant variables and rename them

varSelect <- c("ppid", "day_new", "beep", "angry" , "depr", "anx", "sad", "relax", "happy")
newVarLabels <- c("ID", "day", "beep","Angry" , "Depressed", "Anxiety", "Sad", "Relaxed", "Happy")


#Get clean data
data <- dataIn %>% 
  split(dataIn$ppid) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
rownames(data) <- 1:nrow(data)

#Write data as csv
write.csv(data, file = "Data/Clean/Brans_2013_ds2.csv")

#Vector of emotions in study
Emotions <- c("Angry" , "Depressed", "Anxiety", "Sad", "Relaxed", "Happy")

#createMetaDataFile(vars = Emotions, obsIDvars = newVarLabels[1:3], scale = "1-100", 
#                   author = "Brans", year = "2013_DS2", data = data)
