
#Read data
dataIn <- haven::read_sav("Data/FromOriginalPapers/Dejonckheere_2019_raw.sav")

#Split into separate df for every individual 
dataIn_split <- split(dataIn, dataIn$ppid)


dateToDayBeep <- function(x){
  #Create empty df with day and beep vectors
  out <- data.frame(day_new = vector(mode = "integer", length = nrow(x)),
                    beep_new = vector(mode = "integer", length = nrow(x)))
  
  #Transform separate year, month and day column to date vector 
  Date <- as.Date(paste0(x$beepyr, "-", x$beepmo, "-", x$beepday), "%Y-%m-%d")
  
  #Get count variable of day starting with 1 by subtracting the lowest date of all the other dates and transforming it to an integer
  out$day_new <- as.integer(Date - min(Date) + 1)
  
  #Get beep 
  out$beep_new <- out$day_new %>%
    table() %>% # Get the frequency of each day
    as.vector() %>% #Transform frequencies of each day to a vector
    sequence() #Use this vector to get sequences ranging from 1 to the frequency of each day
  
  cbind(x,out)
}


varSelect <- c("ppid", "day_new", "beep_new", "missing", "observationID", "angry" ,  "sad", "stress", "relax", "happy")
newVarLabels <- c("ID", "day", "beep", "missing", "observationID", "Angry" , "Sad", "Stressed", "Relaxed", "Happy")



#create day to beep for every pp and rbind them
data <- lapply(dataIn_split, dateToDayBeep) %>%
  rlist::list.rbind() %>%
  mutate(missing = case_when(is.na(angry) ~ 1, !is.na(angry) ~ 0,),
         observationID = rep(1:98, length(unique(ppid)))) %>%
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
rownames(data) <- 1:nrow(data)



#Write data as csv
write.csv(data, file = "Data/Clean/Dejonckheere_2019.csv")


#Vector of emotions in study
Emotions <- c("Angry" , "Sad", "Stressed", "Relaxed", "Happy")

# 
# createMetaDataFile(vars = Emotions , obsIDvars = newVarLabels[1:5], scale = "0-100", 
#                    author = "Dejonckheere", year = "2019", data = data)


