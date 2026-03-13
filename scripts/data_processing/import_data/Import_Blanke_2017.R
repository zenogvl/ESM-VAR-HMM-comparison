


#Get input data
dataIn <- haven::read_sav("Data/FromOriginalPapers/Blanke_2017_raw.sav")


dateToDayBeep <- function(x){
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

varSelect <- c("ID_anonym", "day", "beep","a_ftl_0", "A_miss", "M_nerv", "M_nieder", "M_bekuem", "M_Aff1", "M_entsp", "M_zufr")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Nervous", "Downhearted", "Distressed", "Happy", "Relaxed", "Content")

data <- dataIn %>%
  mutate(date = as.Date(A_time), #Get date as variable
         time =  format(as.POSIXct(A_time), format = "%H:%M:%S") ) %>% # Get time as variable
  split(dataIn$ID_anonym) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  filter(A_miss != 3) %>% #Remove all skipped timepoints
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns


#Get normal row numbering
rownames(data) <- 1:nrow(data)

#Write data as csv
write.csv(data, file = "Data/Clean/Blanke_2017.csv")


#Vector of emotions in study
Emotions <- c("Nervous", "Downhearted", "Distressed", "Happy", "Relaxed", "Content")


#createMetaDataFile(vars = Emotions, recode = varSelect[-(1:5)], obsIDvars = newVarLabels[1:5], 
#                   scale = "1-7", author = "Blanke", year = "2017", data = data)
