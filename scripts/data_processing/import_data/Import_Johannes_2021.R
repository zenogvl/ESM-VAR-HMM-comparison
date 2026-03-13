
#Read data
dataIn <- read.csv("Data/FromOriginalPapers/Johannes_2021_raw.csv")


#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  x <- x[order(x$number_probe_overall), ] #Sort data so the first measurement is the first row 
  x$Day <- as.integer(as.Date(x$day_probe)  - min(as.Date(x$day_probe) ) + 1) #Get integer for number of day 

  x
}


#Old and new Var labels 
varSelect <- c("pp", "Day", "number_probe_day", "number_probe_overall", "missed", 
               "well_being_discontent", "well_being_calm", "well_being_without_energy", "well_being_well", "well_being_tense", "well_being_happy")
newVarLabels <- c("ID", "day", "beep", "obsNumber", "missed", "Discontent", "Calm", "Energetic", "Wellbeing", "Tense", "Happy")

#Get clean data
data <- dataIn %>% 
  mutate(missed = recode(survey_taken, "Yes" = 0, "No" = 1)) %>% #Change the labels for missed survays from yes and no to 0 and 1. 
  split(dataIn$pp) %>%  #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep 
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns


#Write clean data to csv
write.csv(data, "Data/Clean/Johannes_2021.csv")

Emotions <- c("Discontent", "Calm", "Energetic", "Wellbeing", "Tense", "Happy")
