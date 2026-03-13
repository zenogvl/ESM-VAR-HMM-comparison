



#Read data 
dataIn <- haven::read_sav("Data/FromOriginalPapers/Naragon_2023_raw.sav")

#Create external objects for variables  
varSelect <- c("ID", "DAY", "SIG", "Signal", "Active", "Interested", "Excited", "Strong", "Irritable", "Upset", "AfraidAnx", "Sad")
varSelect2 <-  c("ID", "DAY", "SIG", "Signal", "missing","Active", "Interested", "Excited", "Strong", "Irritable", "Upset", "AfraidAnx", "Sad")
newVarLabels <- c("ID", "day", "beep", "observationID", "missing","Active", "Interested", "Excited", "Strong", "Irritable", "Upset", "Anxious", "Sad")
Emotions <- c("Active", "Interested", "Excited", "Strong", "Irritable", "Upset", "Anxious", "Sad")


data <- dataIn %>%
  select(all_of(varSelect)) %>% # Select the columns with relevant data
  mutate(count_NA = rowSums(is.na(.))) %>% #Get the sum of NA in every row
  #Create a variable to code missing data
  mutate(missing = case_when(count_NA == 0  ~ 0,
                             count_NA == 8 ~ 1,
                             count_NA > 0 & count_NA < 8 ~ 2
  )) %>%
  select(all_of(varSelect2)) %>%
  rename_with(~newVarLabels, all_of(varSelect2)) 

#Remove pp with less than 30% of responses 
IDppRemove <- data %>% 
  filter(missing == 1) %>%
  group_by(ID) %>%
  summarise(missingCount = sum(missing)) %>%
  filter(missingCount > 29)

data <- data %>%
  filter(!(ID  %in% IDppRemove$ID))



#Write clean data
write.csv(data, "Data/Clean/Naragon_2023.csv")

