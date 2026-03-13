

#Read data
dataIn <- haven::read_sav("Data/FromOriginalPapers/Hensel_2022_raw.sav")

#Select relevant columns and change names 
varSelect <- c("ID", "Day", "angry","happy", "friendly", "irritable")
newVarLabels <- c("ID", "day", "Angry","Happy", "Friendly", "Irritable")

data <- dataIn %>%
  select(all_of(varSelect)) %>%
  rename_with(~newVarLabels, all_of(varSelect))

#Write data
write.csv(data, "Data/Clean/Hensel_2022.csv")


Emotions <- c( "Angry","Happy", "Friendly", "Irritable")