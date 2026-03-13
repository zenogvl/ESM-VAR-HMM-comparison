

#Read data
dataInPart1 <- readxl::read_xlsx("Data/FromOriginalPapers/Forbes_2022_raw/Forbes_2022_part1.xlsx")
dataInPart2 <- readxl::read_xlsx("Data/FromOriginalPapers/Forbes_2022_raw/Forbes_2022_part2.xlsx")

#Select relevant column form both datasets and combine them
varSelect <- c("Participant_e",  "day","MZP","count", "stress", "valence_mood", "calmness_mood", "energy_mood")
dataIn <- rbind(dataInPart1[,varSelect], dataInPart2[,varSelect])

#Change ID from string with random letters and numbers to integers 
dataIn$Participant_e <- match(dataIn$Participant_e, unique(dataIn$Participant_e), nomatch = NA) %>%
  as.integer()

#Change labels
varLabels <- c("ID", "day", "beep", "numberOfMeasurment", "Stressed", "Valance", "Calmness", "Energetic" )
Emotions <-  c("Stressed", "Valance", "Calmness", "Energetic" )


out <- dataIn
colnames(out) <- varLabels



#Save dataset
write.csv(out, "Data/Clean/Forbes_2022.csv")



