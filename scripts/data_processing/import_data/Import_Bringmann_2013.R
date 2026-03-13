
#Read data
dataIn <- read.table("Data/FromOriginalPapers/Bringmann_2013_raw.txt", sep = ",", header = TRUE)

#Select and rename vars
varSelect <- c("subjno", "dayno", "beepno", "opgewkt_", "pieker", "angstig_", "somber__", "ontspann")
varLabels <- c("ID", "day", "beep", "Excited", "Worried", "Anxiety", "Sad", "Relaxed")
Emotions <- c("Excited", "Worried", "Anxiety", "Sad", "Relaxed")


data <- dataIn[,varSelect]
colnames(data) <- varLabels

#Print data
write.csv(data, "Data/Clean/Bringmann_2013.csv")


#createMetaDataFile(vars = varLabels[4:8], recode = varSelect[4:8], obsIDvars = varLabels[1:3], scale = "1-7", 
#                   author = "Bringmann", year = "2013", data = data)

