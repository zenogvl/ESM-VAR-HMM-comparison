
#Import both datasets
DataInS2 <- read.csv("Data/FromOriginalPapers/Bar_2020_S2_raw.csv")

#Add 1 to day to make the day's start at 1. 
DataInS2$DIARYDAY <- DataInS2$DIARYDAY + 1

#Separate the data of the couples 
DataInS2_men <- DataInS2[,c(1:6)]
DataInS2_women <- DataInS2[,c(1:2,7:10)]

#Create new ID number for women by multipling ID with 2. 
DataInS2_women$ID <- DataInS2_women$ID  * 2


#Rename columns 
varLabels <- c("ID", "day", "Anxiety", "Sad", "Vigor", "Contentment")

colnames(DataInS2_men) <- varLabels
colnames(DataInS2_women) <- varLabels

#Rbind man and women objects
DataS2 <- rbind(DataInS2_men, DataInS2_women)

#Read csv
write.csv(DataS2, "Data/Clean/Bar_2020_ds2.csv")

Emotions <- c("Anxiety", "Sad", "Vigor", "Contentment")
