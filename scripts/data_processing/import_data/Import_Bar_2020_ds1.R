
#Import both datasets
DataInS1 <- read.csv("Data/FromOriginalPapers/Bar_2020_S1_raw.csv")

#Add 1 to day to make the day's start at 1. 
DataInS1$DIARYDAY <- DataInS1$DIARYDAY + 1

#Separate the data of the couples 
DataInS1_men <- DataInS1[,c(1:6)]
DataInS1_women <- DataInS1[,c(1:2,7:10)]

#Create new ID number for women by multipling ID with 2. 
DataInS1_women$ID <- DataInS1_women$ID  * 2

#Rename columns 
varLabels <- c("ID", "day", "Anxiety", "Sad", "Vigor", "Contentment")

colnames(DataInS1_men) <- varLabels
colnames(DataInS1_women) <- varLabels

#Rbind man and women objects
DataS1 <- rbind(DataInS1_men, DataInS1_women)

#Read csv
write.csv(DataS1, "Data/Clean/Bar_2020_ds1.csv")

Emotions <- c("Anxiety", "Sad", "Vigor", "Contentment")

#createMetaDataFile(vars = Emotions, recode = c("Anx", "Sad","Vig", "Con"), obsIDvars = c("ID", "day"), scale = "1-5", 
#                   author = "Bar Kalifa", authorShort = "Bar", year = "2019_DS2", data = DataS2)


