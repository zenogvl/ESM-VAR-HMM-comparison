

#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/deVos_2017_raw.csv", sep = "")

#Function to get the day and/or beep 
dateToDayBeep <- function(x){
  if(nrow(x) != 90){
    error("not 90 obs")
  }
  x$day <- rep(1:30, each = 3)
  x$beep <- rep(1:3, 30)
  x$observationID <- 1:90
  x$missing <- x %>% 
    is.na() %>%
    apply(1, any) %>%
    as.integer()
  x
}

#Selected variables and new labels 
varSelect <- c("subjno", "day", "beep", "observationID", "missing","b1_talkative", "b2_energetic", 
               "b3_tense", "b4_anxious", "b5_enthousiastic", "b6_confident", "b7_distracted" , 
               "b8_restless", "b9_irritated", "b10_satisfied","b11_happy", "b12_depressed", 
               "b13_cheerful", "b14_guilty" )
newVarLabels <- c("ID", "day", "beep", "observationID", "missing", "Talkative", "Energetic", "Tense", 
                  "Anxiety", "Enthousiastic", "Confident", "Distracted" , "Restless", "Irritated", 
                  "Satisfied", "Happy", "Depressed", "Cheerful", "Guilty" )
Emotions <- c("Talkative", "Energetic", "Tense", "Anxiety", "Enthousiastic", "Confident", "Distracted" , "Restless", "Irritated", 
              "Satisfied", "Happy", "Depressed", "Cheerful", "Guilty" )


#Get clean data
data <- dataIn %>% 
  split(dataIn$subjno) %>% #Split into separate df for every individual 
  lapply(dateToDayBeep) %>% #Get beep and day
  rlist::list.rbind() %>% #Combine list of df's into a single one
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Get normal row numbering
rownames(data) <- 1:nrow(data)


#Write clean data
write.csv(data, "Data/Clean/deVos_2017.csv")


