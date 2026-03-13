

#Read data 
dataIn <- read.csv("Data/FromOriginalPapers/Wright_2017_raw.csv", sep = ";")


#Some steps taken from the code from Wendt et al. (2020)

dat <- dataIn %>% 
  filter(!(tmstmp1 == 2)) %>% #Remove missing data
  unique() %>% #Remove duplicate rows
  mutate(count_NA = rowSums(is.na(.))) %>% # Count the number of NA in every row 
  filter(count_NA < 5)

###################################################################
######### Completely copied from https://osf.io/zj9uc #############
###################################################################

dat1 <- dat[
  with(dat, order(UsrID, Date, tmstmp1)),
]

trial <- 0
id <- 0
i <- 2
date <- 0

for (i in 1:nrow(dat1)){
  if (dat1$UsrID[i]==id) {
    trial <- trial+1
  } else {    
    trial <- 0
  }
  id <- dat1$UsrID[i]
  dat1$trial[i] <- trial+1
}

date <- 1160
id <- 80000
trial <- 1

for (i in 1:nrow(dat1)){
  
  if (dat1$UsrID[i]==id) { 
    
    if (dat1$Date[i]==date) {trial <- trial+1} else {trial <- trial+2}
    
    dat1$trial[i] <- trial
    date <- dat1$Date[i]
    
  } else 
    
  {
    trial <- 1
    date <- dat1$Date[i]
    id <- dat1$UsrID[i]
  }
  
}

dat2 <- dat1[with(dat1, order(UsrID, trial)),]

dat <- dat2

daily_trial <- 0
Date <- 0
i <- 1

for (i in 1:nrow(dat)){
  if (dat$Date[i]==Date) {
    daily_trial <- daily_trial+1
  } else {
    daily_trial <- 0
  }
  Date <- dat$Date[i]
  dat$daily_trial[i] <- daily_trial+1
}


day <- 1
Date <- 0
i <- 1
UsrID <- 0

for (i in 1:nrow(dat)){
  
  if (dat$UsrID[i]==UsrID) {
  } else {
    day <- 0
  }
  
  if (dat$Date[i]==Date) {
  } else {
    day <- day+1
  }
  Date <- dat$Date[i]
  UsrID <- dat$UsrID[i]
  dat$day[i] <- day
}



##################################################################
###################  End of copied code  #########################
##################################################################


data <- dat[c(2, 40, 39, 38, 6:36)]

colnames(data) <- c("ID", "day", "beep", "observationID",
                   "Afraid","Ashamed", "Distressed",
                   "Guilty" ,"Hostile","Irritable","Jittery","Nervous","Scared",
                   "Upset","Frightened","Shaky","Angry","Scornful", "Disgusted","Loathing","Sad","Blue",
                   "Downhearted", "Alone", "Lonely", "Active", "Alert","Attentive" ,
                   "Determined", "Enthusiastic", "Excited","Inspired","Interested",
                   "Proud", "Strong")

Emotions <- c("Afraid","Ashamed", "Distressed",
              "Guilty" ,"Hostile","Irritable","Jittery","Nervous","Scared",
              "Upset","Frightened","Shaky","Angry","Scornful", "Disgusted","Loathing","Sad","Blue",
              "Downhearted", "Alone", "Lonely", "Active", "Alert","Attentive" ,
              "Determined", "Enthusiastic", "Excited","Inspired","Interested",
              "Proud", "Strong")

#Write clean data
write.csv(data, "Data/Clean/Wright_2017.csv")



