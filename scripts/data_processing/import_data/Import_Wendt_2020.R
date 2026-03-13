
#Read data 
dataIn <- haven::read_sav("Data/FromOriginalPapers/Wendt_2020_raw.sav")

#Change ID label
dataIn$ID <- match(dataIn$vp_id2, unique(dataIn$vp_id2), nomatch = NA)

#Selected variables and new labels 
varSelect <- c("ID", "Day", 
               "emo_1", "emo_2", "emo_3", "emo_4", "emo_5", "emo_6", "emo_7", 
               "emo_8", "emo_9", "emo_10", "emo_11", "emo_12", "emo_13", "emo_14", 
               "emo_15", "emo_16", "emo_17", "emo_18", "emo_19", "emo_20", "emo_21", 
               "emo_22", "emo_23", "emo_24", "emo_25", "emo_26", "emo_27", "emo_28", 
               "emo_29", "emo_30", "emo_31", "emo_32", "emo_33", "emo_34", "emo_35", 
               "emo_36", "emo_37", "emo_38")
newVarLabels <- c("ID", "day", 
                  "Atease", "Calm", "Relaxed", "Drowsy", "Sluggish", "Tired", "Content",
                  "Pleased", "Happy", "Alert", "Aroused", "Hyperactivated", "Miserable", "Troubled",
                  "Unhappy", "Sleepy", "Quiet", "Still", "Afraid", "Ashamed", "Distressed", "Guilty",
                  "Hostile", "Irritable", "Jittery", "Nervous", "Scared", "Upset", "Active", "Alone",
                  "Attentive", "Determined", "Enthusiastic", "Excited", "Inspired", "Interested", 
                  "Proud", "Strong")
EmotionsAll <- c("Atease", "Calm", "Relaxed", "Drowsy", "Sluggish", "Tired", "Content",
              "Pleased", "Happy", "Alert", "Aroused", "Hyperactivated", "Miserable", "Troubled",
              "Unhappy", "Sleepy", "Quiet", "Still", "Afraid", "Ashamed", "Distressed", "Guilty",
              "Hostile", "Irritable", "Jittery", "Nervous", "Scared", "Upset", "Active", "Alone",
              "Attentive", "Determined", "Enthusiastic", "Excited", "Inspired", "Interested", 
              "Proud", "Strong")


Emotions <- c( "Calm", "Relaxed", "Tired", "Content",
                 "Pleased", "Happy", "Alert", 
                 "Unhappy", "Afraid", "Ashamed", "Distressed", "Guilty",
                 "Hostile", "Irritable", "Jittery", "Nervous", "Scared", "Upset", "Active", "Alone",
                 "Attentive", "Determined", "Enthusiastic", "Excited", "Inspired", "Interested", 
                 "Proud", "Strong")


#Get clean data
data <- dataIn %>% 
  filter(k_days == 25) %>% 
  mutate(Day = as.integer(tbday_nrs)) %>%
  select(all_of(varSelect))  %>% #Select relevant columns 
  rename_with(~newVarLabels, all_of(varSelect)) #Change the names of these columns

#Write clean data
write.csv(data, "Data/Clean/Wendt_2020.csv")
