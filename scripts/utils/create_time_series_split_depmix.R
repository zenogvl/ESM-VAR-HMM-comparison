

create_time_series_split_depmix <-  function(data, variables, method = c("both", "days", "missing")){
  
  #Create a vector that indicates missing data 
  missing_vector <- data[,variables] %>%
    is.na() %>%
    apply(1, any) %>%
    as.vector()
  
  if("beep" %in% colnames(data)){
    if(method == "days"){
      time_serie_split <- data$day %>% 
        table() %>%
        as.vector()
      data_out <- data
    }
    if(method == "missing") {
      #Creates a vector that shows for every row to what seperate time series it belongs. If there is missing data, a new ts starts. 
      seperate_series <- rep(0, length(missing_vector))
      ts <- 1 #Start with 1 for the first ts
      for(i in 2:length(missing_vector)){
        seperate_series[i-1] <- ts 
        if(missing_vector[i] & !missing_vector[i-1]){ #If there is a missing value, a 1 is added to the ts counter to indicate the next time series start. 
          ts <- ts + 1
        }
      }
      #Add ts to the last row (if this is missing, a new ts might start but since it then will be removed this doesn't matter)
      seperate_series[length(missing_vector)] <- ts
      
      data$seperate_series <- seperate_series
      
      #Remove all the missing data 
      data_out <- data[!missing_vector,]
      
      #Now get the vector that indicates how long every time serie is
      time_serie_split <- data_out$seperate_series %>% 
        table() %>%
        as.vector()
    }
    if(method == "both"){
      #The same logic as for the missing data is used, but now a new ts also starts when a new day starts. 
      seperate_series <- rep(0, length(missing_vector))
      ts <- 1 
      data$day
      for(i in 2:length(missing_vector)){
        seperate_series[i-1] <- ts 
        if((missing_vector[i] & !missing_vector[i-1]) | (data$day[i-1] != data$day[i])) { #Also update the ts counter when a new day starts
          ts <- ts + 1
        }
      }
      seperate_series[length(missing_vector)] <- ts
      
      data$seperate_series <- seperate_series
      data_out <- data[!missing_vector,]
      time_serie_split <- data_out$seperate_series %>% 
        table() %>%
        as.vector()
    }
  } else {
    #if(method != "missing") warning("data has no beep variable, days")
    if(method == "days"){ #When obs/day=1, no split takes places 
      time_serie_split <- nrow(data)
      data_out <- data
    } else { #Method both and missing are the same 
      seperate_series <- rep(0, length(missing_vector))
      ts <- 1 #Start with 1 for the first ts
      for(i in 2:length(missing_vector)){
        seperate_series[i-1] <- ts 
        if(missing_vector[i] & !missing_vector[i-1]){ #If there is a missing value, a 1 is added to the ts counter to indicate the next time series start. 
          ts <- ts + 1
        }
      }
      seperate_series[length(missing_vector)] <- ts 
      data$seperate_series <- seperate_series
      data_out <- data[!missing_vector,]
      time_serie_split <- data_out$seperate_series %>% 
        table() %>%
        as.vector()
      
    }
  }
  
  output <- list(data = data_out, 
                 variables = variables, 
                 length_time_serie = time_serie_split)
 output 
}
