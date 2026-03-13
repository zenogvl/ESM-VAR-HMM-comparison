
# 
# data <- read.csv("Data/CleanRescaled/Rescaled_Kroencke_2020.csv")
# variables <- variables_list$Kroencke_2020
# oos_validation_method = c("percentage", "number")[1]
# oos_validation_size = 0.2
# study_name = ds
# output_location = ""
# return_results = TRUE
# min_sample_size = 20
# na_handling_method = c("listwise_delition", "kalman_filter")[2]


add_night <- function(x){
  emptyRow <- rep(NA,ncol(x))
  i <- 1
  while(i < nrow(x)){
    if(x$day[i] != x$day[i + 1]){
      x <- rbind(x[1:i,], emptyRow, x[(i + 1):nrow(x),])
      i <- i + 1
      x[i, which(colnames(x) == "ID")] <- x$ID[i - 1]
    }
    
    i <- i + 1
  }
  x
}

fit_var_models <- function(data, # data of a single study
                           variables, #vector with variable names 
                           oos_validation_method = c("percentage", "number"), #Use a percentage or a specific value to create test data for out of sample validataion 
                           oos_validation_size, #The percentage of data or number of observation for the test dataset. 
                           study_name,
                           output_location,
                           return_results = TRUE,
                           min_sample_size = 20,
                           na_handling_method = c("listwise_delition", "kalman_filter")
) {
  
  ID_enough_data <- data %>%
    na.omit() %>%
    count(ID) %>%
    filter(n > min_sample_size ) %>%
    pull(ID)
  data <- data[data$ID %in% ID_enough_data,]
  
  if(na_handling_method == "kalman_filter") {
    for(pp in unique(data$ID)){
      data[data$ID == pp, variables] <- kalman_filter(data[data$ID == pp, variables]) 
    }
  }
  
  if("beep" %in% colnames(data)){
    data <- data %>%  
      add_night() %>% #add a empty row after each day 
      mutate(across(all_of(variables), ~lag(.), .names =  "{.col}_lag"))  #Create a lag 
  } else {
    data <- data %>%  
      mutate(across(all_of(variables), ~lag(.), .names =  "{.col}_lag"))  #Create a lag 
    
  }
  
  pp <- unique(data$ID)[1]
  #Create output 
  n <- length(unique(data$ID))
  results <- as.data.frame(matrix(NA, nrow = n, ncol = length(variables) + 1))
  colnames(results) <- c("ID", variables)
  results$ID <- unique(data$ID)
  for(pp in unique(data$ID)){
    #print(pp)
    x <- data[data$ID == pp,]
    #Remove all the rows with missing values 
    x <- na.omit(x)
    
    if(nrow(x) < min_sample_size) next #Remove pp with less than 10 observations 
    
    #Split the data in train and test data set
    if(oos_validation_method == "percentage"){
      if(!(oos_validation_size > 0 & oos_validation_size < 1)){
        stop("for method percentage, oosValidationSize should be between 0 and 1")
      } else {
        split_index <- nrow(x) - round(nrow(x) * oos_validation_size, 0)
        x_train <- x[1:(split_index-1),]
        x_test <- x[split_index:nrow(x),]
      }
    } else if (oos_validation_method == "number"){
      if(!(oos_validation_size > 1 & is.integer(oos_validation_size))){
        stop("for method number, oosValidationSize should be an integer")
      } else {
        x_train <- x[1:(oos_validation_size-1),]
        x_test <- x[oos_validation_size:nrow(x),]
      }
    }
    
    
    #Create a fromula where all the lagged variables are the predictors 
    ff_predictors <-  paste0(variables, "_lag", collapse = "+")
    
    predictions <-  as.data.frame(matrix(NA, nrow = nrow(x_test), ncol = length(variables)))
    colnames(predictions) <- c(variables)
    
    v <- variables[1]
    for(v in variables){
      ff <- paste0(v,  "~", ff_predictors)
      #Fit the VAR model on the training data
      fit <- stats::lm(ff , x_train)
      
      #Predict the values 
      pred <- stats::predict(fit, x_test)
      RMSE <- sqrt(1 / nrow(x_test) *
                     sum((x_test[, v] - pred)^2))
      
      results[results$ID == pp, v] <- RMSE
    }
    
  }
  
  write.csv(results, paste0(output_location, "RMSE_", study_name, ".csv"))
  if(return_results) {
    results
  }
}


for(ds in dataset_names_all){
  print(ds)
  fit_var_models(data = read.csv(paste0("Data/CleanRescaled/Rescaled_", ds, ".csv")),
                 variables = variables_list[[ds]],
                 oos_validation_method = "percentage", 
                 oos_validation_size = 0.2, 
                 study_name = ds,
                 output_location = "Output/N1_Results/Rescaled_results/VAR_kalman_filter_RMSE/",
                 return_results = FALSE,
                 min_sample_size = 20,
                 na_handling_method = "kalman_filter"
                 )
}



for(ds in dataset_names_all){
  print(ds)
  fit_var_models(data = read.csv(paste0("Data/CleanRescaled/Rescaled_", ds, ".csv")),
                 variables = variables_list[[ds]],
                 oos_validation_method = "percentage", 
                 oos_validation_size = 0.2, 
                 study_name = ds,
                 output_location = "Output/N1_Results/Rescaled_results/VAR_listwise_deletion_RMSE/",
                 return_results = FALSE,
                 min_sample_size = 20,
                 na_handling_method = "listwise_delition"
  )
}










varN1 <- function(data, # data of a single study
                  variables, #vector with variable names 
                  oosValidationMethods = c("percentage", "number"), #Use a percentage or a specific value to create test data for out of sample validataion 
                  oosValidationSize, #The percentage of data or number of observation for the test dataset. 
                  studyName,
                  outputLocation,
                  returnResults = TRUE,
                  minDataPoints = 10,
                  naMethod = c("listwiseDelition", "kalmanFilter")
                  ) 
  {
  
  

  #Get the ID's of all the PP that have more than min data points
  
  minDataID <- data %>%
    na.omit() %>%
    count(ID) %>%
    filter(n > minDataPoints ) %>%
    pull(ID)
  
  data <- data[data$ID %in% minDataID,]
  
  if(naMethod == "kalmanFilter") {
    for(pp in unique(data$ID)){
      data[data$ID == pp, variables] <- kalmanFilter(data[data$ID == pp, variables]) 
    }
  }
  
  
  
  if("beep" %in% colnames(data)){
    data <- data %>%  
      addNight() %>% #add a empty row after each day 
      mutate_at(all_of(variables), list(lag = ~lag(.)))  #Create a lag 
  } else {
    data <- data %>%  
      mutate_at(all_of(variables), list(lag = ~lag(.)))  #Create a lag 
    
  }
  
    
  #Create output 
  n <- length(unique(data$ID))
  results <- as.data.frame(matrix(NA, nrow = n, ncol = length(variables) + 1))
  colnames(results) <- c("ID", variables)
  results$ID <- unique(data$ID)
  
  pp <- 1091     
  strangepred <- vector()
  for(pp in unique(data$ID)){
    #print(pp)
    x <- data[data$ID == pp,]
    #Remove all the rows with missing values 
    x <- na.omit(x)
    
    if(nrow(x) < minDataPoints) next #Remove pp with less than 10 observations 

    #Split the data in train and test data set
    if(oosValidationMethods == "percentage"){
      if(!(oosValidationSize > 0 & oosValidationSize < 1)){
        stop("for method percentage, oosValidationSize should be between 0 and 1")
      } else {
        splitIndex <- nrow(x) - round(nrow(x) * oosValidationSize, 0)
        x_train <- x[1:(splitIndex-1),]
        x_test <- x[splitIndex:nrow(x),]
      }
    } else if (oosValidationMethods == "number"){
      if(!(oosValidationSize > 1 & is.integer(oosValidationSize))){
        stop("for method number, oosValidationSize should be an integer")
      } else {
        x_train <- x[1:(oosValidationSize-1),]
        x_test <- x[oosValidationSize:nrow(x),]
      }
    }
    
    
    #Create a fromula where all the lagged variables are the predictors 
    ff_predictors <-  paste0(variables, "_lag", collapse = "+")
    # 
    # 2.635953e-36
    # 
    # x_train$Sad
    # 
    v <- variables[11]
    
    predictions <-  as.data.frame(matrix(NA, nrow = nrow(x_test), ncol = length(variables)))
    colnames(predictions) <- c(variables)
    

    
    for(v in variables){
      ff <- paste0(v,  "~", ff_predictors)
      #Fit the VAR model on the training data
      fit <- stats::lm(ff , x_train)
      
      # 
      # fit$coefficients %>% 
      #   broom::tidy() %>%
      #   filter(names == "Sad_lag") %>% 
      #   select(x) %>%
      #   round(1000)
      # 
      # 
      #   filter( x > 10 | x < 10)
      
      #Predict the values 
      pred <- stats::predict(fit, x_test)
      
      predictions[,v] <- pred
      
      Metrics::rmse(x_test[,v], pred )
      
      RMSE <- sqrt(1/nrow(x_test)*
                     sum((x_test[, v] - pred)^2))
      
      results[results$ID == pp, v] <- RMSE
    }
    strangepred <- c(strangepred, sum(predictions > 2 | predictions < -2)/(nrow(predictions)*ncol(predictions)))
    
    
    
  }
  

}






for(ds in ds_names){
  print(ds)
  varN1(data = read.csv(paste0("Data/CleanRescaled/Rescaled_", ds, ".csv")), 
        variables = variables_list[[ds]], 
        oosValidationMethods = "percentage", 
        oosValidationSize = 0.2, 
        studyName = ds, 
        outputLocation = "Output/N1_Results/RescaledResults/VAR_kalman_filter_RMSE/",
        minDataPoints = 20,
        naMethod = "kalmanFilter",
        returnResults = FALSE
  )
}

for(ds in ds_names){
  print(ds)
  varN1(data = read.csv(paste0("Data/CleanRescaled/Rescaled_", ds, ".csv")), 
        variables = variables_list[[ds]], 
        oosValidationMethods = "percentage", 
        oosValidationSize = 0.2, 
        studyName = ds, 
        outputLocation = "Output/N1_Results/RescaledResults/VAR_listwise_deletion_RMSE/",
        minDataPoints = 20,
        naMethod = "listwiseDelition",
        returnResults = FALSE
  )
}