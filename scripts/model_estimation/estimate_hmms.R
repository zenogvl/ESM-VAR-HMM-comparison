depmix_wrapper_pp_level  <- function(input){
  x <- input$data
  ff <- input$ff
  variables <- input$variables
  max_numb_states <- input$max_numb_states
 
  #Split the dataset in train and test dataset 
  if(input$oos_validation_methods == "percentage"){
    if(!(input$oos_validation_size > 0 & input$oos_validation_size < 1)){
      stop("for method percentage, oos_validation_size should be between 0 and 1")
    } else {
      split_index <- nrow(x) - round(nrow(x) * input$oos_validation_size, 0)
      x_train <- x[1:(split_index-1),]
      x_test <- x[split_index:nrow(x),]
    }
  } else if (input$oos_validation_methods == "number"){
    if(!(input$oos_validation_size > 1 & is.integer(input$oos_validation_size))){
      stop("for method number, oos_validation_size should be an integer")
    } else {
      x_train <- x[1:(input$oos_validation_size-1),]
      x_test <- x[input$oos_validation_size:nrow(x),]
    }
  }
  
  #Get a vector that indicates how long all the seperate time series are 
  obj_train <- create_time_series_split_depmix(x_train, variables=variables, method = input$series_split_method)
  x_train <- obj_train$data
  length_time_serie <- obj_train$length_time_serie
  
  
  depmix_object <- vector("list", (max_numb_states-1))
  depmix_fit <- vector("list", (max_numb_states-1))
  tries <- rep(NA, (max_numb_states-1))
  IC <- data.frame(nstate = 2:max_numb_states,
                   AIC = rep(NA, max_numb_states-1),
                   BIC = rep(NA, max_numb_states-1))
  for(s in 2:max_numb_states){
    counter <- 0
    #Repeat loop that continues until a converged model is found. 
    repeat{
      counter <- counter + 1
      depmix_object[[s-1]] <-  depmixS4::depmix(ff, 
                                               data = x_train, 
                                               nstates = s, 
                                               family = rep(list(gaussian()),length(variables)),
                                               ntimes = length_time_serie)
      invisible(capture.output(try(depmix_fit[[s-1]] <- depmixS4::fit(depmix_object[[s-1]]), silent = TRUE))) #Fit model and continu even with an error. 
      
      #If model is empty or did not converge, repeat, if converged, stop. 
      if(is.null(depmix_fit[[s-1]])){
        if(counter == input$max_estimation_tries) break 
      } else if(grepl("without convergence", depmix_fit[[s-1]]@message)){
        if(counter == input$max_estimation_tries) break 
      } else {
        break
      }
    }
    if(counter < input$max_estimation_tries){
      IC$AIC[s-1] <- AIC(depmix_fit[[s-1]])
      IC$BIC[s-1] <- BIC(depmix_fit[[s-1]])
    }
    tries[s-1] <- counter

  }
  #depmix_fit
  #only continues when at least one model is fitted 
  if(any(!is.na(IC[,2:3]))){
    #Select best model based on a IC
    if(input$ic_method == "AIC"){
      best_s_state <- IC$nstate[which.min(IC$AIC)]
    } else if (input$ic_method == "BIC"){
      best_s_state <- IC$nstate[which.min(IC$BIC)]
    }
    
    #Create the TS split for the the test data
    obj_test <- create_time_series_split_depmix(x_test, variables=variables, method = "both")
    x_test <- obj_test$data
    length_time_serie_test <- obj_test$length_time_serie
   
    predict_object <- create_depmix_predict_object(depmix_fit[[best_s_state-1]], data_new = x_test, length_time_serie = length_time_serie_test)
    if(exists("posterior_predict")) rm("posterior_predict")
  
    #Get the most likely state sequence 
    try(posterior_predict <- depmixS4::viterbi(predict_object))
    if(exists("posterior_predict")){
      predict_object@posterior <- posterior_predict
      
      prediction <- predicting_hmms(predict_object)
      RMSE <- vector(mode = "numeric",length = length(variables))
      names(RMSE) <- variables
      
      for(v in variables){
        RMSE[v] <- sqrt(1/nrow(x_test)*
                          sum((x_test[,v] - prediction$prediction[,v])^2)
        )
      }
      return(
        list(depmix_objects = depmix_object,
             depmix_fit = depmix_fit, 
             IC = IC,
             nstateUsed = best_s_state,
             posterior = posterior_predict,
             depmixPredictionObject = predict_object,
             RMSE = RMSE,
             Ntries = tries)
      )
    } else {
      return(NA)
    }
  } 
}

fit_hmms <-  function(data, # data of a single study
                  variables, #vector with variable names 
                  max_numb_states, #The maximam number of states for the HMM to fit. 
                  oos_validation_methods = c("percentage", "number"), #Use a percentage or a specific value to create test data for out of sample validataion 
                  oos_validation_size, #The percentage of data or number of observation for the test dataset. 
                  study_name,
                  output_location,
                  return_results = TRUE,
                  minimal_sample_size = 10,
                  multivar = c(TRUE, FALSE),
                  ic_method = c("AIC", "BIC"),
                  max_estimation_tries = 100,
                  series_split_method = c("both", "days", "missing"),
                  na_method = c("none", "kalman_filter")
) {
  oos_validation_methods <- match.arg(oos_validation_methods)
  #multivar <- match.arg(multivar)
  ic_method <- match.arg(ic_method)
  na_method <- match.arg(na_method)
  #series_split_method <- mathc.arg(series_split_method)

  #Get the ID's of all the PP that have more than min data points
  id_minimal_sample_size <- table(na.omit(data)$ID)[table(na.omit(data)$ID) > minimal_sample_size] %>% #Index the table for above this minimum
    names() %>% #Take the names of this vector to get the IDs
    as.integer() #Transfrom to numeric 
  
  data_small_sample_removed <- data[data$ID %in% id_minimal_sample_size,]
  
  if(na_method == "kalman_filter") {
    for(pp in unique(data$ID)){
      #print(pp)
      data_small_sample_removed[data_small_sample_removed$ID == pp, variables] <- kalman_filter(data_small_sample_removed[data_small_sample_removed$ID == pp, variables]) 
    }
  }
  
  data_list <- split(data_small_sample_removed, data_small_sample_removed$ID)
  
  #Create the formula list:
  ff <- list()
  if(multivar){
    for(v in variables){
      ff[[v]] <-  paste0(v, "~", "1+", paste0(variables[!(variables == v)], collapse  = "+")) %>% 
        stats::as.formula()
    }
  } else {
    for(v in variables){
      ff[[v]] <-  paste0(v, "~", "1") %>% 
        stats::as.formula()
    }
  }
  
  input_list <- list()
  for(pp in names(data_list)){
    input_list[[pp]] <-  list(data = data_list[[pp]],
                             variables = variables,
                             ff = ff,
                             oos_validation_methods = oos_validation_methods,
                             oos_validation_size = oos_validation_size,
                             max_numb_states = max_numb_states,
                             max_estimation_tries = max_estimation_tries,
                             ic_method = ic_method,
                             series_split_method = series_split_method)
    
  }
 
  cores <- parallel::detectCores()
  clus <- parallel::makeCluster(cores)
  parallel::clusterExport(clus, c("create_time_series_split_depmix", 
                                  "depmix_wrapper_pp_level", 
                                  "create_depmix_predict_object",
                                  "densManual",
                                  "predictManual", 
                                  "predicting_hmms"))
  parallel::clusterEvalQ(clus, {
    library(tidyverse)
    library(depmixS4)
  })
  res <- parallel::parLapply(clus, input_list, depmix_wrapper_pp_level)
  parallel::stopCluster(clus)
 

  #Create output matrix RMSE
  RMSE <- matrix(NA, length(res), (length(variables)+1))
  colnames(RMSE) <- c("ID", variables)
  RMSE[,1] <- as.integer(names(data_list))
  
  #Create output matrix for number of tries 
  Ntries <- matrix(NA, length(res), (max_numb_states))
  colnames(Ntries) <- c("ID", 2:max_numb_states)
  Ntries[,1] <- as.integer(names(data_list))
  
  #Create output matrix for number of states
  nstateUsed <- matrix(NA, length(res), 2)
  colnames(nstateUsed) <- c("ID", "nstate")
  nstateUsed[,1] <- as.integer(names(data_list))
  
  
  if(length(res)>0){
    for(pp in 1:length(res)){
      #print(pp)
      if(!is.null(res[[pp]]) & !is.na(res[pp])){
        RMSE[pp,-1] <- res[[pp]]$RMSE
        Ntries[pp,-1] <- res[[pp]]$Ntries 
        nstateUsed[pp,-1] <- res[[pp]]$nstateUsed 
      }
    }
    
    #write.csv(RMSE, "Output/N1_Results/RescaledResults/")
    write.csv(RMSE, paste0(output_location, "RMSE/RMSE_", study_name,".csv"))
    write.csv(Ntries, paste0(output_location, "Ntries/Ntries_", study_name,".csv"))
    write.csv(nstateUsed, paste0(output_location, "Nstates/Nstates_", study_name,".csv"))
    
    
  } 

  saveRDS(res, paste0(output_location, "HMMres_", study_name, ".rds"))
  if(return_results) {
    res
  }
}  

ds = dataset_names[1]
data = read.csv(paste0("Data/CleanRescaled/", "Rescaled_", ds, ".csv"))
variables = variables_list[[ds]]
max_numb_states = 5
oos_validation_methods = "percentage" #Use a percentage or a specific value to create test data for out of sample validataion 
oos_validation_size = .2 #The percentage of data or number of observation for the test dataset. 
study_name = ds
output_location = "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results_2/"
return_results = TRUE
minimal_sample_size = 20
multivar = F
ic_method = "AIC"
max_estimation_tries = 100
series_split_method = "both"
na_method = "kalman_filter"


data.frame(n = 1:length(dataset_names), 
           ds = dataset_names)


for(ds in  dataset_names){
  print(ds)
  fit_hmms(data = read.csv(paste0("Data/CleanRescaled/", "Rescaled_", ds, ".csv")), 
           variables = variables_list[[ds]],
           max_numb_states = 5, 
           oos_validation_methods = "percentage", #Use a percentage or a specific value to create test data for out of sample validataion 
           oos_validation_size = .2, #The percentage of data or number of observation for the test dataset. 
           study_name = ds,
           output_location = "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/",
           return_results = FALSE,
           minimal_sample_size = 20,
           multivar = FALSE,
           ic_method = "AIC",
           max_estimation_tries = 100,
           series_split_method = "days", 
           na_method = "kalman_filter"
  )
}








x <- readRDS("Output/N1_Results/Rescaled_results/HMM_kalman_filter_results_2/HMMres_Bernstein_2019.rds")



variables_list <- readRDS("Data/variables_list.rds")
f <- list.files("Data/CleanRescaled/")
ds_names <- gsub("Rescaled_|.csv", "", f)

ds_names <- ds_names[!(ds_names %in% c("Kroencke_2020", "Wendt_2020"))]


timeKF <- readRDS("Output/HMMfitTimes/")

ds <- ds_names[37]

timeKF <- data.frame(ds = ds_names, time = rep(NA, length(ds_names)))
#Kalman Filter
for(ds in ds_names[36:37]){
  print(ds)
  startTime <- Sys.time()
  fit_hmms(data = read.csv(paste0("Data/CleanRescaled/", "Rescaled_", ds, ".csv")), 
           variables = variables_list[[ds]],
           max_numb_states = 5, 
           oos_validation_methods = "percentage", #Use a percentage or a specific value to create test data for out of sample validataion 
           oos_validation_size = .2, #The percentage of data or number of observation for the test dataset. 
           study_name = ds,
           output_location = "Output/N1_Results/RescaledResults/HMM_kalman_filter_results/",
           return_results = FALSE,
           minimal_sample_size = 20,
           multivar = FALSE,
           ic_method = "AIC",
           max_estimation_tries = 100,
           series_split_method = "days", 
           na_method = "kalman_filter"
  )
  endTime  <- Sys.time()
  timeKF[timeKF$ds == ds, 2] <- difftime(endTime, startTime, units = "mins")
  print(endTime - startTime)
}
saveRDS(timeKF, "Output/HMMfitTimes/timeKF.rds")

timeOD <- data.frame(ds = ds_names, time = rep(NA, length(ds_names)))
#Only split on days
for(ds in ds_names){
  print(ds)
  startTime <- Sys.time()
  fit_hmms(data = read.csv(paste0("Data/CleanRescaled/", "Rescaled_", ds, ".csv")), 
           variables = variables_list[[ds]],
           max_numb_states = 5, 
           oos_validation_methods = "percentage", #Use a percentage or a specific value to create test data for out of sample validataion 
           oos_validation_size = .2, #The percentage of data or number of observation for the test dataset. 
           study_name = ds,
           output_location = "Output/N1_Results/RescaledResults/HMM_onlyday_split_results/",
           return_results = FALSE,
           minimal_sample_size = 20,
           multivar = FALSE,
           ic_method = "AIC",
           max_estimation_tries = 100,
           series_split_method = "days"
  )
  endTime  <- Sys.time()
  timeOD[timeOD$ds == ds, 2] <- difftime(endTime, startTime, units = "mins")
  print(endTime - startTime)
  
}
saveRDS(timeOD, "Output/HMMfitTimes/timeOD.rds")


#Split both
timeDM <- data.frame(ds = ds_names, time = rep(NA, length(ds_names)))
for(ds in ds_names){
  print(ds)
  startTime <- Sys.time()
  fit_hmms(data = read.csv(paste0("Data/CleanRescaled/", "Rescaled_", ds, ".csv")), 
           variables = variables_list[[ds]],
           max_numb_states = 5, 
           oos_validation_methods = "percentage", #Use a percentage or a specific value to create test data for out of sample validataion 
           oos_validation_size = .2, #The percentage of data or number of observation for the test dataset. 
           study_name = ds,
           output_location = "Output/N1_Results/RescaledResults/HMM_daymissing_split_results/",
           return_results = FALSE,
           minimal_sample_size = 20,
           multivar = FALSE,
           ic_method = "AIC",
           max_estimation_tries = 100,
           series_split_method = "both"
  )
  endTime  <- Sys.time()
  timeDM[timeDM$ds == ds, 2] <- difftime(endTime, startTime, units = "mins")
  print(endTime - startTime)
  
}
saveRDS(timeDM, "Output/HMMfitTimes/timeDM.rds")


data = read.csv(paste0("Data/CleanRescaled/Rescaled_", ds, ".csv"))
variables = variables_list[[ds]]
max_numb_states = 5
oos_validation_methods = "percentage" #Use a percentage or a specific value to create test data for out of sample validataion 
oos_validation_size = .2 #The percentage of data or number of observation for the test dataset. 
study_name = ds
output_location = "Output/N1_HMM_Results_V2/"
return_results = TRUE
minimal_sample_size = 20
multivar = F
ic_method = "AIC"
max_estimation_tries = 100
series_split_method = "both"
