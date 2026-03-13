
# 
# for(ds in dataset_names){
#   for(object_name in model_names){
#     
#   }
# }



collect_rmse_data <- function(path_list, 
                              dataset_names,
                              model_names){
  for(object_name in model_names){
    data_input <- list()
    for(ds in dataset_names){
      specific_path <- paste0(path_list[[which(model_names == object_name)]], "RMSE_", ds, ".csv")
      data_input[[ds]] <- read.csv(specific_path) %>% 
        na.omit()
    }
    assign(paste0("rmse_", object_name), data_input, envir = .GlobalEnv)
  } 
}

get_intersect_fitted_models <- function(dataset_names, model_names){
  pp_all_models_fitted <- list()
  for(ds in dataset_names){
    pp_fitted_models <- list()
    for(object_name in model_names){
      pp_fitted_models[[object_name]] <- get(paste0("rmse_", object_name))[[ds]][,"ID"]
    }
    pp_all_models_fitted[[ds]] <- Reduce(intersect ,pp_fitted_models )
  }
  return(pp_all_models_fitted)
}


find_lowest_rmse <- function(path_list, 
                             dataset_names,
                             model_names = list(),
                             show_equal_results = FALSE,
                             format_output = "long"){
  
  if(length(path_list) != length(model_names)){
    warning("No names for models")
    model_names <- as.character(1:length(path_list))
  }
  if(show_equal_results & length(path_list) > 2) {
    warning("Equal results only possible with two models")
    show_equal_results <- FALSE
  }
  
  collect_rmse_data(path_list, dataset_names, model_names) #Assigns to global enviroment 
  
  #get the individuals where the model is fitted everywhere. 
  pp_all_models_fitted <- get_intersect_fitted_models(dataset_names, model_names)
  
  lowest_rmse_list <- list()
  for(ds in dataset_names){
    if(length(pp_all_models_fitted[[ds]]) < 1) next
    for(object_name in model_names){
      assign(paste0("mean_rmse_", object_name), 
             get(paste0("rmse_", object_name))[[ds]]  %>%
               select(-X) %>%
               filter(ID %in% pp_all_models_fitted[[ds]]) %>%
               rowwise() %>%
               mutate(!!paste0("mean_", object_name) := mean(c_across(!ID))) %>%
               select(ID, paste0("mean_", object_name))
      ) 
    }
    if(show_equal_results){
      lowest_rmse_list[[ds]] <- mget(paste0("mean_rmse_", model_names)) %>%
        purrr::reduce(inner_join, by = "ID") %>%
        mutate(
          rmse_values = list(c_across(-ID)),
          unique_rmse = length(unique(rmse_values)),
          lowest_rmse = if_else(unique_rmse == 1, NA_integer_, which.min(rmse_values)),
          is_tie = unique_rmse < length(model_names) &
            unique_rmse > 1
        ) %>%
        mutate(lowest_rmse_factor = case_when(is.na(lowest_rmse) ~ "all_equal", TRUE ~ model_names[lowest_rmse])) %>%
        mutate(lowest_rmse_factor = factor(lowest_rmse_factor, levels = c(model_names, "all_equal"))) %>%
        ungroup() %>%
        count(lowest_rmse_factor) %>%
        complete(lowest_rmse_factor = factor(
          c(model_names, "all_equal"),
          levels = c(model_names, "all_equal")
        ),
        fill = list(n = 0)) %>%
        mutate(dataset = ds, lowest_rmse_rate = n / sum(n))
      
    } else {
      lowest_rmse_list[[ds]] <- mget(paste0("mean_rmse_", model_names)) %>%
        purrr::reduce(inner_join, by = "ID")  %>%
        mutate(ndistinct = n_distinct(c_across(!ID))) %>%
        filter(ndistinct == length(model_names)) %>%
        select(-ndistinct)  %>%
        mutate(lowest_rmse_factor = factor(
          which.min(c_across(!ID)),
          levels = 1:length(model_names),
          labels = model_names
        )) %>%
        ungroup() %>%
        count(lowest_rmse_factor) %>%
        complete(
          lowest_rmse_factor = factor(model_names, levels = model_names),
          fill = list(n = 0)
        ) %>%
        mutate(dataset = rep(ds, n_distinct(lowest_rmse_factor))) %>%
        mutate(lowest_rmse_rate = n / sum(n)) 
      
    }
  }
  if(format_output == "long"){
    return(do.call("rbind", lowest_rmse_list))
  } else if(format_output == "matrix"){
    do.call("rbind", lowest_rmse_list) %>%
      select(!n) %>%
      pivot_wider(names_from = dataset, values_from = lowest_rmse_rate)  %>% 
      column_to_rownames("lowest_rmse_factor") %>%
      as.matrix() %>%
      return()
  } else {
    lowest_rmse_list
  }

}



calc_delta_rmse <- function(path_list,
                            dataset_names,
                            model_names = list()){
  if(length(path_list) != length(model_names)){
    warning("No names for models")
    model_names <- as.character(1:length(path_list))
  }
  
  cat(paste0("Delta RMSE is calculated by substracting ", model_names[2], " from ", model_names[1], 
             ": \n so a positive Delta RMSE indicates that ", model_names[2], 
             " results in beter model \n  while a negative Delta RMSE indicates that ", model_names[1], 
             "is the better model. "
  ))
  
  if(any(c(length(path_list) != 2, length(model_names) != 2))) stop("Delta rmse can only be caculated between two values.")
  
  
  collect_rmse_data(path_list, dataset_names, model_names) #Assigns to global enviroment 
  
  #get the individuals where the model is fitted everywhere. 
  pp_all_models_fitted <- get_intersect_fitted_models(dataset_names, model_names)
  
  delta_rmse_list <- list()
  for(ds in dataset_names){
    if(length(pp_all_models_fitted[[ds]]) < 1) next
    for(object_name in model_names){
      assign(paste0("mean_rmse_", object_name), 
             get(paste0("rmse_", object_name))[[ds]]  %>%
               select(-X) %>%
               filter(ID %in% pp_all_models_fitted[[ds]]) %>%
               rowwise() %>%
               mutate(!!paste0("mean_", object_name) := mean(c_across(!ID))) %>%
               select(ID, paste0("mean_", object_name)) 
      ) 
    }
    delta_rmse_list[[ds]] <- mget(paste0("mean_rmse_", model_names)) %>%
      purrr::reduce(inner_join, by = "ID") %>%
      rowwise() %>%
      mutate(delta_rmse = !!sym(paste0("mean_", model_names[1])) -  !!sym(paste0("mean_", model_names[2])), 
             dataset = ds) %>%
      ungroup()
    
  }
  return(do.call("rbind", delta_rmse_list))
}



















# 
# 
# ds <- dataset_names[3]
# 
# path_list <- list("Output/N1_Results/Rescaled_results/VAR_kalman_filter_RMSE/",
#                   "Output/N1_Results/Rescaled_results/VAR_listwise_deletion_RMSE/")
# 
# model_names <- c("kalman_filter", "listwise_deletion")
# 
# show_equal_results = TRUE
# 
# path_list = list("Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/RMSE/",
#                  "Output/N1_Results/Rescaled_results/HMM_no_update_results/RMSE/",
#                  "Output/N1_Results/Rescaled_results/HMM_seperate_series_results/RMSE/")
# 
# model_names <- c("Kalman", "no_upadata", "sep_series")
# 
# 
# dataset_names <- gsub("RMSE_|.csv", "", list.files("Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/RMSE/"))
# 



