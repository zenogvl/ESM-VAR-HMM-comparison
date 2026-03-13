

predicting_hmms <- function(object){
  
  prediction_distribution <- lapply(1:object@nstates, matrix, data = NA, nrow = sum(object@ntimes), ncol = object@nresp)
  variable_check <- vector(length = object@nresp)
  
  for(s in 1:object@nstates){
    for(v in 1:object@nresp){
      prediction_distribution[[s]][,v] <-  predictManual(object@response[[s]][[v]])
      if(s == 1) {
        variable_check[v] <- as.character(object@response[[1]][[v]]@formula[[2]]) 
      }
    }
  }
  
  prediction <- matrix(NA, nrow = sum(object@ntimes), ncol = object@nresp)
  colnames(prediction) <- variable_check
  for(i in 1:sum(object@ntimes)){
    prediction[i,] <- prediction_distribution[[object@posterior$state[i]]][i,]
  }

  output <- list(prediction = prediction,
                 variable_check = variable_check, 
                 prediction_distributiontribution = prediction_distribution)
  output
}


