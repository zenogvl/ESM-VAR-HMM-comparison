
create_depmix_predict_object <- function(object, data_new, length_time_serie){
  
  object@ntimes <- length_time_serie
  
  #Change the dens array of the object 
  new_densities <- array(NA, dim = c(nrow(data_new), object@nresp, object@nstates))
 
  for(s in 1:object@nstates){
    for(v in 1:object@nresp){
      
      #Create new response object
      response <- object@response[[s]][[v]]
      
      #Change outcome data
      response@y <- as.matrix(data_new[,as.character(response@formula[[2]])]) 
      
      #Change desing matrix
      newX <- matrix(NA, nrow = nrow(data_new), ncol = ncol(response@x))
      colnames(newX) <- colnames(response@x)  
      newX[,1] <- rep(1, nrow(data_new)) #Intercept 
      
      p <- colnames(response@x)[-1][1]
      #For multivariate change the data for the other variables
      if(ncol(response@x) > 1){
        for(p in colnames(response@x)[-1]){
          newX[,p] <- data_new[,p]
        }
      }
      response@x <- newX
      
      object@response[[s]][[v]] <- response
      
      #Use new response object to get the dens
      new_densities[,v,s] <-  densManual(object@response[[s]][[v]])
    }
  }
  object@dens <- new_densities
  
  object
}

densManual <- function(object,log=FALSE) {
  dnorm(x=object@y,mean=predictManual(object),sd=object@parameters$sd,log=log)
}

predictManual <- function(object) {
  object@x%*%object@parameters$coefficients
}

