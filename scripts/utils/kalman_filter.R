
kalman_filter <- function(x){
  #Stop if there are no missing values 
  if(!any(is.na(x))) return(x)
  
  emptyRow <- 0
  mis <- any(is.na(x[1,]))
  while(mis){
    emptyRow <- emptyRow + 1
    x <- x[-1,]
    mis <- any(is.na(x[1,]))
  }
  i <- 1
  for(i in 1:ncol(x)){
    if(length(unique(na.omit(x[,i]))) == 1){
      x[is.na(x[,i]),i] <- unique(na.omit(x[,i]))
      next
    }
    x[,i] <- as.numeric(x[,i]) 
    model <- stats::StructTS(x[,i], type = "level")$model0 #Get the in
    kalman <- stats::KalmanSmooth(x[,i], model, -1)
    x[is.na(x[,i]),i] <- kalman$smooth[is.na(x[,i]), , drop = FALSE]  %*% as.matrix(model$Z)
  }
  if(emptyRow > 0){
    missingMatrix <- matrix(NA, nrow = emptyRow, ncol = ncol(x))
    colnames(missingMatrix) <- colnames(x)
    x <- rbind(missingMatrix, x)
  }
  
  return(x)
}
