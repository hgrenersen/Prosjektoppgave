LOOCV <- function(data){
  # Function to perform LOOCV given a dataset
  # data, matrix/dataframe with the first two columns
  # representing the x and y coordinates, and the last
  # representing the actual data value in the point
  num_of_data_points <- nrow(data)
  
  rownames(data) <- NULL
  
  num_per_fold <- 1
  
  pred_mat <- matrix(nrow = num_per_fold, ncol = num_of_data_points)
  
  pred_var_mat <- matrix(nrow = num_per_fold, ncol = num_of_data_points)
  
  data_cov_mat <- calc_cov_mat_data(data[ , 1:2])
  
  results<-foreach(foldId = 1:num_of_data_points, .export = c("make_predictions")) %dopar% {
    
    indices_in_fold <- (num_per_fold * (foldId - 1) + 1):((num_per_fold * (foldId)))
    
    remaining_data <- data[-indices_in_fold, ]
    
    sub_cross_cov_mat <- as.matrix(data_cov_mat[-indices_in_fold,
                                              indices_in_fold])
    
    sub_data_cov_mat <- data_cov_mat[-indices_in_fold, 
                                     -indices_in_fold]
    
    result <- make_predictions(sub_data_cov_mat, sub_cross_cov_mat, remaining_data)
    
    preds <- result[[1]]
    
    pred_var <- result[[2]]
    
    list(preds = preds, pred_var = pred_var)
  }
  
  for (foldId in 1:num_of_data_points) {
    pred_mat[, foldId] <- results[[foldId]]$preds
    pred_var_mat[, foldId] <- results[[foldId]]$pred_var
  }
  
  data$pred <- c(pred_mat)
  data$pred_var <- c(pred_var_mat)
  
  return(data)
}

LOOCV_downdate <- function(data, sigma = 1, nugget = 0.1){
  # Function to perform LOOCV given a dataset
  # data, matrix/dataframe with the first two columns
  # representing the x and y coordinates, and the last
  # representing the actual data value in the point.
  # This function makes use of relations regarding the 
  # cholesky factor when removing a row and column
  # from matrix which you know the Cholesky factor of
  # through choldrop() from mgcv()
  num_of_data_points <- nrow(data)
  
  rownames(data) <- NULL 
  
  num_per_fold <- 1
  
  pred_mat <- matrix(nrow = num_per_fold, ncol = num_of_data_points)
  
  pred_var_mat <- matrix(nrow = num_per_fold, ncol = num_of_data_points)
  
  data_cov_mat <- calc_cov_mat_data(data[, 1:2])
  
  data_L <- chol(data_cov_mat)
  
  results <- foreach(foldId = 1:num_of_data_points, .export = c("make_predictions"),
                   .packages = c("mgcv")) %dopar% { 
       indices_in_fold <- foldId
       
       remaining_data <- data[-indices_in_fold,]
       
       sub_cross_cov_mat <- as.matrix(data_cov_mat[-indices_in_fold,
                                                 indices_in_fold])
       
       L_sub <- choldrop(data_L, foldId)
       
       W <- backsolve(L_sub, forwardsolve(t(L_sub), sub_cross_cov_mat))
       
       preds <- t(W) %*% remaining_data$z
       
       pred_var <- rep(sigma^2, ncol(sub_cross_cov_mat)) + nugget - colSums(sub_cross_cov_mat * W)
       
       list(preds = preds, pred_var = pred_var)
     }
  
  for (foldId in 1:num_of_data_points) {
    pred_mat[, foldId] <- results[[foldId]]$preds
    pred_var_mat[, foldId] <- results[[foldId]]$pred_var
  }
  
  data$pred <- c(pred_mat)
  data$pred_var <- c(pred_var_mat)
  
  return(data)
}
