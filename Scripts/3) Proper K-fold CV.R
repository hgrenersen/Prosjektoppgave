K_fold_CV <- function(data, nFolds){
  # Function to perform K-fold CV given a dataset.
  #   - data, matrix/dataframe with the first two columns
  # representing the x and y coordinates, and the last representing the actual
  # response in the point
  #   -nFolds, integer giving the number of folds. If the number of data points
  # is not divisible by the number of folds, general_K is called to account for this
  
  # returns: The input data ordered according to the folds, with the additional columns
  # preds and pred_var, which contain the predictions and prediction variance 
  # calculated for the point. If general_K is called, the column foldId is also
  # added to show which fold a point belongs to. 
  num_of_data_points <- nrow(data)
  
  if (num_of_data_points %% nFolds == 0) { # The number of data points is divisible
  # by the number of folds. Set up, ordering the data according to folds
  # and calculating the covariance matrix.
    
    indices_in_fold <- sample(1:num_of_data_points, num_of_data_points)
    
    ordered_data <- data[indices_in_fold, ]
    
    ordered_data$originalId <- rownames(ordered_data)
    
    num_per_fold <- num_of_data_points / nFolds
    
    pred_mat <- matrix(nrow = num_per_fold, ncol = nFolds)
    
    pred_var_mat <- matrix(nrow = num_per_fold, ncol = nFolds)
    
    data_cov_mat <- calc_cov_mat_data(ordered_data[, 1:2])
    
    #The actual K-fold CV
    for (foldId in 1:nFolds) {
      
      indices_in_fold <- (num_per_fold * (foldId - 1) + 1):((num_per_fold * (foldId)))
      
      remaining_data <-ordered_data[-indices_in_fold,]
      
      sub_cross_cov_mat<-data_cov_mat[-indices_in_fold,
                                      indices_in_fold]
      
      sub_data_cov_mat <- data_cov_mat[-indices_in_fold,
                                       -indices_in_fold]
      
      result <- make_predictions(sub_data_cov_mat, sub_cross_cov_mat, remaining_data)
      
      preds <- result[[1]]
      
      pred_var <- result[[2]]
      
      pred_mat[, foldId] <- preds
      
      pred_var_mat[, foldId] <- pred_var
    }
    
    ordered_data$pred <- c(pred_mat)
    ordered_data$pred_var <- c(pred_var_mat)
    
    return(ordered_data[order(as.integer(rownames(ordered_data))),])
  }
  else{
    return(general_K(data, nFolds))
  }
}

general_K <- function(data, nFolds){
  # Virtually the same as K_fold_CV, but accounting for the case where
  # the number of data points isn't divisible by nFolds.
  # Here we'll first distribute the n//k=n%/%k points along the k folds,
  # and then the remainder will be distributed along the k folds
  # randomly.
  num_of_data_points <- nrow(data)
  
  indices_in_fold <- sample(1:num_of_data_points, num_of_data_points %/% nFolds * nFolds)
  
  ordered_data <- data[indices_in_fold,]
  
  ordered_data$foldId <-  sort(rep(1:nFolds,num_of_data_points %/% nFolds))
  
  surplus_points <- data[-indices_in_fold, ]
  if (nrow(surplus_points) != 0) {
    larger_folds <- sample(1:nFolds, num_of_data_points %% nFolds) #maximally has length
    #nFolds - 1
    
    for (i in 1:nrow(surplus_points)) {
      ordered_data[nrow(ordered_data)+1, ] <- surplus_points[i, ]
      ordered_data[nrow(ordered_data), ]$foldId <- larger_folds[i]
    }
    
    ordered_data <- ordered_data[order(ordered_data$foldId),] 
  }
  else{
    larger_folds <- c()
  }
  
  ordered_data$originalId <- rownames(ordered_data)
  
  num_per_fold <- num_of_data_points%/%nFolds #for the smaller folds
  
  preds_list <- list()
  
  pred_vars_list <- list()
  
  data_cov_mat <- calc_cov_mat_data(ordered_data[,1:2])
  
  first_data_ind_in_fold <- 1
    
  last_data_ind_in_fold <- first_data_ind_in_fold + num_per_fold - 1
  
  for (foldId in 1:nFolds) {
    if (foldId %in% larger_folds) {
      last_data_ind_in_fold <- last_data_ind_in_fold + 1
    }
    indices_in_fold <- first_data_ind_in_fold:last_data_ind_in_fold
    
    remaining_data <- ordered_data[-indices_in_fold,]
    
    sub_cross_cov_mat <- data_cov_mat[-indices_in_fold,
                                    indices_in_fold]
    
    sub_data_cov_mat <- data_cov_mat[-indices_in_fold, 
                                     -indices_in_fold]
    
    result <- make_predictions(sub_data_cov_mat, sub_cross_cov_mat, remaining_data)
    
    preds <- result[[1]]
    
    pred_var <- result[[2]] 
    
    preds_list[[foldId]] <- preds
    
    pred_vars_list[[foldId]] <- pred_var  
    
    first_data_ind_in_fold <- last_data_ind_in_fold + 1
    
    last_data_ind_in_fold <- last_data_ind_in_fold + num_per_fold 
  }
  
  ordered_data$pred <- unlist(preds_list)
  ordered_data$pred_var <- unlist(pred_vars_list)
  
  return(ordered_data[order(as.integer(rownames(ordered_data))),])
}

K_fold_CV_parallell <- function(data, nFolds){
  # Parallellized version of K_fold_CV using foreach-
  num_of_data_points <- nrow(data)
  
  if (num_of_data_points %% nFolds == 0) {
    
    indices_in_fold <- sample(1:num_of_data_points, num_of_data_points)
    
    ordered_data <- data[indices_in_fold, ]
    
    num_per_fold <- num_of_data_points / nFolds
    
    preds_vec <- c()
    
    pred_vars_vec <- c()
    
    data_cov_mat <- calc_cov_mat_data(ordered_data[, 1:2])
    
    results<-foreach(foldId = 1:nFolds, .export = c("make_predictions")) %dopar% {
      
      indices_in_fold <- (num_per_fold * (foldId - 1) + 1):((num_per_fold * (foldId)))
      fold <- ordered_data[indices_in_fold,]
      
      remaining_data <-ordered_data[-indices_in_fold,]
      
      sub_cross_cov_mat<-data_cov_mat[-indices_in_fold,
                                      indices_in_fold]
      
      sub_data_cov_mat <- data_cov_mat[-indices_in_fold, 
                                       -indices_in_fold]
      
      result <- make_predictions(sub_data_cov_mat, sub_cross_cov_mat, remaining_data)
      
      preds <- result[[1]]
      
      pred_var <- result[[2]] 
      
      list(preds = preds, pred_var = pred_var) 
    }
    
    for (foldId in 1:nFolds) {
      preds_vec <- c(preds_vec, results[[foldId]]$preds)
      pred_vars_vec <- c(pred_vars_vec, results[[foldId]]$pred_var)
    }
    
    ordered_data$pred <- c(preds_vec)
    ordered_data$pred_var <- c(pred_vars_vec)
    
    return(ordered_data[order(as.integer(rownames(ordered_data))),])
  }
  else{
    return(general_K_parallell(data, nFolds))
  }
}

general_K_parallell<-function(data, nFolds){
  # Parallellized version of general_K
  num_of_data_points <- nrow(data)
  
  indices_in_fold <- sample(1:num_of_data_points, num_of_data_points%/%nFolds*nFolds)
  
  ordered_data <- data[indices_in_fold,]
  
  ordered_data$foldId <-  sort(rep(1:nFolds,num_of_data_points%/%nFolds))
  
  surplus_points <- data[-indices_in_fold, ]
  if(nrow(surplus_points) != 0){
    larger_folds <- sample(1:nFolds, num_of_data_points %% nFolds) #maximally has length
    #nFolds - 1
    
    for (i in 1:nrow(surplus_points)) {
      ordered_data[nrow(ordered_data) + 1, ] <- surplus_points[i, ]
      ordered_data[nrow(ordered_data), ]$foldId <- larger_folds[i]
    }
    
    ordered_data <- ordered_data[order(ordered_data$foldId),] 
  }
  else{
    larger_folds <- c()
  }
  
  preds_vec <- c()
  
  pred_vars_vec <- c()
  
  data_cov_mat <- calc_cov_mat_data(ordered_data[ , 1:2])
  
  obs_per_fold <- as.numeric(table(factor(ordered_data$foldId, levels = 1:nFolds)))
  
  cum_sums <- cumsum(obs_per_fold)
  
  start_index_per_fold <- c(1, cum_sums[-length(cum_sums)] + 1)
  
  results <- foreach(foldId = 1:nFolds, .export = c("make_predictions")) %dopar% {
    
    if (foldId == nFolds) {
      indices_in_fold <- start_index_per_fold[foldId]:nrow(data)
    }
    else {
      indices_in_fold<- start_index_per_fold[foldId]:(start_index_per_fold[foldId+1]-1)
    }
    
    remaining_data <- ordered_data[-indices_in_fold,]
    
    sub_cross_cov_mat <- data_cov_mat[-indices_in_fold,
                                    indices_in_fold]
    
    sub_data_cov_mat <- data_cov_mat[-indices_in_fold, 
                                     -indices_in_fold]
    
    result <- make_predictions(sub_data_cov_mat, sub_cross_cov_mat, remaining_data)
    
    preds <- result[[1]]
    
    pred_var <- result[[2]]
    
    list(preds = preds, pred_var = pred_var)
  }
  
  for (foldId in 1:nFolds) {
    preds_vec <- c(preds_vec, results[[foldId]]$preds)
    pred_vars_vec <- c(pred_vars_vec, results[[foldId]]$pred_var)
  }
  
  ordered_data$pred <- c(preds_vec)
  ordered_data$pred_var <- c(pred_vars_vec)
  
  return(ordered_data[order(as.integer(rownames(ordered_data))),])
}


calculate_score <- function(pred_mat, weights=1, func = "CRPS"){
  # Funtion to calculate the scores for a matrix of predictions
  # pred_mat: location(x,y), actual value, prediction, prediction variance
  # weights: weights for the individual scores when calculating
  # the aggregated score
  # func: score function
  if (func == "CRPS") {
    scores <- CRPS(pred_mat$pred, pred_mat$pred_var, pred_mat$z)
  }
  if (func == "SQE") {
    scores <- (pred_mat$z - pred_mat$pred)^2
  }
  
  return(c(weights %*% scores))
}

make_parallell_predictions<-function(cov_mat_data, cross_cov_mat, data,
                                     sigma = 1, nugget = 0.1, chunk_size=1000){
  # Parallelized version to make predictions at the locations
  # given by the cross_cov_mat
  # cov_mat_data: covariance matrix of all the data locations
  # cross_cov_mat: cross covariance matrix between the data 
  # and prediction locations.
  # data: Matrix containing two columns of locations, along
  # with the response
  # sigma: Marginal variace
  # nugget: Nugget effect
  # chunk_size: Size of chunks of the data points processed
  # in parallell
  L <- chol(cov_mat_data)
  
  n_data <- nrow(cov_mat_data)
  n_pred <- ncol(cross_cov_mat)
  
  num_chunks <- ceiling(n_pred / chunk_size)
  
  W_mat_list <- foreach(chunk = 1:num_chunks) %dopar% {
    
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_pred)
    
    chunk_cross_cov <- cross_cov_mat[, start_idx:end_idx]
    
    chunk_W <- W<-backsolve(L, forwardsolve(t(L), chunk_cross_cov))
    
    return(chunk_W)  # Return the chunk of the covariance matrix
  }
  
  W <- do.call(cbind, W_mat_list)
  
  preds <- t(W) %*% data$z
  
  pred_var <- rep(sigma^2, ncol(cross_cov_mat)) + nugget - colSums(cross_cov_mat * W)
  
  result <- list(preds = preds, pred_var = pred_var)
  
  return(result)
}

pred_grid <- function(data, grid=exp_grid){
  # Function to make predictions over an entire grid
  # data, matrix/dataframe with three columns, the first two
  # representing the x and y coordinates of the data point,
  # and the third representing the data value in the point
  
  # grid, grid to predict over
  
  data_cov_mat <- calc_cov_mat_data(data[ , c(1,2)])
  
  cross_cov_mat <- calc_cross_cov_mat_parallel(data[ , c(1,2)],
                                               grid)
  return(make_parallell_predictions(data_cov_mat, cross_cov_mat, data))
}