findBlocks<-function(data,numBlocks){
  #Find the block each point in data belongs to,
  #assuming that the unit square can be partitioned into numBlocks blocks
  #and adds a column of indices of the block the point belongs to
  i <- floor(data$x * sqrt(numBlocks)) + 1
  j <- floor(data$y * sqrt(numBlocks)) + 1
  
  data$block <- i + sqrt(numBlocks) * (j - 1)
  
  return(data)
}

block_CV_parallell <- function(data, numBlocks = 25){
  # Function that calculates predictions and prediction variance for data points
  # according to a block CV with number of blocks equal to numBlocks. The function
  # also reorders data according to the block they belong to, and adds a column to
  # the dataset (together with pred and pred_var), which identifies the block 
  # observations belong to
  
  data <- findBlocks(data, numBlocks)
  
  ordered_data <- data[order(data$block), ]
  
  ordered_data$originalId <- rownames(ordered_data)
  
  obs_per_fold <- as.numeric(table(factor(ordered_data$block, levels = 1:numBlocks))) 
  #Counts how many observations there are in each block
  
  cum_sums <- cumsum(obs_per_fold)
  #Needed to get the starting index of the first point in each block
  
  start_index_per_fold <- c(1, cum_sums[-length(cum_sums)] + 1)
  
  preds_vec <- c()
  
  pred_vars_vec <- c()
  
  data_cov_mat <- calc_cov_mat_data(ordered_data[, 1:2])
  
  results<-foreach(foldId = 1:numBlocks, .export = c("make_predictions")) %dopar% {
    
    if(obs_per_fold[foldId] == 0){
      return #next, as there are no observations in this block
    } else{
      if(foldId == numBlocks){
        # First block
        indices_in_fold <- start_index_per_fold[foldId]:nrow(data)
      } else{
        indices_in_fold <- start_index_per_fold[foldId]:(start_index_per_fold[foldId + 1] - 1)  
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
  }
  #Unpacking results
  for (foldId in 1:numBlocks) {
    if(obs_per_fold[foldId] != 0){ #only have results for blocks that 
      #contain observations
      preds_vec <- c(preds_vec, results[[foldId]]$preds)
      pred_vars_vec <- c(pred_vars_vec, results[[foldId]]$pred_var)  
    }
    
  }
  
  ordered_data$pred <- preds_vec
  ordered_data$pred_var <- pred_vars_vec
  
  return(ordered_data[order(as.integer(rownames(ordered_data))),])
}
