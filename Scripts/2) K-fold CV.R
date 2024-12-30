library(geoR)

calc_cov_mat_data <- function(data_locs, sigma=1, aRange=0.05, kappa=1, nugget=0.1){
  # Function to calculate the covariance matrix of a dataset with locations
  # given by data_locv, according to a Matern covariance function with parameters
  # as above, and nugget variance
  
  # This is needed in case anyone (read, me...) should be so dumb as to 
  # include the response when calculating the data covariance matrix
  if(ncol(data_locs)!=2){
    return("You have probably included the response...")
  }
  dist_mat <- as.matrix(dist(data_locs))
  
  dist_mat <- cov.spatial(dist_mat, cov.pars = c(sigma^2, aRange), kappa = kappa)
  
  diag(dist_mat) <- sigma^2 + nugget
  
  cov_mat_data <- dist_mat
  
  return(cov_mat_data)
}

calc_cross_cov_mat <-function(data_locs, pred_locs, aRange=0.05, kappa=1, sigma=1){
  # Function that calculates the cross covariance matrix between data 
  # and prediction locations. Note that this is equivalent to 
  # calc_cov_mat_data(data_locs, data_locs), but then without nugget variance
  if(ncol(data_locs)!=2 | ncol(pred_locs)!=2){
    return("You've included more than needed somewhere")
  }
  dist_mat <- as.matrix(cdist(data_locs, pred_locs))
  
  cross_cov_mat <- cov.spatial(dist_mat, cov.pars = c(sigma^2, aRange), kappa = kappa) #This takes approximately
  #6 times as long time as calculating cdist()
  
  return(cross_cov_mat)
}

calc_cross_cov_mat_parallel <- function(data_locs, pred_locs, aRange = 0.05, sigma = 1, kappa = 1, chunk_size = 1500) {
  # Function that calculates the, often large, cross_covariance matrix between
  # the data and prediction locations, but in parallell using foreach.
  n_data <- nrow(data_locs)
  n_pred <- nrow(pred_locs)
  
  num_chunks <- ceiling(n_pred / chunk_size)
  
  cross_cov_mat_list <- foreach(chunk = 1:num_chunks, .packages = c('fields', 'geoR', "rdist")) %dopar% {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_pred)
    
    chunk_pred_locs <- pred_locs[start_idx:end_idx, ]
    
    dist_mat <- as.matrix(cdist(data_locs, chunk_pred_locs))
    
    chunk_cov_mat <- cov.spatial(dist_mat, cov.pars = c(sigma^2, aRange), kappa = kappa)
    
    return(chunk_cov_mat)  # Return the chunk of the covariance matrix
  }
  
  # Combine the chunks into the full cross-covariance matrix
  cross_cov_mat <- do.call(cbind, cross_cov_mat_list)
  
  return(cross_cov_mat)
}



make_predictions<-function(cov_mat_data, cross_cov_mat, data,
                           sigma=1, nugget=0.1){
  # Function that makes predictions given a data covariance and a cross
  # covariance matrix, and some data. Sigma denotes the marginal variance.
  L<-chol(cov_mat_data)
  
  W<-backsolve(L, forwardsolve(t(L), cross_cov_mat))
  
  preds <- t(W)%*%data$z
  
  if(is.null(ncol(cross_cov_mat))){ #This is a vector when nrow=1, and
    #then we must reshape to a whatever matrix by n matrix
    cross_cov_mat <- matrix(cross_cov_mat, ncol=1, nrow=length(cross_cov_mat))
  }
  
  pred_var <- rep(sigma^2, ncol(cross_cov_mat))+nugget - colSums(cross_cov_mat * W) 
  
  result <- list(preds = preds, pred_var = pred_var)
  
  return(result)
}

krige <- function(data, pred_locations, sigma=1){
  #Function to do kriging for the prediction locations using the data
  #which contains coordinates and responses (z)
  
  data_locations <- data[, c(1,2)]
  
  cov_mat_data <- calc_cov_mat_data(data_locations)
  
  cross_cov_mat <- calc_cross_cov_mat_parallel(data_locations, pred_locations)
  
  preds<-make_predictions(cov_mat_data, cross_cov_mat, data)
  
  return(preds)
}

CRPS<-function(mu, sigma2, x){
  #Function to calculate the CRPS given a prediction (mean)
  #mu, a prediction variance, sigma2, and the response, x
  sigma <- sqrt(sigma2)
  standard <- (x-mu) / sigma
  -(sigma * (1 / sqrt(pi)-
           2 * dnorm(standard)-
           standard * (2 *pnorm(standard)-1)))
}
