# Libraries and loading required data/functions

library(progress)
library(doParallel)
library(rdist)
library(mgcv)

source("Scripts/2) K-fold CV.R")
source("Scripts/3) Proper K-fold CV.R")
source("Scripts/4) Large-scale K-fold CV.R")

gridded_data <- readRDS("Data/Gridded sets.rds")

realizations <- readRDS("Data/realizations.rds")

create_grid <- function(num_of_points){
  xseq <- seq( 0,1,,num_of_points+1)[-(num_of_points+1)]
  dx <- diff(xseq)[1]
  xseq<-xseq + 1/2*dx
  
  yseq <- seq( 0,1,,num_of_points+1)[-(num_of_points+1)]
  dy <- diff(yseq)[1]
  yseq<-yseq + 1/2*dx
  
  grid<- list( x= xseq, y= yseq)
  
  exp_grid <- expand.grid(grid)  
}

# Function to generate gridded data given some realization

generate_gridded_data <- function(realizations, realization_grid, data_grid,
                                  realization_index = 1, nugget = 0.1){
  # Realizations: List of realizations
  # realization_grid: List of coordinates on the grid that the 
  # realization is sampled on
  # data_grid: grid the data should we want data to be on
  # realization_index: which realization to use
  # nugget: sigma_N^2 (nugget effect)

  cross_distance <- cdist(data_grid, realization_grid)
  
  closest_realization_centroid <- apply(cross_distance, 1, which.min)
  
  data_grid$z <- c(realizations[[realization_index]])[closest_realization_centroid] + rnorm(nrow(data_grid), sd=sqrt(nugget)) 
  
  return(data_grid)

}

expectedCRPS <- function(truth, truth.var, est, est.var, getAverage=TRUE, na.rm=FALSE, weights=rep(1, length(truth))) {
  # Function to calculate the ECRPS, code by John L. Paige with minor adjustments
  # truth: True values
  # truth.var: Nugget effect
  # est: Prediction
  # est.var: Prediction variance
  # getAvearge: Bool to decide if an average should be calculated
  # na.rm: Bool to decide if NAns should be removed
  # weights: vector of scalars that might be used for a weighted average
  erf <- function(x) {2 * (pnorm(sqrt(2) * x) - 0.5)}
  
  SSS <- truth.var + est.var # sigStarSq
  b <- est-truth # bias
  
  pt1<-sqrt(2*SSS/pi) * exp(-b^2/(2*SSS)) + b*erf(b/sqrt(2*SSS))
  pt2<-sqrt(est.var)/sqrt(pi)
  
  res<-pt1 - pt2
  
  if(getAverage) {
    weights<-weights*(1/sum(weights, na.rm=TRUE))
    sum(res*weights, na.rm=na.rm)
  } else {
    res
  }
}

num_cores <- 10 
cl <- makeCluster(num_cores)
registerDoParallel(cl)

LOOCV_and_ECRPS_for_grid <- function(m, realization_grid = create_grid(200),
num_of_datasets = 500){
  # Function to calculate the LOOCV estimate and the average ECRPS over a grid
  #m : m by m grid
  # Things needed for both LOOCV and true error:
  # datasets, data_grid
  # data covaraince matrix,
  # cholesky factor of the above
  # variance parameters
  data_grid <- create_grid(m)

  num_of_datapoints <- nrow(data_grid)

  grid_datasets <- list()

  # Start by generating data
  set.seed(4)
  print("Generating data")
  pb = progress_bar$new(
      format = paste0("Generating data for grid", m, ": [:bar] :percent eta: :eta"),
      total = num_of_datasets, clear = FALSE, width= 60)
  for(i in 1:num_of_datasets){
    grid_datasets[[i]] <- generate_gridded_data(realizations, realization_grid,
                                    data_grid, realization_index = i)
                                    pb$tick()
  }
  saveRDS(grid_datasets, paste0("Data/Gridded sets/Grid ", m, " by ", m, ".rds"))

  # Preparing things needed to calculating the LOOCV estimate.
  # Much of this code is just the LOOCV function, which is arguably better
  # documented, but with some minor changes to make use of the fact that
  # the weights are equal for all datasets as the sampling scheme is deterministic

  data_cov_mat <- calc_cov_mat_data(data_grid[, 1:2])

  data_L<-chol(data_cov_mat)

  sigma <- 1
  nugget <- 0.1

  LOOCV_result_list <- list()

  for(i in 1:num_of_datasets){
    LOOCV_result_list[[i]] <- matrix(rep(0, num_of_datapoints*2), ncol = 2)
  }

  pb = progress_bar$new(
      format = paste0("LOOCV over grid", m, ": [:bar] :percent eta: :eta"),
      total = num_of_datapoints, clear = FALSE, width= 60)
  
  for(i in 1:num_of_datapoints){
    indices_in_fold<- i
    fold <- data_grid[indices_in_fold,]
    
    remaining_data <-data_grid[-indices_in_fold,]
    
    sub_cross_cov_mat<-as.matrix(data_cov_mat[-indices_in_fold,
                                              indices_in_fold])
    
    sub_data_cov_mat <- data_cov_mat[-indices_in_fold, 
                                      -indices_in_fold]
    #Downdating to get the Cholesky decomposition of 
    L_sub<-choldrop(data_L, i)
    
    W<-backsolve(L_sub, forwardsolve(t(L_sub), sub_cross_cov_mat))
    
    pred_var <- rep(sigma^2, ncol(sub_cross_cov_mat))+nugget - colSums(sub_cross_cov_mat * W)
    WT <- t(W)
    remove(W)
    
    # Doing the actual prediction

    pred_list<-foreach(foldId = 1:num_of_datasets, .export = c("make_predictions"),
                  .packages = c("mgcv")) %dopar% {
      
      pred <- WT%*%grid_datasets[[foldId]][-indices_in_fold,]$z
      
      pred
    }
    pred_for_point <- unlist(pred_list)
    for(j in 1:num_of_datasets){
      LOOCV_result_list[[j]][i, 1] <- pred_for_point[j]
      LOOCV_result_list[[j]][i, 2] <- pred_var
    }
    pb$tick()
  }
  saveRDS(LOOCV_result_list, paste0("Results/Gridded/LOOCV Results on ", m, " by ", m, " grid.rds"))

  # Using the LOOCV predictions to calculate the average CRPS and thus
  # the actual LOOCV estimate

  CRPS_vals <- rep(0, num_of_datasets)

  for(i in 1:num_of_datasets){
    pred <- LOOCV_result_list[[i]][,1]
    pred_var <- LOOCV_result_list[[i]][, 2]
    CRPS_vals[i] <- mean(CRPS(pred, pred_var, c(grid_datasets[[i]]$z)))  
  }
  saveRDS(CRPS_vals, paste0("Results/Gridded/LOOCV CRPS on ", m, " by ", m, " grid.rds"))

  # True error over the domain, resuses things calculated above,
  # to use the data to make predictions over a finer grid and
  # then calculating the ECRPS and averaging
  ECRPS_vals <- rep(0, num_of_datasets)

  cross_cov_mat <- calc_cross_cov_mat_parallel(data_grid,
                                             realization_grid)

  n_data <- nrow(data_cov_mat)
  n_pred <- ncol(cross_cov_mat)

  chunk_size <- 1000

  num_chunks <- ceiling(n_pred / chunk_size)

  W_mat_list <- foreach(chunk = 1:num_chunks) %dopar% {
      
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, n_pred)
      
      chunk_cross_cov <- cross_cov_mat[, start_idx:end_idx]
      
      chunk_W <- backsolve(data_L, forwardsolve(t(data_L), chunk_cross_cov))
      
      return(chunk_W)  # Return the chunk of the weights matrix
    }

  W <- do.call(cbind, W_mat_list)
  print("Weights calculated")
  remove(W_mat_list)
  remove(data_L)

  pred_var_true_error <- rep(1, ncol(cross_cov_mat))+nugget - colSums(cross_cov_mat * W)

  WT <- t(W)
  remove(W) #Don't need this as I have WT
  remove(cross_cov_mat)
  remove(data_cov_mat)

  for(i in 1:num_of_datasets){
    pred <- WT%*%grid_datasets[[i]]$z

    ECRPS_vals[i] <- expectedCRPS(c(realizations[[i]]), nugget,
                            pred, pred_var_true_error)
  }
  saveRDS(ECRPS_vals, paste0("Results/Gridded/ECRPS on ", m, " by ", m, " grid.rds"))
}

# Actually calculating the estimates for different grid sizes

ms <- c(10, 30, 50, 70)

for(m in ms){
  LOOCV_and_ECRPS_for_grid(m)
}