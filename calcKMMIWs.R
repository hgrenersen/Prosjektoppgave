source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/7.5) KMM_checking.R")
source("/home/shomed/h/henriag/R-scripts/Markov-folder/calcAllFuncs.R")
#install.packages("doBy")

library(progress)
setwd("/home/shomed/h/henriag/R-scripts/Markov-folder")
#Needs library(reticulate) too some places. Then, these should be run:
#py_install("numpy")
#py_install("cvxopt")

beta <- 2

print(paste0("Beta er ", beta))

folder_prefix = paste0("Data/Beta=", beta, "/")

datasets <- readRDS(paste0(folder_prefix, "/Data/Datasets.rds"))

num_of_datasets <- 500

test_locs <- as.matrix(create_grid(20))

num_of_nearest_neighbors <- 10

optimal_bandwidths <- rep(0, num_of_datasets)

KMM_IWs <- matrix(rep(0, 2000*num_of_datasets), nrow=2000, ncol=num_of_datasets)

pb <- progress_bar$new(
    format = "  Estimating optimal bandwidth and KMM IW [:bar] :percent eta: :eta",
    total = num_of_datasets, clear = FALSE, width= 60)

for(i in 1:num_of_datasets){
  data <- datasets[[i]]

  data_locs <- as.matrix(data[, 1:2])
  #Tuning
  bandwidths <- median(rdist(data_locs)) * c(seq(0.1, 2, by = 0.1))

  tune_res <- bandwidth_tuning(data_locs, test_locs, bandwidths)
  
  opt_ind <- which.min(tune_res$J_vals)

  optimal_bandwidths[i] <- bandwidths[opt_ind]

  KMM_IWs[, i] <- tune_res$KMM_weights_mat[, opt_ind]


  # KMM_IWs using median distance

  #dist_mat <- rdist(data_locs)

  #optimal_bandwidths[i] <- median(dist_mat)

  #KMM_IWs[, i] <- pythonized_KMM_weights(data_locs, test_locs, optimal_bandwidths[i])$weights

  pb$tick()
}
#saveRDS(KMM_IWs, paste0(folder_prefix, "IWs/Median KMM IWs.rds"))
saveRDS(KMM_IWs, paste0(folder_prefix,"IWs/KMM IWs.rds"))
saveRDS(optimal_bandwidths, paste0(folder_prefix, "IWs/KMM Bandwidths.rds"))