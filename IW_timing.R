library(profvis)
#install.packages("glmnet")

beta <- 0
setwd("/home/shomed/h/henriag/R-scripts/Markov-folder")
folder_prefix <- paste0("Data/Beta=", beta, "/IWs/")

source("calcAllFuncs.R")

KLIEP_fit <- function(dataset, bwd, test_locs = create_grid(20)){
     alphas <- KLIEP_alphas(dataset[, 1:2],
            test_locs, sigma=bwd)$alphas

    KLIEP_weights(alphas, dataset[, 1:2],
    test_locs, sigma=bwd)
}

vor_fit <- function(dataset){
  # Assumes that dataset is a nx2 matrix
    dir <- deldir(dataset[, 1], dataset[, 2], rw=c(0,1,0,1))

    tiles <- tile.list(dir)
    areas <- rep(1, nrow(dataset))
    for(j in 1:nrow(dataset)){
        areas[j] <- tiles[[j]]$area
    }
    areas # Not normalized
}

save_immediate_results <- function(KMM_tuned_IWs,
KMM_median_IWs, KLIEP_IWs, Voronoi_IWs, folder_prefix=paste0("Data/Beta=", beta, "/IWs/")){
  
  print("Saving underways")
  total_KMM_tuned_IWs <- readRDS(paste0(folder_prefix, "Tuned KMM IWs.rds"))
  saveRDS(cbind(total_KMM_tuned_IWs,KMM_tuned_IWs), paste0(folder_prefix, "Tuned KMM IWs.rds"))
  remove(total_KMM_tuned_IWs)
  gc()

  total_KMM_median_IWs <- readRDS(paste0(folder_prefix, "Median KMM IWs.rds"))
  saveRDS(cbind(total_KMM_median_IWs,KMM_median_IWs), paste0(folder_prefix, "Median KMM IWs.rds"))
  remove(total_KMM_median_IWs)
  gc()

  total_KLIEP_IWs <- readRDS(paste0(folder_prefix, "KLIEP IWs.rds"))
  saveRDS(cbind(total_KLIEP_IWs, KLIEP_IWs), paste0(folder_prefix, "KLIEP IWs.rds"))
  remove(total_KLIEP_IWs)
  gc()


  total_Voronoi_IWs <- readRDS(paste0(folder_prefix, "Voronoi IWs.rds"))
  saveRDS(cbind(total_Voronoi_IWs, Voronoi_IWs), paste0(folder_prefix, "Voronoi IWs.rds"))
  remove(total_Voronoi_IWs)
  gc()
}

print(paste0("Beta er", beta))

datasets <- readRDS(paste0("/home/shomed/h/henriag/R-scripts/Markov-folder/Data/Beta=", beta,"/Data/Datasets.rds"))

test_locs <- as.matrix(create_grid(20))

num_of_timings <- 500

save_per <- 50

timing_results <- matrix(rep(0, num_of_timings * 4), ncol = 4)

KMM_tune_bandwidths <- rep(0, num_of_timings)
KMM_median_bandwidths <- rep(0, num_of_timings)
KLIEP_bandwidths <- rep(0, num_of_timings)

for(i in 1:(num_of_timings/save_per)){

  KMM_tuned_IWs <- matrix(rep(0, 2000 * save_per), ncol = save_per)

  KMM_median_IWs <- matrix(rep(0, 2000 * save_per), ncol = save_per)

  KLIEP_IWs <- matrix(rep(0, 2000 * save_per), ncol = save_per)

  Voronoi_IWs <- matrix(rep(0, 2000 * save_per), ncol = save_per)

  for(j in 1:save_per){
    
    global_ind <- (i-1) * save_per + j
    print(global_ind)
    data <- as.matrix(datasets[[global_ind]][, 1:2])

    data_locs <- data

    start_tuned_KMM <- Sys.time()

    dist_mat <- as.matrix(rdist(data))

    median_dist <- median(dist_mat[upper.tri(dist_mat)])

    # KMM, tuned bandwidth
    bandwidths <- median_dist * c(seq(0.1, 2, by = 0.1))
    tune_res <- bandwidth_tuning(data_locs, test_locs, bandwidths)

    opt_ind <- which.min(tune_res$J_vals)

    KMM_tune_bandwidths[global_ind] <- bandwidths[opt_ind]

    KMM_tuned_IWs[, j] <- tune_res$KMM_weights_mat[, opt_ind]

    end_tuned_KMM <- Sys.time()

    timing_results[global_ind, 1] <- end_tuned_KMM[3] - start_tuned_KMM[3]

    # Just the fitting of KMM
    start_med_KMM <- Sys.time()
    dist_mat <- as.matrix(rdist(data))

    median_dist <- median(dist_mat[upper.tri(dist_mat)])
    KMM_median_bandwidths[global_ind] <- median_dist
    KMM_median_IWs[, j] <- pythonized_KMM_weights(data_locs, test_locs, median_dist)$weights
    end_med_KMM <- Sys.time()

    timing_results[global_ind, 2] <- end_med_KMM - start_med_KMM
    # KLIEP
    # Selecting the bandwidth
    
    start_KLIEP <- Sys.time()

    optim_obj <- optim(par = 0.1, fn = bandwidth_objective, 
                        Dist_mat = dist_mat, Data = data,
    correction="None")

    KLIEP_bandwidths[global_ind] <- optim_obj$par

    # Fitting the weights
    KLIEP_IWs[, j] <- KLIEP_fit(data, bwd= optim_obj$par, test_locs = test_locs)
    end_KLIEP <- Sys.time()

    timing_results[global_ind, 3] <- end_KLIEP[3] - start_KLIEP[3]

    # Voronoi
    start_Vor <- Sys.time()

    Voronoi_IWs[, j] <- vor_fit(data)

    end_Vor <- Sys.time()

    timing_results[global_ind, 4] <- end_Vor[3] - start_KLIEP[3]
  }
  if(i == 1){
    saveRDS(KMM_tuned_IWs, paste0(folder_prefix, "Tuned KMM IWs.rds"))
    
    saveRDS(KMM_median_IWs, paste0(folder_prefix, "Median KMM IWs.rds"))
    
    saveRDS(KLIEP_IWs, paste0(folder_prefix, "KLIEP IWs.rds"))
    
    saveRDS(Voronoi_IWs, paste0(folder_prefix, "Voronoi IWs.rds"))
  }
  else{
    save_immediate_results(KMM_tuned_IWs,
KMM_median_IWs, KLIEP_IWs, Voronoi_IWs)
  }
}


saveRDS(KMM_tune_bandwidths, paste0(folder_prefix, "KMM tuned bandwidths.rds"))
saveRDS(KMM_median_bandwidths, paste0(folder_prefix, "KMM median bandwidths.rds"))
saveRDS(KLIEP_bandwidths, paste0(folder_prefix, "KLIEP bandwidths.rds"))

saveRDS(timing_results, paste0(folder_prefix, "Timing results.rds"))