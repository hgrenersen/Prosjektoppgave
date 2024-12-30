library(pracma)
library(doParallel)
library(foreach)

source("/home/shomed/h/henriag/R-scripts/Markov-folder/calcAllFuncs.R")

#setwd("/home/shomed/h/henriag/R-scripts/Markov-folder/Data")
beta = 0
print(paste0("Beta er ", beta))
folder_prefix = paste0("Data/Beta=", beta, "/")

num_of_estimates = 500
num_of_data_pts = 2000

#print(getwd())

realizations = readRDS("Data/realizations.rds")

datasets = readRDS(paste0(folder_prefix, "Data/Datasets.rds"))

intensities = readRDS(paste0(folder_prefix, "Data/Intensities.rds"))

print("Bandwidth optimization starting")

pb <- progress_bar$new(
    format = "Bandwidths [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

bandwidths <- rep(0, num_of_estimates)

calc_bandwidths <- FALSE

if(calc_bandwidths){
    for(i in 1:num_of_estimates){
        data <- as.matrix(datasets[[i]][, 1:2])
        dist_mat <- dist(data)
        optim_obj <- optim(par=0.1, fn=bandwidth_objective, Dist_mat = dist_mat, Data=data,
        correction="None")

        bandwidths[i] <- optim_obj$par
        pb$tick()
    }
    saveRDS(bandwidths, paste0(folder_prefix, "IWs/Bandwidths.rds"))
} else{
    bandwidhts <- readRDS(paste0(folder_prefix, "IWs/Bandwidths.rds"))
    print("Bandwidths loaded")
}

print("Bandwidth optimization done")

### KMM

#KMM_IWs <- calc_all_IWs(datasets, "KMM", bandwidths)
#saveRDS(KMM_IWs, paste0(folder_prefix, "\\IWs\\KMM IWs.rds"))
#remove(KMM_IWs)

### KLIEP

#KLIEP_IWs <- calc_all_IWs(datasets, "KLIEP", bandwidths)
#saveRDS(KLIEP_IWs, paste0(folder_prefix, "\\IWs\\KLIEP IWs.rds"))
#remove(KLIEP_IWs)

## Voronoi

pb <- progress_bar$new(
    format = "Voronoi [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

Vor_IWs <- calc_all_IWs(datasets, "Voronoi", pb)
saveRDS(Vor_IWs, paste0(folder_prefix, "IWs/Voronoi IWs.rds"))
remove(Vor_IWs)

# True error
pb <- progress_bar$new(
    format = "True error [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)
exp_grid <- create_grid(200)

num_cores <- 6 #detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

preds_over_grid_mat <- matrix(rep(0, 200^2*num_of_estimates),
ncol=num_of_estimates)

pred_vars_over_grid_mat <- matrix(rep(0, 200^2*num_of_estimates),
ncol=num_of_estimates)

#preds_over_grid_mat <- readRDS(paste0(folder_prefix, "True/Predictions.rds"))

#pred_vars_over_grid_mat <- readRDS(paste0(folder_prefix, "True/Prediction variances.rds"))

ECRPS <- rep(0, num_of_estimates)

nugget <- 0.1

for(i in 1:num_of_estimates){
    pred_over_grid_list <- pred_grid(datasets[[i]], exp_grid)
    preds_over_grid_mat[, i] <- pred_over_grid_list[[1]]
    pred_vars_over_grid_mat[, i] <- pred_over_grid_list[[2]]

    ECRPS[i] <- expectedCRPS(c(realizations[[i]]), nugget,
    preds_over_grid_mat[, i], pred_vars_over_grid_mat[, i])
    pb$tick()
}

saveRDS(preds_over_grid_mat, paste0(folder_prefix, "True/Predictions.rds"))
saveRDS(pred_vars_over_grid_mat, paste0(folder_prefix, "True/Prediction variances.rds"))
saveRDS(ECRPS, paste0(folder_prefix, "True/ECRPS.rds"))

# True importance weights
library(rdist)
pb <- progress_bar$new(
    format = "True weights [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

true_IWs <- matrix(rep(0, num_of_data_pts*num_of_estimates), ncol=num_of_estimates)

for(i in 1:num_of_estimates){
  data_locs <- datasets[[i]][, 1:2]
  
  cross_dist <- cdist(data_locs, exp_grid)
  
  min_dist_ind <- apply(cross_dist, MARGIN=1, which.min)
  
  true_IWs[, i] <- 1/(c(intensities[[i]])[min_dist_ind])
  
  pb$tick()
}
saveRDS(true_IWs, paste0(folder_prefix, "True/Importance weights.rds"))