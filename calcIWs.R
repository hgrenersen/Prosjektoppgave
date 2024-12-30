source("/home/shomed/h/henriag/R-scripts/Markov-folder/calcAllFuncs.R")
source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/9) KLIEP.R")
library(CVXR)
library(progress)
#install.packages("clarabel")
#print(CVXR::installed_solvers())
beta = 2
print(paste0("Beta er ", beta))
folder_prefix = paste0("Data/Beta=", beta, "/")

num_of_estimates = 500
num_of_data_pts = 2000

datasets = readRDS(paste0(folder_prefix, "Data/Datasets.rds"))

bandwidths = readRDS(paste0(folder_prefix, "IWs/Bandwidths.rds"))
### KLIEP

pb <- progress_bar$new(
        format = paste0("Calculating KLIEP, for beta=", beta, " [:bar] :percent eta: :eta"),
        total = num_of_estimates, clear = FALSE, width= 60)

KLIEP_IWs <- calc_all_IWs(datasets, "KLIEP", pb=pb, bandwidths=bandwidths)
saveRDS(KLIEP_IWs, paste0(folder_prefix, "IWs/KLIEP IWs.rds"))
remove(KLIEP_IWs)