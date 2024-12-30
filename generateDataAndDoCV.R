source("/home/shomed/h/henriag/R-scripts/Markov-folder/calcAllFuncs.R")

beta = 0
print(paste0("Beta er ", beta))
folder_prefix <- paste0("Data/Beta=", beta, "/")

num_of_estimates <- 500

realizations <- readRDS("Data/realizations.rds")

generate_data <- FALSE

# Generating data
if(generate_data){
    set.seed(10)

    pb <- progress_bar$new(
        format = "Generate data [:bar] :percent eta: :eta",
        total = num_of_estimates, clear = FALSE, width= 60)

    simulated_data <- generate_data_with_intensity(realizations, beta, pb=pb)

    intensities <- simulated_data[[1]]
    datasets <- simulated_data[[2]]

    saveRDS(intensities, paste0(folder_prefix, "Data/Intensities.rds"))
    saveRDS(datasets, paste0(folder_prefix, "Data/Datasets.rds"))
    print("Data sucessfully GENERATED")
} else{
    datasets <- readRDS(paste0(folder_prefix, "Data/Datasets.rds"))
    intensities <- readRDS(paste0(folder_prefix, "Data/Intensities.rds"))

    if(nrow(datasets[[1]])==2000 & ncol(datasets[[1]])==3){
        print("Data loaded succesfully")
    } else{
        print("yuck")
    }
}


# CV

num_cores <- 6 #detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Block

## 9-block

set.seed(1)

pb <- progress_bar$new(
    format = "  9-block [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(datasets, method="BLOCK", numBlocks = 9,  Pb=pb,
prefix=paste0(folder_prefix, "CV/")) 


## 25-block

set.seed(1)

pb <- progress_bar$new(
    format = "  25-block  [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)
 
do_CV(datasets, method="BLOCK", numBlocks = 25,  Pb=pb, prefix=paste0(folder_prefix, "CV/")) 


# K-fold

## 9-fold

set.seed(1)

pb <- progress_bar$new(
    format = "  9-fold [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)


do_CV(datasets, nFolds = 9,  method="K",
                   nFolds = 9, Pb=pb, prefix=paste0(folder_prefix, "CV/")) 


## 25-fold

set.seed(1)

pb <- progress_bar$new(
    format = "  25-fold [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(datasets,  nFolds = 25,  method="K",
                   nFolds = 9, Pb=pb, prefix=paste0(folder_prefix, "CV/")) 


# LOOCV

set.seed(1)

pb <- progress_bar$new(
    format = "  LOOCV [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(datasets, Pb=pb, prefix=paste0(folder_prefix, "CV/"))
