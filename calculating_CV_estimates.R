library(doParallel)
library(progress)


source("R-scripts/Markov-folder/Scripts/2) K-fold CV.R")
source("R-scripts/Markov-folder/Scripts/3) Proper K-fold CV.R")
source("R-scripts/Markov-folder/Scripts/4) Large-scale K-fold CV.R")
source("R-scripts/Markov-folder/Scripts/5) Block CV.R")
source("R-scripts/Markov-folder/Scripts/11) CV - predictions and prediction variance.R") 


rep_datasets <-  readRDS("R-scripts/Markov-folder/Data/Noisy representative data.rds")

non_rep_datasets <- readRDS("R-scripts/Markov-folder/Data/Noisy non-representative data.rds")

setwd("R-scripts/Markov-folder")

num_of_estimates <- 500

exp_grid <- create_grid(200)

# LOOCV

num_cores <- 6 #detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

set.seed(1)

pb <- progress_bar$new(
    format = "  LOOCV representative [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(rep_datasets, Pb=pb, prefix="Rep")


## Non-representative

set.seed(1)

pb <- progress_bar$new(
    format = "  LOOCV non-representative [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(non_rep_datasets, prefix = "Non-rep", Pb=pb) 

stopCluster(cl)
# K-fold

## 9-fold

## Representative
#set.seed(1)

#pb <- progress_bar$new(
#    format = "  9-fold, representative [:bar] :percent eta: :eta",
#    total = num_of_estimates, clear = FALSE, width= 60)


#do_CV(rep_datasets, nFolds = 9,  method="K",
#                   nFolds = 9, Pb=pb, prefix="Rep") 

### Non-representative

set.seed(1)

pb <- progress_bar$new(
    format = "  9-fold [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)


do_CV(non_rep_datasets, nFolds = 9,  method="K",
                   nFolds = 9, Pb=pb, prefix = "Non-rep") 


## 25-fold

### Representative

set.seed(1)

pb <- progress_bar$new(
    format = "  25-fold representative [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(rep_datasets,  nFolds = 25,  method="K",
                   nFolds = 9, Pb=pb, prefix="Rep") 


### Non-representative


set.seed(1)

pb <- progress_bar$new(
    format = "  25-fold non-representative [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(non_rep_datasets,  nFolds = 25,  method="K",
                   nFolds = 9, Pb=pb, prefix = "Non-rep") 

# Block

## 9-block

### Representative

set.seed(1)

pb <- progress_bar$new(
    format = "  9-block [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(rep_datasets, method="BLOCK", numBlocks = 9,  Pb=pb,
prefix="Rep") 


### Non-representative


set.seed(1)

pb <- progress_bar$new(
    format = "  9-block non-representative [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(non_rep_datasets, method="BLOCK", numBlocks = 9,  Pb=pb,
                   prefix = "Non-rep") 


## 25-block

### Representative


set.seed(1)

pb <- progress_bar$new(
    format = "  25-block representative [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)
 
do_CV(rep_datasets, method="BLOCK", numBlocks = 25,  Pb=pb, prefix="Rep") 

### Non-representative

set.seed(1)

pb <- progress_bar$new(
    format = "  25-block non-representative [:bar] :percent eta: :eta",
    total = num_of_estimates, clear = FALSE, width= 60)

do_CV(non_rep_datasets, method="BLOCK", numBlocks = 25,  Pb=pb,
                    prefix = "Non-rep") 
