create_grid <- function(num_of_points){
  # Function to create a grid on the unit square with 
  # num_of_points grid points
  xseq <- seq(0, 1, , num_of_points + 1)[-(num_of_points + 1)]
  dx <- diff(xseq)[1]
  xseq <- xseq + 1 / 2 * dx
  
  yseq <- seq(0, 1, , num_of_points + 1)[-(num_of_points + 1)]
  dy <- diff(yseq)[1]
  yseq <- yseq + 1 / 2 * dy
  
  grid <- list(x = xseq, y = yseq) #Yes, I see now that
  # this could just've been grid <- list(x = xseq, y = xseq)
  # and dropped yseq
  
  exp_grid <- expand.grid(grid)
  return(exp_grid)
}

save_underways <- function(i, prefix, foldername, preds, pred_var, foldIds, originalIds,
                            save_per = 50, num_of_estimates = 500){
  # Function used to store results underways from a CV run over all datasets,
  # in case something happens
  # i: starting index of the predictions
  # prefix: folder to store results in, with subfolders for CV methods
  # foldername: name of the CV method used
  # preds: predictions to store
  # pred_var: prediction variances to store
  # foldIds: ids to store so that it can be identified which fold/block
  # a point belongs to
  # originalIds: datapoints original index before shuffling in a CV method

  #Storing the predictions
  filename <- paste0(prefix, foldername,"/", "Preds ", (i - 1) * save_per + 1,
  "_", i * save_per," of ", num_of_estimates,  ".rds", sep="")
  saveRDS(preds[,((i - 1) * save_per):(i * save_per)], file = filename)
  
  #Storing the prediction variances
  filename <- paste0(prefix, foldername,"/", "Pred vars ", (i - 1) * save_per + 1,
  "_", i * save_per," of ", num_of_estimates,  ".rds", sep="")
  saveRDS(pred_var[,((i - 1)*save_per):(i * save_per)], file = filename)
  print("Saved succesfully")
  
  if(foldername != "LOOCV"){ # Storing foldIds for other methods than LOOCV
    filename <- paste0(prefix, foldername, "/", "FoldIds ", (i - 1) * save_per + 1,
    "_", i * save_per," of ", num_of_estimates,  ".rds", sep="")
    saveRDS(foldIds[,((i - 1) * save_per):(i * save_per)], file = filename)
  }
}

final_save <- function(prefix, foldername, preds, pred_var, foldIds, originalIds){
  # Code to store the results of a completed CV run over all datasets.
  # prefix: folder to store results in, with subfolders for CV methods
  # foldername: name of the CV method used
  # preds: predictions to store
  # pred_var: prediction variances to store
  # foldIds: ids to store so that it can be identified which fold/block
  # a point belongs to
  # originalIds: datapoints original index before shuffling in a CV method

  filename <- paste0(prefix, foldername, "/", "Final Preds.rds", sep = "")
  saveRDS(preds, file = filename)

  filename <- paste0(prefix, foldername, "/", "Final Pred vars.rds", sep = "")
  saveRDS(pred_var, file = filename)
  
  if (foldername != "LOOCV") { # Again storing the foldids for methods that aren't LOOCV
    filename <- paste0(prefix, foldername, "/", "Final Fold indices.rds", sep = "")
    
    saveRDS(foldIds, filename)
  }  
}

library(progress)

do_CV <- function(datasets, method = "LOOCV", num_of_estimates = 500, 
                  save_per = 50, save = TRUE, save_underways = FALSE,
                  prefix = "Representative", Pb = FALSE, ...){
  # Function to actually calculate CV estimates over many datasets.
  # datasets: list of datasets, each element is a dataframe with
  # columns for x and y coordinates + response
  # method: CV method to be used, either "LOOCV", "K" or "BLOCK"
  # num_of_estimates: number of datasets CV estimates should be calculated for
  # prefix: Folder to store results to, depending on beta parameter in LGCP
  # Pb : Progress bar, the progress bar object must be passed if this should be used
  # Possible extra arguments: 
  # nFolds = Number of folds to use in K-fold
  # numBlocks = Number of folds to use in Block CV
  
  #Set-up of matrices for storing results and for storing files
  preds <- matrix(rep(0, nrow(datasets[[1]]) * num_of_estimates), ncol = num_of_estimates,
                  nrow = nrow(datasets[[1]]))
  pred_var <- preds
  
  foldIds <- preds
  
  originalIds <-  preds
  
  foldername = "LOOCV"
  
  if (method == "K") {
    foldername = paste0(list(...)[["nFolds"]], "-fold")
  }
  if (method == "BLOCK") {
    foldername = paste0(list(...)[["numBlocks"]], "-block")
  }
  
  num_cores <- detectCores() - 1
  
  for (i in 1:(num_of_estimates/save_per)) {
    for (j in 1:save_per) {
      k <- (i - 1) * save_per + j
      res <- 0
      
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      
      if(method == "LOOCV"){
        
        res <- LOOCV_downdate(datasets[[k]])
      }
      if(method=="K"){
        args <- list(...)
        res <- K_fold_CV_parallell(datasets[[k]], args[["nFolds"]])
        if(nrow(datasets[[1]]) %% args[["nFolds"]] != 0){ #Storing the foldids
        #if the number of datapoints isn't divisible by K, otherwise the
        # foldid is not needed
          foldIds[, k] <- res$foldId
        }
        
      }
      if(method == "BLOCK"){
        args <- list(...)
        res <- block_CV_parallell(datasets[[k]], args[["numBlocks"]])
        foldIds[, k] <- res$block
        
      }
      
      preds[, k] <- res$pred
      
      pred_var[, k] <- res$pred_var
      
      Pb$tick()
      
      stopCluster(cl)
    }
    if (save_underways) {
      save_underways(i, prefix, foldername, preds, pred_var, foldIds,
                     save_per = save_per,
                     num_of_estimates = num_of_estimates)
    }
  }
  if (save) {
    final_save(prefix, foldername, preds, pred_var, foldIds)
  }
  return(list(preds = preds, pred_var = pred_var,
              foldIds = foldIds))
}

loadData <- function(method = "LOOCV", nFolds = 0, prefix = "Representative"){
  # Function to load data from files, provided all results are stored after
  # a complete CV run has been completed
  prefix <- paste0("..\\", prefix, "\\")
  if (method == "LOOCV" ){
    preds <- readRDS(paste0(prefix, "LOOCV\\Final Preds.rds"))
    pred_var <- readRDS(paste0(prefix, "LOOCV\\Final Pred vars.rds"))
    return(list(preds=preds, pred_var=pred_var))
  }
  if(method=="K"){
    prefix <- paste0(prefix, nFolds, "-fold\\")
    
    preds <- readRDS(paste0(prefix, "Final Preds.rds"))
    pred_var <- readRDS(paste0(prefix, "Final Pred vars.rds"))
    foldIds <- readRDS(paste0(prefix, "Final Fold indices.rds"))
    
    return(list(preds=preds, pred_var=pred_var, foldIds=foldIds))
  }
  
  if(method=="BLOCK"){
    prefix <- paste0(prefix, nFolds, "-block\\")
    
    preds <- readRDS(paste0(prefix, "Final Preds.rds"))
    pred_var <- readRDS(paste0(prefix, "Final Pred vars.rds"))
    foldIds <- readRDS(paste0(prefix, "Final Fold indices.rds"))
    
    return(list(preds=preds, pred_var=pred_var, foldIds=foldIds))
    
  }
  
}

# The functions below aren't used!
do_IW <- function(datasets, method="KMM", num_of_estimates=500, save_per=50,
                  Exp_grid = exp_grid, 
                  save=TRUE, save_prefix="..\\Estimates2\\",
                  Sigma_opts=sigma_opts, Pb=FALSE,...){
  # Possible extra arguments: 
  # nFolds = Number of folds to use in K-fold
  # numBlocks = Number of folds to use in Block CV
  
  IWs <- matrix(rep(0, nrow(datasets[[1]])*num_of_estimates), ncol=num_of_estimates,
                nrow=nrow(datasets[[1]]))
  
  foldername <- method
  
  for(i in 1:(num_of_estimates/save_per)){
    for(j in 1:save_per){
      k <- (i-1)*save_per + j
      #print(k/num_of_estimates*100)
      weights <- 0
      
      train <- datasets[[k]][,1:2]
      
      sigma <- Sigma_opts[k]
      
      if(method=="KMM"){
        weights = KMM_weights(train,
                              Exp_grid,
                              sigma)$weights
      }
      if(method=="KLIEP"){
        alphas = KLIEP_alphas(train,
                              Exp_grid,
                              sigma)$alphas
        weights = KLIEP_weights(alphas, 
                                train,
                                Exp_grid,
                                sigma)
        
      }
      IWs[, k]<-weights
      Pb$tick()
      
    }
    if(save){
      filename <- paste0(save_prefix, foldername,"\\", "IWs ", (i-1)*save_per,"_", (i*save_per-1)," of ", num_of_estimates,  ".rds", sep="")
      saveRDS(IWs[,((i-1)*save_per):(i*save_per-1)], file=filename)
      print("Saved succesfully")  
    }
  }
  if(save){
    filename <- paste0(save_prefix, foldername,"\\", "Final IWs.rds", sep="")
    saveRDS(IWs, file=filename) 
  }
  return(IWs)
}

gather_IWs_from_files <- function(method="KMM", num_of_estimates = 500, save_per = 50, 
                                  save_prefix="..\\Estimates2\\") {
  i<-1
  filename <- paste0(save_prefix, method,"\\", "IWs ", (i-1)*save_per,"_", (i*save_per-1)," of ", num_of_estimates,  ".rds", sep="")
  loaded_IWs <-readRDS(filename)
  print(dim(loaded_IWs))
  for(i in 2:(num_of_estimates/save_per)){
    filename <- paste0(save_prefix, method,"\\", "IWs ", (i-1)*save_per,"_", (i*save_per-1)," of ", num_of_estimates,  ".rds", sep="")
    loaded_IWs<-cbind(loaded_IWs, readRDS(filename))
    
    print(all.equal(loaded_IWs[,i*save_per],
                    readRDS(filename)[, 1]))
  }
  return(loaded_IWs)
}

