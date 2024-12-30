library(doParallel)
library(progress)
library(deldir)

source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/1) Simulating data.R")

# Functions to generate data

create_grid <- function(num_of_points){
    xseq <- seq( 0,1,,num_of_points+1)[-(num_of_points+1)]
    dx <- diff(xseq)[1]
    xseq<-xseq + 1/2*dx

    yseq <- seq( 0,1,,num_of_points+1)[-(num_of_points+1)]
    dy <- diff(yseq)[1]
    yseq<-yseq + 1/2*dx

    grid<- list( x= xseq, y= yseq)

    exp_grid <- expand.grid(grid)
    return(exp_grid)
}
library(fields)
generate_data_with_intensity <- function(realizations, beta, num_of_estimates=500, finer_size = 200,
                                         dataset_size = 2000, pb="",
                                         nugget = 0.1){
  num_of_points <- 200
  
  xseq <- seq( 0,1,,num_of_points+1)[-(num_of_points+1)]
  dx <- diff(xseq)[1]
  xseq<-xseq + 1/2*dx
  
  yseq <- seq( 0,1,,num_of_points+1)[-(num_of_points+1)]
  dy <- diff(yseq)[1]
  yseq<-yseq + 1/2*dx
  
  grid<- list( x= xseq, y= yseq)
  obj<- circulantEmbeddingSetup(grid, Covariance="Matern", aRange=.05, smoothness=1)
  
  intensities <- list()
  updated_datasets <- list()
  
  for(i in 1:num_of_estimates){
    
    underlying_realization<- circulantEmbedding(obj)
    
    intensity <- exp(beta*c(underlying_realization))/sum(exp(beta*c(underlying_realization)))
    
    intensities[[i]] <- intensity
    
    updated_datasets[[i]] <-  generate_data(realizations[[i]],
                                            intensity,
                                            create_grid(finer_size), dataset_size)
    #Add noise to the response
    updated_datasets[[i]]$z <- updated_datasets[[i]]$z + rnorm(num_of_points, sd=sqrt(nugget))
    
    pb$tick()
  }
  return(list(intensities=intensities,
              updated_datasets=updated_datasets))
}

# Functions for CV

source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/2) K-fold CV.R")
source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/3) Proper K-fold CV.R")
source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/4) Large-scale K-fold CV.R")
source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/5) Block CV.R")
source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/9) KLIEP.R")
source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/11) CV - predictions and prediction variance.R") 

# Functions for IWs

source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/6) Automatic Bandwidth Selection.R")


calc_all_IWs <- function(datasets, type, pb, bandwidths = 0.01){
    num_of_datasets <- length(datasets)

    num_of_data_points <- nrow(datasets[[1]])

    IWs_mat <- matrix(rep(0, num_of_data_points * num_of_datasets),
    ncol=num_of_datasets)

    if(type=="KMM"){
        for(i in 1:num_of_datasets){
            IWs_mat[, i] <- KMM_weights(datasets[[i]][, 1:2],
            create_grid(200), sigma=bandwidths[i])
            pb$tick()
        }
    }
    if(type=="KLIEP"){
        for(i in 1:num_of_datasets){
            alphas <- KLIEP_alphas(datasets[[i]][, 1:2],
            create_grid(20), sigma=bandwidths[i])$alphas

            IWs_mat[, i] <- KLIEP_weights(alphas, datasets[[i]][, 1:2],
            create_grid(20), sigma=bandwidths[i])
            pb$tick()
        }
    }
    if(type=="Voronoi"){
        for(i in 1:num_of_datasets){
          dir <- deldir(datasets[[i]]$x, datasets[[i]]$y, rw=c(0,1,0,1))

          tiles <- tile.list(dir)
          areas <- rep(1, num_of_data_points)
          for(j in 1:num_of_data_points){
            areas[j] <- tiles[[j]]$area
          }
          IWs_mat[, i] <- areas
          pb$tick()
        }
        
        
    }
    return(IWs_mat)
}

#ECRPS 

expectedCRPS <- function(truth, truth.var, est, est.var, getAverage=TRUE, na.rm=FALSE, weights=rep(1, length(truth))) {
  erf<-function(x) {2*(pnorm(sqrt(2)*x) - 0.5)}
  
  SSS<-truth.var + est.var # sigStarSq
  b<-est-truth # bias
  
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

# Tuning bandwidth in KMM
library(rdist)
library(CVXR)
library(glmnet)
library(doBy)
source("/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/7.5) KMM_checking.R")

create_des_mat <- function(data_locs, IWs,test_locs=NULL, n = 5){
  if(ncol(data_locs)!=2){
      print("WHAT ARE YOU DOING")
      return("YUCK")
  }
  if(is.null(test_locs)){
    dist_mat <- rdist(data_locs)
    
    dist_mat <- as.matrix(dist_mat)
    diag(dist_mat) <- 10
    
    nobs <- nrow(data_locs)
    
  }
  else{
    dist_mat <- as.matrix(cdist(test_locs, data_locs))
    
    nobs <- nrow(test_locs)
  }
  
  min_indices <- apply(dist_mat, MARGIN=1, FUN=which.minn, n=n) # Column i gives
    # the indices of the points closest to point i, sorted decreasingly. I.e.
    # we have n rows and number of columns=number of data locations
  X <- matrix(IWs[t(min_indices)], nrow=nobs, ncol=n)
  X <- cbind(X, matrix(dist_mat[t(min_indices)], nrow=nobs, ncol=n))
  return(X)
}

bandwidth_tuning <- function(data_locs, test_locs, bandwidths,
                              num_of_nearest_neighbors=10, num_cv = 50, use_optim=FALSE){
  n <- nrow(data_locs)
  m <- length(bandwidths)
  J_vals <- rep(0, m)

  num_cv <- 10
  
  
  if(use_optim){
    opt <- optim(par=median(bandwidths), fn =tuning_obj , data_locs = data_locs, test_locs = test_locs)
  }
  else{
    KMM_weights_mat <- matrix(rep(0, n*m), nrow=n, ncol=m)

    for(i in 1:m){
    #print(i/m*100)
      bandwidth <- bandwidths[i]
      # Fitting the weights and storing them
      
      KMM_weights_mat[, i] <- pythonized_KMM_weights(data_locs, test_locs, bandwidth)$weights
      
      normalization <- sum(KMM_weights_mat[, i])
      
      KMM_weights_mat[, i] <- KMM_weights_mat[, i]/normalization
      
      # RLS

      data_des_mat <- create_des_mat(data_locs,  log10(KMM_weights_mat[, i]),
                    n=num_of_nearest_neighbors)
      
      lambdas <- rep(0, num_cv)
      errors <- rep(0, num_cv)
      set.seed(1)
      for(j in 1:num_cv){
        cv_model <- cv.glmnet(data_des_mat, log10(KMM_weights_mat[, i]), alpha=0,
        parallel = FALSE)  
        lambdas[j] <- cv_model$lambda.min
        errors[j] <-min(cv_model$cvm)
      }
    
      model <- glmnet(x = data_des_mat, y= log10(KMM_weights_mat[, i]),
                    lambda=lambdas[which.min(errors)],
                    alpha=0)
      
      preds_test <- predict(model, newx= create_des_mat(data_locs, log10(KMM_weights_mat[, i]), test_locs=test_locs,n=num_of_nearest_neighbors))
      
      preds_test <- exp(preds_test) #Invert the log scaling
      
      #print(ggplot(test_locs, aes(x=x, y=y))+geom_point(aes(color=preds_test))+coord_equal()+geom_point(data=data_locs, color="red"))

      # J score
      J_vals[i] <- mean(KMM_weights_mat[, i]^2) - 2*mean(preds_test)
      
      KMM_weights_mat[, i] <- KMM_weights_mat[, i]*normalization
      
    }
    
  return(list(J_vals=J_vals, KMM_weights_mat = KMM_weights_mat))
  }
}

tuning_obj <- function(bandwidth, data_locs, test_locs){
  num_of_nearest_neighbors <- 10
  num_cv <- 20
  
  KMM_weights <- pythonized_KMM_weights(data_locs, test_locs, bandwidth)$weights
  
  normalization <- sum(KMM_weights)
  
  KMM_weights <- KMM_weights/normalization
  
  # RLS

  data_des_mat <- create_des_mat(data_locs,  log10(KMM_weights),
                n=num_of_nearest_neighbors)
  
  lambdas <- rep(0, num_cv)
  errors <- rep(0, num_cv)
  set.seed(1)
  for(j in 1:num_cv){
    cv_model <- cv.glmnet(data_des_mat, as.matrix(log10(KMM_weights), ncol=1), alpha=0)
    lambdas[j] <- cv_model$lambda.min
    errors[j] <- min(cv_model$cvm)
  }


  model <- glmnet(x = data_des_mat, y= log10(KMM_weights),
                lambda=lambdas[which.min(errors)],
                alpha=0)
  
  preds_test <- predict(model, newx= create_des_mat(data_locs, log10(KMM_weights), test_locs=test_locs,n=num_of_nearest_neighbors))
  
  preds_test <- exp(preds_test) #Invert the log scaling

  # J score
  return(mean(KMM_weights^2) - 2 * mean(preds_test))
}
