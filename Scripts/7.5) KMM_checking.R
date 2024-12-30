GRBF<-function(distance, sigma){
  #Function that evaluates the gaussian kernel given a distance
  # and some bandwidth. Not normalized!
  exp(- 1 / 2 * (1 / sigma * distance)^2)
}

KMM_weights_opti <- function(train, test, sigma, posdef = 0){
  # Not used, as the library CVXR had better performance, see next definition
  # Function that calculates the IWs for training observations according to the 
  # KMM method using the library(optiSolve). Will not be used except for in
  # comparisons.
  #   -train: Training dataset, very important that this does not include a 
  # response column 
  #   -test: Test dataset, also only two columns 
  #   -sigma: Bandwidth to be used in the kernel 
  #   -posdef: Regularization parameter that can be used to make the problem
  # strongly convex, meaning that the estimated weights are unique
  # 
  # returns: A list, containing
  #   -The estimated IWs
  #   -The value of the objective function at the optima
  #   -The Gram matrix
  #   -The kappa vector
  
  if (ncol(train) != 2 | ncol(test) != 2) {
    return("You've probably included the response")
  }
  
  n <- nrow(train)
  m <- nrow(test)
  
  dist_mat <- dist(train)
  
  K <- as.matrix(GRBF(dist_mat, sigma))
  diag(K) <- 1+posdef
  
  kappa <- n / m * as.matrix(GRBF(cdist(train, test), sigma)) %*% rep(1, m)
  
  objective <- quadfun(1 / 2 * K, a = - kappa)
  
  lower_constraints <- lbcon(0, 1:n)
  
  B<-1000
  
  upper_constraints <- ubcon(B, 1:n)
  
  epsilon <- 1 - 1 / sqrt(n)
  
  constraint_matrix <- 1 / n * matrix(data = rep(1, 2 * n), nrow = 2)
  
  eps_constraints <- lincon(constraint_matrix,
   val = c(epsilon + 1, 1 - epsilon), dir = c("<=", ">="), name = c("upper", "lower"))
  
  prob <- cop(objective,
              lb = lower_constraints,
              ub = upper_constraints,
              lc = eps_constraints)
  solution <- solvecop(prob)
  return(list(weights = solution$x,
              objective = t(solution$x) %*% K %*% solution$x - t(kappa) %*% solution$x,
              K = K,
              kappa = kappa))
}

KMM_weights <- function(train, test, sigma, epsilon = 1 - 1 / sqrt(nrow(train)),
                                   solver = "CLARABEL", posdef = 0){
  # Function that calculates the IWs for training observations according to the 
  # KMM method using the library(CVXR). 
  #   -train: Training dataset, very important that this does not include a 
  # response column 
  #   -test: Test dataset, also only two columns 
  #   -sigma: Bandwidth to be used in the kernel 
  #   -epsilon: Parameter for the constraint on the mean
  #   -solver: Name of solver to be used, CLARABEL has worked best 
  # with respect to speed and similarity to the Python implementation (adapt)
  #   -posdef: Regularization parameter that can be used to make the problem
  # strongly convex, meaning that the estimated weights are unique
  # 
  # returns: A list, containing
  #   -The estimated IWs
  #   -Result of the optimization procedure
  #   -The Gram matrix
  #   -The kappa vector
  
  n <-nrow(train)
  m <- nrow(test)
  
  dist_mat <- dist(train)
  
  K <- as.matrix(GRBF(dist_mat, sigma))
  diag(K) <- 1 + posdef
  
  kappa <- apply(n / m * as.matrix(GRBF(cdist(train, test), sigma)),
                 MARGIN = 1, FUN = sum)
  
  B<-1000
  
  w <- Variable(n, pos = TRUE)
  
  obj <- Minimize(1 / 2 * quad_form(w, K) - t(kappa) %*% w)
  
  constraint <- list(1 / n * t(w) %*% rep(1, n) <= epsilon + 1,
                     1 / n * t(w) %*% rep(1, n) >= 1 - epsilon,
                     w <= B)
  
  prob <- Problem(obj, constraint)
  
  result <- solve(prob,  reltol = 1e-6,
                  abstol = 1e-7, warm_start = TRUE,
                  feastol = 1e-7,
                  solver = solver, 
                  num_iter = 100
                  )
  
  return(list(weights = result$getValue(w), result = result,
         K = K, kappa = kappa))
}

pythonized_KMM_weights <- function(train, test, sigma, epsilon = 1 - 1 / sqrt(nrow(train)),
                        solver = "CLARABEL", posdef = 0){
  # Function that calculates the IWs for training observations according to the 
  # KMM method using the library(CVXR), but where the actual optimization is
  # done in Python, yielding faster results than the above function. 
  #   -train: Training dataset, very important that this does not include a 
  # response column 
  #   -test: Test dataset, also only two columns 
  #   -sigma: Bandwidth to be used in the kernel 
  #   -epsilon: Parameter for the constraint on the mean
  #   -solver: Name of solver to be used, CLARABEL has worked best 
  # with respect to speed and similarity to the Python implementation (adapt)
  #   -posdef: Regularization parameter that can be used to make the problem
  # strongly convex, meaning that the estimated weights are unique
  # 
  # returns: A list, containing
  #   -The estimated IWs
  #   -Result of the optimization procedure
  #   -The Gram matrix
  #   -The kappa vector
  
  n <- nrow(train)
  m <- nrow(test)
  
  dist_mat <- dist(train)
  
  K <- as.matrix(GRBF(dist_mat, sigma))
  diag(K) <- 1 + posdef
  
  kappa <- apply(n / m * as.matrix(GRBF(cdist(train, test), sigma)),
                 MARGIN = 1, FUN = sum)
  reticulate::source_python('/home/shomed/h/henriag/R-scripts/Markov-folder/Scripts/Pythonized_KMM.py')
  weights <- py_KMM_optimization(K, kappa)
  
  return(list(weights = weights,
              K = K, kappa = kappa))
}

KMM_objective<-function(weights, K, kappa){
  #Function to evaluate the objective function given IWs,
  #the Gram matrix and the kappa vector
  1 / 2 * t(weights) %*% K %*% weights - t(kappa) %*% weights
}

KMM_constraint <- function(weights, epsilon = 1 - 1 / sqrt(length(weights))) {
  #Checking if the constraint is satisfied for the estimated weights
  abs(mean(weights) - 1) <= epsilon
}