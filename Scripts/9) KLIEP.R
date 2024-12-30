library(rdist)

GRBF <- function(distance, sigma){
  exp(-1 / 2 * (1 / sigma * distance)^2)
}

KLIEP_alphas <- function(train, test, sigma,
                         solver = "CLARABEL"){
  # Function to calculate the coefficients of the basis functions in KLIEP.
  # Per now this assumes that basis functions are centered at all test locations,
  # so please don't use a very fine grid/many test locations for your own good.
  
  #returns: A list containing
  # -alphas: The estimated coefficients of the basis functions
  # -result: The result object of solving the optimization problem.
  # -Phi_test: The matrix used in the optimization that contains basis functions
  # evaluated at the test locations.
  # -Phi_train: The matrix used in the optimization that contains basis functions
  # evaluated at the train locations.
  
  n <- nrow(train)
  m <- nrow(test)
  
  test_dists <- dist(test)
  
  Phi_test <- as.matrix(GRBF(test_dists, sigma))
  diag(Phi_test) <- 1
  
  alpha <- Variable(m, pos = TRUE)
  
  obj <- Maximize(
    sum(log(Phi_test %*% alpha))
  )
  
  cross_dists <- cdist(train, test)
  
  Phi_train <- as.matrix(GRBF(cross_dists, sigma))
  
  constraint <- list(sum(Phi_train %*% alpha) == n,
                     alpha >= 0)
  
  prob <- Problem(obj, constraint)
  
  result <- solve(prob, warm_start = TRUE, reltol = 1e-12,
                  solver = solver,
                  verbose = FALSE)
  
  return(list(alphas = result$getValue(alpha), result = result, 
              Phi_test = Phi_test, Phi_train = Phi_train))
}

KLIEP_alphas_subset <- function(train, test, B, sigma,
                                solver = "CLARABEL"){
  #Same as KLIEP_alphas, except that only a subset, B,
  # of the test locations are used. 
  n <- nrow(train)
  m <- length(B)
  
  test_subset <- sample(1:nrow(test), B)

  test <- test[test_subset, ]
  
  return(KLIEP_alphas(train, test, sigma = sigma, solver = solver))
}

KLIEP_weights <- function(alphas, train, test, sigma){
  #Function to calculate the actual weights at some locations given
  #by train given by the estimated coefficients of the basis functions, 
  #centered at the locations given by test
  cross_dists <- cdist(train, test)
  
  Phi_train <- as.matrix(GRBF(cross_dists, sigma))
  
  return(Phi_train %*% alphas)
}

KLIEP_objective <- function(alpha, Phi_test){
  #Function to calculate the KLIEP objective function
  sum(log(Phi_test %*% alpha))
}

KLIEP_constraint <- function(alpha, Phi_train){
  #Function to check if the constraint in KLIEP is satisfied.
  print(sum(Phi_train %*% alpha))
  sum(Phi_train %*% alpha) == nrow(Phi_train)
}