# Note that this function was written before deciding to incorporate the scaling 
# with the bandwidth into the kernel, so the factor h^(-d) is present in many
# places as in Cronie
calculate_Gram <- function(h, dist_mat){
  # Function to calculate the gram matrix for a given
  # bandwidth h (same as sigma) given a matrix
  # of distances between objects
  upper_K <- 1 / (2 * pi) * exp(-1 / 2 * 1 / h^2 * dist_mat^2)
  
  K <- as.matrix(upper_K)
  
  diag(K) <- 1 / (2 * pi)
  
  return(K)
}

calculate_edge_correction <- function(h, data,
                                      a = 0, b = 1, c = 0, d = 1){
  # Function to calculate a vector of edge correction factors
  # as presented in the report
  edge_correction <- foreach(dataId = 1 : nrow(data),
                             .combine = "c",
                             .packages=c("pracma")) % dopar % {
                               data_point <- data[dataId, ]
                               x <- data_point[1]
                               
                               y <- data_point[2]
                               
                               (pnorm((b - x) / h) - pnorm((a - x) / h)) * (pnorm((d - y)/ h) - pnorm((c - y) / h))
                             }
  return(edge_correction)
}

calc_edge_correction <- function(x, y, h = 1, a = 0, b = 1, c = 0, d = 1){
  # Functoin to calculate edge correction at a given location for a bandwidth h
  (pnorm((b - x) / h) - pnorm((a - x) / h)) * (pnorm((d - y) / h) - pnorm((c - y) / h))
}

calc_lhat_global <- function(h, w, Gram, d = 2){
  #Calculating the intensity estimate using global edge
  #correction
  h^(-d) * 1 / w * (Gram %*% rep(1, ncol(Gram)))
}

calc_lhat_local <- function(h, w, Gram, d = 2){
  #Calculating the intensity estimate using local
  #edge correction
  h^(-d) * Gram %*% (1/w)
}


T_k <- function(h, data, dist_mat, correction = "data"){
  #Function to calculate the estimate of the area
  #of the domain
  Gram <- calculate_Gram(h, dist_mat)
  if (correction == "global") {
    w <- calculate_edge_correction(h, data)
    return(sum(1 / calc_lhat_global(h, w, Gram)))
  }
  else if (correction == "local") {
    w <- calculate_edge_correction(h, data)
    return(sum(1 / calc_lhat_local(h, w, Gram)))
  }
  else { #no edge correction
    return(sum(1 / calc_lhat_global(h, rep(1, nrow(data)), Gram)))
  }
}

bandwidth_objective <- function(h, A = 1, Dist_mat = dist_mat, Data=data, correction="data"){
  #Objective function to minimize for the bandwidth
  return((T_k(h, Data, Dist_mat, correction) - A)^2)
}
