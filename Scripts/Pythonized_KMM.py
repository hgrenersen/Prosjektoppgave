# Python code written to speed up the calculation of KMM IWs
# compared to the approach implemented in R.

import numpy as np
from cvxopt import matrix, solvers
import math

def py_KMM_optimization(K, kappa, B=1000):
  # K: gram matrix
  # kappa: kappa-vector from KMM
  # B: upper bound on weights

  
  K_cop = np.copy(np.array(K))

  kappa_cop = np.copy(np.array(kappa))
  
  n = K_cop.shape[0]
  epsilon=1-1/math.sqrt(n)
  
  # Upper limit:
  h = B * np.ones(n)
  G = np.eye(n)
  
  #Non-negativity constraints, on the form Gx <= h (elementwise)
  h = np.append(h, np.zeros(n))
  G = np.vstack((G, - np.eye(n)))
  
  # Soft constraints on the mean
  h = np.append(h, np.array([n*(epsilon-1), n*(epsilon +1)]))
  G = np.vstack((G, np.array([-np.ones(n), np.ones(n)])))
  
  G = matrix(G)
  h = matrix(h)
  
  K_cop = matrix(K_cop) #No 1/2 here, since both KMM and this library
  #solves the problem 1/2x^T K x + kappa^T x
  kappa_cop= -matrix(kappa_cop) #Need a minus here since there's a 
  # plus in the problem solved by the solver
  
  solvers.options['max_iter']=100
  
  solvers.options['abstol']=1e-7
  solvers.options['reltol']=1e-6
  solvers.options['feastol']=1e-7
  solvers.options['show_progress']=False
  weights = np.array(solvers.qp(K_cop, kappa_cop, G, h)['x']).ravel()
  return weights
