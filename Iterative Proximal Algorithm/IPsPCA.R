library(CVXR)

proximal <- function(T_t, G, tau){
  p <-  dim(T_t)[1]  ## row size
  k <-  dim(T_t)[2]  ## col size
  
  X <-  Variable(rows = p, cols = k)
  prob <-  Problem(Minimize(norm(X - T_t, "F") + tau * mixed_norm(G %*% X, 2, 1)))
  CVXR_result <-  solve(prob, solver = "ECOS")
  
  return(CVXR_result$getValue(X))}
  
IPsPCA <- function(S, v, G, tau, max_it = 20){
  
  ### S: sample "covariance" matrix; leading eigenvectors to be estimated
  ### v: initialization of orthonormal Q with size p*k (or n*k)
  ### G: incidence matrix of given graph
  ### tau: regularization parameter
  Q = v   
  
  t = 0
  
  while (t < max_it){
    T_t = S %*% Q  
    T_t = proximal(T_t, G, tau)
    Q = qr.Q(qr(T_t))  
    
    t = t+ 1}
    
    
    return (list("Q" = Q))  
  }