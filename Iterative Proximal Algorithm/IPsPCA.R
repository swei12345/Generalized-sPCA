library(CVXR)

proximal <- function(T_t, G, tau){
  p <-  dim(T_t)[1]  ## row size
  k <-  dim(T_t)[2]  ## col size
  X <-  Variable(rows = p, cols = k)
  prob <-  Problem(Minimize(norm(X - T_t, "F") + tau * mixed_norm(G %*% X, 2, 1)))
  CVXR_result <-  solve(prob, solver = "ECOS")
  
  return(CVXR_result$getValue(X))}
  
IPsPCA <- function(S, v, G, tau, max_it = 20){
  Q = v   ### need to check if Q full-colum rank
  T_t = S %*% Q
  t = 0
  Q_ts = rep(0, max_it)
  while (t < max_it){
      
    T_t = proximal(T_t, G, tau)
      
    Q = qr.Q(qr(T_t))  
      
    T_t = S %*% Q
      
      
    t = t+ 1
      Q_ts[t] = Q}
    
    Pi = Q %*% t(Q)
    return (list("Pi" = Pi, "Q" = Q, "Q_ts" = Q_ts))  
  }