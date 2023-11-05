library(CVXR)

proximal <- function(T_t, G, tau){
  p <-  dim(T_t)[1]  ## row size
  k <-  dim(T_t)[2]  ## col size
  
  X <-  Variable(rows = p, cols = k)
  prob <-  Problem(Minimize(norm(X - T_t, "F") + tau * mixed_norm(G %*% X, 2, 1)))
  CVXR_result <-  solve(prob, solver = "ECOS")
  
  return(CVXR_result$getValue(X))}
  
BiIPsPCA <- function(X, v, G, tau, max_it = 20){
  
  ### X: data matrix of size n*p with structures on p
  ### v: initialization of orthonormal Q with size p*k (or n*k)
  ### G: incidence matrix of given graph
  ### tau: regularization parameter
  Qa = v   
  t = 0
  off_diags = rep(0, max_it)
  while (t < max_it){
    Tb_t = X %*% Qa  ### space over n
    Qb = qr.Q(qr(Tb_t))  
    
    Ta_t = t(X) %*% Qb ### space over p
    Ta_t = proximal(Ta_t, G, tau)  ### incorporate structure over p
    Qa = qr.Q(qr(Ta_t))
    
    
    
    t = t+ 1
    R_t = t(Qb) %*% X %*% Qa
    
    
    off_diags[t] = norm(as.matrix(R_t[row(R_t)!=col(R_t)]), type = "F")
  }
  
  
  return (list("U" = Qb, "V" = Qa, "off_diags" = off_diags))  
}