### Bi-Iterative Proximal Algorithm
Similar to the orthogonal iteration method, the classical Bi-iterative method extract both the left and the right leading eigenspace through QR factorization. Based on this idea, we can incorporate structures on both $n$ and $p$ dimensions by adding proximal projection steps. That is, compared to the `IPsPCA` method, we instead have the followings in each iteration (t):
```r
Tb_t = X %*% Qa  ### space over n (U)
### We can incorporate possible structures over n here ####
Qb = qr.Q(qr(Tb_t))  
    
Ta_t = t(X) %*% Qb ### space over p (V)
Ta_t = proximal(Ta_t, G, tau)  ### incorporate structures over p
Qa = qr.Q(qr(Ta_t))
```
Eventually we obtain `Qb` as the estimate for $U$ and `Qa` for $V$ in $X_k = U\Sigma_kV^T$. At convergence, we should expect that $Qb^T X Qa$ is a diagonal matrix approximating $\Sigma_k$. As a metric for such convergence, we check the norm of all off-diagonal entries of $Qb_t^T X Qa_t$ at each iteration. The below plot shows the log of norms when `tau = 0.25` as iterations $t$ increase. The details can be found in `BiIPsPCA_experiment1.Rmd`
<p align="center">
<img src="https://github.com/swei12345/Generalized-sPCA/assets/114754235/77d93cee-8972-4a55-8d4d-a67493e02ea3" width="650" height="350"> 
</p>

#### Currently, the iterative algorithms are sometimes very slow. The most time-consuming step is the proximal projection step implemented naively with `CVXR`. Such implementations use the large, but very sparse incidence matrix of a k-NN graph, as a dense matrix in the regularization. This leads to issues on storage and runtime.
