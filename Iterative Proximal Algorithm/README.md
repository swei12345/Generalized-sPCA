### Iterative Proximal Algorithm
To incorporate structure over eigenvectors, we propose this iterative proximal algorithm. This algorithm is modified from the usual orthogonal iteration method to find the leading eigenvectors of a symmetric matrix S. 
We add an additional proximal projection step before the standard orthogonal iteration step. That is, at the t-th iteration, we add this step below:
```math
T_t^{(1)} = \arg\min_X \|X - T_t^{(0)}\|_F^2 + \tau\|\Gamma X\|_{2,1}$
```
### R implementation
A naive implementation of this Iterative Proximal algorithm for structured PCA `IPsPCA` is written in R. The proximal step is computed by the `CVXR` R package. 
Currently, codes for raising exceptions and errors are missing. E.g., need to check if $\tau$ is too large so that $T_t^{(1)}$ is not of full column rank, and QR factorization fails. 

### Example: spatial clustering
For details, see the Rmarkdown file `NAME` (R version 4.3.2) 
Suppose we have a spatial genetic data matrix $X$ of size $n \times p$ and we want to estimate its k leading left eigenvectors $U_k$, which are aligned with respect to a graph $\Gamma$. We can apply `IPsPCA` on $S = X X^T/n$. 

In this example, we generate a spatial genetic data matrix from [SRTsim](https://jiaqiangzhu.shinyapps.io/srtsim/) with locations of genes randomly spaced over a 2D grid. We can visualize the clusters as:
<p align="center">
<img src="https://github.com/swei12345/Generalized-sPCA/assets/114754235/3fae59d3-5119-4682-85cf-8785eb4db02c" width="350" height="300">
</p>


Based on the coordinates, we can form a k-NN graph. In the next step, we use the incidence matrix of this k-NN graph for the $\Gamma$ regularization for the algorithm `IPsPCA`. The tunings and max iterations are arbitrarily chosen. 
<p align="center">
<img src="https://github.com/swei12345/Generalized-sPCA/assets/114754235/cda52700-83cd-4297-9a09-5c2a3f3981af" width="250" height="250"> 
</p>

After we extract the structured $U$, we can apply kmeans clustering on each of its $n$ rows, which are vectors of size $k$. (That is, we treat the k-vectors as features for $n$ locations. This practice is not very meaningful. This is just to provide an example to show what we may do for the downstream analysis). The estimated clusters from the structured $U$ are shown on the **left**; clusters from standard PCA are shown on the **right**. 


<img src="https://github.com/swei12345/Generalized-sPCA/assets/114754235/e55e7650-2d0f-4cf6-8cb4-2f351d72c257" width="300" height="250" align = "left"> 
<img src="https://github.com/swei12345/Generalized-sPCA/assets/114754235/9b7121cf-1f41-41a4-ac3f-b7f3d86b5aa5" width="300" height="250" align = "right"> 

Again, since the parameters, prepocessings, and the clustering techniques are arbitrarily chosen, the numerical results are not that meaningful or interpretable, and its only purpose is to show a possible pipeline of analysis. 
