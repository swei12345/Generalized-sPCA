### Iterative Proximal Algorithm
To incorporate structure over eigenvectors, we propose this iterative proximal algorithm. This algorithm is modified from the usual orthogonal iteration method to find the leading eigenvectors of a symmetric matrix S. 
We add an additional proximal projection step before the standard orthogonal iteration step. That is, at the t-th iteration, we add this step below:
```math
T_t^{(1)} = \arg\min_X \|X - T_t^{(0)}\|_F^2 + \tau\|\Gamma X\|_{2,1}$
```
### R implementation
A naive implementation of this Iterative Proximal algorithm for sparse PCA `IPsPCA` is written in R. The proximal step is computed by the `CVXR` R package. 
Currently, codes for raising exceptions and errors are missing. E.g., need to check if $\tau$ is too large so that $T_t^{(1)}$ is not of full column rank, and QR factorization fails. 

### Example: spatial clustering
For deatils, see the Rmarkdown file `NAME` (R version 4.3.2) 
Suppose we have a spatial genetic data matrix $X$ of size $n \times p$ and we want to estimate its k leading left eigenvectors $U_k$, which are aligned with respect to a graph $\Gamma$. We can apply `IPsPCA` on $S = X X^T/n$. 

In this example, we generate a spatial genetic data matrix from [SRTsim](https://jiaqiangzhu.shinyapps.io/srtsim/)

