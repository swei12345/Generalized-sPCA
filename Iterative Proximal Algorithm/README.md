### Iterative Proximal Algorithm
To incorporate structure over eigenvectors, we propose this iterative proximal algorithm. This algorithm is modified from the usual orthogonal iteration method to find the leading eigenvectors of a symmetric matrix S. 
We add an additional proximal projection step before the standard orthogonal iteration step. That is, at the t-th iteration, we add this step below:
```math
T_t^{(1)} = \arg\min_X \|X - T_t^{(0)}\|_F^2 + \tau\|\Gamma X\|_{2,1}$
```
