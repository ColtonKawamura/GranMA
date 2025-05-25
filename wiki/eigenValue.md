## Eigenvalue and Eigen Vector Calulation

The intent of this workflow is to take a granular packing and get the eigenvalues and eigenvectors of these packings in damped conditions. We solve 

$$
\mathbf{M} \ddot{\mathbf{x}} + \mathbf{\Gamma} \dot{\mathbf{x}} + \mathbf{K} \mathbf{x} = \mathbf{0}
$$

 where $\mathbf{M}$ is the damping matrix and $\mathbf{\Gamma}$ is the damping matrix, and $\mathbf{K}$ is the stiffness matrix.


These matrices are generator from `matSpringDampMass()` within `processEigenModesDampedPara()`.

MATLAB's `polyeig()` solver assumes a solution of the form $\mathbf{x} = \mathbf{v} e^{\lambda t}$ and for eigenvalues can rewrite the equation as:

$$
\lambda = -\frac{\mathbf{\Gamma}}{2\mathbf{M}} \pm \sqrt{ \frac{\mathbf{\Gamma}^2}{4\mathbf{M}^2} - \frac{\mathbf{K}}{\mathbf{M}} }.
$$

For the underdamped case, we have:
$$
\lambda = -\frac{\mathbf{\Gamma}}{2\mathbf{M}} \pm i \sqrt{ \frac{\mathbf{K}}{\mathbf{M}} - \left(\frac{\mathbf{\Gamma}}{2\mathbf{M}}\right)^2 }.
$$

Finally we define $\omega_0 = \sqrt{\frac{\mathbf{K}}{\mathbf{M}}}$ and $b = \frac{\mathbf{\Gamma}}{2\sqrt{\mathbf{K}\mathbf{M}}}$, which gives us the eigenvalues in the form:

$$
\lambda = -b \omega_0 \pm i \omega_0 \sqrt{1 - b^2}
$$





**Relevant Files:**
- [processEignModesDampredPara.m](https://github.com/ColtonKawamura/GranMA/src/matlab_functions/processEignModesDampredPara.m)
- [processEigenModesDampedPara.m](https://github.com/ColtonKawamura/GranMA/src/matlab_functions/processEignModesDampredPara.m)
- [generateQEP.m](https://github.com/ColtonKawamura/GranMA/blob/main/src/eigen/generateQEP.m)

**To run:**
```matlab
addpath('src/eigen')
eigSolver('inputMatrix.mat')
```