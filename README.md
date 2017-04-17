# Eigen-QP #

A C++ template library for solving small quadradic programs very quickly.  Based on the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) Linear algebra template library.  The primary design goal is to be small, fast, and simple.

Currently supports the following features:
*   Quadradic programs with only equality constraints
*   Quadradic programs with only inequality constraints. This is implemented using primal/dual interior point methods (specifically, the Mehrotra predictor/corrector method).  The algorithm follows the exposition in [this Masters thesis](http://etd.dtu.dk/thesis/220437/ep08_19.pdf)
*   Dynamic and fixed array sizes, as well as both float/double support.  Fixed size arrays can improve performance for very small problems (e.g., 4 dimensions) by close to an order of magnitude
*   Experimental [Julia](http://julialang.org) implementation written entirely in Julia (to avoid needing external solvers).  Tested with Julia 0.3.11

Things missing/TODO:
- [ ]   Combined equality/inequality constraints
- [ ]   Specialized code for box constraints (e.g., $l_i \le x_i \le u_i$)
- [ ]   Non-convergence/failure testing (currently assumes problem is feasible and convex)
- [ ]   (Stretch goal) Quadradically constrained quadradic programs (QCQPs)
- [ ]   Currently, the <float> template specialization seems to not converge.  This could just be a sign that the algorithm requires double precision, but this should be investigated.

## Installation and Usage ##

Keeping with the Eigen header-only philosophy, the only installation needed is to copy `eigen-qp.hpp` into your include path and then:

```C++
#include "eigen-qp.hpp"
```

This also requires a C++11 compiliant compiler, which in the case of GCC and Intel requires the `--std=c++11` flag.

The main interface to the code can be accessed through an API that mimics MATLAB's `quadprog` function:

```C++
namespace EigenQP {
template<typename Scalar, int NVars, int NIneq>
void quadprog(Eigen::Matrix<Scalar,NVars,NVars> &Q, Eigen::Matrix<Scalar,NVars,1> &c, 
              Eigen::Matrix<Scalar,NIneq,NVars> &A, Eigen::Matrix<Scalar,NIneq,1> &b,
              Eigen::Matrix<Scalar,NVars,1> &x);
}
```

In this notation, quadprog solves the following inequality constrained quadradic program:

![Quadprog equation](https://latex.codecogs.com/png.latex?%5CLARGE%20%5Cmin_%7Bx%7D%20x%5EtQx%20&plus;%20b%5Etx%20%5C%5C%20%5Ctext%7BSuch%20that%7D%7E%20Ax%20%5Cpreceq%20c)

The template parameters `NVars` and `NIneq` can be set at compile time to allow the compiler to generate faster code for solving small problems; if set to -1 (or `Eigen::Dynamic`), arbitrary (real) matrices are allowed (although the user is responsible for assuring compatible dimensions).

For equality constrained problems, the class `EigenQP::QPEqSolver` can be instantiated.  This is provided for completeness, but the implementation is fairly trivial (simply using Eigen to do a block matrix inversion).
