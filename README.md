# Eigen-QP #

A C++ template library for solving small quadradic programs very quickly.  Based on the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) Linear algebra template library.  The primary design goal is to be small, fast, and simple.

Currently supports the following features:
*   Quadradic programs with only equality constraints
*   Quadradic programs with only inequality constraints
    This is implemented using primal/dual interior point methods (specifically, the Mehrotra predictor/corrector method).  The algorithm follows the exposition in [this Masters thesis](http://etd.dtu.dk/thesis/220437/ep08_19.pdf)
*   Dynamic and fixed array sizes, as well as both float/double support
    Fixed size arrays can improve performance for very small problems (e.g., 4 dimensions) by close to an order of magnitude
*   Experimental [Julia](http://julialang.org) implementation written entirely in Julia (to avoid needing external solvers).  Tested with Julia 0.3.11

Things missing/TODO:
*   Combined equality/inequality constraints
*   Specialized code for box constraints (e.g., $l_i \le x_i \le u_i$)
*   Non-convergence/failure testing (currently assumes problem is feasible and convex)
*   (Stretch goal) Quadradically constrained quadradic programs (QCQPs)
