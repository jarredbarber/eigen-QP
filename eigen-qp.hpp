#ifndef _EIGEN_QP_H_
#define _EIGEN_QP_H_

#include <memory>
#include <iostream>
#include <Eigen/Core>

using namespace Eigen;

template<typename Scalar, int NVars, int NIneq>
class qp_solver
{
    typedef Matrix<Scalar,NVars,1> PVec;
    typedef Matrix<Scalar,NIneq,1> DVec; // Dual (i.e., Lagrange multiplier) vector
    typedef Matrix<Scalar,NVars,NVars> PMat;
private:
    // Problem size
    const int n;
    const int m;

    const Scalar eps = 1E-9;
    const Scalar eta = 0.95;

    // Work buffers
    DVec s;
    DVec z;

    PVec rd;
    DVec rp;
    DVec rs;

    PVec dx;
    DVec ds;
    DVec dz;

public:
    qp_solver(int n_vars=NVars, int n_const=NIneq) : n(n_vars),m(n_const), s(m), z(m), rd(n), rp(m), rs(m), dx(n), ds(m), dz(m)
        {
        }

    ~qp_solver() {}

    int solve(Matrix<Scalar,NVars,NVars> &Q, Matrix<Scalar,NVars,1> &c, 
              Matrix<Scalar,NIneq,NVars> &A, Matrix<Scalar,NIneq,1> &b,
              Matrix<Scalar,NVars,1> &x)
    {
        // Initialization
        s.setOnes();
        z.setOnes();
        x.setZero();

        // Initial residuals. Uses fact that x=0 here.
        rd = c - A.adjoint()*z;
        rp = s + b;
        rs = (s.array()*z.array());

        Scalar ms = 1/(Scalar)m;
        Scalar mu = (Scalar)n*ms; // Initial mu based on knowing that s,z are ones.
        Scalar alpha;

        int iter;

        for (iter=0; iter < 250; iter++)
        {
            // Precompute decompositions for this iteration
            LLT<PMat> Gbar = (Q + A.adjoint()*((z.array()/s.array()).matrix().asDiagonal())*A).llt();

            for (int ii=0; ii < 2; ii++)
            {   
                // Prediction/correction step
                {
                    auto tmp = (rs.array() - z.array()*rp.array())/s.array();
                    dx = -Gbar.solve(rd + A.adjoint()*tmp.matrix());
                    ds = A*dx - rp;
                    dz.array() = -(rs.array() - z.array()*ds.array())/s.array();
                }

                // Compute alph,mu 
                alpha = 1.0;
                for (int jj=0; jj < m; jj++)
                {
                    Scalar a = -z(jj)/dz(jj);
                    alpha    = (a < alpha) && (a > 0) ? a : alpha;
                    a        = -s(jj)/ds(jj);
                    alpha    = (a < alpha) && (a > 0) ? a : alpha;
                }

                if (ii)
                    break; // Don't need to compute any more

                // Centering    
                Scalar mu_aff = (s + alpha*ds).dot(z+alpha*dz)*ms;
                Scalar sigma  = (mu_aff/mu); sigma *= sigma*sigma;

                // Corrector residual
                rs.array() += ds.array()*dz.array() - sigma*mu;
            }

            // Step
            alpha *= eta; // rescale step size
            x += alpha*dx;
            s += alpha*ds;
            z += alpha*dz;

            // Update residuals
            rd = Q*x + c - A.adjoint()*z;
            rp = s - A*x + b;
            rs = (s.array()*z.array());

            mu = s.dot(z)*ms;

            // Convergence test
            if ( (mu < eps) && 
                 (rd.norm() < eps) && 
                 (rs.norm() < eps) )
            {
                break;
            }
        }
        return iter;
    }
};


template<typename Scalar, int NVars, int NIneq>
void quadprog(Matrix<Scalar,NVars,NVars> &Q, Matrix<Scalar,NVars,1> &c, 
              Matrix<Scalar,NIneq,NVars> &A, Matrix<Scalar,NIneq,1> &b,
              Matrix<Scalar,NVars,1> &x)
{
    qp_solver<Scalar,NVars,NIneq> qp(c.size(),b.size());
    qp.solve(Q,c,A,b,x);
}
#endif