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

    const Scalar eps = 1E-8;
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

    void solve(Matrix<Scalar,NVars,NVars> &Q, Matrix<Scalar,NVars,1> &c, 
              Matrix<Scalar,NIneq,NVars> &A, Matrix<Scalar,NIneq,1> &b,
              Matrix<Scalar,NVars,1> &x)
    {
        // Initialization
        s.setOnes();
        z.setOnes();
        x.setZero();

        rd = Q*x + c - A.transpose()*z;
        rp = s - A*x + b;
        rs = (s.array()*z.array()).matrix();

        Scalar ms = (Scalar)m;
        Scalar mu = (Scalar)n/ms; // Initial mu based on knowing that s,z are ones.
        Scalar alpha;

        int iter;

        for (iter=0; iter < 250; iter++)
        {
            // Precompute decompositions for this iteration
            LLT<PMat> Gbar = (Q + A.transpose()*((z.array()/s.array()).matrix().asDiagonal())*A).llt();

            for (int ii=0; ii < 2; ii++)
            {   
                // Prediction step
                {
                    auto tmp = (rs.array() - z.array()*rp.array())/s.array();
                    dx = -Gbar.solve(rd + A.transpose()*tmp.matrix());
                    ds = A*dx - rp;
                    dz.array() = -(rs.array() - z.array()*ds.array())/s.array();
                }

                // Compute alph,mu 
                alpha = 1.0;
                for (int jj=0; jj < m; jj++)
                {
                    Scalar a = -z(jj)/dz(jj);
                    alpha = (a < alpha) && (a > 0) ? a : alpha;
                    a = -s(jj)/ds(jj);
                    alpha = (a < alpha) && (a > 0) ? a : alpha;
                }

                if (ii)
                    break; // Don't need to compute any more

                Scalar mu_aff = (s + alpha*ds).dot(z+alpha*dz)/ms;

                // Centering    
                Scalar sigma = (mu_aff/mu); sigma *= sigma*sigma;

                // Corrector
                rs.array() += ds.array()*dz.array() - sigma*mu;
            }

            // Step
            alpha *= eta;
            x += alpha*dx;
            s += alpha*ds;
            z += alpha*dz;

            // Update residuals
            rd = Q*x + c - A.transpose()*z;
            rp = s - A*x + b;
            rs.array() = (s.array()*z.array());

            mu = s.dot(z)/ms;

            // Convergence test
            if ( (mu < eps) && 
                 (rd.norm() < eps) && 
                 (rs.norm() < eps) )
            {
                break;
            }
        }
    }
};

template<typename Scalar, int NVars, int NIneq>
void quadprog(Matrix<Scalar,NVars,NVars> &Q, Matrix<Scalar,NVars,1> &c, 
              Matrix<Scalar,NIneq,NVars> &A, Matrix<Scalar,NIneq,1> &b,
              Matrix<Scalar,NVars,1> &x)
{
    typedef Matrix<Scalar,NVars,1> PVec;
    typedef Matrix<Scalar,NIneq,1> DVec; // Dual (i.e., Lagrange multiplier) vector
    typedef Matrix<Scalar,NVars,NVars> PMat;

    const int n = c.size();
    const int m = b.size();
    const Scalar eps = 1E-8;
    const Scalar eta = 0.95;

    // Slack variables
    DVec s(m);
    DVec z(m);

    PVec rd(n);
    DVec rp(m);
    DVec rs(m);

    PVec dx(n);
    DVec ds(m);
    DVec dz(m);

    // Initialization
    s.setOnes();
    z.setOnes();
    x.setZero();

    rd = Q*x + c - A.transpose()*z;
    rp = s - A*x + b;
    rs = (s.array()*z.array()).matrix();

    Scalar ms = (Scalar)m;
    Scalar mu = (Scalar)n/ms; // Initial mu based on knowing that s,z are ones.
    Scalar alpha;

    int iter;

    for (iter=0; iter < 250; iter++)
    {
        // Precompute decompositions for this iteration
        LDLT<PMat> Gbar = (Q + A.transpose()*((z.array()/s.array()).matrix().asDiagonal())*A).ldlt();

        for (int ii=0; ii < 2; ii++)
        {   
            // Prediction step
            {
                auto tmp = (rs.array() - z.array()*rp.array())/s.array();
                dx = -Gbar.solve(rd + A.transpose()*tmp.matrix());
                ds = A*dx - rp;
                dz.array() = -(rs.array() - z.array()*ds.array())/s.array();
            }

            // Compute alph,mu 
            alpha = 1.0;
            for (int jj=0; jj < m; jj++)
            {
                Scalar a = -z(jj)/dz(jj);
                alpha = (a < alpha) && (a > 0) ? a : alpha;
                a = -s(jj)/ds(jj);
                alpha = (a < alpha) && (a > 0) ? a : alpha;
            }

            if (ii)
                break; // Don't need to compute any more

            Scalar mu_aff = (s + alpha*ds).dot(z+alpha*dz)/ms;

            // Centering    
            Scalar sigma = (mu_aff/mu); sigma *= sigma*sigma;

            // Corrector
            rs.array() += ds.array()*dz.array() - sigma*mu;
        }

        // Step
        alpha *= eta;
        x += alpha*dx;
        s += alpha*ds;
        z += alpha*dz;

        // Update residuals
        rd = Q*x + c - A.transpose()*z;
        rp = s - A*x + b;
        rs.array() = (s.array()*z.array());

        mu = s.dot(z)/ms;

        // Convergence test
        if ( (mu < eps) && 
             (rd.norm() < eps) && 
             (rs.norm() < eps) )
        {
            break;
        }
    }
}

#endif