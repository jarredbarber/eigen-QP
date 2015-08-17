#include <Eigen/Eigenvalues>
#include <boost/timer/timer.hpp>
#include <iostream>

#include "eigen-qp.hpp"

using namespace std;
using namespace Eigen;
using namespace EigenQP;
/*
 *  min 0.5 x.Q.x + c.x
 *
 *  Ax <= b
 *  Ex = d
 *
 * See: http://etd.dtu.dk/thesis/220437/ep08_19.pdf
 */
// #define NV_FIXED 8
// #define NC_FIXED 16
// #define NE_FIXED 3

#define PROB_VARS 8
#define PROB_INEQ 16
#define PROB_EQ 3

#define N_TEST   2048

template<typename Scalar, int NV_FIXED, int NC_FIXED, int NE_FIXED>
void test()
{

    typedef Eigen::Matrix<Scalar,-1,-1> MatrixXs;
    typedef Eigen::Matrix<Scalar,-1,1> VectorXs;

    // Make a random problem
    int num_vars = NV_FIXED;
    int num_ineq = NC_FIXED;
    int num_eq   = NE_FIXED;

    // Random matrices
    MatrixXs Q = MatrixXs::Random(num_vars,num_vars);
    Q *= Q.adjoint()/sqrt(num_vars); // Make it pos def

    VectorXs c = VectorXs::Random(num_vars);

    MatrixXs A = MatrixXs::Random(num_ineq,num_vars);
    VectorXs b = VectorXs::Random(num_ineq);

    MatrixXs E = MatrixXs::Random(num_eq,num_vars);
    VectorXs f = VectorXs::Random(num_eq);

    VectorXs x_unc;
    // Solve unconstrainted system
    cout << "Unconstrained..." << endl;
    {
        boost::timer::auto_cpu_timer t;
        for (int ii=0; ii < N_TEST; ii++)
          x_unc = -Q.ldlt().solve(c);
    }
    VectorXs x(num_vars);
    // Generate inequality constraints
    b.array() = (A*x_unc).array() - 0.15;
    // Inequality constrained problem
    cout << "quadprog, dynamic code" << endl;
    {
        boost::timer::auto_cpu_timer t;
        for (int ii=0; ii < N_TEST; ii++)
            quadprog(Q,c,A,b,x);
    }
    cout << "    error: " << (x - x_unc).norm() << endl;
    {
        QPIneqSolver<Scalar,-1,-1> *solver;
        cout <<"QPIneqSolver, dynamic, obj creation" << endl;
        {
            boost::timer::auto_cpu_timer t;
            for (int ii=0; ii < N_TEST; ii++)
            {
                solver = new QPIneqSolver<Scalar,-1,-1>(num_vars,num_ineq);
                solver->solve(Q,c,A,b,x);
            }
        }
        cout << "    error: " << (x - x_unc).norm() << endl;
        cout << "QPIneqSolver, dynamic, object reuse" << endl;
        {
            boost::timer::auto_cpu_timer t;
            for (int ii=0; ii < N_TEST; ii++)
                solver->solve(Q,c,A,b,x);
        }
        cout << "    error: " << (x - x_unc).norm() << endl;
    }

    // Fixed size
    Matrix<Scalar, NV_FIXED,NV_FIXED> Q_fixed(Q);
    Matrix<Scalar, NV_FIXED, 1> c_fixed(c);
    Matrix<Scalar, NC_FIXED,NV_FIXED> A_fixed(A);
    Matrix<Scalar, NC_FIXED,1> b_fixed(b);

    // Q_fixed.setRandom();
    // c_fixed.setRandom();
    // A_fixed.setRandom();

    x_unc = -Q_fixed.ldlt().solve(c_fixed);

    b_fixed.array() = (A_fixed*x_unc).array() - 0.12;

    Matrix<Scalar, NV_FIXED, 1> x_fixed;
    cout << "quadprog, fixed code" << endl;
    {
        boost::timer::auto_cpu_timer t;
        for (int ii=0; ii < N_TEST; ii++)
            quadprog(Q_fixed,c_fixed,A_fixed,b_fixed,x_fixed);
    }
    cout << "    error: " << (x_fixed - x_unc).norm() << endl;
    {
        QPIneqSolver<Scalar,NV_FIXED,NC_FIXED> *solver;
        cout << "QPIneqSolver, fixed" << endl;
        {
            boost::timer::auto_cpu_timer t;
            for (int ii=0; ii < N_TEST; ii++)
            {
                solver = new QPIneqSolver<Scalar,NV_FIXED,NC_FIXED>(num_vars,num_ineq);
                solver->solve(Q_fixed,c_fixed,A_fixed,b_fixed,x_fixed);
            }
        }
        cout << "    error: " << (x_fixed - x_unc).norm() << endl;
        cout << "QPIneqSolver, fixed, object reuse" << endl;
        {
            boost::timer::auto_cpu_timer t;
            for (int ii=0; ii < N_TEST; ii++)
                solver->solve(Q_fixed,c_fixed,A_fixed,b_fixed,x_fixed);       
        }
        cout << "    error: " << (x_fixed - x_unc).norm() << endl;
    }

    cout << "quadprog, equality constraints, dynamic" << endl;
    {
        QPEqSolver<Scalar> solver(num_vars,num_eq);
        {
            boost::timer::auto_cpu_timer t;
            for (int ii=0; ii < N_TEST; ii++)
                solver.solve(Q,c,E,f,x);
        }
    }

    cout << "quadprog, ineq/eq constraints, dynamic" << endl;
    {
        //b = A*x - 0.12;
        QPGenSolver<Scalar> solver(num_vars,num_ineq,num_eq);
        {
            boost::timer::auto_cpu_timer t;
            for (int ii=0; ii < N_TEST; ii++)
                solver.solve(Q,c,A,b,E,f,x);
        }
    }
}

int main(int argc, char ** argv)
{
    srand((unsigned int) time(0));
    cout << "TESTING DOUBLE" << endl;
    test<double,PROB_VARS,PROB_INEQ,PROB_EQ>();
    cout << "TESTING FLOAT" << endl;
    test<float,PROB_VARS,PROB_INEQ,PROB_EQ>();
    return 0;
}