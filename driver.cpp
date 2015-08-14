#include "eigen-qp.hpp"

#include <Eigen/Eigenvalues>
#include <boost/timer/timer.hpp>
#include <iostream>

using namespace std;
using namespace Eigen;
/*
 *  min 0.5 x.Q.x + c.x
 *
 *  Ax <= b
 *  Ex = d
 *
 * See: http://etd.dtu.dk/thesis/220437/ep08_19.pdf
 */
#define NV_FIXED 16
#define NC_FIXED 32

int main(int argc, char **argv)
{

	// Make a random problem
	int num_vars = NV_FIXED;
	int num_ineq = NC_FIXED;

	// Random matrices
	MatrixXd Q = MatrixXd::Random(num_vars,num_vars);
	Q *= Q.adjoint()/sqrt(num_vars); // Make it pos def

	VectorXd c = VectorXd::Random(num_vars);

	MatrixXd A = MatrixXd::Random(num_ineq,num_vars);
	VectorXd b = VectorXd::Random(num_ineq);

    VectorXd x_unc;
	// Solve unconstrainted system
	{
		boost::timer::auto_cpu_timer t;
		x_unc = -Q.ldlt().solve(c);
	}
    VectorXd x(num_vars);
    // Generate inequality constraints
    b.array() = (A*x_unc).array() - 0.5;
	// Inequality constrained problem
    {
        boost::timer::auto_cpu_timer t;

        quadprog(Q,c,A,b,x);
    }
    cout << "quadprog error: " << (x - x_unc).norm() << endl;

    // Fixed size
    Matrix<double, NV_FIXED,NV_FIXED> Q_fixed(Q);
    Matrix<double, NV_FIXED, 1> c_fixed(c);
    Matrix<double, NC_FIXED,NV_FIXED> A_fixed(A);
    Matrix<double, NC_FIXED,1> b_fixed(b);

    // Q_fixed.setRandom();
    // c_fixed.setRandom();
    // A_fixed.setRandom();

    x_unc = -Q_fixed.ldlt().solve(c_fixed);

    b_fixed.array() = (A_fixed*x_unc).array() - 0.12;

    Matrix<double, NV_FIXED, 1> x_fixed;
    {
        boost::timer::auto_cpu_timer t;
        for (int ii=0; ii < 100; ii++)
            quadprog(Q_fixed,c_fixed,A_fixed,b_fixed,x_fixed);
    }
    cout << "quadprog error: " << (x_fixed - x_unc).norm() << endl;

    cout << "Testing precaching" << endl;
    qp_solver<double,NV_FIXED,NC_FIXED> *solver;
    {
        boost::timer::auto_cpu_timer t;
        solver = new qp_solver<double,NV_FIXED,NC_FIXED>();
        solver->solve(Q_fixed,c_fixed,A_fixed,b_fixed,x_fixed);
    }
    cout << "qp_solver error: " << (x_fixed - x_unc).norm() << endl;
    {
        boost::timer::auto_cpu_timer t;
        solver->solve(Q_fixed,c_fixed,A_fixed,b_fixed,x_fixed);       
    }
    cout << "qp_solver error: " << (x_fixed - x_unc).norm() << endl;
    delete solver;
	return 0;	
}