#include "eigen-qp.hpp"

#include <Eigen/Eigenvalues>

#include <iostream>

using namespace std;
using namespace Eigen;
/*
 *  min 0.5 x.Q.x + c.x
 *
 *  Ax <= b
 *  Ex = d
 */

int main(int argc, char **argv)
{

	// Make a random problem
	int num_vars = 5;
	int num_ineq = 3;
	int num_eq   = 1;

	// Random matrices
	MatrixXd Q = MatrixXd::Random(num_vars,num_vars);
	Q *= Q.adjoint()/sqrt(num_vars); // Make it pos def

	VectorXd c = VectorXd::Random(num_vars);

	MatrixXd A = MatrixXd::Random(num_ineq,num_vars);
	VectorXd b = VectorXd::Random(num_ineq);

	MatrixXd E = MatrixXd::Random(num_eq,num_vars);
	VectorXd d = VectorXd::Random(num_eq);

	// Solve unconstrainted system
	{
		VectorXd x_unconstrained = Q.ldlt().solve(c);
		cout << "Unconstrained solution: " << endl << "  " << x_unconstrained << endl;
	}
	// Solve equality constrained system: easy
	{
		MatrixXd aug_A(num_vars+num_eq,num_vars+num_eq);
		VectorXd aug_b(num_vars+num_eq);

		aug_A.block(0,0,num_vars,num_vars) = Q;
		aug_A.block(num_vars,0,num_eq,num_vars) = E;
		aug_A.block(0,num_vars,num_vars,num_eq) = E.transpose();
		aug_A.block(num_vars,num_vars,num_eq,num_eq).setZero();

		aug_b.head(num_vars) = c;
		aug_b.tail(num_eq)   = d;

		VectorXd x_equality = aug_A.ldlt().solve(aug_b).head(num_vars);
		cout << "Equality constrained solution: " << endl << x_equality << endl;
	}
	// OK, now the big guns: the full problem
	{
		
	}
	return 0;	
}