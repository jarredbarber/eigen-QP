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

void solve_qp(MatrixXd &Q, VectorXd &c, MatrixXd &A, VectorXd &b, VectorXd &x)
{
	int N = c.size();
	int m = b.size();

	// Slack vars
	VectorXd s = VectorXd::Random(m).array().abs().matrix();
	VectorXd z = VectorXd::Random(m).array().abs().matrix();
	x.setRandom();

	// Jacobian
	MatrixXd J(N + 2*m, N + 2*m );

	J.setZero();
	J.block(0,0,N,N) = Q;
	J.block(N,N+m,m,m).setIdentity();
	J.block(N,0,m,N) = -1.0*A;
	J.block(0,N,N,m) = -1.0*A.adjoint();

	// Residual thing?
	bool converged=false;

	VectorXd r(N+2*m);
	VectorXd delt_aff(N+2*m);

	r.head(N) = Q*x + c - A.transpose()*z;
	r.segment(N,m) = s - A*x + b;
	r.tail(m) = (s.array()*z.array()).matrix();

	double mu = s.dot(z)/m;
	double alpha = 1.0;


	cout << "mu = " << mu << endl;
	for (int iter=0; iter < 46; iter++)
	{
		// Build Jacobian
		J.block(N+m,N,m,m).diagonal() = s;
		J.block(N+m,N+m,m,m).diagonal() = z;
		auto Jinv = J.colPivHouseholderQr();

		delt_aff = Jinv.solve(-1.0*r);

		// Compute predictor step length
		alpha = 1.0;
		for (int jj=0; jj < m; jj++)
		{
			double a = -z(jj)/delt_aff(N+jj);
			alpha = (a < alpha) && (a >= 0) ? a : alpha;
			a = -s(jj)/delt_aff(N+m+jj);
			alpha = (a < alpha) && (a >= 0) ? a : alpha;
		}

		double mu_aff = (s + alpha*delt_aff.tail(m)).dot(z+alpha*delt_aff.segment(N,m))/m;
		double sigma  = (mu_aff/mu);
		sigma *= sigma*sigma; // sigma^3

		r.tail(m) +=  (delt_aff.tail(m).array()*delt_aff.segment(N,m).array() - sigma*mu).matrix();

		delt_aff = Jinv.solve(-1.0*r);

		// Compute alpha again
		alpha = 1.0;
		for (int jj=0; jj < m; jj++)
		{
			double a = -z(jj)/delt_aff(N+jj);
			alpha = (a < alpha) && (a > 0) ? a : alpha;
			a = -s(jj)/delt_aff(N+m+jj);
			alpha = (a < alpha) && (a > 0) ? a : alpha;
		}

		// Step
		x += alpha*delt_aff.head(N);
		z += alpha*delt_aff.segment(N,m);
		s += alpha*delt_aff.tail(m);

		// Convergence check
		mu = s.dot(z)/m;
		double step = alpha*delt_aff.norm();

		if (mu < 1E-6 && step < 1E-6)
		{
			cout << "Finished in " << iter << " iterations." << endl;
			break;
		}
		// cout << "iter: " << iter << " step: " << alpha*delt_aff.norm() << " mu: " << mu << endl;

		// Update residuals
		r.head(N) = Q*x + c - A.transpose()*z;
		r.segment(N,m) = s - A*x + b;
		r.tail(m) = (s.array()*z.array()).matrix();
	}

}
int main(int argc, char **argv)
{

	// Make a random problem
	int num_vars = 1000;
	int num_ineq = 30;
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
		boost::timer::auto_cpu_timer t;
		VectorXd x_unconstrained = Q.ldlt().solve(c);
	}
	// Solve equality constrained system: easy
	{
		boost::timer::auto_cpu_timer t;

		MatrixXd aug_A(num_vars+num_eq,num_vars+num_eq);
		VectorXd aug_b(num_vars+num_eq);

		aug_A.block(0,0,num_vars,num_vars) = Q;
		aug_A.block(num_vars,0,num_eq,num_vars) = E;
		aug_A.block(0,num_vars,num_vars,num_eq) = E.transpose();
		aug_A.block(num_vars,num_vars,num_eq,num_eq).setZero();

		aug_b.head(num_vars) = c;
		aug_b.tail(num_eq)   = d;

		VectorXd x_equality = aug_A.ldlt().solve(aug_b).head(num_vars);
	}
	// Inequality constrained problem
	{
		boost::timer::auto_cpu_timer t;

		VectorXd x(num_vars);
		solve_qp(Q,c,A,b,x);
	}
	return 0;	
}