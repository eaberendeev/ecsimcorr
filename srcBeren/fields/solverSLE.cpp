#include "solverSLE.h"

void solve_SLE(const Operator &A, const Field& rhs, Field& x, const Field& x0)
{
	//prepare GMRES solver
	bicgstab solver(A);	
	//gmres solver(A);	
	solver.setTolerance(1e-9);
	solver.setMaxIterations(300);

	x = solver.solveWithGuess(rhs, x0);

	if(solver.info() != Eigen::Success) {
		std::cout << "Field GMRES failed!" << std::endl;
		std::cout << solver.error() << std::endl;
	}
}
