#include "solverSLE.h"

void solve_SLE(const Operator &A, const Field& rhs, Field& x, const Field& x0)
{
	Eigen::setNbThreads(2);
	//prepare GMRES solver
	bicgstab solver(A);	
	//gmres solver(A);	
	solver.setTolerance(1e-9);
	solver.setMaxIterations(300);
	//gmres.set_restart(6);

	//Solve equation, distribute border values, solve again...
	//Here, a fixed number of iterations is used
	//int i = 0;
	//do {
		x = solver.solveWithGuess(rhs, x0);
		//x = solver.solve(rhs);
	//	comm.move_fields(E, 3);
	//	i++;
	//} while(gmres.info() == Eigen::Success && i < SCHWARZ_ITERS);

	if(solver.info() != Eigen::Success) {
		std::cout << "Field GMRES failed!" << std::endl;
		std::cout << solver.error() << std::endl;
	}
	Eigen::setNbThreads(128);
}
