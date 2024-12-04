#include "solverSLE.h"

void solve_SLE(const Operator &A, const Field& rhs, Field& x, const Field& x0)
{
	//prepare GMRES solver
	BicgstabSolver solver(A);	
	//gmres solver(A);	
	solver.setTolerance(1e-16);
	solver.setMaxIterations(300);

	x = solver.solveWithGuess(rhs, x0);

	if(solver.info() != true) {
		std::cout << "Field GMRES failed!" << std::endl;
		std::cout << solver.error() << std::endl;
	}
}

bool BicgstabSolver::solve(const Operator &mat, const Eigen::VectorXd &rhs,
                           Eigen::VectorXd &x, DiagonalPreconditioner &precond,
                           size_t &iters, double &tol_error) {
    using std::abs;
    using std::sqrt;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorType;
    double tol = tol_error;
    int maxIters = iters;
    int n = mat.cols();
    VectorType r = rhs - mat * x;
    VectorType r0 = r;
    double r0_sqnorm = r0.squaredNorm();
    double rhs_sqnorm = rhs.squaredNorm();
    if (rhs_sqnorm == 0) {
        x.setZero();
        return true;
    }
    double rho = 1;
    double alpha = 1;
    double w = 1;
    VectorType v = VectorType::Zero(n), p = VectorType::Zero(n);
    VectorType y(n), z(n);
    VectorType s(n), t(n);
    double tol2 = tol * tol * rhs_sqnorm;
    double eps2 = Eigen::NumTraits<double>::epsilon() * Eigen::NumTraits<double>::epsilon();
    int i = 0;
    int restarts = 0;
    while (r.squaredNorm() > tol2 && i < maxIters) {
        double rho_old = rho;
        rho = r0.dot(r);
        if (abs(rho) < eps2 * r0_sqnorm) {
            r = rhs - mat * x;
            r0 = r;
            rho = r0_sqnorm = r.squaredNorm();
            if (restarts++ == 0)
                i = 0;
        }
        double beta = (rho / rho_old) * (alpha / w);
        p = r + beta * (p - w * v);
        y = precond.solve(p);   // Применение предобуславливателя
        v.noalias() = mat * y;
        alpha = rho / r0.dot(v);
        s = r - alpha * v;
        z = precond.solve(s);   // Применение предобуславливателя
        t.noalias() = mat * z;
        double tmp = t.squaredNorm();
        if (tmp > 0)
            w = t.dot(s) / tmp;
        else
            w = 0;
        x += alpha * y + w * z;
        r = s - w * t;
        ++i;
    }
    tol_error = sqrt(r.squaredNorm() / rhs_sqnorm);
    iters = i;
	m_error = tol_error;
    return true;
}
