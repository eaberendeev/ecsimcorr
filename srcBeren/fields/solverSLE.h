#ifndef SOLVERSLE_H_
#define SOLVERSLE_H_

#include "util.h"
#include "containers.h"
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/IterativeLinearSolvers>

#define SOLVER_TOLERANCE 1.e-6
// System of linear equations solvers
typedef Eigen::GMRES<Eigen::SparseMatrix<double, MAJOR> > gmres;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<double, MAJOR> > bicgstab;

void solve_SLE(const Operator &A, const Field& rhs, Field& x, const Field& x0);

class DiagonalPreconditioner {
   public:
    DiagonalPreconditioner(const Eigen::VectorXd &diagonal)
        : m_diagonal(diagonal) {}
    Eigen::VectorXd solve(const Eigen::VectorXd &p) const {
        return p.array() / m_diagonal.array();
    }

   private:
    Eigen::VectorXd m_diagonal;
};

class BicgstabSolver {
   public:
    BicgstabSolver(const Operator &A) : m_A(A), m_precond(A.diagonal()) { 
        max_iterations = 1000;
        m_tolerance = 1e-16;
    }
    void setTolerance(double tolerance) { m_tolerance = tolerance; }
    void setMaxIterations(size_t max_iterations) {
        max_iterations = max_iterations;
    }
    Eigen::VectorXd solveWithGuess(const Eigen::VectorXd &rhs, const Eigen::VectorXd &x0) {
        Eigen::VectorXd x = x0;
        m_success = solve(m_A, rhs, x, m_precond, max_iterations, m_tolerance);
        return x;
    }
    bool info() const { return m_success; }
    double error() const { return m_error; }
   private:
    const Operator &m_A;
    DiagonalPreconditioner m_precond;
    bool solve(const Operator &mat, const Eigen::VectorXd &rhs,
               Eigen::VectorXd &x, DiagonalPreconditioner &precond, size_t &iters,
               double &tol_error);
    size_t max_iterations;
    double m_tolerance;
    double m_error;
    bool m_success;
};

#endif 	
