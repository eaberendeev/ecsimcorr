#pragma once

#include <type_traits>
#include "bmatrix.h"
#include "containers.h"
#include "util.h"

template <typename T>
struct is_valid_matrix_type : std::false_type {};

template <>
struct is_valid_matrix_type<BMatrix> : std::true_type {};

template <>
struct is_valid_matrix_type<Operator> : std::true_type {};

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

template <typename MatrixType>
class BicgstabSolver {
    static_assert(is_valid_matrix_type<MatrixType>::value,
                  "MatrixType must be either BMatrix or Operator");
   public:
    BicgstabSolver(const MatrixType &A) : m_A(A), m_precond(A.diagonal()) {
        max_iterations = 1000;
        m_tolerance = 1e-16;
    }
    void setTolerance(double tolerance) { m_tolerance = tolerance; }
    void setMaxIterations(size_t max_iters) {
        max_iterations = max_iters;
    }
    Eigen::VectorXd solveWithGuess(const Eigen::VectorXd &rhs, const Eigen::VectorXd &x0) {
        Eigen::VectorXd x = x0;
        m_iterations = max_iterations;
        m_error = m_tolerance;
        m_success = solve(m_A, rhs, x, m_precond, m_iterations, m_error);
        return x;
    }
    bool info() const { return m_success; }
    double error() const { return m_error; }
    size_t iterations() const { return m_iterations; }

   private:
    const MatrixType &m_A;
    DiagonalPreconditioner m_precond;
    bool solve(const MatrixType &mat, const Eigen::VectorXd &rhs,
               Eigen::VectorXd &x, DiagonalPreconditioner &precond,
               size_t &iters, double &tol_error);
    size_t max_iterations;
    size_t m_iterations;
    double m_tolerance;
    double m_error;
    bool m_success;
};

template <typename MatrixType>
bool BicgstabSolver<MatrixType>::solve(const MatrixType &mat,
                                       const Eigen::VectorXd &rhs,
                                       Eigen::VectorXd &x,
                                       DiagonalPreconditioner &precond,
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
    double eps2 = Eigen::NumTraits<double>::epsilon() *
                  Eigen::NumTraits<double>::epsilon();
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
    return true;
}

template <typename MatrixType, typename SolverType>
void solve_linear_system(const MatrixType &A, const Field &rhs, Field &x,
                         const Field &x0) {
    static_assert(is_valid_matrix_type<MatrixType>::value,
                  "MatrixType must be either BMatrix or Operator");
    // Create solver instance
    SolverType solver(A);
    solver.setTolerance(SLE_SOLVER_TOLERANCE);
    solver.setMaxIterations(SLE_SOLVER_MAX_ITERATIONS);

    x = solver.solveWithGuess(rhs, x0);

    if (solver.iterations() >= SLE_SOLVER_MAX_ITERATIONS) {
        std::cout << "Field solver failed!" << std::endl;
        std::cout << solver.error() << std::endl;
    }
}
