#pragma once

#include <cmath>
#include <iterator>
#include <stdexcept>
#include <type_traits>

#include "utils/bmatrix.h"
#include "include/config.h"
#include "utils/containers.h"
#include "runtime/env_options.h"
#include "utils/sparse.h"
#include "runtime/profiler.h"
#include "utils/util.h"

#define DEFAULT_MAX_ITERATIONS 1000
#define DEFAULT_TOLERANCE      1.e-9

template <typename VectorType>
bool bicgstab_iteration(const Operator& A, const VectorType& rhs, VectorType& x, const VectorType& diagonal,
                        size_t& iters, double& tol_error, double& divergenceNorm);

template <typename VectorType>
class BicgstabSolverBase {
   public:
    BicgstabSolverBase() : max_iterations(DEFAULT_MAX_ITERATIONS), m_tolerance(DEFAULT_TOLERANCE), m_success(false) {
    }

    void setTolerance(double tolerance) {
        if (!std::isfinite(tolerance) || tolerance <= 0.0)
            throw std::invalid_argument("solver tolerance must be finite and positive");
        m_tolerance = tolerance;
    }
    void setMaxIterations(size_t max_iters) {
        max_iterations = max_iters;
    }
    bool info() const {
        return m_success;
    }
    double error() const {
        return m_error;
    }
    size_t iterations() const {
        return m_iterations;
    }

    virtual VectorType solveWithGuess(const VectorType& rhs, const VectorType& x0) = 0;

   protected:
    virtual ~BicgstabSolverBase() = default;

    VectorType m_diagonal;
    size_t max_iterations = -1;
    size_t m_iterations = -1;
    double m_tolerance = -1.0;
    double m_error = -1.0;
    double divergenceNorm = -1.0;
    bool m_success = false;
};

template <typename VectorType>
class BicgstabSolver : public BicgstabSolverBase<VectorType> {
   public:
    using Base = BicgstabSolverBase<VectorType>;
    using Base::divergenceNorm;
    using Base::m_diagonal;
    using Base::m_error;
    using Base::m_iterations;
    using Base::m_success;
    using Base::m_tolerance;
    using Base::max_iterations;
    BicgstabSolver(const Operator& A) : m_A(A) {
        initializePreconditioner(A.rows());
    }
    BicgstabSolver(const Operator& A, const VectorType& diagonal) : m_A(A) {
        computeDiagonalPreconditioner(diagonal);
    }
    VectorType solveWithGuess(const VectorType& rhs, const VectorType& x0) {
        VectorType x = x0;
        m_iterations = max_iterations;
        m_error = m_tolerance;

        m_success = bicgstab_iteration(m_A, rhs, x, m_diagonal, m_iterations, m_error, divergenceNorm);
        return x;
    }

   private:
    void computeDiagonalPreconditioner(const VectorType& diag) {
        const int rows = m_A.rows();
        if (std::ssize(diag) != rows)
            throw std::invalid_argument("preconditioner dimension differs from matrix");
        RECORD_TIMER_PARAMS(sizeof(diag[0]) * std::ssize(diag), timer::MeasureUnit::byte);
        m_diagonal.resize(rows);
        bool valid = true;
#pragma omp parallel for schedule(static, 16) reduction(&& : valid)
        for (int i = 0; i < rows; ++i) {
            const double value = diag[i];
            m_diagonal[i] = value;
            if (!std::isfinite(value) || value == 0.0)
                valid = false;
        }
        if (!valid)
            throw std::invalid_argument("diagonal preconditioner must be finite and nonzero");
    }
    void initializePreconditioner(int rows) {
        RECORD_TIMER_PARAMS(sizeof(m_diagonal[0]) * rows, timer::MeasureUnit::byte);
        m_diagonal.resize(rows);
#pragma omp parallel for schedule(static, 16)
        for (int i = 0; i < std::ssize(m_diagonal); ++i) {
            m_diagonal[i] = 1.0;
        }
    }
    const Operator& m_A;
};

template <typename SolverType, typename VectorType>
void solve_linear_system_impl(SolverType& solver, const VectorType& rhs, VectorType& x, const VectorType& x0,
                              double tol) {
    solver.setTolerance(tol);
    solver.setMaxIterations(SLE_SOLVER_MAX_ITERATIONS);

    x = solver.solveWithGuess(rhs, x0);

    if (!solver.info() || !std::isfinite(solver.error())) {
        throw std::runtime_error("field solver failed after " + std::to_string(solver.iterations()) +
                                 " iterations; relative residual=" + std::to_string(solver.error()));
    }
}

template <typename SolverType, typename VectorType>
double solve_linear_system(const Operator& A, const VectorType& rhs, VectorType& x, const VectorType& x0,
                           double tol = SLE_SOLVER_TOLERANCE) {
    RECORD_TIMER;
    SolverType solver(A);
    solve_linear_system_impl(solver, rhs, x, x0, tol);
    return solver.divergenceNorm;
}
