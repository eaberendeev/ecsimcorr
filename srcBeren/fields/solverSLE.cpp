#include "solverSLE.h"

#include "utils/containers.h"
#include "utils/sparse.h"
#include "utils/thread_partitioned_matrix.h"

template <typename OperatorType, typename VectorType>
bool bicgstab_iteration_impl(const OperatorType& A, const VectorType& rhs, VectorType& x, const VectorType& diagonal,
                             size_t& iters, double& tol_error, double& divergenceNorm) {
    RECORD_TIMER;

    using std::abs;
    using std::sqrt;
    const double tol = tol_error;
    const int maxIters = iters;
    const int n = x.size();

    using base_t = std::remove_reference_t<decltype(x[0])>;

    // r = A*x;
    VectorType r(n);
    spmv(A, x, r);
    // r = rhs - r;
    blas::axpby(base_t{1}, rhs, base_t{-1}, r);

    VectorType r0(n);
    blas::copy(r, r0);
    double r0_sqnorm = blas::squaredNorm(r0);
    const double rhs_sqnorm = blas::squaredNorm(rhs);
    if (rhs_sqnorm == 0) {
        blas::fill(x, base_t{0});
        iters = 0;
        tol_error = 0.0;
        divergenceNorm = 0.0;
        return true;
    }
    double rho = 1;
    double alpha = 1;
    double w = 1;
    VectorType v(n);
    VectorType p(n);
    blas::fill(v, base_t{0});
    blas::fill(p, base_t{0});
    VectorType y(n), z(n);
    VectorType s(n), t(n);
    const double tol2 = tol * tol * rhs_sqnorm;
    const double eps2 = Eigen::NumTraits<double>::epsilon() * Eigen::NumTraits<double>::epsilon();
    int i = 0;

    double rSquared = blas::squaredNorm(r);
    while (std::isfinite(rSquared) && rSquared > tol2 && i < maxIters) {
        timer::FlatTimer loopTimer("single iteration", i);

        double rho_old = rho;
        rho = blas::dot(r0, r);
        if (abs(rho) < eps2 * r0_sqnorm) {
            // Restart: rebuild the true residual and start BiCGSTAB over with
            // the first-iteration state (rho_old=1, alpha=1, w=1, p=0, v=0).
            spmv(A, x, r);
            r = rhs - r;
            r0 = r;
            rho = r0_sqnorm = blas::squaredNorm(r);
            rho_old = 1.0;
            alpha = 1.0;
            w = 1.0;
            blas::fill(p, base_t{0});
            blas::fill(v, base_t{0});
        }
        if (!std::isfinite(rho) || rho == 0.0 || rho_old == 0.0 || w == 0.0)
            break;
        const double beta = (rho / rho_old) * (alpha / w);
        if (!std::isfinite(beta))
            break;

        timer::FlatTimer timerOmp1("OMP section 1", n * sizeof(p[0]) * 4, timer::MeasureUnit::byte);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            p(i) = r(i) + beta * (p(i) - w * v(i));
            y(i) = p(i) / diagonal(i);
        }
        timerOmp1.finish();

        // p = r + beta * (p - w * v);
        // y = precond.solve(p);   // Применение предобуславливателя
        // v = Spmv(y);
        spmv(A, y, v);

        const double denominator = blas::dot(r0, v);
        if (!std::isfinite(denominator) || denominator == 0.0)
            break;
        alpha = rho / denominator;
        if (!std::isfinite(alpha))
            break;

        // s = r - alpha * v;
        // z = precond.solve(s);   // Применение предобуславливателя

        timer::FlatTimer timerOmp2("OMP section 2", n * sizeof(s[0]) * 4, timer::MeasureUnit::byte);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            s(i) = r(i) - alpha * v(i);
            z(i) = s(i) / diagonal(i);
        }
        timerOmp2.finish();
        // An exact alpha step needs neither an omega division nor another SpMV.
        const double sSquared = blas::squaredNorm(s);
        if (std::isfinite(sSquared) && sSquared <= tol2) {
#pragma omp parallel for simd
            for (int j = 0; j < n; ++j) x(j) += alpha * y(j);
            ++i;
            break;
        }
        if (!std::isfinite(sSquared))
            break;
        spmv(A, z, t);

        w = blas::normalizedDot(t, s);
        if (!std::isfinite(w) || w == 0.0)
            break;

        // x += alpha * y + w * z;
        // r = s - w * t;
        timer::FlatTimer timerOmp3("OMP section 3", n * sizeof(x[0]) * 5, timer::MeasureUnit::byte);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            x(i) += alpha * y(i) + w * z(i);
            r(i) = s(i) - w * t(i);
        }
        timerOmp3.finish();
        rSquared = blas::squaredNorm(r);
        ++i;
    }

    // A recurrence can underestimate residuals; success is decided from A*x-b.
    spmv(A, x, r);
    r = rhs - r;
    rSquared = blas::squaredNorm(r);
    tol_error = sqrt(rSquared / rhs_sqnorm);
    divergenceNorm = sqrt(rSquared);
    iters = i;
    return std::isfinite(rSquared) && rSquared <= tol2;
}

// Relative accuracy of each float correction; measured on ECSIM predictor matrices (c dt/dx = 3..5), looser
// corrections need more refinement steps but fewer float iterations in total.
constexpr double correctionTolerance = 1e-2;

// Mixed-precision iterative refinement: the residual is formed with the double matrix, each correction solves
// A y = r / |r| with the float matrix and vectors. Corrections get at most half of the iteration budget; one that
// does not halve its float residual is dropped, and one that does not halve the double residual (float accuracy
// reached) ends refinement. Double BiCGSTAB finishes with the rest of the budget.
template <typename DoubleMatrix, typename MatrixArray>
bool refine_in_float(const Operator& A, const Field3d& rhs, Field3d& x, const Field3d& diagonal, size_t& iters,
                     double& tol_error, double& divergenceNorm) {
    RECORD_TIMER;
    if (!envOptions::useMixedPrecision()) {
        const DoubleMatrix full(A);
        return bicgstab_iteration_impl(full, rhs, x, diagonal, iters, tol_error, divergenceNorm);
    }
    const MatrixArray matrices(A);
    const typename MatrixArray::template view_t<double> full = matrices.template get<double>();
    const double rhs_norm = std::sqrt(blas::squaredNorm(rhs));
    if (rhs_norm == 0)
        return bicgstab_iteration_impl(full, rhs, x, diagonal, iters, tol_error, divergenceNorm);

    const typename MatrixArray::template view_t<float> lower = matrices.template get<float>();
    const Field3dFp32 diagonal_lower = diagonal;
    const double target = tol_error * rhs_norm;
    const int n = x.size();
    Field3d residual(n);
    Field3dFp32 correction_rhs(n);
    Field3dFp32 correction(n);
    const size_t float_budget = iters / 2;
    size_t used = 0;
    double previous = std::numeric_limits<double>::infinity();
    while (true) {
        timer::commonTimer timerLoop("refinement iteration");
        spmv(full, x, residual);
        blas::axpby(1.0, rhs, -1.0, residual);
        const double norm = std::sqrt(blas::squaredNorm(residual));
        if (norm <= target) {
            iters = used;
            tol_error = norm / rhs_norm;
            divergenceNorm = norm;
            return true;
        }
        // A nonfinite residual also ends refinement; double BiCGSTAB reports it.
        if (used >= float_budget || !std::isfinite(norm) || norm > 0.5 * previous)
            break;
        previous = norm;
        // A unit right-hand side keeps the float dot products in range as the residual shrinks.
        residual *= 1.0 / norm;
        blas::copy(residual, correction_rhs);
        correction.setZero();
        size_t correction_iters = float_budget - used;
        double correction_tol = correctionTolerance;
        bicgstab_iteration_impl(lower, correction_rhs, correction, diagonal_lower, correction_iters, correction_tol,
                                divergenceNorm);
        used += correction_iters;
        // A breakdown can leave a correction that does not reduce its residual; x keeps the last good iterate.
        if (!(correction_tol <= 0.5))
            break;
#pragma omp parallel for simd
        for (int i = 0; i < n; ++i) x(i) += norm * correction(i);
    }
    size_t rest = iters - std::min(iters, used);
    const bool success = bicgstab_iteration_impl(full, rhs, x, diagonal, rest, tol_error, divergenceNorm);
    iters = used + rest;
    return success;
}

bool bicgstab_iteration_mixed_precision(const Operator& A, const Field3d& rhs, Field3d& x, const Field3d& diagonal,
                                        size_t& iters, double& tol_error, double& divergenceNorm) {
    return refine_in_float<ThreadPartitionedSparseMatrix<double>, ThreadPartitionedSparseMatrixArray<double, float>>(
        A, rhs, x, diagonal, iters, tol_error, divergenceNorm);
}

bool bicgstab_iteration_mixed_precision_greedy(const Operator& A, const Field3d& rhs, Field3d& x,
                                               const Field3d& diagonal, size_t& iters, double& tol_error,
                                               double& divergenceNorm) {
    return refine_in_float<GreedyThreadPartitionedSparseMatrix<double>,
                           GreedyThreadPartitionedSparseMatrixArray<double, float>>(A, rhs, x, diagonal, iters,
                                                                                    tol_error, divergenceNorm);
}

template <typename VectorType>
bool bicgstab_iteration(const Operator& A, const VectorType& rhs, VectorType& x, const VectorType& diagonal,
                        size_t& iters, double& tol_error, double& divergenceNorm) {
    static_assert(std::is_same_v<VectorType, Field3d>);   // implemented only for Field3d yet
    return bicgstab_iteration_mixed_precision(A, rhs, x, diagonal, iters, tol_error, divergenceNorm);
}

template bool bicgstab_iteration<Field3d>(const Operator& A, const Field3d& rhs, Field3d& x, const Field3d& diagonal,
                                          size_t& iters, double& tol_error, double& divergenceNorm);
