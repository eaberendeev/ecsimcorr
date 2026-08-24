#include "solverSLE.h"

#include <atomic>

#include "containers.h"
#include "sparse.h"
#include "thread_partitioned_matrix.h"

template <typename OperatorType, typename VectorType>
bool bicgstab_iteration_impl(const OperatorType &A, const VectorType &rhs, VectorType &x, const VectorType &diagonal,
                             size_t &iters, double &tol_error) {
    RECORD_TIMER;

    using std::abs;
    using std::sqrt;
    double tol = tol_error;
    int maxIters = iters;
    int n = x.size();

    //    VectorType r = rhs - Spmv(x);
    VectorType r(n);
    spmv(A, x, r);
    r = rhs - r;

    VectorType r0 = r;
    double r0_sqnorm = r0.squared();
    double rhs_sqnorm = rhs.squared();
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

    double rSquared = r.squared();
    while (rSquared > tol2 && i < maxIters) {
        timer::flatTimer loopTimer("single iteration", i);

        double rho_old = rho;
        rho = r0.dot(r);
        if (abs(rho) < eps2 * r0_sqnorm) {
            // Restart: rebuild the true residual and start BiCGSTAB over with
            // the first-iteration state (rho_old=1, alpha=1, w=1, p=0, v=0).
            spmv(A, x, r);
            r = rhs - r;
            r0 = r;
            rho = r0_sqnorm = r.squared();
            rho_old = 1.0;
            alpha = 1.0;
            w = 1.0;
            p.setZero();
            v.setZero();
        }
        double beta = (rho / rho_old) * (alpha / w);

        timer::flatTimer timerOmp1("OMP section 1", n * sizeof(p[0]) * 5, timer::MeasureUnit::byte);
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

        alpha = rho / r0.dot(v);

        // s = r - alpha * v;
        // z = precond.solve(s);   // Применение предобуславливателя

        timer::flatTimer timerOmp2("OMP section 2", n * sizeof(s[0]) * 5, timer::MeasureUnit::byte);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            s(i) = r(i) - alpha * v(i);
            z(i) = s(i) / diagonal(i);
        }
        timerOmp2.finish();
        // t = Spmv(z);
        spmv(A, z, t);

        w = t.normalizedDot(s);

        // x += alpha * y + w * z;
        // r = s - w * t;
        timer::flatTimer timerOmp3("OMP section 3", n * sizeof(x[0]) * 6, timer::MeasureUnit::byte);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            x(i) += alpha * y(i) + w * z(i);
            r(i) = s(i) - w * t(i);
        }
        timerOmp3.finish();
        rSquared = r.squared();
        ++i;
    }

    tol_error = sqrt(rSquared / rhs_sqnorm);
    iters = i;
    return true;
}

bool bicgstab_iteration_mixed_precision(const Operator &A, const Field3d &rhs, Field3d &x, const Field3d &diagonal,
                                        size_t &iters, double &tol_error) {
    RECORD_TIMER;
    if (!envOptions::useMixedPrecision()) {
        const ThreadPartitionedSparseMatrix<double> AFull(A);
        return bicgstab_iteration_impl(AFull, rhs, x, diagonal, iters, tol_error);
    }

    ThreadPartitionedSparseMatrixArray<double, float> matrixArray(A);

    timer::commonTimer timer("preparations mixed precision");
    Field3dFp32 rhsLower = rhs;
    Field3dFp32 xLower = x;
    Field3dFp32 diagonalLower = diagonal;
    size_t itersLower = iters;
    double tol_error_lower = std::max(static_cast<float>(tol_error), std::numeric_limits<float>::epsilon() * 100.0f);
    const ThreadPartitionedSparseMatrixView<float> ALower = matrixArray.get<float>();
    timer.finish();
    bicgstab_iteration_impl(ALower, rhsLower, xLower, diagonalLower, itersLower, tol_error_lower);

    blas::copy(xLower.data(), x.data());
    const ThreadPartitionedSparseMatrixView<double> AFull = matrixArray.get<double>();
    const bool res = bicgstab_iteration_impl(AFull, rhs, x, diagonal, iters, tol_error);

    iters += itersLower;
    return res;
}

bool bicgstab_iteration_mixed_precision_greedy(const Operator &A, const Field3d &rhs, Field3d &x,
                                               const Field3d &diagonal, size_t &iters, double &tol_error) {
    RECORD_TIMER;
    if (!envOptions::useMixedPrecision()) {
        const GreedyThreadPartitionedSparseMatrix<double> AFull(A);
        return bicgstab_iteration_impl(AFull, rhs, x, diagonal, iters, tol_error);
    }

    GreedyThreadPartitionedSparseMatrixArray<double, float> matrixArray(A);
    timer::commonTimer timer("preparations mixed precision");
    Field3dFp32 rhsLower = rhs;
    Field3dFp32 xLower = x;
    Field3dFp32 diagonalLower = diagonal;
    size_t itersLower = iters;
    double tol_error_lower = std::max(static_cast<float>(tol_error), std::numeric_limits<float>::epsilon() * 100.0f);
    const GreedyThreadPartitionedSparseMatrixView<float> ALower = matrixArray.get<float>();
    timer.finish();
    bicgstab_iteration_impl(ALower, rhsLower, xLower, diagonalLower, itersLower, tol_error_lower);

    blas::copy(xLower.data(), x.data());
    const GreedyThreadPartitionedSparseMatrixView<double> AFull = matrixArray.get<double>();
    const bool res = bicgstab_iteration_impl(AFull, rhs, x, diagonal, iters, tol_error);

    iters += itersLower;
    return res;
}

template <typename VectorType>
bool bicgstab_iteration(const Operator &A, const VectorType &rhs, VectorType &x, const VectorType &diagonal,
                        size_t &iters, double &tol_error) {
    static const int checkPeriodicity = envOptions::validationPeriodicity();
    static std::atomic<int> counter{0};
    const bool doCheck = counter.fetch_add(1) % checkPeriodicity == 0;

    if (!doCheck) {
        return bicgstab_iteration_mixed_precision(A, rhs, x, diagonal, iters, tol_error);
    }

    const double desiredTol = tol_error;

    VectorType xRef = x;
    size_t itersRef = iters;
    double tol_error_ref = tol_error;

    const bool resRef = bicgstab_iteration_impl(A, rhs, xRef, diagonal, itersRef, tol_error_ref);
    const bool res = bicgstab_iteration_mixed_precision(A, rhs, x, diagonal, iters, tol_error);

    if (res != resRef) {
        std::cerr
            << "Return value of optimized and reference bicgstab_iteration does not coincide: optimized and reference: "
            << res << " " << resRef << std::endl;
    }
    if (iters != itersRef) {
        std::cerr << "Number of iterations of optimized and reference bicgstab_iteration does not coincide: optimized "
                     "and reference: "
                  << iters << " != " << itersRef << std::endl;
    }
    if (tol_error != tol_error_ref) {
        std::cerr << "Error of optimized and reference bicgstab_iteration does not coincide: optimized and reference: "
                  << tol_error << " - " << tol_error_ref << " = " << tol_error - tol_error_ref
                  << "; desired tolerance: " << desiredTol << std::endl;
    }

    const double xRefNorm = xRef.norm();
    const double diffNorm = (x - xRef).norm();
    const double diffNormNormalized = xRefNorm == 0.0 ? diffNorm : diffNorm / xRefNorm;
    const double maxNormalizedDiff = desiredTol * 2.0;
    if (!std::isfinite(diffNorm) || diffNormNormalized >= maxNormalizedDiff) {
        std::cerr << "Normalized difference between solutions of optimized and reference bicgstab_iteration is too "
                     "large: "
                  << diffNormNormalized << " >= " << maxNormalizedDiff << std::endl;
    }

    return res;
}

template bool bicgstab_iteration<Field3d>(const Operator &A, const Field3d &rhs, Field3d &x, const Field3d &diagonal,
                                          size_t &iters, double &tol_error);
