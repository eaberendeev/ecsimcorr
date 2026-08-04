#include "solverSLE.h"

#include "containers.h"
#include "sparse.h"

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
    int restarts = 0;

    while (r.squared() > tol2 && i < maxIters) {
        timer::flatTimer loopTimer("single iteration", i);

        double rho_old = rho;
        rho = r0.dot(r);
        if (abs(rho) < eps2 * r0_sqnorm) {
            // r = rhs - Spmv(x);
            spmv(A, x, r);
            r = rhs - r;
            r0 = r;
            rho = r0_sqnorm = r.squared();
            if (restarts++ == 0)
                i = 0;
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

        double tmp = t.squared();
        if (tmp > 0)
            w = t.dot(s) / tmp;
        else
            w = 0;

        // x += alpha * y + w * z;
        // r = s - w * t;
        timer::flatTimer timerOmp3("OMP section 3", n * sizeof(x[0]) * 6, timer::MeasureUnit::byte);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            x(i) += alpha * y(i) + w * z(i);
            r(i) = s(i) - w * t(i);
        }
        timerOmp3.finish();
        ++i;
    }

    tol_error = sqrt(r.squared() / rhs_sqnorm);
    iters = i;
    return true;
}

bool bicgstab_iteration_mixed_precision(const Operator &A, const Field3d &rhs, Field3d &x, const Field3d &diagonal,
                                        size_t &iters, double &tol_error) {
    RECORD_TIMER;
    size_t itersLower = 0;
    if (envOptions::useMixedPrecision()) {
        timer::commonTimer timer("preparations mixed precision");
        Field3dFp32 rhsLower = rhs;
        Field3dFp32 xLower = x;
        Field3dFp32 diagonalLower = diagonal;
        itersLower = iters;
        double tol_error_lower =
            std::max(static_cast<float>(tol_error), std::numeric_limits<float>::epsilon() * 100.0f);
        {
            ThreadPartitionedSparseMatrix<float> ALower(A);
            timer.finish();
            bicgstab_iteration_impl(ALower, rhsLower, xLower, diagonalLower, itersLower, tol_error_lower);
        }

        blas::copy(xLower.data(), x.data());
    }
    ThreadPartitionedSparseMatrix<double> AFull(A);
    const bool res = bicgstab_iteration_impl(AFull, rhs, x, diagonal, iters, tol_error);
    iters += itersLower;
    return res;
}

template <typename VectorType>
bool bicgstab_iteration(const Operator &A, const VectorType &rhs, VectorType &x, const VectorType &diagonal,
                        size_t &iters, double &tol_error) {
    static const int checkPeriodicity = envOptions::validationPeriodicity();
    static int counter = 0;
    const bool doCheck = counter % checkPeriodicity == 0;
    counter += 1;

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
    if (!std::isfinite(diffNorm) || diffNormNormalized >= 1e-16) {
        std::cerr
            << "Normalized difference between solutions of optimized and reference bicgstab_iteration is too large: "
            << diffNormNormalized << " >= " << 1e-16 << std::endl;
    }

    return res;
}

template bool bicgstab_iteration<Field3d>(const Operator &A, const Field3d &rhs, Field3d &x, const Field3d &diagonal,
                                          size_t &iters, double &tol_error);
