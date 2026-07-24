#include <omp.h>

#include <Eigen/Dense>

#include "timer.h"

namespace blas {

static inline bool useOmp(int size) {
    return size >= 1024 * 16;
}

// Dot of 2 vectors with run-to run reproducibility in the same vector and omp configuration
static inline double dot(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    assert(a.rows() == b.rows());

    const double* x = a.data();
    const double* y = b.data();
    const int size = a.rows();

    double res = 0.0;

    if (useOmp(size)) {
        const int nthr = omp_get_max_threads();
        double resArray[nthr];
        std::fill_n(resArray, nthr, 0.0);
#pragma omp parallel num_threads(nthr)
        {
            double threadLocalRes = 0;
#pragma omp for
            for (int i = 0; i < size; ++i) {
                threadLocalRes += x[i] * y[i];
            }
            resArray[omp_get_thread_num()] = threadLocalRes;
        }

        for (int i = 0; i < nthr; ++i) {
            res += resArray[i];
        }

    } else {
        for (int i = 0; i < size; ++i) {
            res += x[i] * y[i];
        }
    }

    return res;
}

// Squared norm of vector with run-to run reproducibility in the same vector and omp configuration
static inline double squaredNorm(const Eigen::VectorXd& a) {
    const double* x = a.data();
    const int size = a.rows();

    double res = 0.0;

    if (useOmp(size)) {
        const int nthr = omp_get_max_threads();
        double resArray[nthr];
        std::fill_n(resArray, nthr, 0.0);
#pragma omp parallel num_threads(nthr)
        {
            double threadLocalRes = 0;
#pragma omp for
            for (int i = 0; i < size; ++i) {
                threadLocalRes += x[i] * x[i];
            }
            resArray[omp_get_thread_num()] = threadLocalRes;
        }

        for (int i = 0; i < nthr; ++i) {
            res += resArray[i];
        }

    } else {
        for (int i = 0; i < size; ++i) {
            res += x[i] * x[i];
        }
    }

    return res;
}

static inline void fill(Eigen::VectorXd& a, double value) {
    const int size = a.rows();
#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        a[i] = value;
    }
}

// a *= coeff
static inline void scale(Eigen::VectorXd& a, double coeff) {
    const int size = a.rows();
#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        a[i] *= coeff;
    }
}

// y = alpha * x + beta * y
static inline void axpby(double alpha, const Eigen::VectorXd& x, double beta, Eigen::VectorXd& y) {
    assert(x.rows() == y.rows());
    const int size = x.rows();

#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        y[i] = alpha * x[i] + beta * y[i];
    }
}
// z = alpha * x + beta * y
static inline void sum(const double alpha, const Eigen::VectorXd& x, const double beta, const Eigen::VectorXd& y,
                       Eigen::VectorXd& z) {
    RECORD_TIMER;
    assert(x.rows() == y.rows() && x.rows() == z.rows());
    const int size = x.rows();

    if (alpha == 1.0 && beta == 1.0) {
#pragma omp parallel for if (useOmp(size))
        for (int i = 0; i < size; ++i) {
            z[i] = x[i] + y[i];
        }
    } else {
#pragma omp parallel for if (useOmp(size))
        for (int i = 0; i < size; ++i) {
            z[i] = alpha * x[i] + beta * y[i];
        }
    }
}

// y = x
static inline void copy(const Eigen::VectorXd& x, Eigen::VectorXd& y) {
    assert(x.rows() == y.rows());
    const int size = x.rows();

#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        y[i] = x[i];
    }
}

}   // namespace blas
