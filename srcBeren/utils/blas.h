#include <omp.h>

#include <Eigen/Dense>

#include "timer.h"

namespace blas {

static inline bool useOmp(int size) {
    return size >= 1024 * 16;
}

// Dot of 2 vectors with run-to run reproducibility in the same vector and omp configuration
template <typename T, typename inner_t = double>
static inline T dot(const Eigen::VectorX<T>& a, const Eigen::VectorX<T>& b) {
    assert(a.rows() == b.rows());

    const T* x = a.data();
    const T* y = b.data();
    const int size = a.rows();

    inner_t res = 0.0;

    if (useOmp(size)) {
        const int nthr = omp_get_max_threads();
        inner_t resArray[nthr];
        std::fill_n(resArray, nthr, 0.0);
#pragma omp parallel num_threads(nthr)
        {
            double threadLocalRes = 0;
#pragma omp for
            for (int i = 0; i < size; ++i) {
                threadLocalRes += static_cast<inner_t>(x[i]) * static_cast<inner_t>(y[i]);
            }
            resArray[omp_get_thread_num()] = threadLocalRes;
        }

        for (int i = 0; i < nthr; ++i) {
            res += resArray[i];
        }

    } else {
        for (int i = 0; i < size; ++i) {
            res += static_cast<inner_t>(x[i]) * static_cast<inner_t>(y[i]);
        }
    }

    return static_cast<T>(res);
}

// Squared norm of vector with run-to run reproducibility in the same vector and omp configuration
template <typename T, typename inner_t = double>
static inline T squaredNorm(const Eigen::VectorX<T>& a) {
    const T* x = a.data();
    const int size = a.rows();

    inner_t res = 0.0;

    if (useOmp(size)) {
        const int nthr = omp_get_max_threads();
        inner_t resArray[nthr];
        std::fill_n(resArray, nthr, 0.0);
#pragma omp parallel num_threads(nthr)
        {
            inner_t threadLocalRes = 0;
#pragma omp for
            for (int i = 0; i < size; ++i) {
                threadLocalRes += static_cast<inner_t>(x[i]) * static_cast<inner_t>(x[i]);
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

template <typename T>
static inline void fill(Eigen::VectorX<T>& a, T value) {
    const int size = a.rows();
#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        a[i] = value;
    }
}

// a *= coeff
template <typename T>
static inline void scale(Eigen::VectorX<T>& a, T coeff) {
    const int size = a.rows();
#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        a[i] *= coeff;
    }
}

// y = alpha * x + beta * y
template <typename T>
static inline void axpby(T alpha, const Eigen::VectorX<T>& x, T beta, Eigen::VectorX<T>& y) {
    assert(x.rows() == y.rows());
    const int size = x.rows();

#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        y[i] = alpha * x[i] + beta * y[i];
    }
}
// z = alpha * x + beta * y
template <typename T>
static inline void sum(const T alpha, const Eigen::VectorX<T>& x, const T beta, const Eigen::VectorX<T>& y,
                       Eigen::VectorX<T>& z) {
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
template <typename T1, typename T2>
static inline void copy(const Eigen::VectorX<T1>& x, Eigen::VectorX<T2>& y) {
    assert(x.rows() == y.rows());
    const int size = x.rows();

#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        y[i] = static_cast<T2>(x[i]);
    }
}

}   // namespace blas
