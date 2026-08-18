#include <omp.h>

#include <Eigen/Dense>

#include "timer.h"

namespace blas {

static inline bool useOmp(int size) {
    return size >= 1024 * 16;
}

// Dot of 2 vectors with run-to run reproducibility in the same vector and omp configuration
template <typename T, typename inner_t = double>
static inline T dot(const Eigen::VectorX<T>& x, const Eigen::VectorX<T>& y) {
    assert(x.rows() == y.rows());
    const int size = x.rows();

    inner_t res = inner_t{0};

    if (useOmp(size)) {
        const int nthr = omp_get_max_threads();
        // Per-thread partial sums; the final summation runs in a fixed thread
        // order below, which makes the result reproducible for a given vector
        // and thread count (no reduction tree, no uninitialized entries).
        inner_t partialRes[nthr];
        int usedThreads = std::numeric_limits<int>::max();
#pragma omp parallel num_threads(nthr)
        {
#pragma omp master
            usedThreads = omp_get_num_threads();
            inner_t threadLocalRes = 0;
#pragma omp for simd
            for (int i = 0; i < size; ++i) {
                threadLocalRes += static_cast<inner_t>(x[i]) * static_cast<inner_t>(y[i]);
            }
            partialRes[omp_get_thread_num()] = threadLocalRes;
        }

        for (int i = 0; i < usedThreads; ++i) {
            res += partialRes[i];
        }

    } else {
#pragma omp simd
        for (int i = 0; i < size; ++i) {
            res += static_cast<inner_t>(x[i]) * static_cast<inner_t>(y[i]);
        }
    }

    return static_cast<T>(res);
}

// Squared norm of vector with run-to run reproducibility in the same vector and omp configuration
template <typename T, typename inner_t = double>
static inline T squaredNorm(const Eigen::VectorX<T>& x) {
    const int size = x.rows();

    inner_t res = inner_t{0};

    if (useOmp(size)) {
        const int nthr = omp_get_max_threads();
        inner_t partialRes[nthr];
        int usedThreads = std::numeric_limits<int>::max();
#pragma omp parallel num_threads(nthr)
        {
#pragma omp master
            usedThreads = omp_get_num_threads();
            inner_t threadLocalRes = 0;
#pragma omp for simd
            for (int i = 0; i < size; ++i) {
                threadLocalRes += static_cast<inner_t>(x[i]) * static_cast<inner_t>(x[i]);
            }
            partialRes[omp_get_thread_num()] = threadLocalRes;
        }

        for (int i = 0; i < usedThreads; ++i) {
            res += partialRes[i];
        }
    } else {
#pragma omp simd
        for (int i = 0; i < size; ++i) {
            res += static_cast<inner_t>(x[i]) * static_cast<inner_t>(x[i]);
        }
    }

    return static_cast<T>(res);
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
template <typename T, typename inner_t = double>
static inline void axpby(T alpha, const Eigen::VectorX<T>& x, T beta, Eigen::VectorX<T>& y) {
    assert(x.rows() == y.rows());
    const int size = x.rows();

    if (beta == T{1}) {
        const inner_t alpha2 = static_cast<inner_t>(alpha);
#pragma omp parallel for simd if (useOmp(size))
        for (int i = 0; i < size; ++i) {
            const inner_t x2 = static_cast<inner_t>(x[i]);
            const inner_t y2 = static_cast<inner_t>(y[i]);
            y[i] = static_cast<T>(alpha2 * x2 + y2);
        }
    } else {
        const inner_t alpha2 = static_cast<inner_t>(alpha);
        const inner_t beta2 = static_cast<inner_t>(beta);
#pragma omp parallel for simd if (useOmp(size))
        for (int i = 0; i < size; ++i) {
            const inner_t x2 = static_cast<inner_t>(x[i]);
            const inner_t y2 = static_cast<inner_t>(y[i]);
            y[i] = static_cast<T>(alpha2 * x2 + beta2 * y2);
        }
    }
}
// z = alpha * x + beta * y
template <typename T, typename inner_t = double>
static inline void sum(const T alpha, const Eigen::VectorX<T>& x, const T beta, const Eigen::VectorX<T>& y,
                       Eigen::VectorX<T>& z) {
    RECORD_TIMER;
    assert(x.rows() == y.rows() && x.rows() == z.rows());
    const int size = x.rows();

    if (alpha == T{1} && beta == T{1}) {
#pragma omp parallel for simd if (useOmp(size))
        for (int i = 0; i < size; ++i) {
            // don't use inner T, since computation units in CPU use more bits in internal operations
            z[i] = x[i] + y[i];
        }
    } else {
        const inner_t alpha2 = static_cast<inner_t>(alpha);
        const inner_t beta2 = static_cast<inner_t>(beta);
#pragma omp parallel for simd if (useOmp(size))
        for (int i = 0; i < size; ++i) {
            const inner_t x2 = static_cast<inner_t>(x[i]);
            const inner_t y2 = static_cast<inner_t>(y[i]);
            z[i] = static_cast<T>(alpha2 * x2 + beta2 * y2);
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
