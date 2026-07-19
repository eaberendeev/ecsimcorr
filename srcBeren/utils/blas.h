#include <Eigen/Dense>

namespace blas {

static inline bool useOmp(int size) {
    return size >= 1024 * 16;
}

static inline double dot(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    assert(a.rows() == b.rows());

    const double* x = a.data();
    const double* y = b.data();
    const int size = a.rows();

    double res = 0.0;
#pragma omp parallel for reduction(+ : res) if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        res += x[i] * y[i];
    }

    return res;
}

static inline double squaredNorm(const Eigen::VectorXd& a) {
    const double* x = a.data();
    const int size = a.rows();

    double res = 0.0;
#pragma omp parallel for reduction(+ : res) if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        res += x[i] * x[i];
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

// y = alpha * x + beta * y
static inline void axpby(double alpha, const Eigen::VectorXd& x, double beta, Eigen::VectorXd& y) {
    assert(x.rows() == y.rows());
    const int size = x.rows();

#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        y[i] = alpha * x[i] + beta * y[i];
    }
}
// z = x + y
static inline void sum(const Eigen::VectorXd& x, const Eigen::VectorXd& y, Eigen::VectorXd& z) {
    assert(x.rows() == y.rows() && x.rows() == z.rows());
    const int size = x.rows();

#pragma omp parallel for if (useOmp(size))
    for (int i = 0; i < size; ++i) {
        z[i] = x[i] + y[i];
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
