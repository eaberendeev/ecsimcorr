#include <Eigen/Dense>

static inline bool useOmp(int size) {
    return size >= 1024 * 16;
}

static inline double ddot(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
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
