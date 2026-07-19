#pragma once

#define USE_AMGCL

#ifdef USE_AMGCL
#include <amgcl/amg.hpp>
#include <amgcl/backend/eigen.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/preconditioner/dummy.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/iluk.hpp>
#include <amgcl/relaxation/ilup.hpp>
#include <amgcl/relaxation/ilut.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/spai1.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/idrs.hpp>
#include <amgcl/solver/lgmres.hpp>
#endif   // USE_AMGCL

#include <type_traits>

#include "bmatrix.h"
#include "config.h"
#include "containers.h"
#include "timer.h"
#include "util.h"
#define DEFAULT_MAX_ITERATIONS 1000
#define DEFAULT_TOLERANCE      1.e-9

struct PartitionedMatrix {
    void init(const Operator &A) {
        RECORD_TIMER;
        const int numThreads = omp_get_num_threads();
        const int tid = omp_get_thread_num();

        threadOwner = tid;
        rowsGlobal = A.rows();

        const double *val = A.valuePtr();
        const int *inner = A.innerIndexPtr();
        const int *outer = A.outerIndexPtr();

        rowStart = findPost(A, tid, numThreads);
        rowEnd = findPost(A, tid + 1, numThreads);

        outerIndexes.resize(rowEnd - rowStart + 1);
        innerIndexes.resize(outer[rowEnd] - outer[rowStart]);
        data.resize(outer[rowEnd] - outer[rowStart]);

        for (int i = rowStart; i < rowEnd + 1; ++i) {
            outerIndexes[i - rowStart] = outer[i];
        }

        const int dataOffset = outer[rowStart];
        for (int i = rowStart; i < rowEnd; ++i) {
            for (int j = outer[i]; j < outer[i + 1]; ++j) {
                data[j - dataOffset] = val[j];
                innerIndexes[j - dataOffset] = inner[j];
            }
        }
    }

    static int findPost(const Operator &A, int blockId, int numBlocks) {
        const int nnz = A.nonZeros();
        const int bestPos = static_cast<int64_t>(nnz) * blockId / numBlocks;

        const int *outer = A.outerIndexPtr();

        for (int i = 0; i < A.rows(); ++i) {
            if (outer[i] <= bestPos && bestPos <= outer[i + 1]) {
                return i;
            }
        }
        return std::numeric_limits<int>::min();
    }

    int rowsGlobal;
    int rowStart;
    int rowEnd;
    int threadOwner;

    std::vector<int> innerIndexes;
    std::vector<int> outerIndexes;
    std::vector<double> data;
};

struct ThreadPartitionedSparseMatrix {
    ThreadPartitionedSparseMatrix(const Operator &A) : nthr(omp_get_max_threads()), nnz(A.nonZeros()), matrices(nthr) {
        RECORD_TIMER;
#pragma omp parallel
        {
            matrices[omp_get_thread_num()].init(A);
        }

        // int minNnz = 1'000'000'000;
        // int maxNnz = -1;
        // for (int i = 0; i < nthr; ++i) {
        //     const int nnz = matrices[i].outerIndexes.back() - matrices[i].outerIndexes.front();
        //     std::cout << "nnz: " << nnz << std::endl;
        //     minNnz = std::min(minNnz, nnz);
        //     maxNnz = std::max(maxNnz, nnz);
        // }
        // std::cout << "#### min nnz: " << minNnz << std::endl;
        // std::cout << "#### max nnz: " << maxNnz << std::endl;
    }

    template <typename VectorType>
    inline void spmv(const VectorType &v, VectorType &res) {
        RECORD_TIMER_PARAMS(nnz);
#pragma omp parallel
        {
            assert(nthr == omp_get_num_threads());
            const PartitionedMatrix &localMat = matrices[omp_get_thread_num()];

            timer::commonTimer timerOMP("OMP section", localMat.outerIndexes.back() - localMat.outerIndexes.front());

            for (int i = localMat.rowStart; i < localMat.rowEnd; ++i) {
                double sum = 0.0;
                const int localRow = i - localMat.rowStart;
                const int dataOffset = localMat.outerIndexes[0];
                const int start = localMat.outerIndexes[localRow] - dataOffset;
                const int end = localMat.outerIndexes[localRow + 1] - dataOffset;
#pragma omp simd
                for (int j = start; j < end; ++j) {
                    sum += localMat.data[j] * v[localMat.innerIndexes[j]];
                }
                res[i] = sum;
            }
        }
    }

    const int nthr;
    const int nnz;   // for timings only
    std::vector<PartitionedMatrix> matrices;
};

template <typename VectorType>
inline void spmv(const Operator &A, const VectorType &v, VectorType &res) {
    RECORD_TIMER_PARAMS(A.nonZeros());
    int rows = A.rows();

    const double *val = A.valuePtr();
    const int *inner = A.innerIndexPtr();
    const int *outer = A.outerIndexPtr();
#pragma omp parallel
    {
        timer::commonTimer timerOmp("OMP section");
#pragma omp for
        for (int i = 0; i < rows; ++i) {
            double sum = 0;
#pragma omp simd
            for (int j = outer[i]; j < outer[i + 1]; ++j) {
                sum += val[j] * v[inner[j]];
            }
            res[i] = sum;
        }
    }
}

template <typename VectorType>
bool bicgstab_iteration(const Operator &A, const VectorType &rhs, VectorType &x, const VectorType &diagonal,
                        size_t &iters, double &tol_error) {
    RECORD_TIMER;

    ThreadPartitionedSparseMatrix partitionedA(A);

    using std::abs;
    using std::sqrt;
    double tol = tol_error;
    int maxIters = iters;
    int n = x.size();

    //    VectorType r = rhs - Spmv(x);
    VectorType r(n);
    partitionedA.spmv(x, r);
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
        timer::flatTimer loopTimer("single iteration");

        double rho_old = rho;
        rho = r0.dot(r);
        if (abs(rho) < eps2 * r0_sqnorm) {
            // r = rhs - Spmv(x);
            partitionedA.spmv(x, r);
            r = rhs - r;
            r0 = r;
            rho = r0_sqnorm = r.squared();
            if (restarts++ == 0)
                i = 0;
        }
        double beta = (rho / rho_old) * (alpha / w);

        timer::flatTimer timerOmp1("OMP section 1", n);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            p(i) = r(i) + beta * (p(i) - w * v(i));
            y(i) = p(i) / diagonal(i);
        }
        timerOmp1.finish();

        // p = r + beta * (p - w * v);
        // y = precond.solve(p);   // Применение предобуславливателя
        // v = Spmv(y);
        partitionedA.spmv(y, v);
        alpha = rho / r0.dot(v);

        // s = r - alpha * v;
        // z = precond.solve(s);   // Применение предобуславливателя

        timer::flatTimer timerOmp2("OMP section 2", n);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            s(i) = r(i) - alpha * v(i);
            z(i) = s(i) / diagonal(i);
        }
        timerOmp2.finish();
        // t = Spmv(z);
        partitionedA.spmv(z, t);

        const double tmp = t.squared();
        if (tmp > 0)
            w = t.dot(s) / tmp;
        else
            w = 0;

        // x += alpha * y + w * z;
        // r = s - w * t;
        timer::flatTimer timerOmp3("OMP section 3", n);
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

template <typename VectorType>
class BicgstabSolverBase {
   public:
    BicgstabSolverBase() : max_iterations(DEFAULT_MAX_ITERATIONS), m_tolerance(DEFAULT_TOLERANCE), m_success(false) {
    }

    void setTolerance(double tolerance) {
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

    virtual VectorType solveWithGuess(const VectorType &rhs, const VectorType &x0) = 0;

   protected:
    virtual ~BicgstabSolverBase() = default;

    VectorType m_diagonal;
    size_t max_iterations;
    size_t m_iterations;
    double m_tolerance;
    double m_error;
    bool m_success;
};

template <typename VectorType>
class BicgstabSolver : public BicgstabSolverBase<VectorType> {
   public:
    using Base = BicgstabSolverBase<VectorType>;
    using Base::m_diagonal;
    using Base::m_error;
    using Base::m_iterations;
    using Base::m_success;
    using Base::m_tolerance;
    using Base::max_iterations;
    BicgstabSolver(const Operator &A) : m_A(A) {
        initializePreconditioner(A.rows());
    }
    BicgstabSolver(const Operator &A, const VectorType &diagonal) : m_A(A) {
        computeDiagonalPreconditioner(diagonal);
    }
    VectorType solveWithGuess(const VectorType &rhs, const VectorType &x0) {
        VectorType x = x0;
        m_iterations = max_iterations;
        m_error = m_tolerance;

        m_success = bicgstab_iteration(m_A, rhs, x, m_diagonal, m_iterations, m_error);
        return x;
    }

   private:
    void computeDiagonalPreconditioner(const Eigen::VectorXd &diag) {
        RECORD_TIMER;

#pragma omp parallel for schedule(static, 16)
        for (int i = 0; i < std::ssize(m_diagonal); i++) {
            m_diagonal[i] = diag[i];
        }
    }
    void initializePreconditioner(int rows) {
        RECORD_TIMER;
        m_diagonal.resize(rows);
#pragma omp parallel for schedule(static, 16)
        for (int i = 0; i < std::ssize(m_diagonal); ++i) {
            m_diagonal[i] = 1.0;
        }
    }
    const Operator &m_A;
};

template <typename SolverType, typename VectorType>
void solve_linear_system_impl(SolverType &solver, const VectorType &rhs, VectorType &x, const VectorType &x0) {
    solver.setTolerance(SLE_SOLVER_TOLERANCE);
    solver.setMaxIterations(SLE_SOLVER_MAX_ITERATIONS);

    x = solver.solveWithGuess(rhs, x0);

    if (solver.iterations() >= SLE_SOLVER_MAX_ITERATIONS) {
        std::cout << "Field solver failed!" << std::endl;
        std::cout << solver.error() << std::endl;
    }
}

template <typename SolverType, typename VectorType>
void solve_linear_system(const Operator &A, const VectorType &rhs, VectorType &x, const VectorType &x0) {
    RECORD_TIMER;
    SolverType solver(A);
    solve_linear_system_impl(solver, rhs, x, x0);
}

template <typename VectorType>
void solve_linear_system(const Operator &A, const VectorType &diagonal, const VectorType &rhs, VectorType &x,
                         const VectorType &x0) {
    BicgstabSolver solver(A, diagonal);
    solve_linear_system_impl(solver, rhs, x, x0);
}

#ifdef USE_AMGCL

using namespace amgcl;

template <typename MatrixType>
void solve_amgcl(const MatrixType &A, const Field &rhs, Field &x, const Field &x0) {
    typedef backend::eigen<double> Backend;

    typedef amgcl::make_solver<
        // Use AMG as preconditioner:
        amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>,

        // relaxation::as_preconditioner<Backend,
        //                               amgcl::relaxation::damped_jacobi>,
        //  preconditioner::dummy<Backend>,
        //  And BiCGStab as iterative solver:
        amgcl::solver::gmres<Backend> >
        Solver;

    Solver::params prm;
    prm.solver.tol = 1e-10;
    Solver solve(A, prm);

    int iters;
    double error;
    x = x0;
    std::tie(iters, error) = solve(rhs, x);
    std::cout << iters << " " << error << "\n";

    // auto prm = precond.amgcl_params();
}
#endif   // USE_AMGCL

typedef Eigen::SparseMatrix<float, Eigen::RowMajor> OperatorM;
typedef Eigen::GMRES<Eigen::SparseMatrix<float, Eigen::RowMajor> > gmres_m;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<float, Eigen::RowMajor> > bicgstab_m;
typedef Eigen::VectorXf VectorXf;

static Eigen::SparseMatrix<float, Eigen::RowMajor> convertSparseMatrixDoubleToFloat(
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &matDouble) {
    // Создаём матрицу float с теми же размерами
    Eigen::SparseMatrix<float, Eigen::RowMajor> matFloat(matDouble.rows(), matDouble.cols());

    // Копируем структуру матрицы (индексы строк и столбцов)
    matFloat.resizeNonZeros(matDouble.nonZeros());
    std::copy(matDouble.outerIndexPtr(), matDouble.outerIndexPtr() + matDouble.outerSize() + 1,
              matFloat.outerIndexPtr());
    std::copy(matDouble.innerIndexPtr(), matDouble.innerIndexPtr() + matDouble.nonZeros(), matFloat.innerIndexPtr());

    // Преобразуем значения из double в float
    const double *valuesDouble = matDouble.valuePtr();
    float *valuesFloat = matFloat.valuePtr();
#pragma omp parallel for num_threads(8)
    for (int i = 0; i < matDouble.nonZeros(); ++i) {
        valuesFloat[i] = static_cast<float>(valuesDouble[i]);
    }

    return matFloat;
}

inline void solve_linear_system_mix(const Operator &A, const Field &rhs, Field &x, const Field &x0) {
    OperatorM AM = convertSparseMatrixDoubleToFloat(A);
    bicgstab_m solverM(AM);
    VectorXf rhsM(rhs.size());
    VectorXf xM(rhs.size());
    VectorXf x0M(rhs.size());
    Field xp(rhs.size());
    for (int i = 0; i < rhs.size(); i++) {
        rhsM[i] = static_cast<float>(rhs[i]);
        x0M[i] = static_cast<float>(x0[i]);
    }
    solve_linear_system_impl(solverM, rhsM, xM, x0M);
    bicgstab solver(A);
    for (int i = 0; i < rhs.size(); i++) {
        xp[i] = static_cast<double>(x[i]);
    }
    solve_linear_system_impl(solver, rhs, x, xp);
}
