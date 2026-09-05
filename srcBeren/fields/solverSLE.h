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

#include <chrono>
#include <type_traits>

#include <Eigen/SparseLU>

#include "bmatrix.h"
#include "config.h"
#include "containers.h"
#include "env_options.h"
#include "sparse.h"
#include "thread_partitioned_matrix.h"
#include "timer.h"
#include "util.h"

#define DEFAULT_MAX_ITERATIONS 1000
#define DEFAULT_TOLERANCE      1.e-9

template <typename VectorType>
bool bicgstab_iteration(const Operator &A, const VectorType &rhs, VectorType &x, const VectorType &diagonal,
                        size_t &iters, double &tol_error);

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
    void computeDiagonalPreconditioner(const VectorType &diag) {
        RECORD_TIMER;
        m_diagonal.resize(diag.size());

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
void solve_linear_system_impl(SolverType &solver, const VectorType &rhs, VectorType &x, const VectorType &x0,
                              const char *tag = "") {
    solver.setTolerance(SLE_SOLVER_TOLERANCE);
    solver.setMaxIterations(SLE_SOLVER_MAX_ITERATIONS);

    x = solver.solveWithGuess(rhs, x0);

    static const bool solverLog = getenv("SOLVER_LOG") != nullptr;
    if (solverLog) {
        std::cout << "SOLVER [" << tag << "] iters=" << solver.iterations() << " err=" << solver.error() << "\n";
    }

    if (solver.iterations() >= SLE_SOLVER_MAX_ITERATIONS) {
        std::cout << "Field solver failed!" << std::endl;
        std::cout << solver.error() << std::endl;
    }
}

template <typename SolverType, typename VectorType>
void solve_linear_system(const Operator &A, const VectorType &rhs, VectorType &x, const VectorType &x0,
                         const char *tag = "", const VectorType *diag = nullptr) {
    RECORD_TIMER;
    if (diag) {
        SolverType solver(A, *diag);
        solve_linear_system_impl(solver, rhs, x, x0, tag);
    } else {
        SolverType solver(A);
        solve_linear_system_impl(solver, rhs, x, x0, tag);
    }
}

// Диагональ разреженной (row-major CSR) матрицы; отсутствующие диагональные
// элементы — нули. Используется для Jacobi-предобуславливателя.
template <typename VectorType>
inline VectorType sparse_diagonal(const Operator &A, int n) {
    VectorType d(n);
    d.setZero();
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, MAJOR>::InnerIterator it(A, k); it; ++it) {
            if (it.col() == k) {
                d(k) = it.value();
                break;
            }
        }
    }
    return d;
}

// ---------------------------------------------------------------------------
// M3: предобуславливатели для BiCGSTAB (правопредобусловленный вариант ниже).
// Основное применение — фиксированный оператор corrector
// IMmat = I + 0.25*dt^2*curlB*curlB^T, факторизуется один раз на прогон.
// ---------------------------------------------------------------------------
struct IPreconditioner {
    virtual ~IPreconditioner() = default;
    virtual void apply(const Field3d& in, Field3d& out) = 0;
    // Факторизация ли это ИСХОДНОЙ матрицы системы? Для точного P корректен
    // one-shot: x = x0 + P*(b - A*x0). Внимание: predictor-система A = IMmat + L2
    // отличается от факторизованного IMmat — там isExact() означает лишь
    // «сильный P», итерации продолжаются до истинной невязки (стационарная
    // итерация); one-shot корректен только для corrector (A = IMmat).
    virtual bool isExact() const {
        return false;
    }
};

// Jacobi (диагональ)
struct JacobiPreconditioner : IPreconditioner {
    Field3d diag;
    explicit JacobiPreconditioner(const Field3d& d) : diag(d) {
    }
    void apply(const Field3d& in, Field3d& out) override {
#pragma omp parallel for simd
        for (int i = 0; i < in.size(); ++i) {
            // diag==0 (структурный ноль диагонали) — тождественный шаг вместо
            // деления на ноль
            const double d = diag(i);
            out(i) = (d != 0.0) ? in(i) / d : in(i);
        }
    }
};

// Точная LU-факторизация (fp64): один раз на весь прогон, apply = прямой ход.
// Используется SparseLU (не LDLT): на периодических сетках IMmat НЕ симметричен
// (wrap_index в ghost-строках ломает симметрию curlB/curlE), и симметричная
// часть не SPD — LDLT на ней расходится. SparseLU обрабатывает общий случай.
struct SparseLUPreconditioner : IPreconditioner {
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> lu;
    explicit SparseLUPreconditioner(const Operator& A) {
        const auto t0 = std::chrono::steady_clock::now();
        Eigen::SparseMatrix<double, Eigen::ColMajor> Ac = A;
        lu.compute(Ac);
        const double tf = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "SparseLUPreconditioner: factorized " << A.rows() << "x" << A.cols() << " (nnz="
                  << A.nonZeros() << ") in " << tf << " s\n";
    }
    void apply(const Field3d& in, Field3d& out) override {
        out.data() = lu.solve(in.data());
    }
    bool isExact() const override {
        return true;
    }
};

// Тот же SparseLU в fp32: вдвое меньше памяти, быстрее apply; точность
// восстанавливается внешним fp64-BiCGSTAB (низкоточный предобуславливатель)
struct SparseLUPreconditionerFp32 : IPreconditioner {
    Eigen::SparseLU<Eigen::SparseMatrix<float, Eigen::ColMajor>> lu;
    explicit SparseLUPreconditionerFp32(const Operator& A) {
        Eigen::SparseMatrix<float, Eigen::ColMajor> Ac = A.cast<float>();
        lu.compute(Ac);
        std::cout << "SparseLUPreconditionerFp32: factorized " << A.rows() << "x" << A.cols() << "\n";
    }
    void apply(const Field3d& in, Field3d& out) override {
        Eigen::VectorXf xf = in.data().cast<float>();
        xf = lu.solve(xf);
        out.data() = xf.cast<double>();
    }
};

// Правопредобусловленный BiCGSTAB (fp64, без mixed-фаз): мониторится истинная
// невязка (важно для энергетики, MATH.md §7). AP_cache — опциональный заранее
// собранный (thread-partitioned) оператор, чтобы не пересобирать матрицу
// фиксированной системы каждый вызов.
bool bicgstab_iteration_precond(const Operator& A, const Field3d& rhs, Field3d& x, IPreconditioner& P, size_t& iters,
                                double& tol_error, const ThreadPartitionedSparseMatrix<double>* AP_cache = nullptr);

template <typename VectorType>
void solve_linear_system_precond(const Operator& A, const VectorType& rhs, VectorType& x, const VectorType& x0,
                                 IPreconditioner& P, const char* tag = "",
                                 const ThreadPartitionedSparseMatrix<double>* AP_cache = nullptr) {
    RECORD_TIMER;
    size_t iters = SLE_SOLVER_MAX_ITERATIONS;
    double tol = SLE_SOLVER_TOLERANCE;
    x = x0;
    bicgstab_iteration_precond(A, rhs, x, P, iters, tol, AP_cache);
    static const bool solverLog = getenv("SOLVER_LOG") != nullptr;
    if (solverLog) {
        std::cout << "SOLVER [" << tag << "] iters=" << iters << " err=" << tol << "\n";
    }
    if (iters >= SLE_SOLVER_MAX_ITERATIONS) {
        std::cout << "Field solver failed!" << std::endl;
    }
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
