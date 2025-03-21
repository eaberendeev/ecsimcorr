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
#endif // USE_AMGCL

#include <type_traits>

#include "bmatrix.h"
#include "containers.h"
#include "util.h"

#define DEFAULT_MAX_ITERATIONS 1000
#define DEFAULT_TOLERANCE      1.e-9

template <typename VectorType>
inline void spmv(const Operator &A, const VectorType &v,
                       VectorType &res) {
    int rows = A.rows();

    const double* val = A.valuePtr();
    const int* inner = A.innerIndexPtr();
    const int* outer = A.outerIndexPtr();
#pragma omp parallel for
    for (int i = 0; i < rows; ++i) {
        double sum = 0;
#pragma omp simd
	for (int j = outer[i]; j < outer[i+1]; ++j) {
	__builtin_prefetch(&v[inner[j + 4]]); // Предзагрузка через 4 элемента
            sum += val[j] * v[inner[j]];
        }
        res[i] = sum;
    }
}


template <typename VectorType>
bool bicgstab_iteration(const Operator &A , const VectorType &rhs,
                         VectorType &x, const VectorType &diagonal,
                         size_t &iters, double &tol_error) {
    using std::abs;
    using std::sqrt;
    double tol = tol_error;
    int maxIters = iters;
    int n = x.size();

//    VectorType r = rhs - Spmv(x);
    VectorType r(n);
    spmv(A, x,r);
    r = rhs - r;

    VectorType r0 = r;
    double r0_sqnorm = r0.squaredNorm();
    double rhs_sqnorm = rhs.squaredNorm();
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
    double eps2 = Eigen::NumTraits<double>::epsilon() *
                  Eigen::NumTraits<double>::epsilon();
    int i = 0;
    int restarts = 0;
    double time1 = 0;
    double timeall = omp_get_wtime();

    while (r.squaredNorm() > tol2 && i < maxIters) {
        double rho_old = rho;
        rho = r0.dot(r);
        if (abs(rho) < eps2 * r0_sqnorm) {
            double time10 = omp_get_wtime();
            //r = rhs - Spmv(x);
            spmv(A, x,r);
            r = rhs - r;
            time1 += omp_get_wtime() - time10;
            r0 = r;
            rho = r0_sqnorm = r.squaredNorm();
            if (restarts++ == 0)
                i = 0;
        }
        double beta = (rho / rho_old) * (alpha / w);
#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            p(i) = r(i) + beta * (p(i) - w * v(i));
            y(i) = p(i) / diagonal(i);
        }
        // p = r + beta * (p - w * v);
        // y = precond.solve(p);   // Применение предобуславливателя
        double time10 = omp_get_wtime();
        //v = Spmv(y);
        spmv(A, y, v);
        time1 += omp_get_wtime() - time10;

        alpha = rho / r0.dot(v);

        // s = r - alpha * v;
        // z = precond.solve(s);   // Применение предобуславливателя

#pragma omp parallel for simd
        for (int i = 0; i < n; i++) {
            s(i) = r(i) - alpha * v(i);
            z(i) = s(i) / diagonal(i);
        }
        time10 = omp_get_wtime();
        //t = Spmv(z);
        spmv(A, z,t);
        time1 += omp_get_wtime() - time10;

        double tmp = t.squaredNorm();
        if (tmp > 0)
            w = t.dot(s) / tmp;
        else
            w = 0;

            // x += alpha * y + w * z;
            // r = s - w * t;
#pragma omp parallel for simd 
        for (int i = 0; i < n; i++) {
            x(i) += alpha * y(i) + w * z(i);
            r(i) = s(i) - w * t(i);
        }
        ++i;
    }
    std::cout << "Time: " << omp_get_wtime() - timeall << std::endl;
    std::cout << "Time1: " << time1 << std::endl;
    tol_error = sqrt(r.squaredNorm() / rhs_sqnorm);
    iters = i;
    return true;
}

template <typename VectorType>
class BicgstabSolverBase {
   public:
    BicgstabSolverBase()
        : max_iterations(DEFAULT_MAX_ITERATIONS),
          m_tolerance(DEFAULT_TOLERANCE),
          m_success(false) {}

    void setTolerance(double tolerance) { m_tolerance = tolerance; }
    void setMaxIterations(size_t max_iters) { max_iterations = max_iters; }
    bool info() const { return m_success; }
    double error() const { return m_error; }
    size_t iterations() const { return m_iterations; }

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

        m_success = bicgstab_iteration(m_A, rhs, x, m_diagonal,
                                       m_iterations, m_error);
        return x;
    }

   private:
    void computeDiagonalPreconditioner(const Eigen::VectorXd &diag) {
        for (int i = 0; i < m_diagonal.size(); i++) {
            m_diagonal[i] = diag[i];
        }
    }
    void initializePreconditioner(int rows) {
        m_diagonal.resize(rows);
        std::fill(m_diagonal.begin(), m_diagonal.end(), 1.0);
    }
    const Operator &m_A;
};

template <typename SolverType, typename VectorType>
void solve_linear_system_impl(SolverType &solver, const VectorType &rhs,
                              VectorType &x, const VectorType &x0) {
    solver.setTolerance(SLE_SOLVER_TOLERANCE);
    solver.setMaxIterations(SLE_SOLVER_MAX_ITERATIONS);

    x = solver.solveWithGuess(rhs, x0);

    if (solver.iterations() >= SLE_SOLVER_MAX_ITERATIONS) {
        std::cout << "Field solver failed!" << std::endl;
        std::cout << solver.error() << std::endl;
    }
}

template <typename SolverType, typename VectorType>
void solve_linear_system(const Operator& A, const VectorType& rhs,
                     VectorType& x, const VectorType& x0) {
 SolverType solver(A);
 solve_linear_system_impl(solver, rhs, x, x0);
}

template <typename VectorType>
void solve_linear_system(const Operator &A,
                         const VectorType &diagonal, const VectorType &rhs,
                         VectorType &x, const VectorType &x0) {
    BicgstabSolver solver(A, diagonal);
    solve_linear_system_impl(solver, rhs, x, x0);
}

#ifdef USE_AMGCL

using namespace amgcl;

template <typename MatrixType>
void solve_amgcl(const MatrixType &A, const Field &rhs, Field &x,
                 const Field &x0) {
    typedef backend::eigen<double> Backend;

    typedef amgcl::make_solver<
        // Use AMG as preconditioner:
        amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation,
                   amgcl::relaxation::spai0>,

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

    //auto prm = precond.amgcl_params();
}
#endif // USE_AMGCL

typedef Eigen::SparseMatrix<float, Eigen::RowMajor> OperatorM;
typedef Eigen::GMRES<Eigen::SparseMatrix<float, Eigen::RowMajor> > gmres_m;
typedef Eigen::BiCGSTAB<Eigen::SparseMatrix<float, Eigen::RowMajor> > bicgstab_m;
typedef Eigen::VectorXf VectorXf;

static Eigen::SparseMatrix<float, Eigen::RowMajor> convertSparseMatrixDoubleToFloat(
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &matDouble) {
    // Создаём матрицу float с теми же размерами
    Eigen::SparseMatrix<float, Eigen::RowMajor> matFloat(matDouble.rows(),
                                                         matDouble.cols());

    // Копируем структуру матрицы (индексы строк и столбцов)
    matFloat.resizeNonZeros(matDouble.nonZeros());
    std::copy(matDouble.outerIndexPtr(),
              matDouble.outerIndexPtr() + matDouble.outerSize() + 1,
              matFloat.outerIndexPtr());
    std::copy(matDouble.innerIndexPtr(),
              matDouble.innerIndexPtr() + matDouble.nonZeros(),
              matFloat.innerIndexPtr());

    // Преобразуем значения из double в float
    const double *valuesDouble = matDouble.valuePtr();
    float *valuesFloat = matFloat.valuePtr();
#pragma omp parallel for num_threads(8)
    for (int i = 0; i < matDouble.nonZeros(); ++i) {
        valuesFloat[i] = static_cast<float>(valuesDouble[i]);
    }

    return matFloat;
}

inline void solve_linear_system_mix(const Operator &A, const Field &rhs, Field &x,
                             const Field &x0) {
    OperatorM AM = convertSparseMatrixDoubleToFloat(A);
     bicgstab_m solverM(AM);
     VectorXf rhsM(rhs.size());
     VectorXf xM(rhs.size());
     VectorXf x0M(rhs.size());
     Field xp(rhs.size());
     for (int i = 0; i < rhs.size(); i++) {
         rhsM[i] = static_cast<float> (rhs[i]);
         x0M[i] = static_cast<float> (x0[i]);
     }
    solve_linear_system_impl(solverM, rhsM, xM, x0M);
    bicgstab solver(A);
    for (int i = 0; i < rhs.size(); i++) {
        xp[i] = static_cast<double>(x[i]);
    }
    solve_linear_system_impl(solver, rhs, x, xp);
}
