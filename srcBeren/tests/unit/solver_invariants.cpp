#include "solverSLE.h"
#include "test_support.h"

bool bicgstab_iteration_mixed_precision(const Operator& A, const Field3d& rhs, Field3d& x, const Field3d& diagonal,
                                        size_t& iters, double& tolerance, double& residual);
bool bicgstab_iteration_mixed_precision_greedy(const Operator& A, const Field3d& rhs, Field3d& x,
                                               const Field3d& diagonal, size_t& iters, double& tolerance,
                                               double& residual);

void test_partial_float_solve() {
    constexpr int n = 12;
    Operator matrix(n, n);
    Field3d rhs(n), guess(n), diagonal(n), expected(n);
    rhs.data().setOnes();
    guess.setZero();
    diagonal.data().setOnes();
    for (int i = 0; i < n; ++i) {
        matrix.insert(i, i) = i + 1;
        expected[i] = 1.0 / (i + 1);
    }
    matrix.makeCompressed();

    // Seven float iterations improve the guess without reaching the float
    // tolerance. Both layouts must carry that progress into the double stage.
    for (const auto solve : {bicgstab_iteration_mixed_precision, bicgstab_iteration_mixed_precision_greedy}) {
        Field3d actual = guess;
        size_t iterations = 14;
        double tolerance = 1e-9;
        double reported_residual = 0;
        const bool success = solve(matrix, rhs, actual, diagonal, iterations, tolerance, reported_residual);
        const double residual = (matrix * actual.data() - rhs.data()).norm();
        check::require(success, "partial float solve discarded progress needed within the iteration budget");
        check::require(iterations <= 14, "mixed solve exceeded the total iteration budget");
        check::require(std::isfinite(residual) && residual / rhs.norm() <= 1e-9,
                       "mixed solve failed the independent residual criterion");
        check::near((actual - expected).norm(), 0, 1e-8, "diagonal system has a known analytic solution");
        check::near(reported_residual, residual, 1e-12, "mixed solve reported an incorrect residual");
    }
}

void test_solver() {
    for (int kind = 0; kind < 4; ++kind) {
        const int n = 6;
        Operator matrix(n, n);
        std::vector<Eigen::Triplet<double>> entries;
        for (int i = 0; i < n; ++i) {
            entries.emplace_back(i, i, kind == 0 ? 1 : 4 + i);
            if (kind >= 2 && i + 1 < n)
                entries.emplace_back(i, i + 1, -0.3);
            if (kind >= 2 && i > 0)
                entries.emplace_back(i, i - 1, 0.2);
        }
        matrix.setFromTriplets(entries.begin(), entries.end());
        Field3d expected(Vector3I(2, 1, 1), 3), rhs(Vector3I(2, 1, 1), 3), guess(Vector3I(2, 1, 1), 3), diagonal(n);
        for (int i = 0; i < n; ++i) {
            expected[i] = 0.25 + i * 0.4;
            guess[i] = kind == 3 ? expected[i] : -0.1;
            diagonal[i] = matrix.coeff(i, i);
        }
        rhs.data() = matrix * expected.data();
        BicgstabSolver<Field3d> solver(matrix, diagonal);
        solver.setTolerance(1e-11);
        solver.setMaxIterations(100);
        const auto actual = solver.solveWithGuess(rhs, guess);
        check::require(solver.info(), "well-conditioned system reported failure");
        const double residual = (matrix * actual.data() - rhs.data()).norm();
        check::require(std::isfinite(residual) && residual < 1e-9,
                       "solver returned an incorrect or non-finite solution");
        check::near(solver.error(), residual / rhs.norm(), 1e-10, "reported residual differs from actual");
        if (kind == 3)
            check::require(solver.iterations() == 0, "exact initial guess should require no iterations");
        if (kind == 0)
            check::require(solver.iterations() <= 2, "identity should converge in one alpha step per precision stage");
        rhs.setZero();
        const auto zero = solver.solveWithGuess(rhs, expected);
        check::require(solver.info() && zero.norm() == 0, "zero RHS must give a finite zero solution");
    }
    Operator zero(3, 3);
    zero.makeCompressed();
    Field3d rhs(3), x(3);
    rhs.data().setOnes();
    x.setZero();
    BicgstabSolver<Field3d> singular(zero);
    singular.setMaxIterations(3);
    singular.solveWithGuess(rhs, x);
    check::require(!singular.info(), "inconsistent singular system must report failure");
    check::require(std::isfinite(singular.error()) && singular.error() >= 0.9,
                   "breakdown must preserve an honest residual");
    bool rejected = false;
    try {
        solve_linear_system<BicgstabSolver<Field3d>>(zero, rhs, x, x);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    check::require(rejected, "simulation wrapper must reject solver failure");

    Operator diagonal_matrix(3, 3);
    diagonal_matrix.setIdentity();
    Operator preconditioner_matrix(64, 64);
    preconditioner_matrix.setIdentity();
    for (double invalid : {0.0, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        Field3d diagonal(64);
        diagonal.data().setOnes();
        // With two threads and chunks of 16, thread 1 must report this error.
        diagonal[63] = invalid;
        bool rejected_diagonal = false;
        try {
            BicgstabSolver<Field3d> invalid_solver(preconditioner_matrix, diagonal);
        } catch (const std::invalid_argument&) {
            rejected_diagonal = true;
        }
        check::require(rejected_diagonal, "invalid preconditioner diagonal must be rejected");
    }
    BicgstabSolver<Field3d> limited(diagonal_matrix);
    for (double tolerance : {0.0, -1e-9, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
        bool rejected_tolerance = false;
        try {
            limited.setTolerance(tolerance);
        } catch (const std::invalid_argument&) {
            rejected_tolerance = true;
        }
        check::require(rejected_tolerance, "invalid tolerance must not permit false solver success");
    }
    limited.setMaxIterations(0);
    limited.solveWithGuess(rhs, x);
    check::require(!limited.info(), "iteration budget exhausted without convergence must fail");
}
int main() {
    return check::run([] {
        test_solver();
        test_partial_float_solve();
    });
}
