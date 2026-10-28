#include "solverSLE.h"
#include "test_support.h"

bool bicgstab_iteration_mixed_precision(const Operator& A, const Field3d& rhs, Field3d& x, const Field3d& diagonal,
                                        size_t& iters, double& tolerance, double& residual);
bool bicgstab_iteration_mixed_precision_greedy(const Operator& A, const Field3d& rhs, Field3d& x,
                                               const Field3d& diagonal, size_t& iters, double& tolerance,
                                               double& residual);

void mixed_precision_keeps_progress_before_double_refinement() {
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

// A known vector x_exact defines b=A*x_exact. The test independently measures
// ||A*x-b|| after the solve, rather than accepting the solver's success flag.
enum class KnownSystem { Identity, Diagonal, Nonsymmetric, ExactInitialGuess };

void solve_known_system(KnownSystem problem) {
    constexpr int size = 6;
    const bool identity = problem == KnownSystem::Identity;
    const bool coupled = problem == KnownSystem::Nonsymmetric || problem == KnownSystem::ExactInitialGuess;
    const bool exact_guess = problem == KnownSystem::ExactInitialGuess;
    Operator matrix(size, size);
    std::vector<Eigen::Triplet<double>> entries;
    for (int row = 0; row < size; ++row) {
        entries.emplace_back(row, row, identity ? 1 : 4 + row);
        if (coupled && row + 1 < size)
            entries.emplace_back(row, row + 1, -0.3);
        if (coupled && row > 0)
            entries.emplace_back(row, row - 1, 0.2);
    }
    matrix.setFromTriplets(entries.begin(), entries.end());
    Field3d expected(Vector3I(2, 1, 1), 3), rhs(expected.sizes(), 3), guess(expected.sizes(), 3), diagonal(size);
    for (int i = 0; i < size; ++i) {
        expected[i] = 0.25 + i * 0.4;
        guess[i] = exact_guess ? expected[i] : -0.1;
        diagonal[i] = matrix.coeff(i, i);
    }
    rhs.data() = matrix * expected.data();

    BicgstabSolver<Field3d> solver(matrix, diagonal);
    solver.setTolerance(1e-11);
    solver.setMaxIterations(100);
    const auto actual = solver.solveWithGuess(rhs, guess);
    const double residual = (matrix * actual.data() - rhs.data()).norm();
    std::cout << "measured ||A*x-b||=" << std::scientific << residual
              << " (limit 1e-9), iterations=" << solver.iterations() << '\n';
    check::require(solver.info(), "well-conditioned system reported failure");
    check::require(std::isfinite(residual) && residual < 1e-9, "solver returned an incorrect or non-finite solution");
    check::near(solver.error(), residual / rhs.norm(), 1e-10, "reported residual differs from actual");
    if (exact_guess)
        check::require(solver.iterations() == 0, "exact initial guess should require no iterations");
    if (identity)
        check::require(solver.iterations() <= 2, "identity should converge in one alpha step per precision stage");

    rhs.setZero();
    const auto zero = solver.solveWithGuess(rhs, expected);
    check::require(solver.info() && zero.norm() == 0, "zero RHS must give a finite zero solution");
}

void inconsistent_singular_system_reports_failure() {
    // 0*x=1 has no solution. Neither the solver nor its simulation wrapper may
    // report success just because iteration stopped.
    Operator matrix(3, 3);
    matrix.makeCompressed();
    Field3d rhs(3), guess(3);
    rhs.data().setOnes();
    guess.setZero();
    BicgstabSolver<Field3d> singular(matrix);
    singular.setMaxIterations(3);
    singular.solveWithGuess(rhs, guess);
    check::require(!singular.info(), "inconsistent singular system must report failure");
    check::require(std::isfinite(singular.error()) && singular.error() >= 0.9,
                   "breakdown must preserve an honest residual");
    bool rejected = false;
    try {
        solve_linear_system<BicgstabSolver<Field3d>>(matrix, rhs, guess, guess);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    check::require(rejected, "simulation wrapper must reject solver failure");
}

void invalid_preconditioner_and_stopping_parameters_are_rejected() {
    Operator preconditioner_matrix(64, 64);
    preconditioner_matrix.setIdentity();
    for (double invalid : {0.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        Field3d diagonal(64);
        diagonal.data().setOnes();
        // With two threads and chunks of 16, thread 1 must report this error.
        diagonal[63] = invalid;
        bool rejected = false;
        try {
            BicgstabSolver<Field3d> invalid_solver(preconditioner_matrix, diagonal);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check::require(rejected, "invalid preconditioner diagonal must be rejected");
    }
    Operator matrix(3, 3);
    matrix.setIdentity();
    BicgstabSolver<Field3d> limited(matrix);
    for (double tolerance :
         {0.0, -1e-9, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        bool rejected = false;
        try {
            limited.setTolerance(tolerance);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check::require(rejected, "invalid tolerance must not permit false solver success");
    }
    Field3d rhs(3), guess(3);
    rhs.data().setOnes();
    guess.setZero();
    limited.setMaxIterations(0);
    limited.solveWithGuess(rhs, guess);
    check::require(!limited.info(), "iteration budget exhausted without convergence must fail");
}

int main() {
    return check::run([] {
        check::scenario("Identity matrix: known solution and zero RHS",
                        [] { solve_known_system(KnownSystem::Identity); });
        check::scenario("Diagonal matrix: known solution and zero RHS",
                        [] { solve_known_system(KnownSystem::Diagonal); });
        check::scenario("Nonsymmetric matrix: independent residual and zero RHS",
                        [] { solve_known_system(KnownSystem::Nonsymmetric); });
        check::scenario("An exact initial guess requires no iterations",
                        [] { solve_known_system(KnownSystem::ExactInitialGuess); });
        check::scenario("An inconsistent singular system must fail", inconsistent_singular_system_reports_failure);
        check::scenario("Invalid preconditioner/tolerance and zero iteration budget must fail",
                        invalid_preconditioner_and_stopping_parameters_are_rejected);
        check::scenario("Mixed precision preserves float progress for the final double solve",
                        mixed_precision_keeps_progress_before_double_refinement);
    });
}
