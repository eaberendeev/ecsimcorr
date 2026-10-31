#include <omp.h>

#include <random>

#include "Damping.h"
#include "field_system.h"
#include "maxwell_operators.h"
#include "sparse.h"
#include "test_support.h"
#include "thread_partitioned_matrix.h"

void sparse_sum_and_products_match_eigen() {
    for (int scenario = 0; scenario < 8; ++scenario) {
        const int rows = scenario == 0 ? 0 : 19, cols = 23;
        Operator a(rows, cols), b(rows, cols);
        std::vector<Eigen::Triplet<double>> ta, tb;
        if (scenario > 1)
            for (int r = 0; r < rows; ++r)
                for (int c = 0; c < cols; ++c) {
                    if ((r * 7 + c * 3 + scenario) % 11 == 0)
                        ta.emplace_back(r, c, 0.25 * (r - c));
                    if ((r + c + scenario) % 7 == 0)
                        tb.emplace_back(r, c, -0.25 * (r - c));
                }
        a.setFromTriplets(ta.begin(), ta.end());
        b.setFromTriplets(tb.begin(), tb.end());
        const Operator expected = a + b, actual = parallelSparseSum(a, b);
        if (rows > 0)
            check::near((expected - actual).norm(), 0, 1e-13, "sparse union: empty rows, overlap, cancellation");
        else
            check::require(actual.rows() == 0 && actual.cols() == cols && actual.nonZeros() == 0,
                           "empty sparse matrix shape/portrait changed");
        Field3d x(cols), result(rows);
        for (int i = 0; i < cols; ++i) x[i] = std::cos(i * 0.7);
        spmv(actual, x, result);
        check::near((result.data() - expected * x.data()).norm(), 0, 1e-12, "CSR multiply vs independent Eigen path");
        if (rows > 0) {
            ThreadPartitionedSparseMatrix<double> partitioned(actual);
            spmv(partitioned, x, result);
            check::near((result.data() - expected * x.data()).norm(), 0, 1e-12, "partitioned CSR multiply");
        }
    }
}

void parallel_reductions_match_serial_sums() {
    for (int size : {0, 1, 31, 16383, 32769}) {
        Eigen::VectorXd x(size), y(size);
        long double dot = 0, square = 0;
        for (int i = 0; i < size; ++i) {
            x[i] = std::sin(i * 0.3);
            y[i] = std::cos(i * 0.7);
            dot += static_cast<long double>(x[i]) * y[i];
            square += static_cast<long double>(x[i]) * x[i];
        }
        for (int threads : {1, 2, 4}) {
            omp_set_num_threads(threads);
            check::near(blas::dot(x, y), double(dot), 1e-9, "SIMD/OpenMP dot");
            check::near(blas::squaredNorm(x), double(square), 1e-8, "SIMD/OpenMP norm");
            check::near(blas::normalizedDot(x, y), square ? double(dot / square) : 0, 1e-12,
                        "SIMD/OpenMP normalized dot");
            const auto first = blas::dot(x, y);
            check::near(blas::dot(x, y), first, 0, "fixed-thread reductions must repeat exactly");
        }
    }
    omp_set_num_threads(2);
}

// E=(y,z,x) and the analogous Yee B field have curl=(-1,-1,-1).
// The second check uses an arbitrary field: the discrete identity div(curl)=0
// must hold independently of a particular analytic solution.
void yee_curl_and_divergence_match_known_fields() {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(0.7, 1.1, 0.9));
    BoundaryConditionHandler bc;
    const int n = 3 * domain.total_size();
    Operator curl_e(n, n), curl_b(n, n), div(domain.total_size(), n);
    auto stencils = assemble_maxwell_stencils(domain, bc);
    curl_e = stencils.curlE;
    curl_b = stencils.curlB;
    div = stencils.divE;
    Field3d e(domain.size(), 3), b(domain.size(), 3);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k)
                for (int d = 0; d < 3; ++d) {
                    e(i, j, k, d) = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d)[(d + 1) % 3];
                    b(i, j, k, d) = domain.get_node_position(i, j, k, FieldType::MAGNETIC, d)[(d + 1) % 3];
                }
    const Eigen::VectorXd ce = curl_e * e.data(), cb = curl_b * b.data();
    for (int i = 2; i < domain.size().x() - 2; ++i)
        for (int j = 2; j < domain.size().y() - 2; ++j)
            for (int k = 2; k < domain.size().z() - 2; ++k)
                for (int d = 0; d < 3; ++d) {
                    check::near(ce[domain.vind(i, j, k, d)], -1, 1e-13, "curl E on an affine field");
                    check::near(cb[domain.vind(i, j, k, d)], -1, 1e-13, "curl B on an affine field");
                }
    std::mt19937 engine(17);
    std::uniform_real_distribution<double> uniform(-1, 1);
    for (int i = 0; i < e.capacity(); ++i) e[i] = uniform(engine);
    const Eigen::VectorXd identity = div * (curl_b * e.data());
    for (int i = 2; i < domain.size().x() - 2; ++i)
        for (int j = 2; j < domain.size().y() - 2; ++j)
            for (int k = 2; k < domain.size().z() - 2; ++k)
                check::near(identity[domain.sind(i, j, k)], 0, 1e-12, "discrete div(curl B) must vanish");
}

void damping_accounts_for_removed_field_energy() {
    Domain domain;
    domain.init(Vector3I(6, 6, 6), Vector3R(1, 1, 1));
    for (const std::string type : {"Rectangle", "CircleXY"}) {
        Field3d e(domain.size(), 3), b(domain.size(), 3);
        e.data().setConstant(0.7);
        b.data().setConstant(-0.3);
        const double before = 0.5 * (e.squaredNorm() + b.squaredNorm());
        nlohmann::json cfg = {{"DampingType", type},
                              {"DampingCoeff", 0.8},
                              {"DampCellsX_glob", {2, 1}},
                              {"DampCellsY_glob", {1, 2}},
                              {"DampCellsZ_glob", {2, 2}}};
        const double removed = damping_fields(e, b, domain, cfg);
        check::near(removed, before - 0.5 * (e.squaredNorm() + b.squaredNorm()), 1e-10,
                    "damping must account for all removed energy");
        check::require(removed > 0, "nonempty damping layer did nothing");
    }
}

// An affine E has curl=(-1,-1,-1), so Faraday increases each B component
// by dt when both electric endpoints are equal. A constant E cannot change B.
void faraday_uses_both_electric_endpoints() {
    Domain domain;
    domain.init(Vector3I(4, 5, 6), Vector3R(0.7, 1.1, 0.9));
    BoundaryConditionHandler boundaries;
    MaxwellOperators maxwell;
    maxwell.init(assemble_maxwell_stencils(domain, boundaries));
    Field3d E(domain.size(), 3), En(domain.size(), 3), B(domain.size(), 3), Bn(domain.size(), 3);
    B.data().setConstant(0.3);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k)
                for (int d = 0; d < 3; ++d) {
                    E(i, j, k, d) = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d)[(d + 1) % 3];
                }
    // Distinct endpoints detect using E or En alone in place of their average.
    En = 3.0 * E;
    maxwell.advance_magnetic_field(Bn, B, E, En, 0.2);
    for (int i = 2; i < domain.size().x() - 2; ++i)
        for (int j = 2; j < domain.size().y() - 2; ++j)
            for (int k = 2; k < domain.size().z() - 2; ++k)
                for (int d = 0; d < 3; ++d)
                    check::near(Bn(i, j, k, d), 0.7, 1e-13, "Faraday sign, dt and endpoint average");
    E.data().setConstant(0.6);
    En = E;
    maxwell.advance_magnetic_field(Bn, B, E, En, 0.2);
    check::near((Bn - B).norm(), 0, 1e-13, "constant E leaves B unchanged");
}

void field_system_rebuilds_for_changed_timestep() {
    Domain domain;
    domain.init(Vector3I(3, 4, 5), Vector3R(0.7, 1.1, 0.9));
    BoundaryConditionHandler boundaries;
    auto stencils = assemble_maxwell_stencils(domain, boundaries);
    EcsimFieldSystem system;
    system.init(stencils, domain, 0.1, boundaries);
    const Operator first_M = system.M;
    system.L.setIdentity();
    system.init(stencils, domain, 0.2, boundaries);
    check::require(first_M.norm() > 0, "nontrivial Maxwell operator required");
    check::near((system.M - 4.0 * first_M).norm(), 0, 1e-13, "Maxwell combination scales as dt squared");
    Operator identity(system.IM.rows(), system.IM.cols());
    identity.setIdentity();
    check::near((system.IM + system.M - identity).norm(), 0, 1e-13, "unconstrained IM=I-M");
    check::near(system.L.norm(), 0, 0, "reinitialization clears the previous plasma response");
}

int main() {
    return check::run([] {
        check::scenario("Sparse matrix sum and products agree with independent Eigen operations",
                        sparse_sum_and_products_match_eigen);
        check::scenario("Parallel scalar products agree with long-double serial sums",
                        parallel_reductions_match_serial_sums);
        check::scenario("Yee curl reproduces an affine field and div(curl B)=0",
                        yee_curl_and_divergence_match_known_fields);
        check::scenario("Faraday uses the mean electric field with the correct sign and dt",
                        faraday_uses_both_electric_endpoints);
        check::scenario("Scheme matrices rebuild with dt squared and clear the plasma response",
                        field_system_rebuilds_for_changed_timestep);
        check::scenario("Damping reports the actual decrease of E/B energy", damping_accounts_for_removed_field_energy);
    });
}
