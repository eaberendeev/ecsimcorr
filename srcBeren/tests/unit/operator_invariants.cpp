#include <omp.h>

#include <Eigen/LU>
#include <Eigen/SparseLU>
#include <random>

#include "config/boundary_config.h"
#include "config/field_config.h"
#include "fields/Damping.h"
#include "fields/maxwell_operators.h"
#include "simulation/boundary_conditions.h"
#include "simulation/ecsim/field_system.h"
#include "test_support.h"
#include "utils/sparse.h"
#include "utils/thread_partitioned_matrix.h"

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
    {
        constexpr int ghost = GHOST_CELLS;
        Domain domain;
        domain.init(Vector3I(6, 6, 6), Vector3R(0.7, 1.1, 0.9));
        const int first = std::max(2, ghost + 1);
        const int end = std::min(domain.size().x() - 2, ghost + domain.num_cells().x() - 1);
        FieldBoundaryConditions bc;
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
                        e(i, j, k, d) = domain.storage_node_position(i, j, k, FieldType::ELECTRIC, d)[(d + 1) % 3];
                        b(i, j, k, d) = domain.storage_node_position(i, j, k, FieldType::MAGNETIC, d)[(d + 1) % 3];
                    }
        const Eigen::VectorXd ce = curl_e * e.data(), cb = curl_b * b.data();
        for (int i = first; i < end; ++i)
            for (int j = first; j < end; ++j)
                for (int k = first; k < end; ++k)
                    for (int d = 0; d < 3; ++d) {
                        check::near(ce[domain.vind(i, j, k, d)], -1, 1e-13, "curl E on an affine field");
                        check::near(cb[domain.vind(i, j, k, d)], -1, 1e-13, "curl B on an affine field");
                    }
        std::mt19937 engine(17);
        std::uniform_real_distribution<double> uniform(-1, 1);
        for (int i = 0; i < e.capacity(); ++i) e[i] = uniform(engine);
        const Eigen::VectorXd identity = div * (curl_b * e.data());
        for (int i = first; i < end; ++i)
            for (int j = first; j < end; ++j)
                for (int k = first; k < end; ++k)
                    check::near(identity[domain.sind(i, j, k)], 0, 1e-12, "discrete div(curl B) must vanish");
    }
}

// Compare periodic storage aliases with an independent sum over every raw
// source in the equivalence class. Include periods shorter than a halo band,
// where pairwise folding would reuse already accumulated values.
void periodic_overlap_counts_each_deposit_once() {
    const std::array<Face, 3> faces{Face::XMIN, Face::YMIN, Face::ZMIN};
    constexpr int ghost = GHOST_CELLS;
    for (int period : {1, 2, 3, 7})
        for (int axis = 0; axis < 3; ++axis) {
            Domain domain;
            Vector3I cells(4, 5, 6);
            cells[axis] = period;
            domain.init(cells, Vector3R(0.7, 1.1, 0.9));
            const auto sizes = domain.size();
            PeriodicBoundaryCondition boundary(faces[axis]);
            for (int components : {1, 3}) {
                Field3d field(sizes, components);
                for (int i = 0; i < field.capacity(); ++i) field[i] = 0.2 + std::sin(0.37 * i);
                const Field3d before = field;
                boundary.apply_to_fields(field, components == 1 ? FieldType::DENSITY : FieldType::CURRENT, domain);
                double owner_sum = 0, original_sum = before.data().sum();
                for (int i = 0; i < sizes.x(); ++i)
                    for (int j = 0; j < sizes.y(); ++j)
                        for (int k = 0; k < sizes.z(); ++k) {
                            Vector3I node(i, j, k);
                            const int alias = node[axis];
                            for (int d = 0; d < components; ++d) {
                                long double expected = 0;
                                for (int source = alias % period; source < sizes[axis]; source += period) {
                                    node[axis] = source;
                                    expected += before(node.x(), node.y(), node.z(), d);
                                }
                                check::near(field(i, j, k, d), double(expected), 1e-13,
                                            "periodic density/current aliases sum original deposits once");
                                if (alias >= ghost && alias < ghost + period)
                                    owner_sum += field(i, j, k, d);
                            }
                        }
                check::near(owner_sum, original_sum, 2e-11, "one periodic copy preserves total deposited quantity");
            }

            const int n = 3 * domain.total_size();
            Operator response(n, n);
            std::vector<Trip> trips;
            Field3d electric(sizes, 3);
            for (int row = 0; row < n; ++row) {
                const auto [i, j, k, d] = domain.grid.pos_vind_range(row);
                Vector3I node(i, j, k);
                node[axis] = (node[axis] - ghost + 2 * sizes[axis]) % period;
                electric[row] = std::sin(0.3 * node.x() + 0.5 * node.y() - 0.2 * node.z() + d);
                trips.emplace_back(row, row, 0.25 + 0.01 * (row % 11));
                trips.emplace_back(row, (7 * row + 3) % n, -0.125);
            }
            response.setFromTriplets(trips.begin(), trips.end());
            const Eigen::VectorXd original_action = response * electric.data();
            boundary.apply_to_operator(response, domain);
            const Eigen::VectorXd actual_action = response * electric.data();
            for (int row = 0; row < n; ++row) {
                const auto [i, j, k, d] = domain.grid.pos_vind_range(row);
                Vector3I node(i, j, k);
                long double expected = 0;
                for (int source = node[axis] % period; source < sizes[axis]; source += period) {
                    node[axis] = source;
                    expected += original_action[domain.vind(node.x(), node.y(), node.z(), d)];
                }
                check::near(actual_action[row], double(expected), 1e-12,
                            "periodic operator action equals the sum of original response rows");
            }
        }
}

// Three short periodic axes meet at halo corners: each deposited contribution
// must still enter exactly one physical copy after all three folds.
void periodic_corner_overlap_preserves_total() {
    Domain domain;
    domain.init(Vector3I(1, 2, 3), Vector3R(0.7, 1.1, 0.9));
    std::vector<std::unique_ptr<FieldBoundaryCondition>> conditions;
    for (Face face : {Face::XMIN, Face::YMIN, Face::ZMIN})
        conditions.push_back(std::make_unique<PeriodicBoundaryCondition>(face));
    FieldBoundaryConditions boundaries(std::move(conditions), {true, true, true});
    Field3d density(domain.size(), 1);
    for (int i = 0; i < density.capacity(); ++i) density[i] = 0.2 + std::sin(0.37 * i);
    const Field3d original = density;
    boundaries.apply_to_fields(density, FieldType::DENSITY, domain);
    const auto active = boundaries.active_range(domain.grid);
    double total = 0;
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                long double expected = 0;
                for (int x = i % 1; x < domain.size().x(); x += 1)
                    for (int y = j % 2; y < domain.size().y(); y += 2)
                        for (int z = k % 3; z < domain.size().z(); z += 3) expected += original(x, y, z, 0);
                check::near(density(i, j, k, 0), double(expected), 1e-13, "periodic corner alias sum");
                if (i >= active.start.x() && i < active.end.x() && j >= active.start.y() && j < active.end.y() &&
                    k >= active.start.z() && k < active.end.z())
                    total += density(i, j, k, 0);
            }
    check::near(total, original.data().sum(), 1e-12, "three-axis periodic overlap conserves one-copy charge");
}

// The stored aliases are redundant unknowns; raw symmetry is not required.
// Solve the full folded system and compare with an independently assembled
// physical owner system. Also compare E.L.E over one physical copy, which is
// the response-work contraction used by the field/particle energy balance.
void short_periodic_response_solve_matches_owner_work() {
    for (int period : {1, 2, 3}) {
        constexpr int ghost = GHOST_CELLS;
        Domain domain;
        domain.init(Vector3I(period, 1, 1), Vector3R(0.7, 1.1, 0.9));
        const int n = 3 * domain.total_size();
        const auto owner = [period](int x) { return ((x - ghost) % period + period) % period; };
        const auto raw_row = [&](int x, int d) { return domain.vind(x, ghost, ghost, d); };
        std::vector<Trip> trips;
        Eigen::MatrixXd physical = Eigen::MatrixXd::Zero(3 * period, 3 * period);
        const auto add = [&](int x, int d, int y, int e, double value) {
            trips.emplace_back(raw_row(x, d), raw_row(y, e), value);
            physical(3 * owner(x) + d, 3 * owner(y) + e) += value;
        };
        // A symmetric positive response before periodic assembly, with local
        // spatial and component couplings. Unequal alias counts for N=2/3
        // exercise the representation that cannot remain raw-symmetric.
        for (int x = 0; x < domain.size().x(); ++x)
            for (int d = 0; d < 3; ++d) {
                add(x, d, x, d, 0.3 + 0.05 * x + 0.1 * d);
                if (x + 1 < domain.size().x()) {
                    add(x, d, x + 1, d, 0.05);
                    add(x + 1, d, x, d, 0.05);
                }
                if (d + 1 < 3) {
                    add(x, d, x, d + 1, 0.02);
                    add(x, d + 1, x, d, 0.02);
                }
            }
        Operator response(n, n);
        response.setFromTriplets(trips.begin(), trips.end());
        PeriodicBoundaryCondition boundary(Face::XMIN);
        boundary.apply_to_operator(response, domain);
        Eigen::VectorXd physical_rhs(3 * period), rhs = Eigen::VectorXd::Zero(n);
        for (int q = 0; q < period; ++q)
            for (int d = 0; d < 3; ++d) physical_rhs[3 * q + d] = 1 + 0.1 * q - 0.2 * d;
        for (int x = 0; x < domain.size().x(); ++x)
            for (int d = 0; d < 3; ++d) rhs[raw_row(x, d)] = physical_rhs[3 * owner(x) + d];
        Operator system = 0.2 * response;
        for (int row = 0; row < n; ++row) system.coeffRef(row, row) += 1;
        system.makeCompressed();
        Eigen::SparseMatrix<double> column_system = system;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(column_system);
        check::require(solver.info() == Eigen::Success, "factor the full short-period stored-alias system");
        const Eigen::VectorXd solution = solver.solve(rhs);
        check::require(solver.info() == Eigen::Success, "solve the full short-period stored-alias system");
        const Eigen::MatrixXd physical_system = Eigen::MatrixXd::Identity(3 * period, 3 * period) + 0.2 * physical;
        const Eigen::VectorXd physical_solution = physical_system.fullPivLu().solve(physical_rhs);
        check::near((system * solution - rhs).norm(), 0, 1e-12, "full stored-alias solve residual");
        for (int x = 0; x < domain.size().x(); ++x)
            for (int d = 0; d < 3; ++d)
                check::near(solution[raw_row(x, d)], physical_solution[3 * owner(x) + d], 1e-12,
                            "all solved aliases agree with the independent physical owner solution");
        const Eigen::VectorXd response_action = response * solution;
        double owner_work = 0;
        for (int x = ghost; x < ghost + period; ++x)
            for (int d = 0; d < 3; ++d) owner_work += solution[raw_row(x, d)] * response_action[raw_row(x, d)];
        const double physical_work = physical_solution.dot(physical * physical_solution);
        check::require(physical_work > 0 && std::isfinite(physical_work), "positive finite physical response work");
        check::near(owner_work, physical_work, 1e-12, "one-copy response work equals independent physical contraction");
    }
}

void small_periodic_response() {
    for (int axis = 0; axis < 3; ++axis)
        for (int cells : {1, 2})
            for (int ghost : {1, 2}) {
                Domain domain;
                Vector3I count(6, 7, 8);
                count[axis] = cells;
                domain.init(count, {0.7, 1.1, 0.9}, ghost);
                BoundaryConditionHandler bc = boundary_conditions_from_config(
                    {{"Boundary_conditions",
                      nlohmann::json::parse(
                          R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}},
                    domain);
                Field3d raw(domain.size(), 3), field(domain.size(), 3);
                raw.setZero();
                std::vector<std::pair<int, double>> weights;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        for (int k = 0; k < 3; ++k)
                            for (int c = 0; c < 3; ++c) {
                                const int row = domain.vind(ghost - 1 + i, ghost - 1 + j, ghost - 1 + k, c);
                                const double weight = (1 + c) * (1 + i) * (2 + j) * (3 + k) / 1000.0;
                                raw[row] = weight;
                                weights.emplace_back(row, weight);
                            }
                // Independent periodic extension and expected raw particle gather.
                for (int i = 0; i < domain.size().x(); ++i)
                    for (int j = 0; j < domain.size().y(); ++j)
                        for (int k = 0; k < domain.size().z(); ++k) {
                            const int ix = ((i - ghost) % count.x() + count.x()) % count.x();
                            const int iy = ((j - ghost) % count.y() + count.y()) % count.y();
                            const int iz = ((k - ghost) % count.z() + count.z()) % count.z();
                            for (int c = 0; c < 3; ++c)
                                field(i, j, k, c) = std::sin(0.3 * ix + 0.7 * iy - 0.2 * iz + c);
                        }
                double gather = 0;
                for (const auto& [row, weight] : weights) gather += weight * field[row];
                std::vector<Trip> entries;
                for (const auto& [row, a] : weights)
                    for (const auto& [col, b] : weights) entries.emplace_back(row, col, a * b);
                Operator response(field.capacity(), field.capacity());
                response.setFromTriplets(entries.begin(), entries.end());
                Field3d current = raw;
                bc.apply_to_fields(current, FieldType::CURRENT, domain);
                bc.apply_to_operator(response, domain);
                const Field3d result = response * field;
                check::near((result - gather * current).norm(), 0, 2e-12,
                            "periodic response equals scatter(gather(E)) with overlapping ghost layers");
                const auto range = bc.active_range(domain.grid);
                check::near(dot_product_sum(field, result, range), gather * gather, 2e-12,
                            "periodic matrix preserves physical particle response work");
                check::near(dot_product_sum(field, current, range), gather, 2e-13,
                            "each raw ghost deposit contributes once to physical current work");
            }
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
        const double removed = damping_fields(e, b, domain, damping_parameters_from_config(cfg));
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
    FieldBoundaryConditions boundaries;
    MaxwellOperators maxwell;
    maxwell.init(assemble_maxwell_stencils(domain, boundaries));
    Field3d E(domain.size(), 3), En(domain.size(), 3), B(domain.size(), 3), Bn(domain.size(), 3);
    B.data().setConstant(0.3);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k)
                for (int d = 0; d < 3; ++d) {
                    E(i, j, k, d) = domain.storage_node_position(i, j, k, FieldType::ELECTRIC, d)[(d + 1) % 3];
                }
    // Distinct endpoints detect using E or En alone in place of their average.
    En = 3.0 * E;
    maxwell.advance_magnetic_field(Bn, B, E, En, 0.2);
    for (int i = GHOST_CELLS + 1; i < domain.size().x() - GHOST_CELLS - 1; ++i)
        for (int j = GHOST_CELLS + 1; j < domain.size().y() - GHOST_CELLS - 1; ++j)
            for (int k = GHOST_CELLS + 1; k < domain.size().z() - GHOST_CELLS - 1; ++k)
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
    FieldBoundaryConditions boundaries;
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
        check::scenario(
            "Periodic density/current and response preserve alias sums for the configured halo and short periods",
            periodic_overlap_counts_each_deposit_once);
        check::scenario("Three short periodic axes preserve charge and all halo-corner aliases",
                        periodic_corner_overlap_preserves_total);
        check::scenario("Short-period full-storage solve preserves aliases and physical response work",
                        short_periodic_response_solve_matches_owner_work);
        check::scenario("Damping reports the actual decrease of E/B energy", damping_accounts_for_removed_field_energy);
        check::scenario("Periodic response matches scatter(gather(E)) on thin periodic grids",
                        small_periodic_response);
    });
}
