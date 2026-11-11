#include <cmath>
#include <filesystem>
#include <random>
#include <vector>

#include "diagnostics/DiagnosticMeasurements.h"
#include "simulation/drift_kinetic/drift_kinetic_adaptive.h"
#include "simulation/drift_kinetic/drift_kinetic_particles.h"
#include "simulation/drift_kinetic/drift_kinetic_push.h"
#include "simulation/drift_kinetic/drift_kinetic_yee.h"
#include "simulation/drift_kinetic/simulation_drift_kinetic.h"
#include "simulation/eccapfim/algorithms_eccapfim.h"
#include "test_support.h"

namespace {
using driftKinetic::DriftKineticParticles;
using driftKinetic::Fields;
using driftKinetic::GuidingCenter;
using driftKinetic::GuidingCenterPush;
using driftKinetic::MagneticLevels;
using driftKinetic::YeeCoupling;

Fields uniform_fields(const Vector3R& electric, const Vector3R& magnetic) {
    return {electric, magnetic.length(), magnetic.normalized(), Vector3R(0, 0, 0), Vector3R(0, 0, 0)};
}

GuidingCenter particle(const Vector3R& position, double v_parallel, double mu, double charge = 1.0, double mass = 1.0) {
    GuidingCenter result;
    result.position = position;
    result.v_parallel = v_parallel;
    result.mu = mu;
    result.charge = charge;
    result.mass = mass;
    return result;
}

void uniform_parallel_motion() {
    const Vector3R magnetic(0, 0, 1);
    const auto fields = [&](const Vector3R&, const Vector3R&) { return uniform_fields(Vector3R(0, 0, 0), magnetic); };
    GuidingCenterPush push;
    push.set_tolerances(1e-14, 1e-13, 100);
    const GuidingCenter initial = particle(Vector3R(0.1, 0.2, 0.3), 0.4, 0.2);
    const auto result = push.push(initial, 0.3, fields, true);
    check::require(result.converged, "uniform-B guiding-center push converged");
    check::near(result.state.position.x(), initial.position.x(), 1e-14, "uniform-B keeps x");
    check::near(result.state.position.y(), initial.position.y(), 1e-14, "uniform-B keeps y");
    check::near(result.state.position.z(), initial.position.z() + 0.4 * 0.3, 1e-13, "parallel drift");
    check::near(result.state.v_parallel, 0.4, 1e-13, "parallel velocity is constant");
    check::near(result.state.mu, 0.2, 0.0, "magnetic moment is constant");
}

void parallel_electric_acceleration() {
    const Vector3R magnetic(0, 0, 1), electric(0, 0, 0.5);
    const auto fields = [&](const Vector3R&, const Vector3R&) { return uniform_fields(electric, magnetic); };
    GuidingCenterPush push;
    push.set_tolerances(1e-14, 1e-13, 100);
    const GuidingCenter initial = particle(Vector3R(0, 0, 0), 0.0, 0.1);
    const auto result = push.push(initial, 0.3, fields, true);
    check::require(result.converged, "parallel-E guiding-center push converged");
    check::near(result.state.v_parallel, 0.5 * 0.3, 1e-13, "parallel acceleration q E/m");
    check::near(result.state.position.z(), 0.5 * 0.5 * 0.3 * 0.3, 1e-13, "implicit midpoint displacement");
}

void electric_cross_magnetic_drift() {
    const Vector3R magnetic(0, 0, 2), electric(0.5, 0, 0);
    const auto fields = [&](const Vector3R&, const Vector3R&) { return uniform_fields(electric, magnetic); };
    GuidingCenterPush push;
    push.set_tolerances(1e-14, 1e-13, 100);
    const GuidingCenter initial = particle(Vector3R(0, 0, 0), 0.0, 0.0);
    const auto result = push.push(initial, 0.3, fields, false);
    check::require(result.converged, "E-cross-B guiding-center push converged");
    check::near(result.state.position.x(), 0.0, 1e-14, "E-cross-B has no x displacement");
    check::near(result.state.position.y(), -0.25 * 0.3, 1e-13, "E-cross-B drift");
    check::near(result.state.position.z(), 0.0, 1e-14, "E-cross-B has no parallel displacement");
    check::near(result.state.v_parallel, 0.0, 1e-13, "E-cross-B keeps parallel velocity");
}

void mirror_energy_invariant() {
    const auto mirror = [&](const Vector3R& start, const Vector3R& end) {
        const double z = 0.5 * (start.z() + end.z());
        const double length = 1.0 + 0.2 * z;
        return Fields{Vector3R(0, 0, 0), length, Vector3R(0, 0, 1), Vector3R(0, 0, 0.2), Vector3R(0, 0, 0)};
    };
    GuidingCenterPush push;
    push.set_tolerances(1e-14, 1e-13, 100);
    GuidingCenter state = particle(Vector3R(0, 0, 0), 1.0, 0.5);
    const auto energy = [&](const GuidingCenter& p) {
        const double length = 1.0 + 0.2 * p.position.z();
        return 0.5 * p.mass * p.v_parallel * p.v_parallel + p.mu * length;
    };
    const double initial_energy = energy(state);
    for (int step = 0; step < 20; ++step) {
        const auto result = push.push(state, 0.01, mirror, true);
        check::require(result.converged, "mirror guiding-center push converged");
        state = result.state;
        check::near(energy(state), initial_energy, 1e-8, "mirror motion conserves parallel-plus-mu-B energy");
    }
    check::near(state.mu, 0.5, 0.0, "mirror magnetic moment is constant");
}

void adaptive_substep_mirror() {
    // A mirror field whose |B| depends on the interval: the leaves must use
    // their own time levels and together cover the whole step.
    const auto mirror = [&](const Vector3R& start, const Vector3R& end, double, double) {
        const double z = 0.5 * (start.z() + end.z());
        return Fields{Vector3R(0, 0, 0), 1.0 + 0.2 * z, Vector3R(0, 0, 1), Vector3R(0, 0, 0.2), Vector3R(0, 0, 0)};
    };
    GuidingCenterPush push;
    push.set_tolerances(1e-14, 1e-13, 1);   // too few iterations: forces bisection
    driftKinetic::AdaptiveSubstepping adaptive;
    adaptive.configure(true, 6);
    std::vector<driftKinetic::Leaf> leaves;
    driftKinetic::PushStatistics stats;
    const GuidingCenter initial = particle(Vector3R(0, 0, 0), 1.0, 0.5);
    check::require(adaptive.advance(push, initial, initial, 0.2, mirror, false, leaves, stats), "valid adaptive push");
    check::require(stats.retries > 0 && leaves.size() > 1, "non-converged root step is bisected");
    double covered = 0.0;
    for (size_t i = 0; i < leaves.size(); ++i) {
        covered += leaves[i].end - leaves[i].begin;
        if (i > 0) {
            check::near(leaves[i].begin, leaves[i - 1].end, 0.0, "leaves are contiguous in time");
            check::vector(leaves[i].start.position, leaves[i - 1].finish.position, 0.0, "leaves are chained");
        }
    }
    check::near(covered, 1.0, 1e-15, "leaves cover the step");
    bool rejected = false;
    try {
        adaptive.configure(true, 31);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "adaptive substep depth is validated");
}

// ---- Yee coupling kernels -------------------------------------------------

Field3d zero_field(const Domain& domain, int components) {
    Field3d field(domain.size(), components);
    field.setZero();
    return field;
}

Domain box(const Vector3I& cells, const Vector3R& h) {
    Domain domain;
    domain.init(cells, h);
    return domain;
}

// Fills component c from f(world position of that Yee component).
template <typename F>
void fill_magnetic(Field3d& field, const Domain& domain, int c, F&& f) {
    const auto shift = Vector3R(c == 0 ? 0.0 : 0.5, c == 1 ? 0.0 : 0.5, c == 2 ? 0.0 : 0.5);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                const Vector3R r((i - GHOST_CELLS + shift.x()) * domain.cell_size(0),
                                 (j - GHOST_CELLS + shift.y()) * domain.cell_size(1),
                                 (k - GHOST_CELLS + shift.z()) * domain.cell_size(2));
                field(i, j, k, c) = f(r);
            }
}

void random_field(Field3d& field, const Domain& domain, double offset, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-1, 1);
    for (int p = 0; p < field.capacity(); ++p) field[p] = u(rng);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) field(i, j, k, 2) += offset;
    algorithmsECCAPFIM::copy_periodic_nodes(field, domain);
}

double physical_dot(const Field3d& a, const Field3d& b, const Domain& domain) {
    const int g = GHOST_CELLS;
    const auto n = domain.num_cells();
    double result = 0;
    for (int i = g; i < g + n.x(); ++i)
        for (int j = g; j < g + n.y(); ++j)
            for (int k = g; k < g + n.z(); ++k)
                for (int c = 0; c < 3; ++c) result += a(i, j, k, c) * b(i, j, k, c);
    return result;
}

void uniform_field_interpolation() {
    const Domain domain = box({6, 6, 6}, {0.7, 1.1, 0.9});
    Field3d E = zero_field(domain, 3), B = zero_field(domain, 3);
    for (int c = 0; c < 3; ++c) {
        fill_magnetic(E, domain, c, [&](const Vector3R&) { return 1.0 + c; });
        fill_magnetic(B, domain, c, [&](const Vector3R&) { return c == 2 ? 2.0 : 0.0; });
    }
    const YeeCoupling coupling(domain, nullptr);
    const MagneticLevels levels{&B, &B, &B};
    const auto fields = coupling.interpolate(E, levels, {1.3, 2.9, 1.7}, {2.6, 1.4, 2.2}, {0, 1});
    check::vector(fields.electric, Vector3R(1, 2, 3), 1e-14, "uniform electric orbit average");
    check::near(fields.lenB, 2.0, 1e-14, "uniform |B|");
    check::vector(fields.b, Vector3R(0, 0, 1), 1e-14, "uniform b");
    check::vector(fields.gradB, Vector3R(0, 0, 0), 1e-14, "uniform grad|B|");
    check::vector(fields.rotB, Vector3R(0, 0, 0), 1e-14, "uniform curl B");
}

void linear_magnetic_derivatives() {
    const Domain domain = box({12, 12, 12}, {0.7, 1.1, 0.9});
    const YeeCoupling coupling(domain, nullptr);
    const Vector3R start(3.9, 5.3, 4.6), end(4.8, 6.9, 5.7);   // away from the periodic seam
    {
        // B = (0, 0, 1 + 0.2 z): grad|B| = (0, 0, 0.2) and |B| linear along the orbit.
        Field3d B = zero_field(domain, 3);
        fill_magnetic(B, domain, 2, [](const Vector3R& r) { return 1.0 + 0.2 * r.z(); });
        const MagneticLevels levels{&B, &B, &B};
        Field3d E = zero_field(domain, 3);
        const auto fields = coupling.interpolate(E, levels, start, end, {0, 1});
        check::vector(fields.gradB, Vector3R(0, 0, 0.2), 1e-13, "linear grad|B| is exact");
        check::near(fields.lenB, 1.0 + 0.1 * (start.z() + end.z()), 1e-13, "|B| at the orbit ends");
    }
    {
        // B = (0.4 z, 0.3 x, 1): curl B = (0, 0.4, 0.3).
        Field3d B = zero_field(domain, 3);
        fill_magnetic(B, domain, 0, [](const Vector3R& r) { return 0.4 * r.z(); });
        fill_magnetic(B, domain, 1, [](const Vector3R& r) { return 0.3 * r.x(); });
        fill_magnetic(B, domain, 2, [](const Vector3R&) { return 1.0; });
        const MagneticLevels levels{&B, &B, &B};
        Field3d E = zero_field(domain, 3);
        const auto fields = coupling.interpolate(E, levels, start, end, {0, 1});
        check::vector(fields.rotB, Vector3R(0, 0.4, 0.3), 1e-13, "linear curl B is exact");
    }
}

// Orbits across cell planes, corners and the periodic seam.
const std::vector<std::pair<Vector3R, Vector3R>> orbits = {
    {{0.2, 0.3, 0.4}, {0.2, 0.3, 0.4}},  {{-0.35, 0.55, 0.2}, {0.35, -0.55, -0.2}},
    {{2.4, 6.5, 2.7}, {2.6, 6.7, 2.95}}, {{0.8, 1.2, 0.7}, {3.4, -2.5, 3.7}},
    {{3.4, -2.5, 3.7}, {0.8, 1.2, 0.7}}, {{0.7, 1.1, 0.45}, {0.7, 1.1, 1.35}}};

void work_identity() {
    // sum_g E_g . J_g = q V . E_p: the current does exactly the work of the
    // field interpolated for the push.
    const Domain domain = box({5, 6, 7}, {0.7, 1.1, 0.9});
    Field3d E = zero_field(domain, 3), B = zero_field(domain, 3);
    random_field(E, domain, 0, 1);
    random_field(B, domain, 4, 2);
    const YeeCoupling coupling(domain, nullptr);
    const MagneticLevels levels{&B, &B, &B};
    for (const auto& [start, end] : orbits) {
        Field3d J = zero_field(domain, 3), M = zero_field(domain, 3);
        const Vector3R velocity(0.3, -0.7, 1.1);
        coupling.deposit(J, M, levels, start, end, velocity, 1.7, 0.0, {0, 1});
        const auto fields = coupling.interpolate(E, levels, start, end, {0, 1});
        check::near(physical_dot(E, J, domain), 1.7 * velocity.dot(fields.electric), 1e-13,
                    "guiding-centre current work equals interpolated field work");
    }
}

void prepared_differences_match_reference() {
    // The optimized kernel reads precomputed d_a B_c; it must reproduce the
    // on-the-fly differences, including the blended level of a substep.
    const Domain domain = box({5, 6, 7}, {0.7, 1.1, 0.9});
    Field3d E = zero_field(domain, 3), B0 = zero_field(domain, 3), Bh = zero_field(domain, 3),
            B1 = zero_field(domain, 3);
    random_field(E, domain, 0, 5);
    random_field(B0, domain, 3, 6);
    random_field(B1, domain, 3, 7);
    Bh.data() = 0.5 * (B0.data() + B1.data());
    YeeCoupling coupling(domain, nullptr);
    const MagneticLevels levels{&B0, &Bh, &B1};
    coupling.prepare(levels);
    for (const auto& interval :
         {driftKinetic::StepInterval{0, 1}, driftKinetic::StepInterval{0.25, 0.5}, driftKinetic::StepInterval{0, 0.5}})
        for (const auto& [start, end] : orbits) {
            const Vector3R B_start = coupling.start_magnetic(levels, start, interval);
            const auto fast = coupling.interpolate(E, levels, start, end, interval, &B_start);
            const auto reference = coupling.interpolate_reference(E, levels, start, end, interval);
            check::vector(fast.electric, reference.electric, 1e-13, "optimized E equals reference");
            check::vector(fast.gradB, reference.gradB, 1e-12, "optimized grad|B| equals reference");
            check::vector(fast.rotB, reference.rotB, 1e-12, "optimized curl B equals reference");
            check::near(fast.lenB, reference.lenB, 1e-14, "optimized |B| equals reference");
        }
}

void charge_continuity() {
    // (rho^{n+1} - rho^n) / dt + div J = 0 with the quadratic nodal charge.
    const Domain domain = box({5, 6, 7}, {0.7, 1.1, 0.9});
    Field3d B = zero_field(domain, 3);
    random_field(B, domain, 4, 3);
    const YeeCoupling coupling(domain, nullptr);
    const MagneticLevels levels{&B, &B, &B};
    const double dt = 0.35, q = -1.3;
    const int g = GHOST_CELLS;
    for (const auto& [start, end] : orbits) {
        Field3d J = zero_field(domain, 3), M = zero_field(domain, 3), before = zero_field(domain, 1),
                after = zero_field(domain, 1);
        coupling.deposit(J, M, levels, start, end, (end - start) / dt, q, 0.0, {0, 1});
        algorithmsECCAPFIM::deposit_charge(before, start, q, domain);
        algorithmsECCAPFIM::deposit_charge(after, end, q, domain);
        const auto n = domain.num_cells();
        for (int i = 0; i < n.x(); ++i)
            for (int j = 0; j < n.y(); ++j)
                for (int k = 0; k < n.z(); ++k) {
                    const Vector3I node(i + g, j + g, k + g);
                    double divergence = 0;
                    for (int d = 0; d < 3; ++d) {
                        auto prev = node;
                        prev[d] = algorithmsECCAPFIM::periodic_node(node[d] - g - 1, d, domain);
                        divergence += (J(node.x(), node.y(), node.z(), d) - J(prev.x(), prev.y(), prev.z(), d)) /
                                      domain.cell_size(d);
                    }
                    check::near(after(node.x(), node.y(), node.z(), 0) - before(node.x(), node.y(), node.z(), 0) +
                                    dt * divergence,
                                0, 1e-13, "guiding-centre charge continuity");
                }
    }
}

void magnetization_identity() {
    // sum_g M_g . B_g = -mu hat_b . (B(r0) + B(r1)) / 2, the moment energy of the orbit ends.
    const Domain domain = box({5, 6, 7}, {0.7, 1.1, 0.9});
    Field3d B = zero_field(domain, 3);
    random_field(B, domain, 4, 4);
    const YeeCoupling coupling(domain, nullptr);
    const MagneticLevels levels{&B, &B, &B};
    for (const auto& [start, end] : orbits) {
        Field3d J = zero_field(domain, 3), M = zero_field(domain, 3);
        coupling.deposit(J, M, levels, start, end, Vector3R(0, 0, 0), 0.0, 0.6, {0, 1});
        const Vector3R B0 = coupling.magnetic_at(B, start), B1 = coupling.magnetic_at(B, end);
        const Vector3R b = 0.5 * (B0.normalized() + B1.normalized());
        const Vector3R hat_b = b / b.squaredNorm();
        check::near(physical_dot(M, B, domain), -0.6 * 0.5 * hat_b.dot(B0 + B1), 1e-13,
                    "magnetization is deposited with the interpolation weights of B");
    }
}

// ---- Self-consistent schemes ---------------------------------------------

nlohmann::json settings(const std::string& scheme) {
    return {{"Scheme", scheme},
            {"Dt", 0.35},
            {"NumCellsX", 4},
            {"NumCellsY", 5},
            {"NumCellsZ", 3},
            {"Dx", 0.7},
            {"Dy", 1.1},
            {"Dz", 0.9},
            {"StartFromTime", 0},
            {"n0", 1e10},
            {"k_particles_reservation", 0},
            {"Collisions", nlohmann::json::array()},
            {"diagnostics", {{"outputs", nlohmann::json::array()}}},
            {"DriftKinetic", {{"magnetic", {0.1, 0, 1}}}},
            {"Eccapfim", {{"preconditioner", "mass_matrix"}, {"coupled_atol", 1e-13}, {"coupled_rtol", 1e-13}}},
            {"Boundary_conditions",
             nlohmann::json::parse(
                 R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"periodic":{"face":"ZMIN"}}])")}};
}

std::unique_ptr<SimulationDriftKinetic> simulation(const nlohmann::json& config) {
    auto electrons = check::species("Electrons", -1, 1, 1, 1);
    electrons["model"] = "drift_kinetic";
    auto ions = check::species("Ions", 1, 9, 1, 1);
    ions["model"] = "kinetic";
    const nlohmann::json particles{{"particles", nlohmann::json::array({electrons, ions})}};
    char name[] = "drift_kinetic_invariants";
    char* argv[] = {name};
    auto base = build_simulation(config, particles, 1, argv);
    auto* concrete = dynamic_cast<SimulationDriftKinetic*>(base.get());
    check::require(concrete != nullptr, "drift-kinetic scheme factory dispatch");
    base.release();
    std::unique_ptr<SimulationDriftKinetic> result(concrete);
    result->init();
    return result;
}

double total_energy(SimulationDriftKinetic& sim) {
    const auto range = sim.bc_handler.active_range(sim.domain.grid);
    double result = 0.5 * dot_product_sum(sim.fieldEn, sim.fieldEn, range) +
                    0.5 * dot_product_sum(sim.fieldBn, sim.fieldBn, range) +
                    sim.guiding_center_energy(sim.fieldBn + sim.fieldBInit);
    for (const auto& [name, sp] : sim.species) result += sp->get_kinetic_energy();
    return result;
}

void guiding_charge(SimulationDriftKinetic& sim, Field3d& density) {
    density.resize(sim.domain.size(), 1);
    density.setZero();
    auto& electrons = sim.guiding_centers("Electrons");
    for (int cell = 0; cell < electrons.size(); ++cell)
        for (const auto& p : electrons(cell))
            algorithmsECCAPFIM::deposit_charge(density, p.position, p.charge * p.weight, sim.domain, &sim.bc_handler);
}

void add_guiding_centers(SimulationDriftKinetic& sim) {
    auto& electrons = sim.guiding_centers("Electrons");
    const std::vector<std::pair<Vector3R, double>> markers = {
        {{2.79, 5.49, 2.69}, 0.6}, {{0.35, 1.1, 0.45}, -0.4}, {{1.5, 2.3, 1.4}, 0.2}, {{2.1, 4.0, 0.1}, -0.7}};
    for (const auto& [r, v] : markers) {
        GuidingCenter p;
        p.position = r;
        p.v_parallel = v;
        p.mu = 0.05 + 0.1 * std::abs(v);
        p.charge = electrons.charge();
        p.mass = electrons.mass();
        p.weight = electrons.weight();
        electrons.add(p);
    }
}

void perturb_fields(SimulationDriftKinetic& sim) {
    // Divergence-free self field: gradients, curvature and induction act on the centres.
    const auto size = sim.domain.size();
    for (int i = 0; i < size.x(); ++i)
        for (int j = 0; j < size.y(); ++j)
            for (int k = 0; k < size.z(); ++k) {
                sim.fieldEn(i, j, k, 1) = 0.02 * std::sin(2 * M_PI * (i - GHOST_CELLS) / 4);
                sim.fieldBn(i, j, k, 2) = 0.08 * std::cos(2 * M_PI * (i - GHOST_CELLS) / 4);
                sim.fieldBn(i, j, k, 1) = 0.05 * std::sin(2 * M_PI * (k - GHOST_CELLS) / 3);
            }
    algorithmsECCAPFIM::copy_periodic_nodes(sim.fieldEn, sim.domain);
    algorithmsECCAPFIM::copy_periodic_nodes(sim.fieldBn, sim.domain);
}

void coupled_energy_and_charge(const std::string& scheme, const std::string& solver, bool substepping) {
    auto config = settings(scheme);
    config["Eccapfim"]["coupled_solver"] = solver;
    if (substepping)
        config["DriftKineticPush"] = {
            {"solver", "picard"}, {"substepping", true}, {"max_substep_depth", 8}, {"max_iterations", 6}};
    auto sim = simulation(config);
    perturb_fields(*sim);
    add_guiding_centers(*sim);
    if (scheme == "hybrid_kinetic")
        for (const auto& r : {Vector3R(2.79, 5.49, 2.69), Vector3R(0.35, 1.1, 0.45), Vector3R(1.5, 2.3, 1.4)})
            sim->species.at("Ions")->add_particle(Particle(r, {-0.1, 0.2, 0.05}));
    const double initial = total_energy(*sim);
    const auto range = sim->bc_handler.active_range(sim->domain.grid);
    long retries = 0;
    for (int step = 1; step <= 6; ++step) {
        Field3d before, after;
        guiding_charge(*sim, before);
        const double energy_before = sim->guiding_center_energy(sim->fieldBn + sim->fieldBInit);
        sim->prepare_step(step);
        sim->make_step(step);
        guiding_charge(*sim, after);
        // Pointwise continuity of the guiding-centre charge with the free current.
        const int g = GHOST_CELLS;
        const auto n = sim->domain.num_cells();
        double worst = 0;
        for (int i = 0; i < n.x(); ++i)
            for (int j = 0; j < n.y(); ++j)
                for (int k = 0; k < n.z(); ++k) {
                    const Vector3I node(i + g, j + g, k + g);
                    double divergence = 0;
                    for (int d = 0; d < 3; ++d) {
                        auto prev = node;
                        prev[d] = algorithmsECCAPFIM::periodic_node(node[d] - g - 1, d, sim->domain);
                        divergence += (sim->fieldJ_guiding(node.x(), node.y(), node.z(), d) -
                                       sim->fieldJ_guiding(prev.x(), prev.y(), prev.z(), d)) /
                                      sim->domain.cell_size(d);
                    }
                    worst = std::max(worst, std::abs(after(node.x(), node.y(), node.z(), 0) -
                                                     before(node.x(), node.y(), node.z(), 0) + sim->dt_ * divergence));
                }
        check::near(worst, 0, 1e-12, "guiding-centre charge continuity in the coupled step");
        // Guiding-centre energy theorem dK = dt E . J_gc - M . dB, independent of the field solve.
        const double work = sim->dt_ * dot_product_sum(sim->fieldEp + sim->fieldE_external, sim->fieldJ_guiding, range);
        const double moment_work = dot_product_sum(sim->fieldM, sim->fieldBn - sim->fieldB, range);
        const double energy_after = sim->guiding_center_energy(sim->fieldBn + sim->fieldBInit);
        check::near(energy_after - energy_before - work + moment_work, 0, 1e-11, "guiding-centre energy theorem");
        check::near(total_energy(*sim), initial, 1e-10, "particle + field energy conservation");
        check::require(sim->nonlinear_residual < 1e-11, "accepted coupled field residual");
        retries += sim->push_statistics("Electrons").retries;
        sim->make_diagnostic(step);
    }
    if (substepping)
        check::require(retries > 0, "substepping was exercised");
    check::require(sim->guiding_centers("Electrons").get_total_num_of_particles() == 4, "guiding centres kept");
    sim->finalize();
}

// Grid moments of guiding centres for the diagnostic outputs, checked by
// conservation sums: total charge, the Maxwell current, and the energy as half
// the trace of the second moment.
void guiding_center_moments() {
    auto sim = simulation(settings("hybrid_kinetic"));
    perturb_fields(*sim);
    add_guiding_centers(*sim);
    for (const auto& r : {Vector3R(2.79, 5.49, 2.69), Vector3R(0.35, 1.1, 0.45)})
        sim->species.at("Ions")->add_particle(Particle(r, {-0.1, 0.2, 0.05}));
    for (int step = 1; step <= 2; ++step) {
        sim->prepare_step(step);
        sim->make_step(step);
    }
    const MomentSources* sources = sim->extra_moment_sources();
    check::require(sources && sources->size() == 1 && sources->contains("Electrons"), "guiding-centre moment source");
    const SpeciesMomentSource& electrons = *sources->at("Electrons");
    const auto& markers = sim->guiding_centers("Electrons");
    double charge = 0;
    for (int cell = 0; cell < markers.size(); ++cell)
        for (const auto& p : markers(cell)) charge += p.charge * p.weight;
    const auto n = sim->domain.num_cells();
    const int g = GHOST_CELLS;

    // Charge: the folded density (marker charge per node, as for kinetic species)
    // sums to the total marker charge.
    for (ShapeType shape : {ShapeType::Linear, ShapeType::Quadratic}) {
        Field3d density(sim->domain.size(), 1);
        electrons.deposit_density(density, shape);
        sim->bc_handler.apply_to_fields(density, FieldType::DENSITY, sim->domain);
        double sum = 0;
        for (int i = 0; i < n.x(); ++i)
            for (int j = 0; j < n.y(); ++j)
                for (int k = 0; k < n.z(); ++k) sum += density(i + g, j + g, k + g, 0);
        check::near(sum, charge, 1e-13, "guiding-centre density sums to the marker charge");
    }

    // Current: kinetic plus guiding-centre species reproduce the Maxwell current.
    DiagnosticMeasurements measured({sim->domain, sim->species, sim->fieldEn, sim->fieldBn, sim->fieldE_external,
                                     sim->fieldBInit, sim->bc_handler, ShapeType::Linear, sources});
    measured.request_density("Electrons");
    measured.calculate();
    double worst = 0, scale = 0;
    for (int i = 0; i < n.x(); ++i)
        for (int j = 0; j < n.y(); ++j)
            for (int k = 0; k < n.z(); ++k)
                for (int c = 0; c < 3; ++c) {
                    const double total = measured.current("Ions")(i + g, j + g, k + g, c) +
                                         measured.current("Electrons")(i + g, j + g, k + g, c);
                    worst = std::max(worst, std::abs(total - sim->fieldJ(i + g, j + g, k + g, c)));
                    scale = std::max(scale, std::abs(sim->fieldJ(i + g, j + g, k + g, c)));
                }
    check::require(scale > 0 && worst <= 1e-13 * scale, "species currents sum to the Maxwell current");
    check::require(measured.density("Electrons").capacity() > 0, "measured guiding-centre density");

    // Second moment: half its trace sums to m v_par^2 / 2 + mu |B| per marker.
    double trace = 0;
    for (int a = 0; a < 3; ++a) {
        Field3d moment(sim->domain.size(), 1);
        electrons.deposit_second_moment(moment, a, a);
        for (int p = 0; p < moment.capacity(); ++p) trace += moment(p);
    }
    const double energy = sim->guiding_center_energy(sim->fieldBn + sim->fieldBInit);
    check::near(0.5 * trace, energy, 1e-13 * energy, "half the second-moment trace is the guiding-centre energy");

    DiagnosticMeasurements spectrum({sim->domain, sim->species, sim->fieldEn, sim->fieldBn, sim->fieldE_external,
                                     sim->fieldBInit, sim->bc_handler, ShapeType::Linear, sources});
    spectrum.request_spectrum("Electrons");
    bool rejected = false;
    try {
        spectrum.calculate();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "energy spectrum is rejected for guiding centres");
    sim->finalize();
}

// Restart of a hybrid run from the recovery files: the next step of the
// resumed run reproduces the uninterrupted one.
void guiding_center_restart() {
    auto config = settings("hybrid_kinetic");
    config["diagnostics"]["outputs"] = nlohmann::json::array({{{"type", "recovery"}, {"interval", 1}}});
    auto original = simulation(config);
    perturb_fields(*original);
    add_guiding_centers(*original);
    original->species.at("Ions")->add_particle(Particle({2.79, 5.49, 2.69}, {-0.1, 0.2, 0.05}));
    original->prepare_step(1);
    original->make_step(1);
    original->make_diagnostic(1);
    check::require(std::filesystem::exists("Recovery/Particles/Electrons/Electrons.gc.backup"),
                   "guiding-centre recovery file written");

    config["StartFromTime"] = 1;
    config["diagnostics"]["outputs"] = nlohmann::json::array();
    const auto directory = std::filesystem::current_path();
    std::filesystem::create_directories("restart");
    std::filesystem::current_path("restart");   // the recovery reader uses ../Recovery
    auto resumed = simulation(config);
    std::filesystem::current_path(directory);
    check::require(resumed->guiding_centers("Electrons").get_total_num_of_particles() == 4,
                   "restart guiding-centre count");
    check::near((resumed->fieldEn - original->fieldEn).norm(), 0, 0, "restart E");
    check::near((resumed->fieldBn - original->fieldBn).norm(), 0, 0, "restart B");

    for (auto* sim : {original.get(), resumed.get()}) {
        sim->prepare_step(2);
        sim->make_step(2);
    }
    check::near((resumed->fieldEn - original->fieldEn).norm(), 0, 1e-12, "restart next electric field");
    check::near((resumed->fieldBn - original->fieldBn).norm(), 0, 1e-12, "restart next magnetic field");
    check::near(total_energy(*resumed), total_energy(*original), 1e-12, "restart total energy");
    const auto& expected = original->guiding_centers("Electrons");
    const auto& actual = resumed->guiding_centers("Electrons");
    for (int cell = 0; cell < expected.size(); ++cell) {
        check::require(actual(cell).size() == expected(cell).size(), "restart guiding centres per cell");
        for (size_t n = 0; n < expected(cell).size(); ++n) {
            check::vector(actual(cell)[n].position, expected(cell)[n].position, 1e-12, "restart next position");
            check::near(actual(cell)[n].v_parallel, expected(cell)[n].v_parallel, 1e-12, "restart next v_par");
            check::near(actual(cell)[n].mu, expected(cell)[n].mu, 0, "restart magnetic moment");
        }
    }
    original->finalize();
    resumed->finalize();
}

void prescribed_drift() {
    auto config = settings("drift_kinetic");
    config["DriftKinetic"] = {{"magnetic", {0, 0, 2}}, {"electric", {0.5, 0, 0}}, {"evolve_fields", false}};
    auto sim = simulation(config);
    auto& electrons = sim->guiding_centers("Electrons");
    GuidingCenter p;
    p.position = Vector3R(1.0, 2.0, 1.0);
    p.v_parallel = 0.3;
    p.charge = electrons.charge();
    p.mass = electrons.mass();
    p.weight = electrons.weight();
    electrons.add(p);
    sim->prepare_step(1);
    sim->make_step(1);
    for (int cell = 0; cell < electrons.size(); ++cell)
        for (const auto& q : sim->guiding_centers("Electrons")(cell)) {
            // E x B / B^2 = (0, -0.25, 0); a uniform E perpendicular to b does no work.
            check::vector(q.position, Vector3R(1.0, 2.0 - 0.25 * 0.35, 1.0 + 0.3 * 0.35), 1e-13,
                          "prescribed-field E x B and parallel motion");
            check::near(q.v_parallel, 0.3, 1e-14, "no parallel acceleration");
        }
    sim->finalize();
}

void particle_container_rebin() {
    const Domain domain = box({4, 4, 4}, {1, 1, 1});
    DriftKineticParticles particles(domain, "electrons", -1.0, 1.0, 2.0);
    GuidingCenter p;
    p.position = Vector3R(0.5, 0.5, 0.5);
    particles.add(p);
    const auto cell = domain.storage_cell_index(p.position);
    particles(cell.x(), cell.y(), cell.z())[0].position = Vector3R(2.5, 1.5, 1.5);
    particles.update_cells();
    const auto target = domain.storage_cell_index(Vector3R(2.5, 1.5, 1.5));
    check::require(particles(cell.x(), cell.y(), cell.z()).empty(), "old cell is emptied");
    check::require(particles(target.x(), target.y(), target.z()).size() == 1, "particle rebins to the new cell");
}
}   // namespace

int main() {
    return check::run([] {
        check::scenario("Guiding-center uniform-B parallel motion", uniform_parallel_motion);
        check::scenario("Guiding-center parallel electric acceleration", parallel_electric_acceleration);
        check::scenario("Guiding-center E-cross-B drift", electric_cross_magnetic_drift);
        check::scenario("Guiding-center mirror energy invariant", mirror_energy_invariant);
        check::scenario("Guiding-center adaptive substeps", adaptive_substep_mirror);
        check::scenario("Guiding-center particle container rebin", particle_container_rebin);
        check::scenario("Yee coupling uniform fields", uniform_field_interpolation);
        check::scenario("Yee coupling linear grad|B| and curl B", linear_magnetic_derivatives);
        check::scenario("Yee coupling current work identity", work_identity);
        check::scenario("Yee coupling charge continuity", charge_continuity);
        check::scenario("Yee coupling prepared differences", prepared_differences_match_reference);
        check::scenario("Yee coupling magnetization identity", magnetization_identity);
        check::scenario("Drift-kinetic prescribed fields", prescribed_drift);
        check::scenario("Guiding-centre grid moments for diagnostics", guiding_center_moments);
        check::scenario("Guiding-centre restart from recovery files", guiding_center_restart);
        check::scenario("Drift-kinetic coupled Anderson",
                        [] { coupled_energy_and_charge("drift_kinetic", "anderson", false); });
        check::scenario("Drift-kinetic coupled NGMRES",
                        [] { coupled_energy_and_charge("drift_kinetic", "ngmres", false); });
        check::scenario("Drift-kinetic coupled Newton-Krylov",
                        [] { coupled_energy_and_charge("drift_kinetic", "newton_krylov", false); });
        check::scenario("Drift-kinetic coupled substeps",
                        [] { coupled_energy_and_charge("drift_kinetic", "ngmres", true); });
        check::scenario("Hybrid kinetic ions and guiding-centre electrons",
                        [] { coupled_energy_and_charge("hybrid_kinetic", "anderson", false); });
    });
}
