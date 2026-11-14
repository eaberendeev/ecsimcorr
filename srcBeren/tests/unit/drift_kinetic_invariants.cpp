#include <cmath>
#include <cstdio>
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

// Guiding-centre trajectories in fields given on the Yee grid, with large
// steps, against analytic drifts. The fields are kept away from the periodic
// seams, where the analytic profiles are not periodic.
GuidingCenter advance_on_grid(YeeCoupling& coupling, const Field3d& electric, const MagneticLevels& levels,
                              const GuidingCenter& start, double dt) {
    GuidingCenterPush push;
    push.set_tolerances(1e-16, 1e-14, 200);
    const Vector3R start_magnetic = coupling.start_magnetic(levels, start.position, {0, 1});
    const auto fields = [&](const Vector3R& a, const Vector3R& b) {
        return coupling.interpolate(electric, levels, a, b, {0, 1}, &start_magnetic);
    };
    const auto result = push.push(start, dt, fields, true);
    check::require(result.converged, "guiding-centre push on the grid converged");
    return result.state;
}

// B = B0 (1 + (x - x0) / L) z: the discrete gradient of the linear |B| is
// exact, so the drift V = mu B0 / (q B L) y is reproduced exactly.
void grad_b_drift_on_grid() {
    const Domain domain = box({12, 40, 4}, {1, 1, 1});
    const double B0 = 1.0, L = 10.0, x0 = 6.0, dt = 10.0;
    Field3d B = zero_field(domain, 3), E = zero_field(domain, 3);
    fill_magnetic(B, domain, 2, [&](const Vector3R& r) { return B0 * (1 + (r.x() - x0) / L); });
    const MagneticLevels levels{&B, &B, &B};
    YeeCoupling coupling(domain, nullptr);
    coupling.prepare(levels);
    GuidingCenter p = particle(Vector3R(x0, 10.0, 2.0), 0.0, 0.1);
    const int steps = 20;
    for (int n = 0; n < steps; ++n) p = advance_on_grid(coupling, E, levels, p, dt);
    const double drift = 0.1 * B0 / (B0 * L);   // mu B0 / (q B L), q = 1, |B| = B0 at x0
    check::near(p.position.y(), 10.0 + drift * steps * dt, 1e-11, "grad-B drift along b x grad|B|");
    check::near(p.position.x(), x0, 1e-12, "no displacement along grad|B|");
    check::near(p.position.z(), 2.0, 1e-12, "no parallel displacement");
    check::near(p.v_parallel, 0.0, 1e-13, "no parallel acceleration across the gradient");
}

// On the axis of B_z = B0 (1 + z'^2 / L^2), B_x,y = -B0 (x', y') z' / L^2
// (z' = z - zc), the mirror force is linear: a harmonic bounce with
// omega^2 = 2 mu B0 / (m L^2). The implicit midpoint step rotates the phase by
// theta, tan(theta / 2) = omega dt / 2, so z'_{n+1} + z'_{n-1} = 2 cos(theta) z'_n.
void mirror_bounce_on_grid() {
    const Domain domain = box({8, 8, 48}, {1, 1, 1});
    const double B0 = 1.0, L = 8.0, xa = 4.0, ya = 4.0, zc = 24.0, mu = 0.5, dt = 12.0;
    Field3d B = zero_field(domain, 3), E = zero_field(domain, 3);
    fill_magnetic(B, domain, 0, [&](const Vector3R& r) { return -B0 * (r.x() - xa) * (r.z() - zc) / (L * L); });
    fill_magnetic(B, domain, 1, [&](const Vector3R& r) { return -B0 * (r.y() - ya) * (r.z() - zc) / (L * L); });
    fill_magnetic(B, domain, 2, [&](const Vector3R& r) { return B0 * (1 + std::pow((r.z() - zc) / L, 2)); });
    const MagneticLevels levels{&B, &B, &B};
    YeeCoupling coupling(domain, nullptr);
    coupling.prepare(levels);
    const double omega = std::sqrt(2 * mu * B0) / L;   // m = 1
    const double theta = 2 * std::atan(0.5 * omega * dt);
    std::vector<double> z{3.0};
    GuidingCenter p = particle(Vector3R(xa, ya, zc + 3.0), 0.0, mu);
    for (int n = 0; n < 24; ++n) {
        p = advance_on_grid(coupling, E, levels, p, dt);
        z.push_back(p.position.z() - zc);
        check::near(p.position.x(), xa, 1e-12, "bounce stays on the axis (x)");
        check::near(p.position.y(), ya, 1e-12, "bounce stays on the axis (y)");
    }
    double worst = 0;
    for (size_t n = 1; n + 1 < z.size(); ++n)
        if (std::abs(z[n]) > 0.5)
            worst = std::max(worst, std::abs(z[n + 1] + z[n - 1] - 2 * std::cos(theta) * z[n]));
    check::near(worst, 0, 1e-10, "discrete bounce phase advance tan(theta/2) = omega dt / 2");
    // The midpoint invariant v^2 + omega^2 z^2 of the harmonic bounce.
    check::near(p.v_parallel * p.v_parallel + omega * omega * z.back() * z.back(), omega * omega * 9.0, 1e-11,
                "bounce invariant");
}

// Vacuum toroidal field B = B0 R0 / R phi around an axis outside the box: the
// guiding centre circles the axis and drifts along it with
// V_z = (m v_par^2 / B + mu) / (q R) (curvature plus grad-B). The field is not
// polynomial, so the error must fall at second order in the cell size.
double toroidal_drift_error(double h) {
    const int n = static_cast<int>(std::lround(32 / h));
    const Domain domain = box({n, n, static_cast<int>(std::lround(8 / h))}, {h, h, h});
    const double B0 = 1.0, R0 = 20.0, xa = 16.0 - R0, ya = 16.0, v = 0.04, mu = 0.01, dt = 10.0;
    const auto radius2 = [&](const Vector3R& r) { return std::pow(r.x() - xa, 2) + std::pow(r.y() - ya, 2); };
    Field3d B = zero_field(domain, 3), E = zero_field(domain, 3);
    fill_magnetic(B, domain, 0, [&](const Vector3R& r) { return -B0 * R0 * (r.y() - ya) / radius2(r); });
    fill_magnetic(B, domain, 1, [&](const Vector3R& r) { return B0 * R0 * (r.x() - xa) / radius2(r); });
    const MagneticLevels levels{&B, &B, &B};
    YeeCoupling coupling(domain, nullptr);
    coupling.prepare(levels);
    GuidingCenter p = particle(Vector3R(16.0, 16.0, 4.0), v, mu);
    const int steps = 20;
    for (int step = 0; step < steps; ++step) p = advance_on_grid(coupling, E, levels, p, dt);
    const double expected = steps * dt * (v * v / B0 + mu) / R0;
    return std::abs((p.position.z() - 4.0) - expected) / expected;
}

void toroidal_drift_on_grid() {
    const double coarse = toroidal_drift_error(1.0), fine = toroidal_drift_error(0.5);
    std::printf("toroidal drift relative error: h=1 %.3e, h=0.5 %.3e\n", coarse, fine);
    check::require(fine < 1e-4, "toroidal drift within 1e-4 at h = 0.5");
    check::require(coarse / fine > 3.0, "toroidal drift converges at second order");
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

std::unique_ptr<SimulationDriftKinetic> simulation(const nlohmann::json& config,
                                                   const std::string& electron_model = "drift_kinetic") {
    auto electrons = check::species("Electrons", -1, 1, 1, 1);
    electrons["model"] = electron_model;
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
    std::filesystem::remove_all("Recovery");   // of other scenarios or of an earlier run
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

// Change of the marker model: the formulas alone. Energy m v^2 / 2 =
// m v_par^2 / 2 + mu |B|, the Larmor radius is the same in both descriptions,
// and the way back restores v_par and |v_perp| (the gyrophase is new).
void model_switch_formulas() {
    const Vector3R B(0.3, -0.4, 1.2);
    const double q = -1.0, m = 0.7, w = 2.5;
    const Particle marker(Vector3R(1.1, 2.2, 0.3), Vector3R(0.21, -0.37, 0.52));
    const GuidingCenter centre = driftKinetic::to_guiding_center(marker, B, q, m, w);
    check::vector(centre.position, marker.coord, 0, "the marker keeps its position");
    check::near(0.5 * m * centre.v_parallel * centre.v_parallel + centre.mu * B.length(),
                0.5 * m * marker.velocity.squaredNorm(), 1e-15, "energy of the guiding centre");
    check::near(centre.v_parallel, marker.velocity.dot(B.normalized()), 1e-15, "parallel velocity");
    const Vector3R no_field(0, 0, 0);
    check::near(driftKinetic::larmor_radius(centre, B), driftKinetic::larmor_radius(marker.velocity, no_field, B, q, m),
                1e-14, "Larmor radius in both descriptions");
    const double v_perp = (marker.velocity - centre.v_parallel * B.normalized()).length();
    check::near(driftKinetic::larmor_radius(marker.velocity, no_field, B, q, m), m * v_perp / B.length(), 1e-15,
                "Larmor radius m v_perp / |q B|");
    check::require(std::isinf(driftKinetic::larmor_radius(marker.velocity, no_field, no_field, q, m)),
                   "no field: unmagnetized");
    // In crossed fields the radius is that of the gyration around the E x B drift:
    // the same at every gyrophase, although the speed across B is not.
    const Vector3R E(0.05, 0.02, -0.01), b = B.normalized(), drift = E.cross(b) / B.length();
    const Vector3R e1 = b.cross(Vector3R(1, 0, 0)).normalized(), e2 = b.cross(e1);
    double slowest = 1e9, fastest = 0;
    for (int n = 0; n < 12; ++n) {
        const double phase = 2 * M_PI * n / 12;
        const Vector3R v = 0.3 * b + drift + 0.2 * (std::cos(phase) * e1 + std::sin(phase) * e2);
        check::near(driftKinetic::larmor_radius(v, E, B, q, m), m * 0.2 / B.length(), 1e-14,
                    "gyrophase-independent Larmor radius");
        const double speed = (v - v.dot(b) * b).length();
        slowest = std::min(slowest, speed);
        fastest = std::max(fastest, speed);
    }
    check::require(fastest - slowest > 0.05, "the speed across B does depend on the gyrophase");

    // Field directions along every axis exercise each choice of the perpendicular pair.
    for (const Vector3R& field : {B, Vector3R(2, 0, 0), Vector3R(0, -3, 0), Vector3R(0, 0, 0.5), Vector3R(1, 1, 1)})
        for (const double phase : {0.0, 1.0, 2.5, 4.0, 6.0}) {
            const GuidingCenter start = driftKinetic::to_guiding_center(marker, field, q, m, w);
            const Particle back = driftKinetic::to_kinetic(start, field, phase);
            check::near(back.velocity.squaredNorm(), marker.velocity.squaredNorm(), 1e-14, "round-trip energy");
            check::near(back.velocity.dot(field.normalized()), start.v_parallel, 1e-14, "round-trip v_par");
            const GuidingCenter again = driftKinetic::to_guiding_center(back, field, q, m, w);
            check::near(again.mu, start.mu, 1e-14, "round-trip magnetic moment");
        }
    // The gyrophase depends on the phase argument: two phases give different v_perp directions.
    const Vector3R v0 = driftKinetic::to_kinetic(centre, B, 0.0).velocity;
    const Vector3R v1 = driftKinetic::to_kinetic(centre, B, 0.5 * M_PI).velocity;
    check::near((v0 - centre.v_parallel * b).dot(v1 - centre.v_parallel * b), 0, 1e-14, "quarter turn of v_perp");

    // Hashed gyrophase: in range, repeatable, and uniform enough over many markers.
    double mean_cos = 0, mean_sin = 0;
    const int samples = 4000;
    for (int n = 0; n < samples; ++n) {
        GuidingCenter g = centre;
        g.position.x() += 1e-3 * n;
        const double phase = driftKinetic::hashed_gyrophase(g, 7);
        check::require(phase >= 0.0 && phase < 2 * M_PI, "gyrophase range");
        check::near(phase, driftKinetic::hashed_gyrophase(g, 7), 0, "gyrophase is a function of the state");
        mean_cos += std::cos(phase) / samples;
        mean_sin += std::sin(phase) / samples;
    }
    // Uniform phases: |<exp(i phase)>| ~ 1 / sqrt(samples) = 0.016.
    check::require(std::hypot(mean_cos, mean_sin) < 0.06, "gyrophases are spread over the circle");
    check::require(driftKinetic::hashed_gyrophase(centre, 7) != driftKinetic::hashed_gyrophase(centre, 8),
                   "gyrophase changes with the step");
}

// Charge of both parts of an adaptive species on the grid, with the common shape.
void adaptive_charge(SimulationDriftKinetic& sim, Field3d& density) {
    density.resize(sim.domain.size(), 1);
    density.setZero();
    const auto& kinetic = *sim.species.at("Electrons");
    for (int cell = 0; cell < kinetic.size(); ++cell)
        for (const auto& p : kinetic(cell))
            algorithmsECCAPFIM::deposit_charge(density, p.coord, kinetic.charge() * kinetic.mpw(), sim.domain,
                                               &sim.bc_handler);
    const auto& guiding = sim.guiding_centers("Electrons_gc");
    for (int cell = 0; cell < guiding.size(); ++cell)
        for (const auto& p : guiding(cell))
            algorithmsECCAPFIM::deposit_charge(density, p.position, p.charge * p.weight, sim.domain, &sim.bc_handler);
}

// Adaptive species in the coupled step, sorted by the Larmor radius alone
// (the omega_c dt threshold is zero). h = 0.7 and |B| = 1.005: with the
// thresholds 0.1 h and 0.2 h a marker of the Larmor radius below 0.07 becomes
// a guiding centre, above 0.14 it is kinetic, and in between it keeps its model.
void model_switch_in_coupled_step() {
    auto config = settings("hybrid_kinetic");
    config["ModelSwitch"] = {{"omega_dt_to_guiding_center", 0.0},
                             {"larmor_radius_to_guiding_center", 0.1},
                             {"larmor_radius_to_kinetic", 0.2},
                             {"interval", 1}};
    auto sim = simulation(config, "adaptive");
    perturb_fields(*sim);
    auto& kinetic = *sim->species.at("Electrons");
    auto& guiding = sim->guiding_centers("Electrons_gc");
    check::require(guiding.get_total_num_of_particles() == 0, "adaptive species starts kinetic");
    // b is close to z: v_x, v_y set the Larmor radius.
    const std::vector<std::pair<Vector3R, Vector3R>> markers = {
        {{2.79, 5.49, 2.69}, {0.01, 0.02, 0.6}},    // rho ~ 0.06 with the tilt of b: guiding centre
        {{0.35, 1.1, 0.45}, {0.02, -0.01, -0.4}},   // guiding centre
        {{1.5, 2.3, 1.4}, {0.3, 0.2, 0.2}},         // rho ~ 0.36: kinetic
        {{2.1, 4.0, 0.1}, {-0.25, 0.3, -0.7}},      // kinetic
        {{1.0, 1.0, 1.0}, {0.0, 0.1, 0.0}}};        // rho ~ 0.1: between the thresholds, stays kinetic
    for (const auto& [r, v] : markers) kinetic.add_particle(Particle(r, v));
    sim->species.at("Ions")->add_particle(Particle({2.79, 5.49, 2.69}, {-0.1, 0.2, 0.05}));
    const int total = static_cast<int>(markers.size());

    Field3d before, after;
    adaptive_charge(*sim, before);
    const double energy_before = total_energy(*sim);
    sim->switch_models(0);
    adaptive_charge(*sim, after);
    check::near((after - before).norm(), 0, 1e-15, "the model change keeps the charge density");
    check::near(total_energy(*sim), energy_before, 1e-14, "the model change keeps the energy");
    const auto& statistics = sim->switch_statistics("Electrons_gc");
    check::require(statistics.to_guiding_center == 2 && statistics.to_kinetic == 0, "two magnetized markers moved");
    check::near(statistics.energy_defect, 0, 1e-15, "reported energy defect");
    check::require(guiding.get_total_num_of_particles() == 2 && kinetic.get_total_num_of_particles() == 3,
                   "marker counts after the first sorting");
    sim->switch_models(0);
    check::require(sim->switch_statistics("Electrons_gc").to_guiding_center == 0 &&
                       sim->switch_statistics("Electrons_gc").to_kinetic == 0,
                   "a second call moves nothing");

    // Coupled steps with both parts present: the energy and the marker count are kept.
    const double initial = total_energy(*sim);
    for (int step = 1; step <= 4; ++step) {
        sim->prepare_step(step);
        sim->make_step(step);
        check::near(total_energy(*sim), initial, 1e-10, "energy with markers of both models");
        check::require(sim->nonlinear_residual < 1e-11, "accepted coupled field residual");
        check::require(sim->guiding_centers("Electrons_gc").get_total_num_of_particles() +
                               sim->species.at("Electrons")->get_total_num_of_particles() ==
                           total,
                       "markers are neither lost nor doubled");
        sim->make_diagnostic(step);
    }

    // Heating a guiding centre across the field takes it over the upper
    // threshold: at the next step it is kinetic again, with the same energy.
    auto& centres = sim->guiding_centers("Electrons_gc");
    GuidingCenter* heated = nullptr;
    for (int cell = 0; cell < centres.size() && !heated; ++cell)
        if (!centres(cell).empty())
            heated = &centres(cell).front();
    check::require(heated != nullptr, "a guiding centre to heat");
    const double lenB = 1.0;   // order of |B|; rho = sqrt(2 mu m / |B|) / |q| = 0.4
    heated->mu = 0.5 * 0.4 * 0.4 * lenB;
    const double v_parallel = heated->v_parallel;
    const int centres_before = centres.get_total_num_of_particles();
    Field3d charge_before, charge_after;
    adaptive_charge(*sim, charge_before);
    const double heated_energy = total_energy(*sim);
    sim->prepare_step(5);
    adaptive_charge(*sim, charge_after);
    check::require(sim->switch_statistics("Electrons_gc").to_kinetic == 1, "the heated centre became kinetic");
    // Other markers may cross a threshold at the same step: the count follows the statistics.
    check::require(sim->guiding_centers("Electrons_gc").get_total_num_of_particles() ==
                       centres_before - 1 + sim->switch_statistics("Electrons_gc").to_guiding_center,
                   "guiding-centre count follows the statistics");
    check::near((charge_after - charge_before).norm(), 0, 1e-15, "charge density on the way back");
    check::near(total_energy(*sim), heated_energy, 1e-14, "energy on the way back");
    // The new kinetic marker carries the parallel velocity of the centre.
    const Field3d magnetic = sim->fieldBn + sim->fieldBInit;
    driftKinetic::YeeCoupling coupling(sim->domain, &sim->bc_handler);
    double closest = 1e9;
    const auto& electrons = *sim->species.at("Electrons");
    for (int cell = 0; cell < electrons.size(); ++cell)
        for (const auto& p : electrons(cell))
            closest = std::min(
                closest, std::abs(p.velocity.dot(coupling.magnetic_at(magnetic, p.coord).normalized()) - v_parallel));
    check::near(closest, 0, 1e-14, "parallel velocity of the returned marker");
    sim->make_step(5);
    check::near(total_energy(*sim), heated_energy, 1e-10, "step after the way back");
    sim->finalize();
}

// Sorting by omega_c dt = |q| |B| dt / m: a rule of the position. B_z =
// 1 + 0.08 cos(2 pi x / Lx) with dt = 0.35 gives omega_c dt from 0.32 to 0.38.
// Markers where it exceeds 0.36 become guiding centres whatever their
// velocity, where it is below 0.34 they are kinetic.
void model_switch_by_gyration_step() {
    auto config = settings("hybrid_kinetic");
    config["DriftKinetic"]["magnetic"] = {0, 0, 1};
    config["ModelSwitch"] = {{"omega_dt_to_guiding_center", 0.36},
                             {"omega_dt_to_kinetic", 0.34},
                             {"larmor_radius_to_guiding_center", 100.0}};
    auto sim = simulation(config, "adaptive");
    const auto size = sim->domain.size();
    for (int i = 0; i < size.x(); ++i)
        for (int j = 0; j < size.y(); ++j)
            for (int k = 0; k < size.z(); ++k)
                sim->fieldBn(i, j, k, 2) = 0.08 * std::cos(2 * M_PI * (i - GHOST_CELLS) / 4);
    algorithmsECCAPFIM::copy_periodic_nodes(sim->fieldBn, sim->domain);
    auto& kinetic = *sim->species.at("Electrons");
    const double lx = 4 * 0.7;
    const int count = 40;
    for (int n = 0; n < count; ++n)   // slow and fast markers alternate along x
        kinetic.add_particle(Particle({lx * (n + 0.5) / count, 2.0, 1.0}, {n % 2 ? 0.5 : 0.01, 0.02, 0.1}));
    sim->switch_models(0);
    const Field3d magnetic = sim->fieldBn + sim->fieldBInit;
    driftKinetic::YeeCoupling coupling(sim->domain, &sim->bc_handler);
    int centres = 0, fast_centres = 0;
    const auto& guiding = sim->guiding_centers("Electrons_gc");
    for (int cell = 0; cell < guiding.size(); ++cell)
        for (const auto& p : guiding(cell)) {
            check::require(sim->dt_ * coupling.magnetic_at(magnetic, p.position).length() > 0.36,
                           "guiding centres only where omega_c dt is above the threshold");
            ++centres;
            fast_centres += p.mu > 0.05;
        }
    for (int cell = 0; cell < kinetic.size(); ++cell)
        for (const auto& p : kinetic(cell))
            check::require(sim->dt_ * coupling.magnetic_at(magnetic, p.coord).length() <= 0.36,
                           "kinetic markers only where omega_c dt is below the threshold");
    check::require(centres > 5 && centres < count - 5, "both models are present");
    check::require(fast_centres > 2, "the rule does not look at the velocity");
    // Lowering the field moves the centres back when omega_c dt falls below 0.34.
    sim->fieldBn *= -1.0;   // B_z = 1 - 0.08 cos: strong and weak regions swap
    sim->switch_models(1);
    check::require(sim->switch_statistics("Electrons_gc").to_kinetic > 0 &&
                       sim->switch_statistics("Electrons_gc").to_guiding_center > 0,
                   "markers follow the field in both directions");
    check::near(sim->switch_statistics("Electrons_gc").energy_defect, 0, 1e-14, "energy of the model changes");
    check::near(driftKinetic::numerical_moment_factor(0.5), 1 / 1.0625, 1e-15, "moment factor at omega_c dt = 0.5");
    sim->finalize();
}

// One electron bouncing in a magnetic mirror B_z = B0 (1 + z^2 / L^2) with
// prescribed fields. The marker is kinetic near the centre (omega_c dt below
// 1.5 B0 dt, |z| < 0.71 L) and a guiding centre near the turning points
// (omega_c dt above 2 B0 dt, |z| > L), so it changes its model four times per
// bounce period. Returns the largest |mu / mu0 - 1| among the guiding-centre
// stages; checks the energy, the turning point and that both directions of
// the change occurred.
double mirror_bounce_with_model_change(double v_perp0, double dt, double periods) {
    auto config = settings("hybrid_kinetic");
    config["NumCellsX"] = 8;
    config["NumCellsY"] = 8;
    config["NumCellsZ"] = 48;
    config["Dx"] = config["Dy"] = config["Dz"] = 1.0;
    config["Dt"] = dt;
    config["DriftKinetic"] = {{"magnetic", {0, 0, 1}}, {"evolve_fields", false}};
    config["ModelSwitch"] = {{"omega_dt_to_guiding_center", 2.0 * dt}, {"omega_dt_to_kinetic", 1.5 * dt}};
    auto sim = simulation(config, "adaptive");
    const double B0 = 1.0, L = 8.0, xa = 4.0, ya = 4.0, zc = 24.0;
    fill_magnetic(sim->fieldBInit, sim->domain, 0,
                  [&](const Vector3R& r) { return -B0 * (r.x() - xa) * (r.z() - zc) / (L * L); });
    fill_magnetic(sim->fieldBInit, sim->domain, 1,
                  [&](const Vector3R& r) { return -B0 * (r.y() - ya) * (r.z() - zc) / (L * L); });
    fill_magnetic(sim->fieldBInit, sim->domain, 2,
                  [&](const Vector3R& r) { return B0 * (1 + std::pow((r.z() - zc) / L, 2)); });
    // The electron (q = -1, m = 1) starts at the centre on a Larmor circle around the axis.
    const double v_par0 = std::sqrt(2.0) * v_perp0, mu0 = 0.5 * v_perp0 * v_perp0 / B0;
    sim->species.at("Electrons")->add_particle(Particle({xa + v_perp0 / B0, ya, zc}, {0, v_perp0, v_par0}));
    const double energy = total_energy(*sim);
    const double period = 2 * M_PI * L / v_perp0;   // harmonic bounce, omega_b = sqrt(2 mu B0 / m) / L
    const int steps = static_cast<int>(periods * period / sim->dt_);
    long to_guiding = 0, to_kinetic = 0;
    double z_max = 0, mu_error = 0, energy_error = 0, radius = 0;
    for (int step = 1; step <= steps; ++step) {
        sim->prepare_step(step);
        const auto& statistics = sim->switch_statistics("Electrons_gc");
        to_guiding += statistics.to_guiding_center;
        to_kinetic += statistics.to_kinetic;
        sim->make_step(step);
        energy_error = std::max(energy_error, std::abs(total_energy(*sim) - energy));
        const auto& centres = sim->guiding_centers("Electrons_gc");
        for (int cell = 0; cell < centres.size(); ++cell)
            for (const auto& p : centres(cell)) {
                z_max = std::max(z_max, std::abs(p.position.z() - zc));
                mu_error = std::max(mu_error, std::abs(p.mu / mu0 - 1));
                radius = std::max(radius, std::hypot(p.position.x() - xa, p.position.y() - ya));
            }
    }
    check::require(sim->guiding_centers("Electrons_gc").get_total_num_of_particles() +
                           sim->species.at("Electrons")->get_total_num_of_particles() ==
                       1,
                   "the electron is kept through the model changes");
    // Each passage of a turning point is one change to a guiding centre and one back.
    check::require(to_guiding >= 2 && to_kinetic >= 2, "the model changed in both directions");
    check::near(energy_error, 0, 1e-12 * energy, "energy through the bounces and the model changes");
    // Turning point B_t = energy / mu: z_t = L v_par0 / v_perp0 for mu0, and
    // d z_t / z_t = (3 / 4) |d mu / mu| at B_t = 3 B0.
    check::near(z_max / (L * std::sqrt(2.0)), 1.0, 0.75 * mu_error + 0.01, "mirror turning point");
    check::require(radius < 4 * v_perp0 / B0, "the centre stays within a few Larmor radii of the axis");
    sim->finalize();
    return mu_error;
}

// The moment mu of the guiding-centre stages is not that of the start: the
// kinetic stages do not keep m v_perp^2 / (2 |B|). The orbit of the midpoint
// rule feels the mirror force of the moment mu / (1 + (omega_c dt / 2)^2)
// (numerical_moment_factor), so the deviation is of the order (omega_c dt / 2)^2
// at the upper threshold and falls as dt^2. The displacement R = x of the
// change itself adds a term of the order rho / L = 1 / 40, which is not seen.
void model_switch_in_mirror() {
    const double coarse = mirror_bounce_with_model_change(0.2, 0.3, 1.5);   // omega_c dt = 0.6 at the change
    const double fine = mirror_bounce_with_model_change(0.2, 0.1, 1.5);     // omega_c dt = 0.2
    std::printf("mu deviation after kinetic stages: dt=0.3 %.3e, dt=0.1 %.3e\n", coarse, fine);
    check::require(coarse < 1.5 * 0.3 * 0.3, "mu deviation bounded by (omega_c dt / 2)^2 at dt = 0.3");
    check::require(fine < 1.5 * 0.1 * 0.1, "mu deviation bounded by (omega_c dt / 2)^2 at dt = 0.1");
    check::require(fine < coarse / 4, "mu deviation falls as dt^2");
}

// Restart of an adaptive species: both parts are stored, and the model
// changes of the continued run repeat those of the uninterrupted one.
void model_switch_restart() {
    std::filesystem::remove_all("Recovery");   // of other scenarios or of an earlier run
    auto config = settings("hybrid_kinetic");
    config["ModelSwitch"] = {{"omega_dt_to_guiding_center", 0.0},
                             {"larmor_radius_to_guiding_center", 0.1},
                             {"larmor_radius_to_kinetic", 0.2}};
    config["diagnostics"]["outputs"] = nlohmann::json::array({{{"type", "recovery"}, {"interval", 1}}});
    auto original = simulation(config, "adaptive");
    perturb_fields(*original);
    for (const auto& [r, v] : std::vector<std::pair<Vector3R, Vector3R>>{{{2.79, 5.49, 2.69}, {0.01, 0.02, 0.6}},
                                                                         {{0.35, 1.1, 0.45}, {0.02, -0.01, -0.4}},
                                                                         {{1.5, 2.3, 1.4}, {0.3, 0.2, 0.2}},
                                                                         {{1.0, 1.0, 1.0}, {0.0, 0.065, 0.3}}})
        original->species.at("Electrons")->add_particle(Particle(r, v));
    original->species.at("Ions")->add_particle(Particle({2.79, 5.49, 2.69}, {-0.1, 0.2, 0.05}));
    original->prepare_step(1);
    original->make_step(1);
    original->make_diagnostic(1);
    check::require(std::filesystem::exists("Recovery/Particles/Electrons_gc/Electrons_gc.gc.backup"),
                   "recovery file of the guiding-centre part");

    config["StartFromTime"] = 1;
    config["diagnostics"]["outputs"] = nlohmann::json::array();
    const auto directory = std::filesystem::current_path();
    std::filesystem::create_directories("restart");
    std::filesystem::current_path("restart");
    auto resumed = simulation(config, "adaptive");
    std::filesystem::current_path(directory);
    check::require(resumed->guiding_centers("Electrons_gc").get_total_num_of_particles() ==
                       original->guiding_centers("Electrons_gc").get_total_num_of_particles(),
                   "restart: guiding-centre part");
    check::require(resumed->species.at("Electrons")->get_total_num_of_particles() ==
                       original->species.at("Electrons")->get_total_num_of_particles(),
                   "restart: kinetic part");
    for (int step = 2; step <= 3; ++step)
        for (auto* sim : {original.get(), resumed.get()}) {
            sim->prepare_step(step);
            sim->make_step(step);
        }
    check::near((resumed->fieldEn - original->fieldEn).norm(), 0, 1e-12, "restart electric field");
    check::near((resumed->fieldBn - original->fieldBn).norm(), 0, 1e-12, "restart magnetic field");
    check::near(total_energy(*resumed), total_energy(*original), 1e-12, "restart total energy");
    check::require(resumed->guiding_centers("Electrons_gc").get_total_num_of_particles() ==
                       original->guiding_centers("Electrons_gc").get_total_num_of_particles(),
                   "restart: same model of every marker");
    original->finalize();
    resumed->finalize();
}

// Velocity distributions: weight sums, bin placement and the spatial box.
void distribution_functions() {
    auto sim = simulation(settings("hybrid_kinetic"));
    add_guiding_centers(*sim);
    sim->species.at("Ions")->add_particle(Particle({2.79, 5.49, 2.69}, {-0.1, 0.2, 0.05}));
    sim->species.at("Ions")->add_particle(Particle({0.35, 1.1, 0.45}, {0.3, -0.25, 0.0}));
    const MomentSources* sources = sim->extra_moment_sources();
    const auto& box = sim->domain.geom;
    DistributionRequest gc{"vpar_mu", {{{-1.0, 1.0, 8}, {0.0, 0.2, 4}}}, box.box_min, box.box_max};
    DistributionRequest kinetic{"vx_vy", {{{-0.5, 0.5, 10}, {-0.5, 0.5, 10}}}, box.box_min, box.box_max};
    DistributionRequest corner = gc;
    corner.region_max = Vector3R(1.0, 2.0, 1.0);   // holds only the marker at (0.35, 1.1, 0.45)
    DiagnosticMeasurements measured({sim->domain, sim->species, sim->fieldEn, sim->fieldBn, sim->fieldE_external,
                                     sim->fieldBInit, sim->bc_handler, ShapeType::Linear, sources});
    measured.request_distribution("Electrons", "all", gc);
    measured.request_distribution("Electrons", "corner", corner);
    measured.request_distribution("Ions", "all", kinetic);
    measured.calculate();
    const auto integral = [](const Distribution& f) {
        const auto& [x, y] = f.request.axes;
        double sum = 0;
        for (double value : f.values) sum += value;
        return sum * (x.max - x.min) / x.bins * (y.max - y.min) / y.bins;
    };
    const double w = sim->guiding_centers("Electrons").weight();
    check::near(integral(measured.distribution("Electrons", "all")), 4 * w, 1e-14,
                "guiding-centre f integrates to the weight");
    check::near(integral(measured.distribution("Electrons", "corner")), w, 1e-14, "spatial box selects markers");
    check::near(integral(measured.distribution("Ions", "all")), 2 * sim->species.at("Ions")->mpw(), 1e-14,
                "kinetic f integrates to the weight");
    // The marker (v_par, mu) = (-0.4, 0.09) of the corner box: bin (2, 1) of the 8 x 4 grid.
    const auto& f = measured.distribution("Electrons", "corner");
    check::near(f.values[2 * 4 + 1] * 0.25 * 0.05, w, 1e-14, "guiding centre in its (v_par, mu) bin");

    DiagnosticMeasurements wrong({sim->domain, sim->species, sim->fieldEn, sim->fieldBn, sim->fieldE_external,
                                  sim->fieldBInit, sim->bc_handler, ShapeType::Linear, sources});
    wrong.request_distribution("Electrons", "wrong", kinetic);
    bool rejected = false;
    try {
        wrong.calculate();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "kinetic projection is rejected for guiding centres");
    sim->finalize();
}

// Momentum columns: E x B of uniform fields on the Yee grid and the
// parallel momentum of the guiding centres.
void momentum_diagnostics() {
    auto sim = simulation(settings("hybrid_kinetic"));
    add_guiding_centers(*sim);
    Field3d electric(sim->domain.size(), 3), magnetic(sim->domain.size(), 3);
    electric.setZero();
    magnetic.setZero();
    for (int p = 0; p < electric.capacity(); p += 3) {
        electric[p + 1] = 0.3;    // E_y
        magnetic[p + 2] = -0.7;   // B_z
        magnetic[p] = 0.2;        // B_x
    }
    const Vector3R expected = 60.0 * Vector3R(0, 0.3, 0).cross(Vector3R(0.2, 0, -0.7));   // 4 x 5 x 3 cells
    check::vector(sim->field_momentum(electric, magnetic), expected, 1e-13, "uniform field momentum");
    Vector3R parallel(0, 0, 0);
    const Vector3R b = Vector3R(0.1, 0, 1).normalized();   // the prescribed B of the settings, B_self = 0
    const auto& markers = sim->guiding_centers("Electrons");
    for (int cell = 0; cell < markers.size(); ++cell)
        for (const auto& p : markers(cell)) parallel += (p.mass * p.weight * p.v_parallel) * b;
    check::vector(sim->extra_momentum(), parallel, 1e-14, "guiding-centre parallel momentum");
    sim->finalize();
}

// Open z walls: a guiding centre streaming into the upper wall is cut at the
// wall within the step and absorbed with its energy at the hit; the others stay.
// Particle + field energy plus the absorbed energy is conserved.
void guiding_center_walls() {
    auto config = settings("drift_kinetic");
    config["Boundary_conditions"] = nlohmann::json::parse(
        R"([{"periodic":{"face":"XMIN"}},{"periodic":{"face":"YMIN"}},{"open":[{"face":"ZMIN"},{"face":"ZMAX"}]}])");
    if constexpr (GHOST_CELLS < 2) {
        bool rejected = false;
        try {
            simulation(config);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check::require(rejected, "walls need two ghost layers");
        return;
    } else {
        auto sim = simulation(config);
        // The container is replaced at every accepted step: take it anew after a step.
        auto& electrons = sim->guiding_centers("Electrons");
        const double weight = electrons.weight();
        const double top = sim->domain.geom.box_max.z();
        const std::vector<std::pair<Vector3R, double>> markers = {
            {{1.5, 2.3, top - 0.3}, 3.0}, {{0.9, 1.1, 1.3}, 0.4}, {{2.1, 4.0, 0.9}, -0.3}};
        for (const auto& [r, v] : markers) {
            GuidingCenter p;
            p.position = r;
            p.v_parallel = v;
            p.mu = 0.05;
            p.charge = electrons.charge();
            p.mass = electrons.mass();
            p.weight = electrons.weight();
            electrons.add(p);
        }
        const double initial = total_energy(*sim);
        double lost = 0;
        const auto range = sim->bc_handler.active_range(sim->domain.grid);
        for (int step = 1; step <= 3; ++step) {
            const double energy_before = sim->guiding_center_energy(sim->fieldBn + sim->fieldBInit);
            sim->prepare_step(step);
            sim->make_step(step);
            lost += sim->extra_lost_energy();
            const double work =
                sim->dt_ * dot_product_sum(sim->fieldEp + sim->fieldE_external, sim->fieldJ_guiding, range);
            const double moment_work = dot_product_sum(sim->fieldM, sim->fieldBn - sim->fieldB, range);
            const double energy_after = sim->guiding_center_energy(sim->fieldBn + sim->fieldBInit);
            check::near(energy_after + sim->extra_lost_energy() - energy_before - work + moment_work, 0, 1e-11,
                        "guiding-centre energy theorem with absorption");
            check::near(total_energy(*sim) + lost, initial, 1e-10, "energy conservation with the absorbed energy");
            for (int cell = 0; cell < sim->guiding_centers("Electrons").size(); ++cell)
                for (const auto& p : sim->guiding_centers("Electrons")(cell)) {
                    check::require(p.position.z() >= 0 && p.position.z() <= top, "guiding centres stay in the box");
                }
        }
        check::require(sim->guiding_centers("Electrons").get_total_num_of_particles() == 2, "one centre absorbed");
        // Energy at the wall: m v^2 / 2 (v_par changes little in the weak fields) plus mu |B|.
        check::require(lost > 0.5 * 9.0 * weight * 0.9, "absorbed energy is the streaming energy");
        sim->finalize();
    }
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
        check::scenario("Grid guiding centre: grad-B drift", grad_b_drift_on_grid);
        check::scenario("Grid guiding centre: mirror bounce phase", mirror_bounce_on_grid);
        check::scenario("Grid guiding centre: toroidal curvature drift", toroidal_drift_on_grid);
        check::scenario("Yee coupling magnetization identity", magnetization_identity);
        check::scenario("Drift-kinetic prescribed fields", prescribed_drift);
        check::scenario("Guiding-centre grid moments for diagnostics", guiding_center_moments);
        check::scenario("Guiding-centre restart from recovery files", guiding_center_restart);
        check::scenario("Velocity distributions of kinetic and guiding-centre species", distribution_functions);
        check::scenario("Momentum of fields and guiding centres", momentum_diagnostics);
        check::scenario("Guiding centres absorbed by open walls", guiding_center_walls);
        check::scenario("Model change of a marker: formulas", model_switch_formulas);
        check::scenario("Model change in the coupled hybrid step", model_switch_in_coupled_step);
        check::scenario("Model change by omega_c dt", model_switch_by_gyration_step);
        check::scenario("Model change in a magnetic mirror", model_switch_in_mirror);
        check::scenario("Model change: restart", model_switch_restart);
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
