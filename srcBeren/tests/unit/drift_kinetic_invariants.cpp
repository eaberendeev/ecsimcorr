#include <cmath>

#include "simulation/drift_kinetic/drift_kinetic_adaptive.h"
#include "simulation/drift_kinetic/drift_kinetic_grid.h"
#include "simulation/drift_kinetic/drift_kinetic_particles.h"
#include "simulation/drift_kinetic/drift_kinetic_push.h"
#include "test_support.h"

namespace {
using driftKinetic::DriftKineticGrid;
using driftKinetic::DriftKineticParticles;
using driftKinetic::Fields;
using driftKinetic::VectorField;
using driftKinetic::GuidingCenter;
using driftKinetic::GuidingCenterPush;

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
    const auto mirror = [&](const Vector3R& start, const Vector3R& end) {
        const double z = 0.5 * (start.z() + end.z());
        const double length = 1.0 + 0.2 * z;
        return Fields{Vector3R(0, 0, 0), length, Vector3R(0, 0, 1), Vector3R(0, 0, 0.2), Vector3R(0, 0, 0)};
    };
    GuidingCenterPush push;
    push.set_tolerances(1e-14, 1e-13, 100);
    driftKinetic::AdaptiveSubstepping adaptive;
    adaptive.configure(true, 4);
    GuidingCenter state = particle(Vector3R(0, 0, 0), 1.0, 0.5);
    const auto energy = [&](const GuidingCenter& p) {
        const double length = 1.0 + 0.2 * p.position.z();
        return 0.5 * p.mass * p.v_parallel * p.v_parallel + p.mu * length;
    };
    const double initial_energy = energy(state);
    for (int step = 0; step < 20; ++step) {
        const auto result = adaptive.push(push, state, 0.01, mirror);
        check::require(result.converged, "adaptive guiding-centre push converged");
        check::require(result.leaves >= 1, "adaptive push has at least one leaf");
        state = result.state;
        check::near(energy(state), initial_energy, 1e-8, "adaptive mirror motion conserves energy");
    }
    bool rejected = false;
    try {
        adaptive.configure(true, 31);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    check::require(rejected, "adaptive substep depth is validated");
}

void uniform_grid_sampling() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    DriftKineticGrid grid(domain);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                grid.electric()(i, j, k) = Vector3R(1, 2, 3);
                grid.magnetic()(i, j, k) = Vector3R(0, 0, 2);
            }
    const auto fields = grid.sample(Vector3R(0.5, 0.5, 0.5));
    check::vector(fields.electric, Vector3R(1, 2, 3), 1e-14, "uniform electric interpolation");
    check::near(fields.lenB, 2.0, 1e-14, "uniform magnetic magnitude");
    check::vector(fields.b, Vector3R(0, 0, 1), 1e-14, "uniform magnetic direction");
    check::vector(fields.gradB, Vector3R(0, 0, 0), 1e-14, "uniform gradient");
    check::vector(fields.rotB, Vector3R(0, 0, 0), 1e-14, "uniform curl");
}

void linear_magnetic_gradient() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    DriftKineticGrid grid(domain);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k)
                grid.magnetic()(i, j, k) = Vector3R(0, 0, 1.0 + 0.2 * ((k - GHOST_CELLS) + 0.5));
    const auto fields = grid.sample(Vector3R(0.5, 0.5, 1.5));
    check::near(fields.lenB, 1.3, 1e-14, "linear magnetic magnitude");
    check::vector(fields.b, Vector3R(0, 0, 1), 1e-14, "linear magnetic direction");
    check::vector(fields.gradB, Vector3R(0, 0, 0.2), 1e-14, "linear magnetic gradient");
    check::vector(fields.rotB, Vector3R(0, 0, 0), 1e-14, "linear magnetic curl");
}

void grid_deposition() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    DriftKineticGrid grid(domain);
    VectorField current(domain.size()), magnetization(domain.size());
    const Vector3R zero(0, 0, 0);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                current(i, j, k) = zero;
                magnetization(i, j, k) = zero;
            }
    grid.deposit_current(current, Vector3R(0.5, 0.5, 0.5), Vector3R(1, 2, 3), 2.0);
    grid.deposit_magnetization(magnetization, Vector3R(0.5, 0.5, 0.5), Vector3R(0, 0, 1), 0.5, 2.0);
    Vector3R current_sum(0, 0, 0), magnetization_sum(0, 0, 0);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                current_sum += current(i, j, k);
                magnetization_sum += magnetization(i, j, k);
            }
    check::vector(current_sum, Vector3R(2, 4, 6), 1e-14, "current deposition preserves total weight");
    check::vector(magnetization_sum, Vector3R(0, 0, -1), 1e-14, "magnetization deposition preserves total moment");
}

void segment_deposition() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    DriftKineticGrid grid(domain);
    VectorField current(domain.size()), magnetization(domain.size());
    const Vector3R zero(0, 0, 0);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                current(i, j, k) = zero;
                magnetization(i, j, k) = zero;
            }
    grid.deposit_segment(current, magnetization, Vector3R(0.5, 0.5, 0.5), Vector3R(2.5, 0.5, 0.5),
                         Vector3R(1, 2, 3), Vector3R(0, 0, 1), 2.0, 0.5, 3.0);
    Vector3R current_sum(0, 0, 0), magnetization_sum(0, 0, 0);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                current_sum += current(i, j, k);
                magnetization_sum += magnetization(i, j, k);
            }
    check::vector(current_sum, Vector3R(6, 12, 18), 1e-13, "segment current preserves charge times velocity");
    check::vector(magnetization_sum, Vector3R(0, 0, -1.5), 1e-13, "segment magnetization preserves moment");
}

void particle_container_rebin() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    DriftKineticParticles particles(domain, "electrons", -1.0, 1.0, 2.0);
    GuidingCenter p;
    p.position = Vector3R(0.5, 0.5, 0.5);
    p.v_parallel = 2.0;
    particles.add(p);
    const auto cell = domain.storage_cell_index(p.position);
    check::require(particles.get_total_num_of_particles() == 1, "container counts one particle");
    check::require(particles(cell.x(), cell.y(), cell.z()).size() == 1, "particle lands in its storage cell");
    check::near(particles.parallel_kinetic_energy(), 2.0, 1e-14, "parallel kinetic energy");
    particles(cell.x(), cell.y(), cell.z())[0].position = Vector3R(2.5, 1.5, 1.5);
    particles.update_cells();
    check::require(particles(cell.x(), cell.y(), cell.z()).empty(), "old cell is emptied");
    const auto target = domain.storage_cell_index(Vector3R(2.5, 1.5, 1.5));
    check::require(particles(target.x(), target.y(), target.z()).size() == 1, "particle rebins to the new cell");
    check::require(particles.weight() == 2.0 && !particles.is_neutral(), "species parameters are preserved");
}
}   // namespace

int main() {
    return check::run([] {
        check::scenario("Guiding-center uniform-B parallel motion", uniform_parallel_motion);
        check::scenario("Guiding-center parallel electric acceleration", parallel_electric_acceleration);
        check::scenario("Guiding-center E-cross-B drift", electric_cross_magnetic_drift);
        check::scenario("Guiding-center mirror energy invariant", mirror_energy_invariant);
        check::scenario("Guiding-center adaptive mirror substeps", adaptive_substep_mirror);
        check::scenario("Guiding-center segment deposition", segment_deposition);
        check::scenario("Guiding-center particle container rebin", particle_container_rebin);
        check::scenario("Drift-kinetic uniform field interpolation", uniform_grid_sampling);
        check::scenario("Drift-kinetic linear magnetic gradient", linear_magnetic_gradient);
        check::scenario("Drift-kinetic current and magnetization deposition", grid_deposition);
    });
}
