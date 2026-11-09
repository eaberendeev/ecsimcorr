#include <cmath>

#include "simulation/drift_kinetic/drift_kinetic_particles.h"
#include "simulation/drift_kinetic/drift_kinetic_push.h"
#include "test_support.h"

namespace {
using driftKinetic::DriftKineticParticles;
using driftKinetic::Fields;
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
        check::scenario("Guiding-center particle container rebin", particle_container_rebin);
    });
}
