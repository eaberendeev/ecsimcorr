#include <cmath>
#include <iostream>
#include <stdexcept>

#include "ParticlesDiagnostic.h"
#include "algorithms_ecsim_corr.h"
#include "boundary_conditions.h"
#include "test_support.h"

namespace {
int checks = 0;
void require(bool condition, const char* message) {
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}

void check_captured_fate_swaps() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    auto cfg = nlohmann::json::parse(R"({"Name":"Electrons","Mass":1,"NumPartPerCell":1,"distribution":[]})");
    ParticlesArray sp(cfg, domain);
    for (double vz : {-0.5, -5.0, -0.4}) sp.add_particle(Particle(Vector3R(2, 2, -0.25), Vector3R(0, 0, vz)));
    BoundaryConditionHandler handler;
    handler.load_from_json(
        nlohmann::json::parse(
            R"({"Boundary_conditions":[{"electron_reflection":{"face":"ZMIN","radius":1,"energy_threshold":100}}]})"),
        domain);
    ParticleBoundaryEvents captured;
    for (int cell = 0; cell < sp.size(); ++cell)
        for (const auto& p : sp.particlesData(cell)) captured.record(p, handler.classify_particle(p, sp, domain));
    require(captured.size() == 3, "snapshot stores only three outgoing particles");
    for (int cell = 0; cell < sp.size(); ++cell)
        for (auto& p : sp.particlesData(cell)) p.velocity *= 100;
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &captured);
    require(sp.get_total_num_of_particles() == 2, "captured fates lost alignment during particle swaps");
    require(sp.diag.boundary.at(Face::ZMIN).reflected_count == 2, "wrong frozen reflection statistics");
    require(sp.diag.boundary.at(Face::ZMIN).lost_count == 1, "wrong frozen absorption statistics");
    cfg["Name"] = "Ions";
    ParticlesArray ions(cfg, domain);
    const Particle trial(Vector3R(2, 2, -0.25), Vector3R(0, 0, -0.01));
    require(handler.classify_particle(trial, ions, domain).fate == ParticleFate::Removed,
            "electron reflection must not reflect ions");
}

template <ShapeType ST>
void exercise(const std::string& condition, Face face, bool reference, bool fused) {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    auto config = nlohmann::json::parse(
        R"({"Name":"Electrons","Charge":-1,"Mass":1,"Density":1,"NumPartPerCell":1,"distribution":[]})");
    ParticlesArray sp(config, domain);
    const bool lower = face == Face::ZMIN;
    const double z = lower ? 0.25 : 3.75;
    const double vz = lower ? -0.5 : 0.5;
    // Hit inside the reflecting aperture, then end outside its radius.
    // The physical wall interaction must use the hit position, not the endpoint.
    sp.add_particle(Particle(Vector3R(2.75, 2.0, z), Vector3R(0.4, 0, vz)));
    BoundaryConditionHandler handler;
    nlohmann::json params = {{"face", lower ? "ZMIN" : "ZMAX"},
                             {"radius", 1.0},
                             {"energy_threshold", 100.0},
                             {"electron_threshold_energy_", 0.2},
                             {"gap", 0.0}};
    handler.load_from_json({{"Boundary_conditions", nlohmann::json::array({{{condition, params}}})}}, domain);
    ParticleBoundaryEvents fates;
    const auto policy = handler.path_policy(sp, domain, &fates);
    sp.density_on_grid_update<ST>();
    const Field3d before = sp.densityOnGrid;
    sp.currentOnGrid.setZero();
    if (fused) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2.0, 1.0, sp.currentOnGrid, work, policy);
    } else if (reference) {
        move_and_calc_current_reference<ST>(sp, 1.0, sp.currentOnGrid, policy);
    } else {
        move_and_calc_current<ST>(sp, 1.0, sp.currentOnGrid, policy);
    }
    require(fates.size() == 1, "deposit must save just the outgoing particle event");
    // Cross the reflection threshold after depositing: the stored fate must win.
    for (int cell = 0; cell < sp.size(); ++cell)
        for (auto& p : sp.particlesData(cell)) p.velocity *= 100;
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &fates);
    sp.update_cells(domain);
    require(sp.get_total_num_of_particles() == 1, "reflected electron disappeared");
    for (int cell = 0; cell < sp.size(); ++cell) {
        for (const auto& p : sp.particlesData(cell)) {
            require(std::abs(p.coord.z() - z) < 1e-14, "wrong reflected endpoint");
            require(p.velocity.z() * vz < 0, "velocity did not reflect");
        }
    }
    sp.density_on_grid_update<ST>();
    double max_residual = 0;
    for (int i = 1; i < domain.size().x() - 1; ++i)
        for (int j = 1; j < domain.size().y() - 1; ++j)
            for (int k = 1; k < domain.size().z() - 1; ++k) {
                const auto& current = sp.currentOnGrid;
                const double div = current(i, j, k, 0) - current(i - 1, j, k, 0) + current(i, j, k, 1) -
                                   current(i, j - 1, k, 1) + current(i, j, k, 2) - current(i, j, k - 1, 2);
                const double residual = sp.densityOnGrid(i, j, k, 0) - before(i, j, k, 0) + div;
                require(std::isfinite(residual), "nonfinite reflected-current continuity residual");
                max_residual = std::max(max_residual, std::abs(residual));
            }
    if (max_residual > 1e-12)
        std::cerr << "residual=" << max_residual << '\n';
    require(max_residual < 1e-12, "reflected current violates charge continuity");

    // Same face and geometry, but energy above threshold: absorption, no return leg.
    Particle trial(Vector3R(2.0, 2.0, lower ? -0.25 : 4.25), Vector3R(0, 0, 10 * vz));
    require(handler.classify_particle(trial, sp, domain).fate == ParticleFate::Removed,
            "energetic electron should be absorbed");
    const auto path = current_particle_path(Vector3R(2, 2, z), trial, domain, policy);
    require(!path.reflected && std::abs(path.end.z() - (lower ? 0.0 : 4.0)) < 1e-14,
            "absorbed current must end at the wall");
    trial.coord.x() = 3.5;
    trial.velocity.z() = vz;
    require(handler.classify_particle(trial, sp, domain).fate == ParticleFate::Removed,
            "electron outside central circle should be absorbed");
}
// Two different physical walls: fate, deposited endpoint charge and loss face
// must all follow the first hit, independently of JSON order.
template <ShapeType ST>
void mixed_walls(bool reverse_order, int implementation, bool both_absorbing = false) {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    ParticlesArray sp(check::species("Electrons", -1, 1, 1, 1), domain);
    const Vector3R start(3.9, 2, 0.1), velocity(0.6, 0, -0.3);
    sp.add_particle(Particle(start, velocity));
    auto conditions = nlohmann::json::parse(
        R"([{"open":[{"face":"XMAX"}]},{"electron_reflection":{"face":"ZMIN","radius":10,"energy_threshold":100}}])");
    if (both_absorbing)
        conditions[1] = nlohmann::json::parse(R"({"open":[{"face":"ZMIN"}]})");
    if (reverse_order)
        std::reverse(conditions.begin(), conditions.end());
    BoundaryConditionHandler handler;
    handler.load_from_json({{"Boundary_conditions", conditions}}, domain);
    ParticleBoundaryEvents events;
    const auto policy = handler.path_policy(sp, domain, &events);
    const Particle trial(start + velocity, velocity);
    const auto path = current_particle_path(start, trial, domain, policy);
    check::require(path.crossing.face == Face::XMAX && path.crossing.fate == ParticleFate::Removed && !path.reflected,
                   "mixed walls must absorb at first XMAX hit");
    check::near(path.crossing.t_hit, 1.0 / 6, 1e-14, "first wall hit fraction");
    check::vector(path.end, Vector3R(4, 2, 0.05), 1e-14, "current ends on XMAX before ZMIN");
    sp.density_on_grid_update<ST>();
    const Field3d before = sp.densityOnGrid;
    sp.currentOnGrid.setZero();
    if (implementation == 2) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2, 1, sp.currentOnGrid, work, policy);
    } else if (implementation == 1) {
        move_and_calc_current_reference<ST>(sp, 1, sp.currentOnGrid, policy);
    } else {
        move_and_calc_current<ST>(sp, 1, sp.currentOnGrid, policy);
    }
    check::near(events.size(), 1, 0, "only the outgoing particle consumes event storage");
    // The deposit must satisfy continuity to the wall endpoint. Removing that
    // endpoint charge afterwards represents the wall sink, not a bulk source.
    ParticlesArray wall_charge(check::species("Electrons", -1, 1, 1, 1), domain);
    wall_charge.add_particle(Particle(path.end, velocity));
    wall_charge.density_on_grid_update<ST>();
    double error = 0;
    for (int i = 1; i < domain.size().x() - 1; ++i)
        for (int j = 1; j < domain.size().y() - 1; ++j)
            for (int k = 1; k < domain.size().z() - 1; ++k) {
                const auto& J = sp.currentOnGrid;
                const double div = J(i, j, k, 0) - J(i - 1, j, k, 0) + J(i, j, k, 1) - J(i, j - 1, k, 1) +
                                   J(i, j, k, 2) - J(i, j, k - 1, 2);
                const double residual = wall_charge.densityOnGrid(i, j, k, 0) - before(i, j, k, 0) + div;
                check::require(std::isfinite(residual), "finite mixed-wall continuity");
                error = std::max(error, std::abs(residual));
            }
    check::near(error, 0, 1e-12, "mixed-wall continuity including wall endpoint charge");
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &events);
    check::near(sp.get_total_num_of_particles(), 0, 0, "particle absorbed at XMAX");
    check::near(sp.diag.boundary.at(Face::XMAX).lost_count, 1, 0, "loss recorded on first wall");
    check::require(!sp.diag.boundary.contains(Face::ZMIN), "no interaction with later ZMIN wall");
}

template <ShapeType ST>
void start_on_upper_wall(int implementation) {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    ParticlesArray sp(check::species("Electrons", -1, 1, 1, 1), domain);
    sp.add_particle(Particle(Vector3R(2, 2, 4), Vector3R(0.1, 0, 0.5)));
    BoundaryConditionHandler handler;
    handler.load_from_json(nlohmann::json::parse(R"({"Boundary_conditions":[{"open":[{"face":"ZMAX"}]}]})"), domain);
    ParticleBoundaryEvents events;
    sp.currentOnGrid.setZero();
    const auto policy = handler.path_policy(sp, domain, &events);
    if (implementation == 2) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2, 1, sp.currentOnGrid, work, policy);
    } else if (implementation == 1) {
        move_and_calc_current_reference<ST>(sp, 1, sp.currentOnGrid, policy);
    } else {
        move_and_calc_current<ST>(sp, 1, sp.currentOnGrid, policy);
    }
    check::near(sp.currentOnGrid.norm(), 0, 1e-15, "outward start on ZMAX has zero internal current path");
    check::near(events.size(), 1, 0, "t=0 hit is saved");
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &events);
    check::near(sp.get_total_num_of_particles(), 0, 0, "outward particle on ZMAX absorbed");
}
}   // namespace

int main() {
    try {
        check_captured_fate_swaps();
        for (int impl = 0; impl < 3; ++impl) {
            start_on_upper_wall<ShapeType::Linear>(impl);
            start_on_upper_wall<ShapeType::Quadratic>(impl);
        }
        for (bool reverse : {false, true})
            for (int impl = 0; impl < 3; ++impl) {
                mixed_walls<ShapeType::Linear>(reverse, impl);
                mixed_walls<ShapeType::Quadratic>(reverse, impl);
                mixed_walls<ShapeType::Linear>(reverse, impl, true);
                mixed_walls<ShapeType::Quadratic>(reverse, impl, true);
            }
        for (const std::string condition : {"electron_reflection", "bphi"})
            for (Face face : {Face::ZMIN, Face::ZMAX})
                for (int impl = 0; impl < 3; ++impl) {
                    exercise<ShapeType::Linear>(condition, face, impl == 1, impl == 2);
                    exercise<ShapeType::Quadratic>(condition, face, impl == 1, impl == 2);
                }
        std::cout << checks << " reflected-current checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
