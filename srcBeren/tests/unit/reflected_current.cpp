#include <cmath>
#include <iostream>
#include <stdexcept>

#include "ParticlesDiagnostic.h"
#include "algorithms_ecsim_corr.h"
#include "boundary_conditions.h"
#include "test_support.h"

namespace {
using check::require;

enum class DepositMethod { Buffered, Reference, BorisAndDeposit };
const char* method_name(DepositMethod method) {
    switch (method) {
        case DepositMethod::Buffered:
            return "buffered deposit";
        case DepositMethod::Reference:
            return "reference deposit";
        case DepositMethod::BorisAndDeposit:
            return "Boris + deposit";
    }
    throw std::runtime_error("unknown deposit method");
}

void frozen_boundary_fates_survive_particle_removal() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    auto cfg = nlohmann::json::parse(R"({"Name":"Electrons","Mass":1,"NumPartPerCell":1,"distribution":[]})");
    ParticlesArray sp(cfg, domain);
    for (double x : {1.0, 2.0, 3.0})
        for (double vz : {-0.5, -5.0, -0.4, -6.0}) sp.add_particle(Particle(Vector3R(x, 2, -0.25), Vector3R(0, 0, vz)));
    BoundaryConditionHandler handler;
    handler.load_from_json(
        nlohmann::json::parse(
            R"({"Boundary_conditions":[{"electron_reflection":{"face":"ZMIN","radius":1,"energy_threshold":100}}]})"),
        domain);
    ParticleBoundaryEvents captured;
    // Each cell has three frozen decisions and one unsaved absorbing particle.
    // Different cells can be recorded and consumed by different OpenMP threads.
#pragma omp parallel for schedule(dynamic, 1)
    for (int cell = 0; cell < sp.size(); ++cell)
        for (size_t i = 0; i + 1 < sp.particlesData(cell).size(); ++i)
            captured.record(sp.particlesData(cell)[i],
                            handler.classify_particle(sp.particlesData(cell)[i], sp, domain));
    require(captured.size() == 9, "snapshot stores three decisions per boundary cell");
    for (int cell = 0; cell < sp.size(); ++cell)
        for (auto& p : sp.particlesData(cell)) p.velocity *= 100;
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &captured);
    require(captured.size() == 9, "boundary handling must not mutate the saved decisions");
    require(sp.get_total_num_of_particles() == 6, "captured fates lost alignment during particle swaps");
    require(sp.diag.boundary.at(Face::ZMIN).reflected_count == 6, "wrong frozen reflection statistics");
    require(sp.diag.boundary.at(Face::ZMIN).lost_count == 6, "wrong frozen absorption statistics");
    cfg["Name"] = "Ions";
    ParticlesArray ions(cfg, domain);
    const Particle trial(Vector3R(2, 2, -0.25), Vector3R(0, 0, -0.01));
    require(handler.classify_particle(trial, ions, domain).fate == ParticleFate::Removed,
            "electron reflection must not reflect ions");
}

template <ShapeType ST>
void reflected_electron_obeys_continuity(const std::string& condition, Face face, DepositMethod method) {
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
    if (method == DepositMethod::BorisAndDeposit) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2.0, 1.0, sp.currentOnGrid, work, policy);
    } else if (method == DepositMethod::Reference) {
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
void first_wall_determines_current_and_absorption(bool reverse_order, DepositMethod method,
                                                  bool both_absorbing = false) {
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
    if (method == DepositMethod::BorisAndDeposit) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2, 1, sp.currentOnGrid, work, policy);
    } else if (method == DepositMethod::Reference) {
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
void outward_wall_start_has_no_internal_current(DepositMethod method) {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    ParticlesArray sp(check::species("Electrons", -1, 1, 1, 1), domain);
    sp.add_particle(Particle(Vector3R(2, 2, 4), Vector3R(0.1, 0, 0.5)));
    BoundaryConditionHandler handler;
    handler.load_from_json(nlohmann::json::parse(R"({"Boundary_conditions":[{"open":[{"face":"ZMAX"}]}]})"), domain);
    ParticleBoundaryEvents events;
    sp.currentOnGrid.setZero();
    const auto policy = handler.path_policy(sp, domain, &events);
    if (method == DepositMethod::BorisAndDeposit) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2, 1, sp.currentOnGrid, work, policy);
    } else if (method == DepositMethod::Reference) {
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
    return check::run([] {
        check::scenario("Saved absorption/reflection fates survive swap-and-pop particle removal",
                        frozen_boundary_fates_survive_particle_removal);
        const DepositMethod methods[] = {DepositMethod::Buffered, DepositMethod::Reference,
                                         DepositMethod::BorisAndDeposit};
        for (auto method : methods) {
            const std::string name = method_name(method);
            check::scenario("Linear outward start on ZMAX has zero internal current: " + name,
                            [&] { outward_wall_start_has_no_internal_current<ShapeType::Linear>(method); });
            check::scenario("Quadratic outward start on ZMAX has zero internal current: " + name,
                            [&] { outward_wall_start_has_no_internal_current<ShapeType::Quadratic>(method); });
        }
        for (bool reverse : {false, true})
            for (auto method : methods) {
                const std::string name =
                    std::string(method_name(method)) + (reverse ? ", reversed JSON order" : ", original JSON order");
                check::scenario("Linear first-wall absorption and continuity: " + name, [&] {
                    first_wall_determines_current_and_absorption<ShapeType::Linear>(reverse, method);
                });
                check::scenario("Quadratic first-wall absorption and continuity: " + name, [&] {
                    first_wall_determines_current_and_absorption<ShapeType::Quadratic>(reverse, method);
                });
                check::scenario("Linear first-wall choice with two absorbing walls: " + name, [&] {
                    first_wall_determines_current_and_absorption<ShapeType::Linear>(reverse, method, true);
                });
                check::scenario("Quadratic first-wall choice with two absorbing walls: " + name, [&] {
                    first_wall_determines_current_and_absorption<ShapeType::Quadratic>(reverse, method, true);
                });
            }
        for (const std::string condition : {"electron_reflection", "bphi"})
            for (Face face : {Face::ZMIN, Face::ZMAX})
                for (auto method : methods) {
                    const std::string name =
                        condition + (face == Face::ZMIN ? " ZMIN: " : " ZMAX: ") + method_name(method);
                    check::scenario("Linear reflected-current continuity and frozen fate, " + name, [&] {
                        reflected_electron_obeys_continuity<ShapeType::Linear>(condition, face, method);
                    });
                    check::scenario("Quadratic reflected-current continuity and frozen fate, " + name, [&] {
                        reflected_electron_obeys_continuity<ShapeType::Quadratic>(condition, face, method);
                    });
                }
    });
}
