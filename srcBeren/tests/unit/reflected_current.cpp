#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "config/boundary_config.h"
#include "coupling/current.h"
#include "coupling/density.h"
#include "coupling/moments.h"
#include "simulation/boundary_conditions.h"
#include "simulation/ecsim_corr/algorithms_ecsim_corr.h"
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

const char* wall_name(Face face) {
    switch (face) {
        case Face::ZMIN:
            return "ZMIN";
        case Face::ZMAX:
            return "ZMAX";
        case Face::XMAX:
            return "XMAX";
        default:
            return "other wall";
    }
}

void frozen_boundary_fates_survive_particle_removal() {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    auto cfg = nlohmann::json::parse(R"({"Name":"Electrons","Mass":1,"NumPartPerCell":1,"distribution":[]})");
    ParticlesArray sp(cfg, domain);
    for (double x : {1.0, 2.0, 3.0})
        for (double vz : {-0.5, -5.0, -0.4, -6.0}) sp.add_particle(Particle(Vector3R(x, 2, -0.25), Vector3R(0, 0, vz)));
    BoundaryConditionHandler handler;
    handler = boundary_conditions_from_config(
        nlohmann::json::parse(
            R"({"Boundary_conditions":[{"electron_reflection":{"face":"ZMIN","radius":1,"energy_threshold":100}}]})"),
        domain);
    ParticleBoundaryEvents captured;
    // Each cell has three frozen decisions and one unsaved absorbing particle.
    // Different cells can be recorded and consumed by different OpenMP threads.
#pragma omp parallel for schedule(dynamic, 1)
    for (int cell = 0; cell < sp.size(); ++cell)
        for (size_t i = 0; i + 1 < sp(cell).size(); ++i)
            captured.record(sp(cell)[i], handler.classify_particle(sp(cell)[i], sp, domain));
    require(captured.size() == 9, "snapshot stores three decisions per boundary cell");
    for (int cell = 0; cell < sp.size(); ++cell)
        for (Particle& p : sp(cell)) p.velocity *= 100;
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &captured);
    require(captured.size() == 9, "boundary handling must not mutate the saved decisions");
    require(sp.get_total_num_of_particles() == 6, "captured fates lost alignment during particle swaps");
    require(sp.step_state.boundary.at(Face::ZMIN).reflected_count == 6, "wrong frozen reflection statistics");
    require(sp.step_state.boundary.at(Face::ZMIN).lost_count == 6, "wrong frozen absorption statistics");
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
    handler = boundary_conditions_from_config({{"Boundary_conditions", nlohmann::json::array({{{condition, params}}})}},
                                              domain);
    ParticleBoundaryEvents fates;
    const auto policy = handler.path_policy(sp, domain, &fates);
    Field3d density(domain.size(), 1);
    coupling::deposit_density<ST>(sp, density);
    const Field3d before = density;
    sp.currentOnGrid.setZero();
    if (method == DepositMethod::BorisAndDeposit) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2.0, 1.0, sp.currentOnGrid, work, policy);
    } else if (method == DepositMethod::Reference) {
        coupling::drift_and_deposit_current_reference<ST>(sp, 1.0, sp.currentOnGrid, policy);
    } else {
        coupling::drift_and_deposit_current<ST>(sp, 1.0, sp.currentOnGrid, policy);
    }
    require(fates.size() == 1, "deposit must save just the outgoing particle event");
    // Cross the reflection threshold after depositing: the stored fate must win.
    for (int cell = 0; cell < sp.size(); ++cell)
        for (Particle& p : sp(cell)) p.velocity *= 100;
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &fates);
    sp.update_cells(domain);
    require(sp.get_total_num_of_particles() == 1, "reflected electron disappeared");
    for (int cell = 0; cell < sp.size(); ++cell) {
        for (const Particle& p : sp(cell)) {
            require(std::abs(p.coord.z() - z) < 1e-14, "wrong reflected endpoint");
            require(p.velocity.z() * vz < 0, "velocity did not reflect");
        }
    }
    coupling::deposit_density<ST>(sp, density);
    double max_residual = 0;
    for (int i = 1; i < domain.size().x() - 1; ++i)
        for (int j = 1; j < domain.size().y() - 1; ++j)
            for (int k = 1; k < domain.size().z() - 1; ++k) {
                const auto& current = sp.currentOnGrid;
                const double div = current(i, j, k, 0) - current(i - 1, j, k, 0) + current(i, j, k, 1) -
                                   current(i, j - 1, k, 1) + current(i, j, k, 2) - current(i, j, k - 1, 2);
                const double residual = density(i, j, k, 0) - before(i, j, k, 0) + div;
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
    require(!path.reflected && std::abs(path.end.z() - (lower ? -0.5 : 4.5)) < 1e-14,
            "absorbed current must end half a cell beyond the wall");
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
    handler = boundary_conditions_from_config({{"Boundary_conditions", conditions}}, domain);
    ParticleBoundaryEvents events;
    const auto policy = handler.path_policy(sp, domain, &events);
    const Particle trial(start + velocity, velocity);
    const auto path = current_particle_path(start, trial, domain, policy);
    check::require(path.crossing.face == Face::XMAX && path.crossing.fate == ParticleFate::Removed && !path.reflected,
                   "mixed walls must absorb at first XMAX hit");
    check::near(path.crossing.t_hit, 1.0 / 6, 1e-14, "first wall hit fraction");
    check::vector(path.end, Vector3R(4.5, 2, 0.05), 1e-14, "current ends beyond XMAX before ZMIN");
    Field3d density(domain.size(), 1);
    coupling::deposit_density<ST>(sp, density);
    const Field3d before = density;
    sp.currentOnGrid.setZero();
    if (method == DepositMethod::BorisAndDeposit) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2, 1, sp.currentOnGrid, work, policy);
    } else if (method == DepositMethod::Reference) {
        coupling::drift_and_deposit_current_reference<ST>(sp, 1, sp.currentOnGrid, policy);
    } else {
        coupling::drift_and_deposit_current<ST>(sp, 1, sp.currentOnGrid, policy);
    }
    check::near(events.size(), 1, 0, "only the outgoing particle consumes event storage");
    // The deposit must satisfy continuity to the wall endpoint. Removing that
    // endpoint charge afterwards represents the wall sink, not a bulk source.
    ParticlesArray wall_charge(check::species("Electrons", -1, 1, 1, 1), domain);
    wall_charge.add_particle(Particle(path.end, velocity));
    coupling::deposit_density<ST>(wall_charge, density);
    double error = 0;
    for (int i = 1; i < domain.size().x() - 1; ++i)
        for (int j = 1; j < domain.size().y() - 1; ++j)
            for (int k = 1; k < domain.size().z() - 1; ++k) {
                const auto& J = sp.currentOnGrid;
                const double div = J(i, j, k, 0) - J(i - 1, j, k, 0) + J(i, j, k, 1) - J(i, j - 1, k, 1) +
                                   J(i, j, k, 2) - J(i, j, k - 1, 2);
                const double residual = density(i, j, k, 0) - before(i, j, k, 0) + div;
                check::require(std::isfinite(residual), "finite mixed-wall continuity");
                error = std::max(error, std::abs(residual));
            }
    check::near(error, 0, 1e-12, "mixed-wall continuity including wall endpoint charge");
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &events);
    check::near(sp.get_total_num_of_particles(), 0, 0, "particle absorbed at XMAX");
    check::near(sp.step_state.boundary.at(Face::XMAX).lost_count, 1, 0, "loss recorded on first wall");
    check::require(!sp.step_state.boundary.contains(Face::ZMIN), "no interaction with later ZMIN wall");
}

// An absorbed particle takes all its charge off the interior nodes: the density
// implied by continuity, rho - div J, vanishes strictly inside the walls. This
// includes a start on the wall itself, where the quadratic shape still reaches
// the first interior node.
template <ShapeType ST>
void absorbed_charge_leaves_interior_nodes(Face face, const Vector3R& start, const Vector3R& velocity,
                                           DepositMethod method) {
    Domain domain;
    domain.init(Vector3I(4, 4, 4), Vector3R(1, 1, 1));
    ParticlesArray sp(check::species("Electrons", -1, 1, 1, 1), domain);
    sp.add_particle(Particle(start, velocity));
    BoundaryConditionHandler handler;
    handler = boundary_conditions_from_config(
        nlohmann::json::parse(R"({"Boundary_conditions":[{"open":[{"face":"XMIN"},{"face":"XMAX"},)"
                              R"({"face":"ZMIN"},{"face":"ZMAX"}]},{"periodic":{"face":"YMIN"}}]})"),
        domain);
    ParticleBoundaryEvents events;
    const PathPolicy policy = handler.path_policy(sp, domain, &events);
    Field3d density(domain.size(), 1);
    coupling::deposit_density<ST>(sp, density);
    sp.currentOnGrid.setZero();
    if (method == DepositMethod::BorisAndDeposit) {
        Field3d zero(domain.size(), 3);
        zero.setZero();
        double work = 0;
        fused_push_and_deposit<ST>(sp, zero, zero, 2, 1, sp.currentOnGrid, work, policy);
    } else if (method == DepositMethod::Reference) {
        coupling::drift_and_deposit_current_reference<ST>(sp, 1, sp.currentOnGrid, policy);
    } else {
        coupling::drift_and_deposit_current<ST>(sp, 1, sp.currentOnGrid, policy);
    }
    check::near(events.size(), 1, 0, "the particle hits the wall");
    const Field3d& J = sp.currentOnGrid;
    double interior = 0;
    for (int i = 1; i < domain.size().x() - 1; ++i)
        for (int j = 1; j < domain.size().y() - 1; ++j)
            for (int k = 1; k < domain.size().z() - 1; ++k) {
                const Vector3R node = domain.storage_node_position(i, j, k, FieldType::DENSITY, 0);
                if (node.x() <= 0 || node.x() >= 4 || node.z() <= 0 || node.z() >= 4)
                    continue;
                const double div = J(i, j, k, 0) - J(i - 1, j, k, 0) + J(i, j, k, 1) - J(i, j - 1, k, 1) +
                                   J(i, j, k, 2) - J(i, j, k - 1, 2);
                interior = std::max(interior, std::abs(density(i, j, k, 0) - div));
            }
    check::near(interior, 0, 1e-12, std::string("charge left on interior nodes at ") + wall_name(face));
    Species empty;
    handler.apply_to_particles(sp, empty, domain, &events);
    check::near(sp.get_total_num_of_particles(), 0, 0, "particle absorbed");
    check::near(sp.step_state.boundary.at(face).lost_count, 1, 0, "loss recorded on the hit wall");
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
            for (Face face : {Face::ZMIN, Face::ZMAX, Face::XMAX}) {
                const bool upper = face != Face::ZMIN;
                const Vector3R velocity = face == Face::XMAX ? Vector3R(0.5, 0, 0.1) : Vector3R(0.1, 0, upper ? 0.5 : -0.5);
                const Vector3R near_wall = face == Face::XMAX ? Vector3R(3.8, 2, 2) : Vector3R(2, 2, upper ? 3.8 : 0.2);
                const Vector3R on_wall = face == Face::XMAX ? Vector3R(4, 2, 2) : Vector3R(2, 2, upper ? 4 : 0);
                const std::string where = std::string(wall_name(face)) + ", " + name;
                check::scenario("Linear absorbed charge leaves interior nodes, " + where, [&] {
                    absorbed_charge_leaves_interior_nodes<ShapeType::Linear>(face, near_wall, velocity, method);
                });
                check::scenario("Quadratic absorbed charge leaves interior nodes, " + where, [&] {
                    absorbed_charge_leaves_interior_nodes<ShapeType::Quadratic>(face, near_wall, velocity, method);
                });
                check::scenario("Quadratic start on the wall leaves interior nodes, " + where, [&] {
                    absorbed_charge_leaves_interior_nodes<ShapeType::Quadratic>(face, on_wall, velocity, method);
                });
            }
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
