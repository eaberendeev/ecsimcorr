#include <cmath>
#include <iostream>
#include <stdexcept>

#include "ParticlesDiagnostic.h"
#include "algorithms_ecsim_corr.h"
#include "boundary_conditions.h"

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
    auto cfg = nlohmann::json::parse(
        R"({"Name":"Electrons","Mass":1,"NumPartPerCell":1,"distribution":[]})");
    ParticlesArray sp(cfg, domain);
    for (double vz : {-0.5, -5.0, -0.4})
        sp.add_particle(Particle(Vector3R(2, 2, -0.25), Vector3R(0, 0, vz)));
    BoundaryConditionHandler handler;
    handler.load_from_json(nlohmann::json::parse(
        R"({"Boundary_conditions":[{"electron_reflection":{"face":"ZMIN","radius":1,"energy_threshold":100}}]})"),
        domain);
    auto captured = handler.capture_particle_fates(sp, domain);
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
    sp.add_particle(Particle(Vector3R(2.1, 2.0, z), Vector3R(0.2, 0.1, vz)));
    BoundaryConditionHandler handler;
    nlohmann::json params = {{"face", lower ? "ZMIN" : "ZMAX"},
                             {"radius", 1.0},
                             {"energy_threshold", 100.0},
                             {"electron_threshold_energy_", 0.2},
                             {"gap", 0.0}};
    handler.load_from_json({{"Boundary_conditions", nlohmann::json::array({{{condition, params}}})}}, domain);
    const auto policy = handler.path_policy(sp, domain);
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
    auto fates = handler.capture_particle_fates(sp, domain);
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
}   // namespace

int main() {
    try {
        check_captured_fate_swaps();
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
