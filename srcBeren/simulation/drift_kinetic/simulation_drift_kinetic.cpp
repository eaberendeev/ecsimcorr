#include "simulation_drift_kinetic.h"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
Vector3R vector_from_json(const nlohmann::json& value, const Vector3R& fallback) {
    if (!value.is_array() || value.size() != 3)
        return fallback;
    return Vector3R(value[0].get<double>(), value[1].get<double>(), value[2].get<double>());
}
}   // namespace

void SimulationDriftKinetic::init_fields() {
    const std::array<bool, 3> periodic{bc_handler.is_periodic(0), bc_handler.is_periodic(1),
                                       bc_handler.is_periodic(2)};
    grid_ = std::make_unique<driftKinetic::DriftKineticGrid>(domain, periodic);
    const auto kinetic = system_config.value("DriftKinetic", nlohmann::json::object());
    const Vector3R electric = vector_from_json(kinetic.value("electric", nlohmann::json::array()), Vector3R(0, 0, 0));
    const Vector3R magnetic = vector_from_json(kinetic.value("magnetic", nlohmann::json::array()), Vector3R(0, 0, 1));
    if (!(magnetic.length() > 0.0) || !std::isfinite(magnetic.length()))
        throw std::invalid_argument("drift_kinetic requires a finite nonzero DriftKinetic.magnetic field");
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                grid_->electric()(i, j, k) = electric;
                grid_->magnetic()(i, j, k) = magnetic;
            }
}

void SimulationDriftKinetic::init() {
    Simulation::init();
    push_.set_tolerances(1e-14, 1e-13, 100);
    convert_initial_particles();
}

void SimulationDriftKinetic::convert_initial_particles() {
    for (auto& [name, sp_ptr] : species) {
        auto& sp = *sp_ptr;
        if (sp.is_neutral())
            throw std::invalid_argument("drift_kinetic requires charged particle species: " + name);
        auto particles = std::make_unique<driftKinetic::DriftKineticParticles>(domain, name, sp.charge(), sp.mass(),
                                                                               sp.mpw());
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const Particle& particle : sp(cell)) {
                const auto fields = grid_->sample(particle.coord);
                const double v_parallel = particle.velocity.dot(fields.b);
                const Vector3R v_perpendicular = particle.velocity - v_parallel * fields.b;
                driftKinetic::GuidingCenter guiding;
                guiding.position = particle.coord;
                guiding.v_parallel = v_parallel;
                guiding.mu = 0.5 * sp.mass() * v_perpendicular.squaredNorm() / fields.lenB;
                guiding.charge = sp.charge();
                guiding.mass = sp.mass();
                guiding.weight = sp.mpw();
                particles->add(guiding);
            }
        particles_.emplace(name, std::move(particles));
    }
    species.clear();
    charged_species.clear();
    particle_sources.clear();
}

void SimulationDriftKinetic::make_step(int timestep) {
    (void) timestep;
    const auto fields = [&](const Vector3R& start, const Vector3R& end) {
        return grid_->sample(0.5 * (start + end));
    };
    for (auto& [name, ptr] : particles_) {
        auto& sp = *ptr;
        for (int cell = 0; cell < sp.size(); ++cell)
            for (auto& particle : sp(cell)) {
                const auto result = push_.push(particle, dt_, fields, true);
                if (!result.converged)
                    throw std::runtime_error("drift_kinetic: guiding-center push did not converge");
                particle = result.state;
            }
        sp.update_cells();
    }
}

void SimulationDriftKinetic::make_diagnostic(int timestep) {
    double energy = 0.0;
    for (const auto& [name, ptr] : particles_) {
        const auto& sp = *ptr;
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const auto& particle : sp(cell)) {
                const double lenB = grid_->sample(particle.position).lenB;
                energy += 0.5 * particle.mass * particle.v_parallel * particle.v_parallel + particle.mu * lenB;
            }
    }
    std::cout << "DriftKinetic step " << timestep << " energy=" << energy << "\n";
}

void SimulationDriftKinetic::finalize() {
    Simulation::finalize();
}
