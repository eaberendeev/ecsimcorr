#include "simulation_drift_kinetic.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "simulation/eccapfim/algorithms_eccapfim.h"
#include "utils/containers.h"

namespace {
void zero_vector_field(driftKinetic::VectorField& field) {
    std::fill(field.data().begin(), field.data().end(), Vector3R(0, 0, 0));
}

bool species_uses_drift_kinetic(const nlohmann::json& particles_config, const std::string& name) {
    if (!particles_config.contains("particles") || !particles_config["particles"].is_array())
        return false;
    for (const auto& entry : particles_config["particles"]) {
        if (!entry.is_object() || entry.value("Name", std::string()) != name)
            continue;
        const auto model = entry.value("model", std::string("kinetic"));
        return model == "drift_kinetic" || model == "drift-kinetic" || model == "guiding_center" ||
               model == "guiding-centre";
    }
    return false;
}

Vector3R vector_from_json(const nlohmann::json& value, const Vector3R& fallback) {
    if (!value.is_array() || value.size() != 3)
        return fallback;
    return Vector3R(value[0].get<double>(), value[1].get<double>(), value[2].get<double>());
}

int yee_node(int storage_index, int axis, const Domain& domain, const BoundaryConditionHandler& bc) {
    if (bc.is_periodic(axis))
        return domain.grid.periodic_node_owner(storage_index, axis);
    return std::clamp(storage_index, 0, domain.size()[axis] - 1);
}

Vector3R cell_field_at(const Field3d& field, const Vector3I& cell, const Domain& domain,
                       const BoundaryConditionHandler& bc) {
    Vector3R value(0, 0, 0);
    for (int di = 0; di < 2; ++di)
        for (int dj = 0; dj < 2; ++dj)
            for (int dk = 0; dk < 2; ++dk) {
                const int i = yee_node(cell.x() + di, 0, domain, bc);
                const int j = yee_node(cell.y() + dj, 1, domain, bc);
                const int k = yee_node(cell.z() + dk, 2, domain, bc);
                for (int c = 0; c < 3; ++c) value[c] += 0.125 * field(i, j, k, c);
            }
    return value;
}
}   // namespace

void SimulationDriftKinetic::init() {
    SimulationEccapfim::init();
    push_.set_tolerances(1e-14, 1e-13, 100);
    convert_initial_particles();
}

void SimulationDriftKinetic::init_fields() {
    SimulationEccapfim::init_fields();

    const auto kinetic = system_config.value("DriftKinetic", nlohmann::json::object());
    if (kinetic.contains("electric") && !system_config.contains("ExternalFieldE")) {
        const Vector3R electric = vector_from_json(kinetic["electric"], Vector3R(0, 0, 0));
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k)
                    for (int c = 0; c < 3; ++c) fieldE_external(i, j, k, c) = electric[c];
    }
    if (kinetic.contains("magnetic") && !system_config.contains("ExternalFieldB")) {
        const Vector3R magnetic = vector_from_json(kinetic["magnetic"], Vector3R(0, 0, 1));
        if (!(magnetic.length() > 0.0) || !std::isfinite(magnetic.length()))
            throw std::invalid_argument("drift_kinetic requires a finite nonzero DriftKinetic.magnetic field");
        for (int i = 0; i < domain.size().x(); ++i)
            for (int j = 0; j < domain.size().y(); ++j)
                for (int k = 0; k < domain.size().z(); ++k)
                    for (int c = 0; c < 3; ++c) fieldBInit(i, j, k, c) = magnetic[c];
    }

    const std::array<bool, 3> periodic{bc_handler.is_periodic(0), bc_handler.is_periodic(1),
                                       bc_handler.is_periodic(2)};
    grid_ = std::make_unique<driftKinetic::DriftKineticGrid>(domain, periodic);
    sync_grid_from_fields(fieldE + fieldE_external, fieldB + fieldBInit);
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k)
                grid_->sample(Vector3R((i - GHOST_CELLS + 0.5) * domain.cell_size(0),
                                       (j - GHOST_CELLS + 0.5) * domain.cell_size(1),
                                       (k - GHOST_CELLS + 0.5) * domain.cell_size(2)));
}

void SimulationDriftKinetic::sync_grid_from_fields(const Field3d& electric, const Field3d& magnetic) {
    zero_vector_field(grid_->electric());
    zero_vector_field(grid_->magnetic());
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                grid_->electric()(i, j, k) = cell_field_at(electric, Vector3I(i, j, k), domain, bc_handler);
                grid_->magnetic()(i, j, k) = cell_field_at(magnetic, Vector3I(i, j, k), domain, bc_handler);
            }
}

void SimulationDriftKinetic::add_cell_field_to_yee(const driftKinetic::VectorField& cell_field, Field3d& yee_field) {
    for (int i = 0; i < domain.size().x(); ++i)
        for (int j = 0; j < domain.size().y(); ++j)
            for (int k = 0; k < domain.size().z(); ++k) {
                const Vector3R value = cell_field(i, j, k);
                for (int di = 0; di < 2; ++di)
                    for (int dj = 0; dj < 2; ++dj)
                        for (int dk = 0; dk < 2; ++dk) {
                            const int ii = yee_node(i + di, 0, domain, bc_handler);
                            const int jj = yee_node(j + dj, 1, domain, bc_handler);
                            const int kk = yee_node(k + dk, 2, domain, bc_handler);
                            for (int c = 0; c < 3; ++c) {
#pragma omp atomic update
                                yee_field(ii, jj, kk, c) += 0.125 * value[c];
                            }
                        }
            }
    algorithmsECCAPFIM::copy_periodic_nodes(yee_field, domain, &bc_handler);
}

Vector3R SimulationDriftKinetic::drift_velocity(const driftKinetic::GuidingCenter& state,
                                               const driftKinetic::Fields& fields) const {
    const double qm = state.charge / state.mass;
    const Vector3R effective_B = fields.lenB * fields.b + (state.v_parallel / qm) * fields.rotB;
    const double effective_len = fields.b.dot(effective_B);
    if (!(std::abs(effective_len) > 0.0) || !std::isfinite(effective_len))
        throw std::runtime_error("drift_kinetic: degenerate effective magnetic field");
    const Vector3R effective_b = effective_B / effective_len;
    const Vector3R effective_force = state.charge * fields.electric - state.mu * fields.gradB;
    return state.v_parallel * effective_b +
           effective_force.cross(fields.b) / (state.charge * effective_len);
}

void SimulationDriftKinetic::convert_initial_particles() {
    const std::string scheme = system_config.value("Scheme", std::string("drift_kinetic"));
    const bool all_species = scheme == "drift_kinetic";
    for (auto it = species.begin(); it != species.end();) {
        const std::string name = it->first;
        auto& sp = *it->second;
        if (!all_species && !species_uses_drift_kinetic(particles_config, name)) {
            ++it;
            continue;
        }
        if (sp.is_neutral())
            throw std::invalid_argument("drift_kinetic requires charged particle species: " + name);
        auto particles = std::make_unique<driftKinetic::DriftKineticParticles>(
            domain, name, sp.charge(), sp.mass(), sp.mpw());
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
        charged_species.erase(name);
        particle_sources.erase(name);
        it = species.erase(it);
    }
}

void SimulationDriftKinetic::begin_nonlinear_step() {
    SimulationEccapfim::begin_nonlinear_step();
    initial_particles_.clear();
    trial_particles_.clear();
    for (const auto& [name, ptr] : particles_) initial_particles_[name] = ptr->clone();
}

void SimulationDriftKinetic::evaluate_current(const Field3d& electric, const Field3d& magnetic) {
    // Kinetic species fill fieldJ through the ECCAPFIM backend; DK species are
    // added below.  With no kinetic species this is simply fieldJ = 0.
    SimulationEccapfim::evaluate_current(electric, magnetic);
    sync_grid_from_fields(electric, magnetic);

    driftKinetic::VectorField current(domain.size()), magnetization(domain.size());
    zero_vector_field(current);
    zero_vector_field(magnetization);

    for (const auto& [name, initial] : initial_particles_) {
        auto trial = initial->clone();
        for (int cell = 0; cell < trial->size(); ++cell)
            for (auto& particle : (*trial)(cell)) {
                const auto start = particle;
                const auto result = push_.push(
                    particle, dt_,
                    [&](const Vector3R& a, const Vector3R& b) { return grid_->sample(0.5 * (a + b)); }, true);
                if (!result.converged)
                    throw std::runtime_error("drift_kinetic: guiding-centre push did not converge");
                particle = result.state;

                const auto fields = grid_->sample(0.5 * (start.position + particle.position));
                const Vector3R velocity = drift_velocity(particle, fields);
                grid_->deposit_current(current, particle.position, particle.charge * velocity, particle.weight);
                grid_->deposit_magnetization(magnetization, particle.position, fields.b, particle.mu, particle.weight);
            }
        trial->update_cells();
        trial_particles_[name] = std::move(trial);
    }

    const auto magnetization_curl = grid_->curl(magnetization);
    for (int cell = 0; cell < current.capacity(); ++cell) current.data()[cell] += magnetization_curl.data()[cell];
    add_cell_field_to_yee(current, fieldJ);
}

void SimulationDriftKinetic::accept_nonlinear_step() {
    SimulationEccapfim::accept_nonlinear_step();
    particles_ = std::move(trial_particles_);
}

void SimulationDriftKinetic::make_diagnostic(int timestep) {
    if (!species.empty()) SimulationEccapfim::make_diagnostic(timestep);
    sync_grid_from_fields(fieldE + fieldE_external, fieldB + fieldBInit);
    double energy = 0.0;
    int count = 0;
    for (const auto& [name, ptr] : particles_) {
        const auto& sp = *ptr;
        for (int cell = 0; cell < sp.size(); ++cell)
            for (const auto& particle : sp(cell)) {
                const double lenB = grid_->sample(particle.position).lenB;
                energy += 0.5 * particle.mass * particle.v_parallel * particle.v_parallel +
                          particle.mu * lenB;
                ++count;
            }
    }
    std::cout << "DriftKinetic step " << timestep << " energy=" << energy << " particles=" << count << "\n";
}

void SimulationDriftKinetic::finalize() {
    SimulationEccapfim::finalize();
}
