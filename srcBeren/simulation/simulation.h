// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H

#include <memory>
#include <nlohmann/json.hpp>

#include "collisions/collision_manager.h"
#include "diagnostics/Diagnostic.h"
#include "fields/maxwell_operators.h"
#include "grid/domain.h"
#include "particles/particle_array.h"
#include "particles/particle_sources.h"
#include "simulation/boundary_conditions.h"
#include "utils/containers.h"

// Main simulation class
class Simulation {
   public:
    Simulation(const nlohmann::json& s_config, const nlohmann::json& p_config, int argc, char** argv);
    Simulation();

    void advance_neutral_particles(double dt);

    virtual ShapeType diagnostic_density_shape() const {
        return SHAPE;
    }

    void collect_current(Field3d& J);
    virtual void calculate();
    virtual void finalize();

    virtual void init();
    virtual void init_particles(const nlohmann::json& config);
    virtual void init_fields() {};

    virtual void prepare_step(const int timestep) {
        std::cout << "Prepare step is not implemented for timestep " << timestep << "\n";
    };

    virtual void collision_step(const int timestep);

    virtual void make_step(const int timestep) {
        std::cout << "Make step is not implemented for timestep " << timestep << "\n";
    };

    virtual void make_diagnostic(const int timestep) {
        std::cout << "Make diagnostic is not implemented for timestep " << timestep << "\n";
    };

    virtual ~Simulation();

    virtual std::unique_ptr<ParticlesArray> make_particles_array(const nlohmann::json& p_config) {
        return std::make_unique<ParticlesArray>(p_config, domain);
    }

    nlohmann::json system_config;
    nlohmann::json particles_config;
    // Declared before species so it outlives their borrowed references during destruction.
    Domain domain;

    // Time step, read once in init() and reused by all per-step code.
    double dt_ = 0.0;

    // Spatial field operators shared by the selected scheme.
    MaxwellOperators maxwell;
    BoundaryConditionHandler bc_handler;

    // Particles
    Species species;
    std::unordered_map<std::string, ParticleSources> particle_sources;
    std::map<std::string, std::reference_wrapper<ParticlesArray>> charged_species;

    // Collisions
    CollisionManager collision_manager;
};

std::unique_ptr<Simulation> build_simulation(const nlohmann::json& system_config,
                                             const nlohmann::json& particles_config, int argc, char** argv);
#endif
