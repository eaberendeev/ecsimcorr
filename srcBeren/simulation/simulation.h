// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H

#include <memory>
#include <nlohmann/json.hpp>

#include "Diagnostic.h"
#include "Mesh.h"
#include "ParticlesArray.h"
#include "Timer.h"
#include "boundary_conditions.h"
#include "collision_manager.h"
#include "containers.h"
#include "grid/domain.h"

class BinaryCollider;
class BinaryColliderWithNeutrals;

// Main simulation class
class Simulation {
   public:
    Simulation(const nlohmann::json& s_config, const nlohmann::json& p_config, int argc, char** argv);
    // Defined with the destructor where the unique_ptr collider types are complete.
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
    Domain domain;

    // Time step, read once in init() and reused by all per-step code.
    double dt_ = 0.0;

    // Fields mesh
    Mesh mesh;
    BoundaryConditionHandler bc_handler;

    // Particles
    Species species;
    std::map<std::string, std::reference_wrapper<ParticlesArray>> charged_species;

    // Collisions
    CollisionManager collision_manager;

    // Legacy "Collider" fallback state. Kept per Simulation instance so random
    // streams do not leak between runs; constructed on first use.
    std::unique_ptr<BinaryColliderWithNeutrals> binary_collider_with_neutrals_;
    std::unique_ptr<BinaryCollider> binary_collider_;

    Timer globalTimer;
};

std::unique_ptr<Simulation> build_simulation(const nlohmann::json& system_config,
                                             const nlohmann::json& particles_config, int argc, char** argv);
#endif
