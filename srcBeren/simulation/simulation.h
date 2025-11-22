// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H

#include <nlohmann/json.hpp>

#include "Diagnostic.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "parameters_map.h"

// Main simulation class
class Simulation {
   public:
    Simulation(const ParametersMap& systemParameters, const nlohmann::json& particles_config,
               const ParametersMap& outputParameters, int argc, char** argv);
    Simulation(){};

    void collect_current(Field3d& J);
    void collect_charge_density(Field3d& field);
    virtual void calculate();
    virtual void finalize(){};

    virtual void init();
    virtual void init_particles(const nlohmann::json& config);
    virtual void init_fields(){};
    virtual void prepare_step(const int timestep){
        std::cout << "Prepare step is not implemented for timestep " << timestep << "\n";
    };
    virtual void collision_step(const int timestep) {
        std::cout << "Collision step is not implemented for timestep " << timestep
                  << "\n";
    };
    virtual void make_step(const int timestep){
        std::cout << "Make step is not implemented for timestep " << timestep
                  << "\n";
    };
    // virtual void diagnostic_energy(Diagnostics &diagnostic, const int timestep){
    //     std::cout << "Diagnostic energy is not implemented\n";
    // };
    virtual void make_diagnostic(const int timestep) {
        std::cout << "Make diagnostic is not implemented for timestep "
                  << timestep << "\n";
    };
    virtual ~Simulation() = default;
    void output_fields2D(
        const int timestep,
        const std::vector<std::pair<Field3d&, std::string>>& fields);

    virtual std::unique_ptr<ParticlesArray> make_particles_array(
        const nlohmann::json& config) {
        return std::make_unique<ParticlesArray>(config,
                                                parameters, domain);
    }
    // Simulation parameters
    ParametersMap parameters;
    ParametersMap outputParameters;
    nlohmann::json particles_config;
    Bounds bounds;
    Domain domain;

    // Fields mesh
    Mesh mesh;
    // Particles
    Species species;
    std::vector<std::reference_wrapper<ParticlesArray>> charged_species;
    // Writer writer;

    // Diagnostics diag;
    Timer globalTimer;

};

std::unique_ptr<Simulation> build_simulation(
    const ParametersMap& systemParameters,
    const nlohmann::json& particles_config,
    const ParametersMap& outputParameters, int argc, char** argv);
#endif
