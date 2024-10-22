// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H

#include "Diagnostic.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "parameters_map.h"

// Main simulation class
class Simulation {
   public:
    Simulation(const ParametersMap& systemParameters,
               const std::vector<ParametersMap>& speciesParameters,
               const ParametersMap& outputParameters, int argc, char** argv);
    Simulation(){};

    void inject_particles(const int timestep);
    void collect_current(Field3d& J);
    void collect_charge_density(Field3d& field);
    virtual void make_all();

    virtual void init();
    virtual void init_particles();
    virtual void init_fields(){};
    virtual void prepare_step(const int timestep){};
    virtual void output_all(const int timestep);
    virtual void make_step(const int timestep){};
    // virtual void diagnostic_energy(Diagnostics &diagnostic, const int timestep){
    //     std::cout << "Diagnostic energy is not implemented\n";
    // };
    virtual void make_diagnostic(const int timestep) {
        std::cout << "Diagnostic is not implemented\n";
    };
    virtual ~Simulation() = default;
    void output_fields2D(
        const int timestep,
        const std::vector<std::pair<Field3d&, std::string>>& fields);

    // Simulation parameters
    ParametersMap parameters;
    std::vector<ParametersMap> speciesParameters;
    ParametersMap outputParameters;
    Bounds bounds;
    Domain domain;

    // Fields mesh
    Mesh mesh;
    // Particles
    std::vector<ParticlesArray> species;
    // Writer writer;

    // Diagnostics diag;
    Timer globalTimer;

};

#endif
