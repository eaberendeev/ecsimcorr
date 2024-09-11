// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H

#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "parameters_map.h"

// Main simulation class
class Simulation {
   public:
    Simulation(int argc, char** argv);
    void make_all();
    void collect_charge_density(Field3d& field,
                                const std::vector<ParticlesArray>& species);

    void collect_current(Field3d& field,
                         const std::vector<ParticlesArray>& species);

    void init(Mesh& mesh);
    void set_particles();
    void set_fields(Mesh& mesh);
    void inject_particles(const int timestep);
    void make_folders() const;

   private:
    void output_all(const Mesh& mesh,
                    const std::vector<ParticlesArray>& species, int timestep);
    void set_init_particles_distribution_rectangle(
        ParticlesArray& sp);
    void set_init_particles_distribution_cilinder(
        ParticlesArray& sp);

    // Simulation parameters
    ParametersMap parameters;
    Bounds bounds;
    Domain domain;

    // Fields mesh
    Mesh mesh;
    // Particles
    std::vector<ParticlesArray> species;
};

#endif
