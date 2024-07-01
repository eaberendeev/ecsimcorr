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
    // Mesh variables 
    // TO DO: move it from Mesh cless 
    // Field3d fieldE;
    // Field3d fieldEn;
    // Field3d fieldEp;
    // Field3d fieldB;
    // Field3d fieldJp; // predict current for EM solver
    // Field3d fieldJp_full; // predict current for EM solver Jp + Lmat(E+E_n);
    // Field3d fieldJe; // Esirkepov current for E correction
    // Field3d fieldB0;
    // Field3d fieldBInit;
    // std::vector<IndexMap> LmatX;

    // //Sources and fields on the grid
    // Field3d chargeDensityOld;
    // Field3d chargeDensity;

    // Fields mesh
    Mesh mesh;
    // Particles
    std::vector<ParticlesArray> species;
};

#endif
