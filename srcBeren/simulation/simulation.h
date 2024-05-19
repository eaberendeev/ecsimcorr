// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_H
#define SIMULATION_H


#include "ParticlesArray.h"
#include "Vec.h"
#include "World.h"
#include "parameters_map.h"

// Main simulation class
class Simulation {
   public:
    Simulation(int argc, char **argv);
    void make_all();
    void collect_charge_density(Field3d& field,
                                const std::vector<ParticlesArray>& species);

    void collect_current(Field3d& field,
                         const std::vector<ParticlesArray>& species);

   private:
   // Simulation parameters
      const ParametersMap parameters;
      Bounds bounds;
      Domain domain;
    // Mesh variables
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
   
    // // Particles
    // std::vector<ParticlesArray> species;

};

#endif
