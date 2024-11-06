// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_IMPLICIT_H
#define SIMULATION_IMPLICIT_H

#include "simulation.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "parameters_map.h"

// Main simulation class
class SimulationImplicit: public Simulation{
    public:
     SimulationImplicit(const ParametersMap& _systemParameters,
                        const std::vector<ParametersMap>& _speciesParameters,
                        const ParametersMap& _outputParameters, int argc,
                        char** argv)
         : Simulation(_systemParameters, _speciesParameters,
                      _outputParameters, argc, argv) {
         std::cout << "Implicit simulation\n";
     }
     void init_fields() override;
     void prepare_step(const int timestep) override;
     void make_step(const int timestep) override;
     //void output_all(const int timestep) override;
     void init_particles() override;
     void diagnostic_energy(Diagnostics& diagnostic);
     void make_diagnostic(const int timestep) override;

     Field3d fieldJ;
     Field3d fieldE;
     Field3d fieldEn;
     Field3d fieldB;
     Field3d fieldBn;
     Field3d fieldE05;
     Field3d fieldB05;
     Field3d fieldBInit;
     Field3d fieldBFull;
};

#endif
