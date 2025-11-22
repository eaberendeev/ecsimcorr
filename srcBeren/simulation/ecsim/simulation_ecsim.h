// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_H
#define SIMULATION_ECSIM_H

#include "simulation.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "parameters_map.h"

// Main simulation class
class SimulationEcsim: public Simulation{
    public:
     SimulationEcsim(const ParametersMap& _systemParameters,
                     const nlohmann::json& particles_config,
                     const ParametersMap& _outputParameters, int argc,
                     char** argv)
         : Simulation(_systemParameters, particles_config, _outputParameters,
                      argc, argv) {}
     void init_fields() override;
     void prepare_step(const int timestep) override;
     void collision_step(const int timestep) override;
     void make_step(const int timestep) override;
     void make_stepNGP(const int timestep);
     virtual void diagnostic_energy(Diagnostics& diagnostic);
     void make_diagnostic(const int timestep) override;
     void prepare_block_matrix(ShapeType type);
     void convert_block_matrix(ShapeType type);
     void first_push();
     void second_push();

     Field3d fieldJp;        // predict current for EM solver
     Field3d fieldJp_full;   // predict current for EM solver Jp + Lmat(E+E_n);
     Field3d fieldJe;        // Esirkepov current for E correction};
     Field3d fieldE;
     Field3d fieldEn;
     Field3d fieldEp;
     Field3d fieldB;
     Field3d fieldBn;
     Field3d fieldBInit;
     Field3d fieldBFull;
};

#endif
