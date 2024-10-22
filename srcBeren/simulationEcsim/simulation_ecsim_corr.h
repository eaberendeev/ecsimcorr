// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_CORR_H
#define SIMULATION_ECSIM_CORR_H

#include "simulation.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "parameters_map.h"

// Main simulation class
class SimulationEcsimCorr: public Simulation{
    public:
     SimulationEcsimCorr(const ParametersMap& _systemParameters,
                         const std::vector<ParametersMap>& _speciesParameters,
                         const ParametersMap& _outputParameters, int argc,
                         char** argv)
         : Simulation(_systemParameters, _speciesParameters, _outputParameters,
                       argc, argv) {}
     void init_fields() override;
     void prepare_step(const int timestep) override;
     void make_step(const int timestep) override;
     //void output_all(const int timestep) override;
     void diagnostic_energy(Diagnostics& diagnostic,
                            const int timestep);
     void make_diagnostic(const int timestep) override;
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
