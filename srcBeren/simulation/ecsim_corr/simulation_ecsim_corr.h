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
#include "simulation_ecsim.h"
// Main simulation class
class SimulationEcsimCorr: public SimulationEcsim{
    public:
     SimulationEcsimCorr(const ParametersMap& _systemParameters,
                         const std::vector<ParametersMap>& _speciesParameters,
                         const ParametersMap& _outputParameters, int argc,
                         char** argv)
         : SimulationEcsim(_systemParameters, _speciesParameters, _outputParameters,
                       argc, argv) {}
     void make_step(const int timestep) override;
   //  void make_stepNGP(const int timestep) override;
     void diagnostic_energy(Diagnostics& diagnostic) override;
     void make_diagnostic(const int timestep) override;
};

#endif
