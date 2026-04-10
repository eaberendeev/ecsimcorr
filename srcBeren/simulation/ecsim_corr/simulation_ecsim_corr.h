// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_CORR_H
#define SIMULATION_ECSIM_CORR_H

#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "simulation.h"
#include "simulation_ecsim.h"
// Main simulation class
class SimulationEcsimCorr : public SimulationEcsim {
   public:
    SimulationEcsimCorr(const nlohmann::json& system_config,
                        const nlohmann::json& particles_config, int argc,
                        char** argv)
        : SimulationEcsim(system_config, particles_config, argc, argv) {}
    void make_step(const int timestep) override;
    //  void make_stepNGP(const int timestep) override;
    void correctv(ParticlesArray& sort, const double dt);
    void diagnostic_energy(Diagnostics& diagnostic) override;
};

#endif
