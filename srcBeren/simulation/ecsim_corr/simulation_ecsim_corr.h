// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_CORR_H
#define SIMULATION_ECSIM_CORR_H

#include <map>
#include <string>

#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "simulation.h"
#include "simulation_ecsim.h"

// Main simulation class
class SimulationEcsimCorr : public SimulationEcsim {
   public:
    SimulationEcsimCorr(const nlohmann::json& system_config, const nlohmann::json& particles_config, int argc,
                        char** argv)
        : SimulationEcsim(system_config, particles_config, argc, argv) {
    }
    void make_step(const int timestep) override;
    void deposit_density(ParticlesArray& particles) override;
    void advance_velocities_and_deposit_second_half_current();
    void rescale_velocities_for_energy_balance(ParticlesArray& sort, const double dt);
    void diagnostic_energy(Diagnostics& diagnostic) override;
    void solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J,
                                                     const double dt);

   private:
    void drift_half_step_and_assemble_response();
    void correct_particle_energy_and_deposit_density();

    // Power in code energy/time units from the Boris update. Overwritten per
    // charged species in the second half-step, consumed in the same step;
    // scratch state is recomputed after restart.
    std::map<std::string, double> predictor_power_by_species_;
    std::map<std::string, ParticleBoundaryEvents> boundary_events_;
};

#endif
