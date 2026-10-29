// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_CORR_H
#define SIMULATION_ECSIM_CORR_H

#include <map>
#include <string>

#include "ParticlesArray.h"
#include "containers.h"
#include "grid/Domain.h"
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
    ShapeType diagnostic_density_shape() const override {
        return SHAPE_CH;
    }
    void rescale_velocities_for_energy_balance(ParticlesArray& sort, const double dt, double predictor_power);
    void diagnostic_energy(Diagnostics& diagnostic) const override;
    void solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J,
                                                     const double dt);

   private:
    std::optional<ChargeContinuityResidual> continuity_;
    void prepare_diagnostic_currents() override;
    // One species' Boris power and wall decisions belong to the same half-step.
    // Particle addresses in boundary_events remain valid until boundary application;
    // do not reorder, inject or redistribute particles between these two stages.
    struct SpeciesSecondHalfStep {
        double predictor_power = 0;   // sum(q*mpw*v_mid . E_predictor), in code energy/time units
        ParticleBoundaryEvents boundary_events;
    };
    using SecondHalfStepResults = std::map<std::string, SpeciesSecondHalfStep>;

    void drift_half_step_and_assemble_response();
    SecondHalfStepResults advance_velocities_and_deposit_second_half_current();
    void correct_particle_energy_and_apply_boundaries(SecondHalfStepResults& second_half_step);
};

#endif
