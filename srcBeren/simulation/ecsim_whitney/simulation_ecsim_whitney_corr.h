#pragma once

#include <map>
#include <optional>
#include <string>

#include "diagnostics/Diagnostic.h"
#include "particles/particle_path.h"
#include "simulation_ecsim_whitney.h"

// Whitney ECSIM predictor followed by a charge-conserving correction on the
// actual two-leg particle path. The correction current is the trajectory
// integral of the same Whitney edge forms; no Esirkepov kernel is used.
// Absorbing walls use clipped paths; reflecting faces are rejected.
class SimulationEcsimWhitneyCorr : public SimulationEcsimWhitney {
   public:
    using SimulationEcsimWhitney::SimulationEcsimWhitney;
    void init() override;
    void make_step(const int timestep) override;
    void diagnostic_energy(Diagnostics& diagnostic) const override;

   private:
    std::optional<ChargeContinuityResidual> continuity_;
    // The second leg owns the boundary decision used after the velocity
    // rescaling. First-leg absorption is already applied before the predictor.
    struct SpeciesSecondHalfStep {
        double predictor_power = 0;
        ParticleBoundaryEvents boundary_events;
    };
    using SecondHalfStepResults = std::map<std::string, SpeciesSecondHalfStep>;

    void drift_half_step_and_assemble_response();
    SecondHalfStepResults advance_velocities_and_deposit_path_current();
    void correct_particle_energy_and_apply_boundaries(SecondHalfStepResults& second_half_step);
    void solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J,
                                                     double dt);
    void rescale_velocities_for_energy_balance(ParticlesArray& particles, double dt, double predictor_power);
};
