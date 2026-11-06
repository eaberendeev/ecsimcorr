#pragma once

#include <map>
#include <optional>
#include <string>

#include "diagnostics/Diagnostic.h"
#include "simulation_ecsim_whitney_frozen.h"

// Ballistic frozen Whitney ECSIM predictor followed by the same trajectory
// integral correction as the point scheme. The correction current is built
// from the actual two-leg path; the predictor uses the orbit-averaged Jp/L.
// The prototype supports periodic boxes only.
class SimulationEcsimWhitneyFrozenCorr : public SimulationEcsimWhitneyFrozen {
   public:
    using SimulationEcsimWhitneyFrozen::SimulationEcsimWhitneyFrozen;
    void init() override;
    void make_step(const int timestep) override;
    void diagnostic_energy(Diagnostics& diagnostic) const override;

   protected:
    bool supports_wall_boundaries() const override {
        return false;
    }

   private:
    std::optional<ChargeContinuityResidual> continuity_;

    void drift_half_step_and_assemble_response();
    std::map<std::string, double> advance_velocities_and_deposit_path_current();
    void correct_particle_energy_and_apply_boundaries(const std::map<std::string, double>& predictor_power);
    void solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J,
                                                     double dt);
    void rescale_velocities_for_energy_balance(ParticlesArray& particles, double dt, double predictor_power);
};
