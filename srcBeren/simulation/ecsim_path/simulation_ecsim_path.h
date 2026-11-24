#pragma once

#include <map>
#include <string>

#include "magnetization.h"
#include "simulation/ecsim_corr/simulation_ecsim_corr.h"

// ECSIM-Path ("Scheme": "ecsim_path"): ECSIM-Corr without the energy
// correction. Each particle is pushed with the mean predicted field on its
// straight step path x^n -> x^{n+1} and deposits that path (path_gather.h), so
// its energy changes by the work of the field on its own current; the
// predictor current is built on the predicted path, and the missing magnetic
// moment of the straight path is added as a magnetization current
// (magnetization.h). One corrector pass is always made ("CorrectorPasses" >= 1):
// without it the step energy error doubles and the fields grow noisy at
// Omega_e dt ~ 10.
// Stability needs a thermal path per step of at most about half a cell,
// v_th dt <~ h/2 (the predictor freezes the path weights); a warning is printed
// at the first step otherwise.
class SimulationEcsimPath : public SimulationEcsimCorr {
   public:
    using SimulationEcsimCorr::SimulationEcsimCorr;
    void init_fields() override;
    void prepare_step(int timestep) override;
    void diagnostic_values(Diagnostics& diagnostic) const override;

   protected:
    void predict_species_current(ParticlesArray& sp) override;
    void after_first_half_step() override;
    void add_extra_current(Field3d& current) const override;
    void push_species_and_deposit(ParticlesArray& sp, const Field3d& electric, const Field3d& magnetic,
                                  SpeciesSecondHalfStep& result) override;
    void after_magnetic_field_advance() override;
    void correct_species_energy(ParticlesArray& sp, double predictor_power) override;

   private:
    // Fixed-point passes of the path gather from the point field.
    static constexpr int gather_passes_ = 2;
    bool path_checked_ = false;
    Field3d magnetization_moment_, magnetization_current_;   // dM (B layout), curl_B dM (E layout)
    std::map<std::string, pathMagnetization::Moments> moments_;
    double magnetization_work_ = 0, magnetization_remainder_ = 0;   // last step
    void warn_about_long_paths() const;
};
