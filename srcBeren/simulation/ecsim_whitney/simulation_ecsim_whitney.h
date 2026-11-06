#pragma once

#include "simulation/ecsim/simulation_ecsim.h"

// ECSIM time staggering and Maxwell solve with compatible spatial forms.
// No orbit iteration or Gauss correction is introduced.
class SimulationEcsimWhitney : public SimulationEcsim {
   public:
    using SimulationEcsim::SimulationEcsim;
    void init() override;
    ShapeType diagnostic_density_shape() const override {
        return ShapeType::Linear;
    }
    void drift_particles_and_assemble_response() override;
    void advance_particle_velocities() override;

   protected:
    // Assemble Jp and Lmat2 at the current particle positions. The caller
    // owns the drift: full step for ECSIM, first half-step for ECSIM-Corr.
    void assemble_point_response();
    // Non-periodic axes require a wall-constrained Maxwell response. Point
    // response drops crossing particles before assembly; the frozen ballistic
    // path instead clips their orbits and keeps their truncated current.
    bool wall_boundaries_ = false;
    virtual bool supports_wall_boundaries() const {
        return true;
    }
    bool parallel_point_ = false;   // sequential reference remains available
    // currentOnGrid already contains q w W (v^n+v^{n+1})/2 from the push.
    void prepare_diagnostic_currents() override {
    }
};
