// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#include "simulation_ecsim_base.h"

// Ordinary ECSIM: x^{n-1/2} -> x^{n+1/2}, then v^n -> v^{n+1}.
class SimulationEcsim : public SimulationEcsimBase {
   public:
    using SimulationEcsimBase::SimulationEcsimBase;
    void make_step(int timestep) override;
    void diagnostic_energy(Diagnostics& diagnostic) const override;
    void drift_particles_and_assemble_response();
    void advance_particle_velocities();

   protected:
    void prepare_diagnostic_currents() override;
};
