#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "drift_kinetic_grid.h"
#include "drift_kinetic_particles.h"
#include "drift_kinetic_push.h"
#include "simulation/simulation.h"

// Phase 2c-1: a selectable drift-kinetic scheme with prescribed uniform fields.
// Initial kinetic particles are converted to guiding centres at init; the step
// advances them with the implicit guiding-centre pusher. Self-consistent
// Maxwell coupling and the hybrid per-species model are later phases.
class SimulationDriftKinetic : public Simulation {
   public:
    using Simulation::Simulation;
    void init() override;
    void init_fields() override;
    void prepare_step(int timestep) override {
        (void) timestep;
    }
    void make_step(int timestep) override;
    void make_diagnostic(int timestep) override;
    void finalize() override;

   private:
    void convert_initial_particles();

    std::unique_ptr<driftKinetic::DriftKineticGrid> grid_;
    std::unordered_map<std::string, std::unique_ptr<driftKinetic::DriftKineticParticles>> particles_;
    driftKinetic::GuidingCenterPush push_;
};
