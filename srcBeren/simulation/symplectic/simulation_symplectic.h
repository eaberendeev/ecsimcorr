#pragma once

#include "DiagnosticOutput.h"
#include "simulation.h"
#include "whitney_pic.h"

class SimulationSymplectic : public Simulation {
   public:
    using Simulation::Simulation;
    void init() override;
    void init_fields() override;
    void deposit_density(ParticlesArray& sp) override;
    void prepare_step(int timestep) override;
    void make_step(int timestep) override;
    void make_diagnostic(int timestep) override;
    void finalize() override;

    // E, B, x and v all belong to the same integer time n, also in recovery.
    // B includes the uniform magnetic flux; no separate imposed field is used.
    Field3d electric, magnetic;

   private:
    void validate_state() const;
    void kinetic_flow(int axis, double dt);
    void electric_flow(double dt);
    std::unique_ptr<Diagnostics> diagnostic_;
    std::vector<std::unique_ptr<IDiagnosticOutput>> outputs_;
    Field3d zero_field_;
    Field3d gauss_initial_;
    double energy_initial_ = 0;
};
