#pragma once

#include "diagnostics/Diagnostic.h"
#include "diagnostics/IDiagnosticOutput.h"
#include "simulation/simulation.h"
#include "whitney_pic.h"

class SimulationSymplectic : public Simulation {
   public:
    using Simulation::Simulation;
    void init() override;
    void init_fields() override;
    void prepare_step(int timestep) override;
    void make_step(int timestep) override;
    void make_diagnostic(int timestep) override;
    void finalize() override;

    // E, B, x and v all belong to the same integer time n, also in recovery.
    // B includes the uniform magnetic flux; no separate imposed field is used.
    Field3d electric, magnetic;

    // Whitney nodal charge deposits. The first adds one species into rho,
    // the second sums every charged species (used by the Gauss checks).
    void deposit_charge_density(const ParticlesArray& sp, Field3d& rho) const;
    void collect_charge_density(Field3d& field) const;

   private:
    void validate_initial_state() const;
    int reset_deposit_buffers(std::vector<Field3d>& buffers, int components, const ParticlesArray& sp) const;
    // Scratch deposits are cleared before EACH physical subflow, reduced in
    // thread-index order, and never checkpointed. Particle/field state owns
    // all information needed to restart at an integer time.
    std::vector<Field3d> current_buffers_;
    mutable std::vector<Field3d> density_buffers_;
    void kinetic_flow(bool reverse, double dt);
    void electric_flow(double dt);
    std::unique_ptr<Diagnostics> diagnostic_;
    std::vector<std::unique_ptr<IDiagnosticOutput>> outputs_;
    Field3d zero_field_;
    Field3d gauss_initial_;
    double energy_initial_ = 0;
};
