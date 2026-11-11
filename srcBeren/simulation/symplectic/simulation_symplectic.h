#pragma once

#include "diagnostics/Diagnostic.h"
#include "diagnostics/IDiagnosticOutput.h"
#include "simulation/simulation.h"

class SimulationSymplectic : public Simulation {
   public:
    using Simulation::Simulation;
    void init() override;
    void init_fields() override;
    void prepare_step(int timestep) override;
    void make_step(int timestep) override;
    void make_diagnostic(int timestep) override;
    void finalize() override;

    // sum(F^2)/2 over one physical copy of the periodic mesh.
    double field_energy(const Field3d& field) const;

    // div E - rho with the Whitney nodal charge of every charged species.
    Field3d gauss_defect() const;

    // E, B, x and v all belong to the same integer time n, also in recovery.
    // B includes the uniform magnetic flux; no separate imposed field is used.
    // Ghost nodes always hold the periodic images of the physical nodes.
    Field3d electric, magnetic;

   private:
    void kinetic_flow(bool reverse, double dt);
    void electric_flow(double dt);
    void magnetic_flow(double dt);

    // A curl or the path-integrated current of one species inside one flow; never checkpointed.
    Field3d scratch_;
    Field3d zero_field_;
    Field3d gauss_initial_;
    double energy_initial_ = 0;
    std::unique_ptr<Diagnostics> diagnostic_;
    std::vector<std::unique_ptr<IDiagnosticOutput>> outputs_;
};
