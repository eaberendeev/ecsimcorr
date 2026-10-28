// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_H
#define SIMULATION_ECSIM_H

#include <memory>
#include <vector>

#include "IDiagnosticOutput.h"
#include "ParticlesArray.h"
#include "World.h"
#include "algorithms_ecsim.h"
#include "containers.h"
#include "memory.h"
#include "simulation.h"

// Main simulation class
class SimulationEcsim : public Simulation {
   public:
    SimulationEcsim(const nlohmann::json& system_config, const nlohmann::json& particles_config, int argc, char** argv)
        : Simulation(system_config, particles_config, argc, argv) {
        const int nthr = omp_get_max_threads();
        rowBlocksGlobal.resizeAndReset(nthr);
        for (int i = 0; i < nthr; ++i) {
            constexpr int reserveSize = 32 * 1024 * 1024;   // 32 MB
            rowBlocksGlobal[i].reserve(reserveSize / sizeof(rowBlocksGlobal[0][0]));
        }
    }
    ~SimulationEcsim();
    void init_fields() override;
    void prepare_step(const int timestep) override;
    void make_step(const int timestep) override;
    virtual void diagnostic_energy(Diagnostics& diagnostic);
    void make_diagnostic(const int timestep) override;
    void finalize() override;
    template <ShapeType ST = SHAPE>
    void prepare_block_matrix();
    void drift_particles_and_assemble_response();
    void predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& B, Field3d& J);
    void predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& E_ex, const Field3d& B, Field3d& J);
    void advance_particle_velocities();

    void collect_per_species_diagnostics(Diagnostics& diagnostic, double& kineticEnergy, double& kineticEnergyNew,
                                         double& totalLostEnergy);
    void compute_field_energy_and_conservation(Diagnostics& diagnostic, const IndexRange& irange, double dt,
                                               double kineticEnergy, double kineticEnergyNew, double totalLostEnergy,
                                               double totalInjectEnergy, double external_field_power,
                                               double dampingEnergy, double filtered_current_power = 0.0);

    void assembleLmat2(double dt);

    // Jp/Je are reset in prepare_step and accumulated during this step.
    Field3d fieldJp;        // Jp(x^{n+1/2}, v^n): field-independent ECSIM response
    Field3d fieldJp_full;   // Jp + L*(E^n+E^{n+1})/dt, reconstructed for diagnostics
    Field3d fieldJe;        // Corr: Esirkepov current averaged over both half-paths
    // prepare_step copies the previous endpoint into E/B after damping;
    // diagnostics read both endpoints before the next prepare_step.
    Field3d fieldE;            // E^n
    Field3d fieldEn;           // E^{n+1}
    Field3d fieldEp;           // predicted E^{n+1/2}
    Field3d fieldB;            // B^n
    Field3d fieldBn;           // B^{n+1}
    Field3d fieldBInit;        // prescribed, time-independent magnetic field
    Field3d fieldBFull;        // scratch B + BInit; diagnostics also use Bn + BInit
    Field3d fieldE_external;   // prescribed, time-independent electric field

    // энергия полей, поглощённая damping-слоями на текущем шаге
    double dampingEnergy_ = 0;

   protected:
    Field3d fieldE0;   // поле-шаблон Дирихле: на кольце er0 стоит предписанное значение, вне кольца нули
    // <E0,Jp_raw> before the current mask, in code energy/time units.
    // Overwritten during response assembly each step; dt times this value
    // enters the discrete energy residual. Initially zero, recomputed on restart.
    double filtered_current_power_ = 0.0;

   private:
    std::unique_ptr<Diagnostics> diagnostic_ptr_;
    std::vector<std::unique_ptr<IDiagnosticOutput>> outputs_;

    // used for optimized branch of function assembleLmat2
    SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>> rowBlocksGlobal;

    // Counter for periodic Lmat2 reference validation (per simulation instance).
    int assembly_check_counter_ = 0;
};

#endif
