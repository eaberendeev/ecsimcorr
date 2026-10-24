// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_H
#define SIMULATION_ECSIM_H

#include <memory>
#include <vector>

#include "DiagnosticOutput.h"
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
    void init_operators() override;
    void init_fields() override;
    void prepare_step(const int timestep) override;
    void make_step(const int timestep) override;
    virtual void diagnostic_energy(Diagnostics& diagnostic);
    void make_diagnostic(const int timestep) override;
    void finalize() override;
    template <ShapeType ST = SHAPE>
    void prepare_block_matrix();
    void first_push();
    void predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& B, Field3d& J);
    void predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& E_ex, const Field3d& B, Field3d& J);
    void second_push();

    void collect_per_species_diagnostics(Diagnostics& diagnostic, double& kineticEnergy, double& kineticEnergyNew,
                                         double& totalLostEnergy);
    void compute_field_energy_and_conservation(Diagnostics& diagnostic, const IndexRange& irange, double dt,
                                               double kineticEnergy, double kineticEnergyNew, double totalLostEnergy,
                                               double totalInjectEnergy, double energyJe_ex, double dampingEnergy,
                                               double er0Work = 0.0);

    void assembleLmat2(double dt);

    Field3d fieldJp;        // predict current for EM solver
    Field3d fieldJp_full;   // predict current for EM solver Jp + Lmat(E+E_n);
    Field3d fieldJe;        // Esirkepov current for E correction};
    Field3d fieldE;
    Field3d fieldEn;
    Field3d fieldEp;
    Field3d fieldB;
    Field3d fieldBn;
    Field3d fieldBInit;
    Field3d fieldBFull;
    Field3d fieldE_external;

    // энергия полей, поглощённая damping-слоями на текущем шаге
    double dampingEnergy_ = 0;

    Operator Mmat;
    Operator IMmat;
    std::vector<IndexMap> LmatX;

   protected:
    Field3d fieldE0;   // поле-шаблон Дирихле: на кольце er0 стоит предписанное значение, вне кольца нули
    // <E0, J_raw> до apply_to_fields(CURRENT): учётная работа отфильтрованного
    // граничного тока, а не физическая мощность порта/Poynting.
    double er0Work_ = 0.0;

   private:
    std::unique_ptr<Diagnostics> diagnostic_ptr_;
    std::vector<std::unique_ptr<IDiagnosticOutput>> outputs_;

    // used for optimized branch of function assembleLmat2
    SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>> rowBlocksGlobal;

    // Counter for periodic Lmat2 reference validation (per simulation instance).
    int assembly_check_counter_ = 0;
};

#endif
