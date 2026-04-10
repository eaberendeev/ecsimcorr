// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_H
#define SIMULATION_ECSIM_H

#include "ParticlesArray.h"
#include "World.h"
#include "algorithms_ecsim.h"
#include "containers.h"
#include "simulation.h"

// Main simulation class
class SimulationEcsim : public Simulation {
   public:
    SimulationEcsim(const nlohmann::json& system_config,
                    const nlohmann::json& particles_config, int argc,
                    char** argv)
        : Simulation(system_config, particles_config, argc, argv) {}
    void init_operators() override;
    void init_fields() override;
    void prepare_step(const int timestep) override;
    void make_step(const int timestep) override;
    void make_stepNGP(const int timestep);
    virtual void diagnostic_energy(Diagnostics& diagnostic);
    void make_diagnostic(const int timestep) override;
    void prepare_block_matrix(ShapeType type);
    void convert_block_matrix(ShapeType type);
    void first_push();
    void predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& B,
                                Field3d& J);
    void predict_electric_field(Field3d& Ep, const Field3d& E,
                                const Field3d& E_ex, const Field3d& B,
                                Field3d& J);
    void calculate_current();
    void second_push();

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

    Operator Mmat;
    Operator IMmat;
};

#endif
