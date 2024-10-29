// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_TEST_TRACK_IN_FIELD_B_H
#define SIMULATION_TEST_TRACK_IN_FIELD_B_H

#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "parameters_map.h"
#include "simulation.h"

class TrackDiagnostic{
    public:
    std::ofstream fEnergy;
    std::unique_ptr<ParticleTracker> tracker;
};

// Main simulation class
class ConstantFieldParticleTrajectorySimulator : public Simulation {
   public:
    ConstantFieldParticleTrajectorySimulator(
        const ParametersMap& _systemParameters,
        const std::vector<ParametersMap>& _speciesParameters,
        const ParametersMap& _outputParameters, int argc, char** argv)
        : Simulation(_systemParameters, _speciesParameters, _outputParameters,
                     argc, argv) {
        std::cout << "Implicit simulation\n";
    }
    void init_fields() override {}
    void init() override;
    void init_fields(const Field3d& E, const Field3d& B);
    void prepare_step(const int timestep) override;
    void make_step(const int timestep) override;
    void move_ecsim();
    void move_implicit();
    // void output_all(const int timestep) override;
    void init_particles() override{};
    void init_particles(
        std::vector<std::pair<std::string, Particle>>& particles);
    void diagnostic_energy(const int timestep);
    void make_diagnostic(const int timestep) override;
    void init_diagnostic(){
        std::string dirname = outputParameters.get_string("TestType") + "_" +
                              parameters.get_string("MoveType");
        std::string type = "Dt_" + parameters.get_string("Dt");
        std::string name = dirname + "//energy_" + type + ".txt";
        diagnostic.fEnergy.open(name);
        diagnostic.tracker =
            std::make_unique<ParticleTracker>(species, 1, dirname, type);
    }

    Field3d fieldJ;
    Field3d fieldE;
    Field3d fieldEn;
    Field3d fieldB;
    Field3d fieldBn;
    Field3d fieldE05;
    Field3d fieldB05;
    Field3d fieldBInit;
    Field3d fieldBFull;
    TrackDiagnostic diagnostic;
};

#endif
