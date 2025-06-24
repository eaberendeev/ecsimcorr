// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_SYMPLECTIC_H
#define SIMULATION_SYMPLECTIC_H

#include "simulation.h"
#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"
#include "parameters_map.h"

// Main simulation class
class SimulationSymplectic: public Simulation{
    public:
     SimulationSymplectic(const ParametersMap& _systemParameters,
                          const std::vector<ParametersMap>& _speciesParameters,
                          const ParametersMap& _outputParameters, int argc,
                          char** argv)
         : Simulation(_systemParameters, _speciesParameters, _outputParameters,
                      argc, argv) {}
     void init_fields() override;
     void prepare_step(const int timestep) override;
     void collision_step(const int timestep){};
     void make_step(const int timestep) override;
     void diagnostic_energy(Diagnostics& diagnostic);
     void make_diagnostic(const int timestep) override;

     void move_particles(const Field3d& fieldB, Field3d& fieldJ) {
         double dt = parameters.get_double("Dt");
         for (auto& sp : species) {
             move_x(sp, fieldB, fieldJ, dt);
             move_y(sp, fieldB, fieldJ, dt);
             move_z(sp, fieldB, fieldJ, dt);
         }
     }
     void move_x(std::unique_ptr<ParticlesArray>& particlesArray,
                 const Field3d& fieldB, Field3d& fieldJ, const double dt);
     void move_y(std::unique_ptr<ParticlesArray>& particlesArray,
                 const Field3d& fieldB, Field3d& fieldJ, const double dt);
     void move_z(std::unique_ptr<ParticlesArray>& particlesArray,
                 const Field3d& fieldB, Field3d& fieldJ, const double dt);
     void move_x_orig(std::unique_ptr<ParticlesArray>& particlesArray,
                 const Field3d& fieldB, Field3d& fieldJ, const double dt);
     void move_y_orig(std::unique_ptr<ParticlesArray>& particlesArray,
                      const Field3d& fieldB, Field3d& fieldJ, const double dt);
     void move_z_orig(std::unique_ptr<ParticlesArray>& particlesArray,
                      const Field3d& fieldB, Field3d& fieldJ, const double dt);
     void update_velocity(std::unique_ptr<ParticlesArray>& particlesArray,
                          const double dt);

     void update_E() {
         const double dt = parameters.get_double("Dt");

         fieldE.data() += dt * mesh.curlB * fieldB.data() - dt * fieldJ.data();

     }
     void update_B() {
         const double dt = parameters.get_double("Dt");
         fieldB.data() -= dt * mesh.curlE * fieldE.data();
     }
     Field3d fieldJ;
     Field3d fieldE;
     Field3d fieldB;
     Field3d fieldBInit;
     Field3d fieldBFull;
};

struct InterpolationWeights {
    alignas(32) int indices[2];
    alignas(32) double weights[2];
};

static inline InterpolationWeights compute_weights(double coord, double shift) {
    InterpolationWeights result;
    const double shifted_coord = coord - shift;
    const int index = static_cast<int>(shifted_coord);
    const double delta = shifted_coord - index;

    result.indices[0] = index;
    result.indices[1] = index + 1;
    result.weights[0] = 1.0 - delta;
    result.weights[1] = delta;

    return result;
}

#endif
