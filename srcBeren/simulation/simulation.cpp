// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation.h"

#include <iostream>
#include <map>
#include <string>

#include "Coil.h"
#include "Damping.h"
#include "Diagnostic.h"
#include "Mesh.h"
#include "ParticlesArray.h"
#include "Read.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "recovery.h"
#include "simulation_ecsim.h"
#include "simulation_ecsim_corr.h"

Simulation::Simulation(const ParametersMap& _systemParameters,
                       const nlohmann::json& s_config,
                       const nlohmann::json& p_config, int argc, char** argv)
    : parameters(_systemParameters),
      system_config(s_config),
      particles_config(p_config) {
    // for skip warning about unused arguments
    (void) argv[argc - 1];

    static char help[] = "Plasma simulation.\n\n";
    std::cout << help;
    parameters.print();
}

void Simulation::init(){
    bounds.setBounds(system_config);
    domain.setDomain(system_config);
    mesh.init(domain, parameters);
    init_operators();
    init_fields();
    init_particles(particles_config);
    globalTimer.init("globalFunctions.time");
    std::cout << "Simulation initialized\n";
}

void Simulation::init_operators() {
    Imat.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlE.resize(domain.total_size() * 3, domain.total_size() * 3);
    curlB.resize(domain.total_size() * 3, domain.total_size() * 3);
    divE.resize(domain.total_size(), domain.total_size() * 3);

    stencil_Imat(Imat, domain);
    stencil_curlE(curlE, domain);
    stencil_curlB(curlB, domain);
    stencil_divE(divE, domain);
}

void Simulation::calculate() {
    RandomGenerator gen;
    const int startTimeStep = parameters.get_int("StartTimeStep");
    const int lastTimestep = parameters.get_int("LastTimestep");

   // Writer writer(mesh, species, domain, parameters);
    std::cout << "Start simulation\n";
    make_diagnostic(0);
    for (auto timestep = startTimeStep + 1; timestep <= lastTimestep;
         ++timestep) {
        prepare_step(timestep);
        make_step(timestep);
        for(auto &sp : species) {
            std::cout << "Moved " << sp->get_total_num_of_particles() << " "
                      << sp->name() << "\n";
        }
        double collision_time = omp_get_wtime();
        if(parameters.get_string("Collider") != "None"){
            collision_step(timestep);
        }
        std::cout << "Collision time: " << omp_get_wtime() - collision_time << "\n";
        make_diagnostic(timestep);
        globalTimer.write(timestep, 1);
        }
}

void Simulation::init_particles(const nlohmann::json& j) {
    if (!j.contains("particles") || !j["particles"].is_array()){
        throw std::runtime_error("mixture requires particles[]");
    }
    
    for (const auto &config : j["particles"]) {
        species.push_back(make_particles_array(config));
    }

    charged_species.reserve(species.size());
    for (auto& sp_up : species) {
        if (!sp_up->is_neutral())
            charged_species.push_back(std::ref(*sp_up));
    }
    for (auto &sp : species) {
        if (parameters.get_int("k_particles_reservation") > 0.) {
            for (auto k = 0; k < sp->size(); ++k) {
                sp->particlesData(k).reserve(
                    parameters.get_int("k_particles_reservation") *
                    sp->NumPartPerCell);
            }
        }
        if (parameters.get_int("StartFromTime") > 0) {
            read_particles_from_recovery(
                sp);   // Only for start simulation from old files!!!
            std::cout << "Upload " + sp->name() + " success!\n";
            continue;
        } else {
            double init_energy = sp->distribute_initial_particles(sp->get_initial_distributions(), domain);
            std::cout << sp->name() << " distibuted with init energy: " << init_energy << "\n";
        }
        sp->density_on_grid_update();
        std::cout << sp->particlesData.size() << " "
                  << sp->particlesData.capacity() << "\n";
    }
}

void Simulation::collect_current(Field3d &J) {
    J.setZero();
    for (const auto &sp : species) {
        J.data() += sp->currentOnGrid.data();
    }
}

void Simulation::collect_charge_density(
    Field3d &field) {
    field.setZero();
    for (const auto &sp : species) {
        field.data() += sp->densityOnGrid.data();
    }
}

std::unique_ptr<Simulation> build_simulation(
    const ParametersMap &systemParameters,
    const nlohmann::json &system_config,
    const nlohmann::json &particles_config, int argc, char **argv) {
    auto scheme_name = get_checked<std::string>(system_config,"Scheme");

    std::unique_ptr<Simulation> simulation = nullptr;

    if (scheme_name == "ecsim") {
        simulation = std::make_unique<SimulationEcsim>(
            systemParameters, system_config, particles_config,
            argc, argv);
    } else if (scheme_name == "ecsim_corr") {
        simulation = std::make_unique<SimulationEcsimCorr>(
            systemParameters, system_config, particles_config,
            argc, argv);
    } else {
        std::cout << "Scheme " << scheme_name << " is not supported\n";
        exit(-1);
    }
    return simulation;
}

void Simulation::collision_step([[maybe_unused]] const int timestep) {
    const double dt = parameters.get_double("Dt");
    const double n0 = parameters.get_double("n0");

    if (parameters.get_string("Collider") == "ColliderWithNeutrals") {
        CollisionScheme scheme = CollisionScheme::PHYSICAL_ONLY;
        CollisionProcessOptions process_opts = CollisionProcessOptions();
        static BinaryColliderWithNeutrals collider(n0, scheme, process_opts);

        collider.collide_with_neutrals_binary(species, domain, dt);
    }
}
