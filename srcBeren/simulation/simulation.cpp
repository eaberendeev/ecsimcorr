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

Simulation::Simulation(const ParametersMap &_systemParameters,
                       const std::vector<ParametersMap> &_speciesParameters,
                       const ParametersMap &_outputParameters, int argc,
                       char **argv)
    : parameters(_systemParameters), speciesParameters(_speciesParameters),
    outputParameters(_outputParameters) {
    static char help[] = "Plasma simulation.\n\n";
    std::cout << help;
    parameters.print();
}

void Simulation::init(){
    bounds.setBounds(parameters);
    domain.setDomain(parameters, bounds);
    mesh.init(domain, parameters);
    init_fields();
    init_particles();
    globalTimer.init("globalFunctions.time");
}


void Simulation::make_all() {
    RandomGenerator gen;
    const int startTimeStep = parameters.get_int("StartTimeStep");
    const int lastTimestep = parameters.get_int("LastTimestep");

   // Writer writer(mesh, species, domain, parameters);

    make_diagnostic(0);
    for (auto timestep = startTimeStep + 1; timestep <= lastTimestep;
         ++timestep) {
        prepare_step(timestep);
        make_step(timestep);
        for(auto &sp : species) {
            std::cout << "Moved " << sp->get_total_num_of_particles() << " "
                      << sp->name() << "\n";
        }

        make_diagnostic(timestep);

        globalTimer.write(timestep, 1);
        }
}

void Simulation::inject_particles(const int timestep) {
    for (auto &sp : species) {
        if (sp->distType != "INJECTION" ) continue;
        sp->distribute_particles(parameters, timestep);
    }
}

void Simulation::init_particles() {

    for (const auto &particlesParameters : speciesParameters) {
        species.push_back(std::make_unique<ParticlesArray>(particlesParameters,
                                                           parameters, domain));
    }

    for (auto &sp : species) {
        if (parameters.get_int("StartFromTime") > 0) {
            read_particles_from_recovery(
                sp);   // Only for start simulation from old files!!!
            std::cout << "Upload " + sp->name() + " success!\n";
            continue;
        } else{
            if (sp->distType == "INITIAL"){
                sp->distribute_particles(parameters, 0);
            }
        }
        if (parameters.get_int("k_particles_reservation") > 0.) {
            for (auto k = 0; k < sp->size(); ++k) {
                sp->particlesData(k).reserve(
                    parameters.get_int("k_particles_reservation") *
                    parameters.get_int("NumPartPerCell"));
            }
        }
        sp->density_on_grid_update();
        std::cout << sp->particlesData.size() << " "
                  << sp->particlesData.capacity() << "\n";
    }
}

void Simulation::collect_current(Field3d &J) {
    J.set_zero();
    for (const auto &sp : species) {
        J.data() += sp->currentOnGrid.data();
    }
}

void Simulation::collect_charge_density(
    Field3d &field) {
    field.set_zero();
    for (const auto &sp : species) {
        field.data() += sp->densityOnGrid.data();
    }
}
