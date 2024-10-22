// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include <fstream>
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
#include "test_track_in_fieldB.h"

void SimulationTestTrackInFiledB::init() {
    bounds.setBounds(parameters);
    domain.setDomain(parameters, bounds);
    mesh.init(domain, parameters);
    globalTimer.init("globalFunctions.time");
    std::string dirname = outputParameters.get_string("TestType");
    create_directory(".//"+dirname);
}

void SimulationTestTrackInFiledB::init_fields(const Field3d &E,
                                              const Field3d &B) {
    fieldE.resize(domain.size(), 3);
    fieldEn.resize(domain.size(), 3);
    fieldB.resize(domain.size(), 3);
    fieldBn.resize(domain.size(), 3);
    fieldB05.resize(domain.size(), 3);
    fieldE05.resize(domain.size(), 3);
    fieldBInit.resize(domain.size(), 3);
    fieldBFull.resize(domain.size(), 3);
    fieldB.set_zero();

    fieldBInit = B;
    fieldE = E;

    //mesh.set_uniform_field(fieldE, 0, 0.001, 0);
    fieldEn = fieldE;
    fieldBn = fieldB;
}

void SimulationTestTrackInFiledB::init_particles(std::vector<std::pair<std::string,Particle>>& pairs) {
    Simulation::init_particles();

    for (auto &pair : pairs) {
        int name = get_num_of_type_particles(species, pair.first);
        species[name].add_particle(pair.second);
    }
}

void move_ecsim(){}
void move_implicit(){}

void SimulationTestTrackInFiledB::make_step(const int timestep) {
    const double dt = parameters.get_double("Dt");

    globalTimer.start("Total");

    globalTimer.start("densityCalc");
    for (auto &sp : species) {
        sp.density_on_grid_update();   // calculate dendity field
    }
    globalTimer.finish("densityCalc");

    //collect_charge_density(mesh.chargeDensityOld);

        fieldE05.data() = 0.5 * (fieldEn.data() + fieldE.data());
        fieldB05.data() = 0.5 * (fieldBn.data() + fieldB.data());
        fieldBFull.data() = fieldB05.data() + fieldBInit.data();

        globalTimer.start("particles1");
        for (auto &sp : species) {
            // first iteration E_n+1 = E_n, B_n+1  = B_n
            sp.push_Chen(mesh, fieldE05, fieldBFull, dt);
        }

        globalTimer.finish("particles1");
}

void SimulationTestTrackInFiledB::prepare_step(const int timestep) {
    fieldE = fieldEn;
    fieldB = fieldBn;

    for (auto &sp : species) {
        sp.prepare();   // save start coord for esirkepov current
    }
}

void SimulationTestTrackInFiledB::make_diagnostic(
                                                  const int timestep) {
    diagnostic_energy(timestep);
#ifdef SET_PARTICLE_IDS
    diagnostic.tracker->track_particles(species, timestep);
#endif
}

void SimulationTestTrackInFiledB::diagnostic_energy(
                                                    const int timestep) {


    for (auto &sp : species) {
        if(timestep == 0){
            diagnostic.fEnergy << sp.name();
        }
        diagnostic.fEnergy << sp.get_kinetic_energy();
    }
    diagnostic.fEnergy << "\n";
}
