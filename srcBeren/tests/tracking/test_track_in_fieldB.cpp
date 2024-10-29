// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "test_track_in_fieldB.h"

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

void ConstantFieldParticleTrajectorySimulator::init() {
    bounds.setBounds(parameters);
    domain.setDomain(parameters, bounds);
    mesh.init(domain, parameters);
    globalTimer.init("globalFunctions.time");
    std::string dirname = outputParameters.get_string("TestType");
    create_directory(".//" + dirname);
}

void ConstantFieldParticleTrajectorySimulator::init_fields(const Field3d &E,
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

    fieldEn = fieldE;
    fieldBn = fieldB;
}

void ConstantFieldParticleTrajectorySimulator::init_particles(
    std::vector<std::pair<std::string, Particle>> &pairs) {
    Simulation::init_particles();

    for (auto &pair : pairs) {
        int name = get_num_of_type_particles(species, pair.first);
        species[name]->add_particle(pair.second);
    }
}

void ConstantFieldParticleTrajectorySimulator::move_ecsim() {
    const double dt = parameters.get_double("Dt");

    fieldBFull.data() = fieldB.data() + fieldBInit.data();

    for (auto &sp: species) {
        sp->move_and_calc_current(0.5 * dt);
        sp->update_cells(domain);
        // +++ get v'_{n+1} from v_{n} and E'_{n+1}
        sp->predict_velocity(fieldE, fieldEn, fieldBFull, domain, dt);
        sp->move_and_calc_current(0.5 * dt);
        sp->update_cells(domain);
    }

}

void ConstantFieldParticleTrajectorySimulator::move_implicit() {
    const double dt = parameters.get_double("Dt");

    fieldE05.data() = 0.5 * (fieldEn.data() + fieldE.data());
    fieldB05.data() = 0.5 * (fieldBn.data() + fieldB.data());
    fieldBFull.data() = fieldB05.data() + fieldBInit.data();

    globalTimer.start("particles");
    for (auto &sp : species) {
        sp->push_Chen(mesh, fieldE05, fieldBFull, dt);
        sp->update_cells(domain);
    }

    globalTimer.finish("particles");
}

void ConstantFieldParticleTrajectorySimulator::make_step(const int timestep) {
    const auto mover = parameters.get_string("MoveType");

    if(mover == "implicit"){
        move_implicit();
    }
    else if(mover == "ecsim"){
        move_ecsim();
    }
}

void ConstantFieldParticleTrajectorySimulator::prepare_step(const int timestep) {
    fieldE = fieldEn;
    fieldB = fieldBn;

    for (auto &sp : species) {
        sp->prepare();   // save start coord for esirkepov current
    }
}

void ConstantFieldParticleTrajectorySimulator::make_diagnostic(const int timestep) {
    diagnostic_energy(timestep);
#ifdef SET_PARTICLE_IDS
    diagnostic.tracker->track_particles(species, timestep);
#endif
}

void ConstantFieldParticleTrajectorySimulator::diagnostic_energy(const int timestep) {
    for (auto &sp : species) {
        if (timestep == 0) {
            diagnostic.fEnergy << sp->name();
        }
        diagnostic.fEnergy << sp->get_kinetic_energy();
    }
    diagnostic.fEnergy << "\n";
}
