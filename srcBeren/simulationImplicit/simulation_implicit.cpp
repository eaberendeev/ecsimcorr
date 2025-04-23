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
#include "simulation_implicit.h"

void SimulationImplicit::init_fields() {

    fieldJ.resize(domain.size(), 3);
    fieldE.resize(domain.size(), 3);
    fieldEn.resize(domain.size(), 3);
    fieldB.resize(domain.size(), 3);
    fieldBn.resize(domain.size(), 3);
    fieldB05.resize(domain.size(), 3);
    fieldE05.resize(domain.size(), 3);
    fieldBInit.resize(domain.size(), 3);
    fieldBFull.resize(domain.size(), 3);

    fieldJ.setZero();

    if (parameters.get_int("StartFromTime") > 0) {
        read_fields_from_recovery(fieldEn, fieldB);
    } else {
        mesh.set_uniform_field(fieldE, 0.0, 0.0, 0.0);
        mesh.set_uniform_field(fieldBInit,
                               parameters.get_double("BUniform", 0),
                               parameters.get_double("BUniform", 1),
                               parameters.get_double("BUniform", 2));
        // mesh.fieldB = mesh.fieldBInit;
        fieldBn = fieldB;
        fieldEn = fieldE;
    }

    fieldE.setZero();
    fieldB.setZero();

    fieldEn = fieldE;
    fieldBn = fieldB;
}

void SimulationImplicit::init_particles() {
    Simulation::init_particles();
 //   int electrons = get_num_of_type_particles(species, "Electrons");
    //Particle test(7.53, 7.5, 7.5, -0.1, -0.1, -0.1);
    //species[electrons].add_particle(test);

 //   RandomGenerator gen;
    // for(int i = 0; i < 1; i++) {
        // Particle test(15+2 * gen.Uniform01(), 15+2 * gen.Uniform01(),
        //               15+2 * gen.Uniform01(), 0.1 * (1 - 2 * gen.Uniform01()),
        //               0.1 * (1 - 2 * gen.Uniform01()),
        //               0.0 * (1 - 2 * gen.Uniform01()));
       // std::cout << test << "\n";
     //       Particle test(2.25,2.25,2.25,1.,-0.0,0.0);
        // Particle test(0.199515, 1.99196, 4.63394, 0.24208, -0.0337293,
        //               -0.0598213);
        // Particle test(0.199515, 1.99196, 4.63394, 0.2, -0.0337293,
        //               -0.0);
    //   species[electrons].add_particle(test);
    //}
   // exit(0);
}

void SimulationImplicit::make_step([[maybe_unused]] const int timestep) {
    const double dt = parameters.get_double("Dt");

    globalTimer.start("Total");

    globalTimer.start("densityCalc");
    for (auto &sp : species) {
        sp->density_on_grid_update();   // calculate dendity field
    }
    globalTimer.finish("densityCalc");

    collect_charge_density(mesh.chargeDensityOld);

    double res = 1;
    int iter_count = 0;
    while (res > 1e-8) {
        iter_count++;
        fieldE05.data() = 0.5 * (fieldEn.data() + fieldE.data());
        fieldB05.data() = 0.5 * (fieldBn.data() + fieldB.data());
        fieldBFull.data() = fieldB05.data() + fieldBInit.data();
        
         globalTimer.start("particles1");
        for (auto &sp : species) {
            // first iteration E_n+1 = E_n, B_n+1  = B_n

            sp->push_Chen(fieldE05, fieldBFull, dt);
            mesh.apply_boundaries(sp->currentOnGrid, domain);
        }
        globalTimer.finish("particles1");
        Field3d J_old(fieldJ);

        collect_current(fieldJ);
        // for (int i = 0; i < fieldJ.size().x(); i++) {
        //     std::cout<< fieldJ(i,10,5,0) <<"\n";
        // }
        std::cout << "current norm " << fieldJ.data().norm() << "  " << (fieldJ.data()-J_old.data()).norm() << "\n";

        globalTimer.start("Fields");
        res = mesh.calculate_residual(fieldEn, fieldE, fieldB,
                                      fieldJ, dt);

        // res = A*E_n+1 - rhs
        // solve_fields(E, B, J); from E_n, B_n, J_n+1/2 => E_n+1, B_n+1
       mesh.impicit_find_fieldE(fieldEn, fieldE, fieldB,
                            fieldJ, dt);
       mesh.compute_fieldB(fieldBn, fieldB, fieldE,
                           fieldEn, dt);

        std::cout << fieldE.data().norm() << " "
                  << fieldEn.data().norm() << " " << res << " res \n";
        globalTimer.finish("Fields");
    }
    std::cout << "iter_count " << iter_count << "\n";

    // todo: bounds for particles
    //            call sp.update_cells(domain);

    for (auto &sp : species) {
        sp->density_on_grid_update();
    }

    // later output data and check conservation layws
    collect_charge_density(mesh.chargeDensity);
    std::cout << mesh.chargeDensity.data().norm()
              << " norm mesh.chargeDensity \n";

    auto divJ = mesh.divE * fieldJ.data();

    auto delta =
        (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (dt) +
        divJ;
    std::cout << delta.norm() << " norm drho / Dt - divJ \n";
}

void SimulationImplicit::prepare_step(const int timestep) {
    damping_fields(fieldEn, fieldB, domain,
                   parameters);
    mesh.prepare();

    fieldE = fieldEn;
    fieldB = fieldBn;
    fieldJ.setZero();

    inject_particles(timestep, domain);
    for (auto &sp : species) {
        sp->prepare();   // save start coord for esirkepov current
    }
}

void SimulationImplicit::make_diagnostic(const int timestep) {
    static Diagnostics diagnostic(outputParameters,
                                  domain, species);
    // output_all(timestep);
    for (auto &sp : species) {
        write_particles_to_recovery(sp, timestep,
                                    parameters.get_int("RecoveryInterval"));
    }
    write_fields_to_recovery(fieldEn, fieldBn, timestep,
                             parameters.get_int("RecoveryInterval"));
    // writer.output(0, parameters, timestep);
    diagnostic_energy(diagnostic);
    diagnostic.write_energy(parameters, timestep);
    const std::string pathToField = ".//Fields//Diag2D//";

    std::vector<std::pair<const Field3d &, std::string>> fields = {
        {fieldEn, pathToField + "FieldE"}, {fieldBn, pathToField + "FieldB"}};
    diagnostic.output_fields2D(timestep,
                               fields);
    for (auto &sp : species) {

        const std::string pathToField =
            ".//Particles//" + sp->name() + "//Diag2D//";

        std::vector<std::pair<const Field3d &, std::string>> fields = {
            {sp->currentOnGrid, pathToField + "Current"},
            {sp->densityOnGrid, pathToField + "Density"}};
        diagnostic.output_fields2D(timestep, fields);
    }
#ifdef SET_PARTICLE_IDS
    static ParticleTracker tracker(species, 1, "Tracking", "");
    tracker.track_particles(species, timestep);
#endif
}

void SimulationImplicit::diagnostic_energy(
    Diagnostics &diagnostic) {
    double kineticEnergy = 0;
    double kineticEnergyNew = 0;
    for (auto &sp : species) {
        diagnostic.addEnergy(sp->name() + "Init",
                             sp->get_init_kinetic_energy());
        diagnostic.addEnergy(sp->name(), sp->get_kinetic_energy());
        diagnostic.addEnergy(sp->name() + "Inject", sp->injectionEnergy);
        diagnostic.addEnergy(sp->name() + "LostZ", sp->lostEnergyZ);
        diagnostic.addEnergy(sp->name() + "LostXY", sp->lostEnergyXY);
        diagnostic.addEnergy(sp->name() + "Z", sp->get_kinetic_energy(Z));
        diagnostic.addEnergy(sp->name() + "XY", sp->get_kinetic_energy(X, Y));
        kineticEnergy += diagnostic.energy[sp->name() + "Init"];
        kineticEnergyNew += diagnostic.energy[sp->name()];
    }

    diagnostic.addEnergy("energyFieldE", mesh.calc_energy_field(fieldEn));
    diagnostic.addEnergy("energyFieldB", mesh.calc_energy_field(fieldBn));
    double energyFieldEold = mesh.calc_energy_field(fieldE);
    double energyFieldBold = mesh.calc_energy_field(fieldB);


    double energyFieldDifference = diagnostic.energy["energyFieldB"] +
                                   diagnostic.energy["energyFieldE"] -
                                   energyFieldBold - energyFieldEold;

    std::cout << "Energy " << kineticEnergyNew - kineticEnergy << " "
              << energyFieldDifference << " "
              << 4 * 0.5 *
                     (calc_JE(fieldE, fieldJ, domain.get_bounds()) +
                      calc_JE(fieldEn, fieldJ, domain.get_bounds()))
              << "\n";

    diagnostic.addEnergy(
        "energyConserve",
        std::abs(kineticEnergyNew - kineticEnergy + energyFieldDifference));
}
