// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_ecsim_corr.h"

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


// Particles have ccordinates and velocities. Mesh have 3D fields in nodes (each field stored in 1D array with 4d index x,y,z,d)
void SimulationEcsimCorr::make_step([[maybe_unused]] const int timestep) {
    const double dt = parameters.get_double("Dt");

    globalTimer.start("Total");

    globalTimer.start("densityCalc");
    for (auto &sp : species) {
        sp->density_on_grid_update(SHAPE_CH);   // calculate dendity field
        }
        globalTimer.finish("densityCalc");

        collect_charge_density(mesh.chargeDensityOld);

        globalTimer.start("particles1");
        fieldBFull.data() = fieldB.data() + fieldBInit.data();

        for (auto &sp : species) {
            sp->move_and_calc_current(0.5 * dt, sp->currentOnGrid,
                                      SHAPE_CH);   //  +++ x_n -> x_{n+1/2}

            sp->update_cells(domain);
            // +++ get J(x_{n+1/2},v_n)_predict

            sp->predict_current(fieldBFull, fieldJp, domain, dt, SHAPE);

            // Stencil LmatX in special format. Very slow function
            // +++ get Lgg'(x_{n+1/2})
            sp->fill_matrixL(mesh, fieldBFull, domain, dt, SHAPE);
        }
        // todo: zeros Lmat + current
        globalTimer.finish("particles1");

        mesh.apply_boundaries(fieldJp);
        mesh.apply_boundaries(mesh.LmatX);
        globalTimer.start("stencilLmat");
        // convert LmatX to Eigen sparce matrix format Lmat
        mesh.stencil_Lmat(domain);
        globalTimer.finish("stencilLmat");


        globalTimer.start("FieldsPredict");
        // --- solve A*E'_{n+1}=f(E_n, B_n, J(x_{n+1/2})). mesh consist En,
        // En+1_predict
        mesh.predictE(fieldEp, fieldE, fieldB, fieldJp,dt); // solve system of linear equations using Lmat for find fieldE

        globalTimer.finish("FieldsPredict");

        globalTimer.start("particles2");

        fieldBFull.data() = fieldB.data() + fieldBInit.data();
        for (auto &sp : species) {
            // +++ get v'_{n+1} from v_{n} and E'_{n+1}
            sp->predict_velocity(fieldE, fieldEp, fieldBFull, domain, dt,
                                 SHAPE);
            // calc new particles velocity using new fieldE
            // +++ x_{n+1/2} -> x_{n+1}
            //sp.move(0.5 * Dt);

            sp->move_and_calc_current(0.5 * dt, sp->currentOnGrid, SHAPE_CH);
            sp->update_cells(domain);

            sp->currentOnGrid.data() *= 0.5;
            mesh.apply_boundaries(sp->currentOnGrid);

        }
        globalTimer.finish("particles2");

        collect_current(fieldJe);
        std::cout << "Current " << fieldJe.data().norm() << "\n";

        globalTimer.start("FieldsCorr");
        // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final
        mesh.correctE(fieldEn, fieldE, fieldB, fieldJe, dt); // solve simple systeomof linear equations for correct fieldE
        globalTimer.finish("FieldsCorr");

        fieldJp_full.data() =
            fieldJp.data() +
            mesh.Lmat2 * (fieldE.data() + fieldEp.data()) / dt;

        globalTimer.start("particles3");
        for (auto &sp : species) {
            // correct particles velocity for save energy
            sp->correctv(fieldE, fieldEp, fieldEn, fieldJp_full, domain, dt);
        }
        for (auto &sp : species) {
            sp->density_on_grid_update(SHAPE_CH);
        }
        globalTimer.finish("particles3");

        globalTimer.start("computeB");
        // calculate fieldB
        mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn,
                            parameters.get_double("Dt"));
        globalTimer.finish("computeB");


        // later output data and check conservation layws
        collect_charge_density(mesh.chargeDensity);
        std::cout << mesh.chargeDensity.data().norm()
                  << " norm mesh.chargeDensity \n";

        auto divJ = mesh.divE * fieldJe.data();

        auto delta =
            (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (dt) +
            divJ;
        std::cout << delta.norm() << " norm drho / Dt - divJ \n";
}

void SimulationEcsimCorr::init_fields(){
    fieldJp.resize(domain.size(), 3);
    fieldJp_full.resize(domain.size(), 3);
    fieldJe.resize(domain.size(), 3);

    fieldE.resize(domain.size(), 3);
    fieldEn.resize(domain.size(), 3);
    fieldEp.resize(domain.size(), 3);
    fieldB.resize(domain.size(), 3);
    fieldBn.resize(domain.size(), 3);
    fieldBInit.resize(domain.size(), 3);
    fieldBFull.resize(domain.size(), 3);

    fieldJp.set_zero();
    fieldJe.set_zero();
    
    fieldEn.set_zero();
    fieldBInit.set_zero();
    set_coils(fieldBInit, domain, parameters);

    if (parameters.get_int("StartFromTime") > 0) {
        read_fields_from_recovery(fieldE, fieldB);
    } else {
        mesh.set_uniform_field(fieldE, 0.0, 0.0, 0.0);
        mesh.set_uniform_field(fieldBInit,
                               parameters.get_double("BUniform", 0),
                               parameters.get_double("BUniform", 1),
                               parameters.get_double("BUniform", 2));
    }
    fieldEn = fieldE;
    fieldBn = fieldB;
}

void SimulationEcsimCorr::prepare_step(const int timestep) {
    inject_particles(timestep);
    damping_fields(fieldEn, fieldBn, domain,
                   parameters);
    fieldE = fieldEn;
    fieldB = fieldBn;
    fieldJp.set_zero();
    fieldJe.set_zero();

#pragma omp parallel for
    for (size_t i = 0; i < mesh.LmatX.size(); i++) {
        for (auto it = mesh.LmatX[i].begin(); it != mesh.LmatX[i].end(); ++it) {
            it->second = 0.;
        }
    }
    for (auto &sp : species) {
        sp->prepare();   // save start coord for esirkepov current
    }
}

void SimulationEcsimCorr::make_diagnostic(const int timestep) {
    static Diagnostics diagnostic(outputParameters, domain, species);

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
    diagnostic.output_fields2D(timestep, fields);
    for (auto &sp : species) {
        sp->get_Pr();
        apply_periodic_border_with_add(sp->Pxx, bounds);
        apply_periodic_border_with_add(sp->Pyy, bounds);
        apply_periodic_border_with_add(sp->Pzz, bounds);

        const std::string pathToField =
            ".//Particles//" + sp->name() + "//Diag2D//";

        std::vector<std::pair<const Field3d &, std::string>> fields = {
            {sp->currentOnGrid, pathToField + "Current"},
            {sp->densityOnGrid, pathToField + "Density"},
            {sp->Pxx, pathToField + "Pxx"},
            {sp->Pyy, pathToField + "Pyy"},
            {sp->Pzz, pathToField + "Pzz"}};
        diagnostic.output_fields2D(timestep, fields);
    }
#ifdef SET_PARTICLE_IDS
    static ParticleTracker tracker(species, 1, "Tracking", "");
    tracker.track_particles(species, timestep);
#endif
}

void SimulationEcsimCorr::diagnostic_energy(
    Diagnostics &diagnostic) {
    double kineticEnergy = 0;
    double kineticEnergyNew = 0;
    for (auto &sp : species) {
        diagnostic.addEnergy(sp->name() + "Init",
                             sp->get_init_kinetic_energy());
        diagnostic.addEnergy(sp->name(), sp->get_kinetic_energy());
        diagnostic.addEnergy(sp->name() + "Inject", sp->injectionEnergy);
        diagnostic.addEnergy(sp->name() + "Lost", sp->lostEnergy);
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
              << fieldJe.data().dot(fieldE.data() + fieldEn.data()) << "\n";

    diagnostic.addEnergy(
        "energyConserve",
        std::abs(kineticEnergyNew - kineticEnergy + energyFieldDifference));
}
