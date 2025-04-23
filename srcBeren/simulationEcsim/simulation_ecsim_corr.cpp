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
#include "ParticlesDiagnostic.h"
#include "Read.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "recovery.h"


// Particles have ccordinates and velocities. Mesh have 3D fields in nodes (each field stored in 1D array with 4d index x,y,z,d)
void SimulationEcsimCorr::make_step([[maybe_unused]] const int timestep) {

    const double dt = parameters.get_double("Dt");
    const bool useCorrection = USE_ECSIM_CORRECTION;
    globalTimer.start("Total");

    globalTimer.start("densityCalc");
    for (auto &sp : species) {
        sp->density_on_grid_update(SHAPE_CH);   // calculate dendity field
    }
    globalTimer.finish("densityCalc");

    collect_charge_density(mesh.chargeDensityOld);
    mesh.apply_density_boundaries(mesh.chargeDensityOld, domain);

    globalTimer.start("particles1");
    fieldBFull.data() = fieldB.data() + fieldBInit.data();

    for (auto &sp : species) {
        sp->move_and_calc_current(0.5 * dt, sp->currentOnGrid, SHAPE_CH);
        //  +++ x_n -> x_{n+1/2}

        sp->update_cells(domain);
        // +++ get J(x_{n+1/2},v_n)_predict

        sp->predict_current(fieldBFull, fieldJp, domain, dt, SHAPE);
    }
    globalTimer.finish("particles1");

    globalTimer.start("particlesLmat2");

    Array3D<int> countInCell(domain.size());
    countInCell.setZero();
    for (auto &sp : species) {
        for(int i = 0; i < sp->countInCell.capacity(); i++) {
            countInCell(i) += sp->countInCell(i);
        }
    }
    mesh.LmatX2.prepare(countInCell);
    for (auto &sp : species) {
        sp->fill_matrixL2(mesh, fieldBFull, domain, dt, SHAPE);
    }
    // todo: zeros Lmat + current
    globalTimer.finish("particlesLmat2");

   mesh.apply_boundaries(fieldJp, domain);

    globalTimer.start("bound1");
   // mesh.apply_boundaries(mesh.LmatX, domain);
    globalTimer.finish("bound1");

    globalTimer.start("stencilLmat2");

    mesh.stencil_Lmat2(mesh.Lmat2, domain);

    globalTimer.finish("stencilLmat2");

    //mesh.print_Lmat(Lmat2);
    globalTimer.start("bound2");
    mesh.apply_boundaries(mesh.Lmat2, domain);
    globalTimer.finish("bound2");

    globalTimer.start("FieldsPredict");
    // --- solve A*E'_{n+1}=f(E_n, B_n, J(x_{n+1/2})). mesh consist En,
    // En+1_predict
    // solve system of linear equations using Lmat for find fieldE
    mesh.predictE(fieldEp, fieldE, fieldB, fieldJp, dt);

    globalTimer.finish("FieldsPredict");

    globalTimer.start("particles2");

    fieldBFull.data() = fieldB.data() + fieldBInit.data();
    for (auto &sp : species) {
        // +++ get v'_{n+1} from v_{n} and E'_{n+1}
        sp->predict_velocity(fieldE, fieldEp, fieldBFull, domain, dt, SHAPE);
        // calc new particles velocity using new fieldE
        // +++ x_{n+1/2} -> x_{n+1}
        sp->move_and_calc_current(0.5 * dt, sp->currentOnGrid, SHAPE_CH);
        sp->update_cells(domain);
        sp->currentOnGrid.data() *= 0.5;
        mesh.apply_boundaries(sp->currentOnGrid, domain);
    }
    globalTimer.finish("particles2");

    collect_current(fieldJe);

    std::cout << "Current " << fieldJe.data().norm() << "\n";

    if (useCorrection) {
        // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final

        globalTimer.start("FieldsCorr");
        // solve simple systeomof linear equations for correct fieldE
        mesh.correctE(fieldEn, fieldE, fieldB, fieldJe, dt);
        globalTimer.finish("FieldsCorr");
    } else {
        fieldEn.data() = fieldEp.data();
    }
    fieldJp_full.data() =
        fieldJp.data() + mesh.Lmat2 * (fieldE.data() + fieldEn.data()) / dt;

    globalTimer.start("particles3");
   // std::cout << "before " << "\n";
    if (useCorrection) {
        for (auto &sp : species) {
            // correct particles velocity for save energy
            sp->correctv(fieldE, fieldEp, fieldEn, fieldJp_full, domain, dt);
        }
    }
    // std::cout << "After " <<"\n";

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
  //  mesh.apply_density_boundaries(mesh.chargeDensity, domain);

    std::cout << mesh.chargeDensity.data().norm()
              << " norm mesh.chargeDensity \n";

    auto divJ = mesh.divE * fieldJe.data();

    auto delta =
        (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (dt) +
        divJ;
    std::cout << delta.norm() << " norm drho / Dt - divJ \n";
    globalTimer.finish("Total");
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

    fieldJp.setZero();
    fieldJe.setZero();

    fieldE.setZero();
    fieldB.setZero();

    if (parameters.get_int("StartFromTime") > 0) {
        read_fields_from_recovery(fieldE, fieldB);
    } else {
        mesh.set_uniform_field(fieldE, 0.0, 0.0, 0.0);
    }
    mesh.set_uniform_field(fieldBInit, parameters.get_double("BUniform", 0),
                           parameters.get_double("BUniform", 1),
                           parameters.get_double("BUniform", 2));
    set_coils(fieldBInit, domain, parameters);

    fieldEn = fieldE;
    fieldBn = fieldB;
}

void SimulationEcsimCorr::prepare_step(const int timestep) {
    inject_particles(timestep, domain);
    damping_fields(fieldEn, fieldBn, domain,
                   parameters);
    fieldE = fieldEn;
    fieldB = fieldBn;
    fieldJp.setZero();
    fieldJe.setZero();

// #pragma omp parallel for
//     for (size_t i = 0; i < mesh.LmatX.size(); i++) {
//         for (auto it = mesh.LmatX[i].begin(); it != mesh.LmatX[i].end(); ++it) {
//             it->second = 0.;
//         }
//     }
    for (auto &sp : species) {
        sp->prepare();   // save start coord for esirkepov current
        std::cout << sp->get_total_num_of_particles() << " size \n";
    }
}

void SimulationEcsimCorr::collision_step([[maybe_unused]] const int timestep) {
    const double dt = parameters.get_double("Dt");
    const double n0 = parameters.get_double("n0");

    static BinaryCollider collider(n0);

    // TODO : make different colliders type
    collider.collide_with_neutrals_binary(species, dt);
    // double start, end;
    // start = 0;

    // for (auto &sp : species) {
    //     start += sp->get_init_kinetic_energy();
    // }

    // collider.collide_same_sort_binary(species, dt);
    // collider.collide_ion_electron_binary(species, dt);
    // end = 0;
    // for (auto &sp : species) {
    //     end += sp->get_init_kinetic_energy();
    // }
    // std::cout << "Kinetic energy " << (end - start)/start << "\n";
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

    fieldBFull.data() = fieldBn.data() + fieldBInit.data();

    std::vector<std::pair<const Field3d &, std::string>> fields = {
        {fieldEn, pathToField + "FieldE"},
        {fieldBFull, pathToField + "FieldB"}};
    diagnostic.output_fields2D(timestep, fields);
    for (auto &sp : species) {
        
        const std::string spectrumPath = ".//Particles//" + sp->name() + "//";
        EnergySpectrum spectrum = sp->calculate_energy_spectrum();
        diagnostic.output_energy_spectrum(spectrum, timestep, spectrumPath);

        const std::string pathToField =
            ".//Particles//" + sp->name() + "//Diag2D//";
        if(sp->is_neutral()){
            std::vector<std::pair<const Field3d &, std::string>> fields = {
                {sp->densityOnGrid, pathToField + "Density"}};
            diagnostic.output_fields2D(timestep, fields);
            continue;
        }
        Field3d pressureRR(domain.size(), 1);
        Field3d pressurePP(domain.size(), 1);
        Field3d pressureZZ(domain.size(), 1);
        Field3d pressureRP(domain.size(), 1);
        Field3d pressureRZ(domain.size(), 1);
        Field3d pressureZP(domain.size(), 1);
        // Использование:
        sp->calculate_pressure_component(pressureRR, RadialVelocity{}, RadialVelocity{});
        sp->calculate_pressure_component(pressurePP, PhiVelocity{}, PhiVelocity{});
        sp->calculate_pressure_component(pressureZZ, ZVelocity{}, ZVelocity{});
        sp->calculate_pressure_component(pressureRP, RadialVelocity{},
                                         PhiVelocity{});
        sp->calculate_pressure_component(pressureRZ, RadialVelocity{}, ZVelocity{});
        sp->calculate_pressure_component(pressureZP, ZVelocity{},
                                         PhiVelocity{});
        std::vector<std::pair<const Field3d &, std::string>> fields = {
            {sp->currentOnGrid, pathToField + "Current"},
            {sp->densityOnGrid, pathToField + "Density"},
            {pressureRR, pathToField + "Prr"},
            {pressurePP, pathToField + "Ppp"},
            {pressureZZ, pathToField + "Pzz"},
            {pressureRP, pathToField + "Prp"},
            {pressureRZ, pathToField + "Prz"},
            {pressureZP, pathToField + "Pzp"}};
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
        diagnostic.addEnergy(sp->name() + "LostZ", sp->lostEnergyZ);
        diagnostic.addEnergy(sp->name() + "LostXY", sp->lostEnergyXY);
        diagnostic.addEnergy(sp->name() + "Z", sp->get_kinetic_energy(Z));
        diagnostic.addEnergy(sp->name() + "XY", sp->get_kinetic_energy(X, Y));
        kineticEnergy += diagnostic.energy[sp->name() + "Init"];
        kineticEnergyNew += diagnostic.energy[sp->name()];
    }

    diagnostic.addEnergy("energyFieldE", mesh.calc_energy_field(fieldEn));
    diagnostic.addEnergy("energyFieldB", mesh.calc_energy_field(fieldBn));
    fieldBFull.data() = fieldBn.data() + fieldBInit.data();
    diagnostic.addEnergy("energyFieldBFull",
                         mesh.calc_energy_field(fieldBFull));
    double energyFieldEold = mesh.calc_energy_field(fieldE);
    double energyFieldBold = mesh.calc_energy_field(fieldB);

    double energyFieldDifference = diagnostic.energy["energyFieldB"] +
                                   diagnostic.energy["energyFieldE"] -
                                   energyFieldBold - energyFieldEold;

    fieldJp_full.data() =
        fieldJp.data() + mesh.Lmat2 * (fieldE.data() + fieldEn.data()) / 1;

    std::cout << "Energy " << kineticEnergyNew - kineticEnergy << " "
              << energyFieldDifference << " "
              << 0.5*fieldJp_full.data().dot(fieldE.data() + fieldEn.data())
              << "\n";

    diagnostic.addEnergy(
        "energyConserve",
        std::abs(kineticEnergyNew - kineticEnergy + energyFieldDifference));
}
