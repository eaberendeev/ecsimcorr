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
#include "cross_section.h"
#include "collisions_with_neutrals.h"

// Particles have ccordinates and velocities. Mesh have 3D fields in nodes (each field stored in 1D array with 4d index x,y,z,d)
void SimulationEcsimCorr::make_step([[maybe_unused]] const int timestep) {
    const double dt = get_checked<double>(system_config, "Dt");
    globalTimer.start("Total");
    std::cout << "timestep CORRECTION"<< "\n";
    globalTimer.start("densityCalc");
    for (auto &kv : species) {
        auto &sp = *kv.second;
        sp.density_on_grid_update(SHAPE_CH);   // calculate dendity field
        bc_handler.apply_to_fields(sp.densityOnGrid, FieldType::DENSITY,
                                   domain);
    }
    globalTimer.finish("densityCalc");

    collect_charge_density(mesh.chargeDensityOld);

    globalTimer.start("particles1");
    fieldBFull.data() = fieldB.data() + fieldBInit.data();

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        sp.move_and_calc_current(0.5 * dt, sp.currentOnGrid, SHAPE_CH);
        //  +++ x_n -> x_{n+1/2}
        double t1 = omp_get_wtime();
        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);

        std::cout << "time update cells " << omp_get_wtime() - t1 << "\n";
        // +++ get J(x_{n+1/2},v_n)_predict

        algorithmsECSIM::predict_current(sp, fieldBFull, fieldJp, dt,
                                         SHAPE);
    }
    globalTimer.finish("particles1");
    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    globalTimer.start("particlesLmat2");


    // for (auto &sp : species) {
    //     sp->fill_matrixL(mesh, fieldBFull, domain, dt, SHAPE);
    // }

    prepare_block_matrix(SHAPE);

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        sp.fill_matrixL2(mesh, fieldBFull, domain, dt, SHAPE);
    }
    // todo: zeros Lmat + current
    globalTimer.finish("particlesLmat2");


    globalTimer.start("bound1");
   // mesh.apply_boundaries(mesh.LmatX, domain);
    globalTimer.finish("bound1");

    globalTimer.start("stencilLmat2");

    mesh.stencil_Lmat2(mesh.Lmat2, domain);
   // convert_block_matrix(SHAPE);


    globalTimer.finish("stencilLmat2");

    //mesh.print_Lmat(Lmat2);
    globalTimer.start("bound2");
    bc_handler.apply_to_operator(mesh.Lmat2, domain);
    globalTimer.finish("bound2");
    //mesh.apply_boundaries(mesh.LmatX, domain);

   // Operator Lmat_compare(domain.total_size() * 3, domain.total_size() * 3); 
   // mesh.stencil_Lmat(Lmat_compare, domain);
  //  std::cout << "norm block convert_matrix " << (Lmat_compare - mesh.Lmat2).norm() << std::endl;

    globalTimer.start("FieldsPredict");
    // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})).
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);

    globalTimer.finish("FieldsPredict");

    globalTimer.start("particles2");

    Field3d fieldE_full = fieldEp + fieldE_external;

    fieldBFull.data() = fieldB.data() + fieldBInit.data();

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        algorithmsECSIM::predict_velocity(sp, fieldE_full, fieldBFull, dt,
                                          SHAPE);

        // calc new particles velocity using new fieldE
        // +++ x_{n+1/2} -> x_{n+1}
        sp.move_and_calc_current(0.5 * dt, sp.currentOnGrid, SHAPE_CH);
        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);

        sp.currentOnGrid.data() *= 0.5;
        //mesh.apply_boundaries(sp.currentOnGrid, domain);
        bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT,
                                   domain);
    }
    globalTimer.finish("particles2");

    collect_current(fieldJe);

    std::cout << "Current " << fieldJe.data().norm() << "\n";
        // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final

    globalTimer.start("FieldsCorr");
    // solve simple systeomof linear equations for correct fieldE
    mesh.correctE(fieldEn, fieldE, fieldB, fieldJe, dt);
    globalTimer.finish("FieldsCorr");

    globalTimer.start("particles3");

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        // correct particles velocity for save energy
        correctv(sp, dt);
    }

    // std::cout << "After " <<"\n";

    for (auto &kv : species) {
        auto &sp = *kv.second;
        sp.density_on_grid_update(SHAPE_CH);
        bc_handler.apply_to_fields(sp.densityOnGrid, FieldType::DENSITY, domain);
    }
    globalTimer.finish("particles3");

    globalTimer.start("computeB");
    // calculate fieldB
    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn,
                        get_checked<double>(system_config, "Dt"));
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
    globalTimer.finish("Total");
}

void SimulationEcsimCorr::diagnostic_energy(Diagnostics &diagnostic) {
    IndexRange irange = bc_handler.active_range(domain.grid);
    Field3d fieldEn12 = 0.5 * (fieldE + fieldEn);
    double kineticEnergy = 0;
    double kineticEnergyNew = 0;
    double energyJe_ex = 0;
    double energyJe = 0;
    for (auto &kv : species) {
        auto &sp = *kv.second;
        diagnostic.addEnergy(sp.name() + "Init", sp.get_init_kinetic_energy());
        diagnostic.addEnergy(sp.name(), sp.get_kinetic_energy());
        diagnostic.addEnergy(sp.name() + "Particles",
                             sp.get_total_num_of_particles());
        diagnostic.addEnergy(sp.name() + "Inject", sp.injectionEnergy);
        diagnostic.addEnergy(sp.name() + "LostEnergyZ", sp.lostEnergyZ);
        diagnostic.addEnergy(sp.name() + "LostEnergyXY", sp.lostEnergyXY);
        diagnostic.addEnergy(sp.name() + "LostParticlesZ", sp.lostParticlesZ);
        diagnostic.addEnergy(sp.name() + "LostParticlesXY", sp.lostParticlesXY);
        diagnostic.addEnergy(sp.name() + "Z", sp.get_kinetic_energy(Z));
        diagnostic.addEnergy(sp.name() + "XY", sp.get_kinetic_energy(X, Y));
        kineticEnergy += diagnostic.energy[sp.name() + "Init"];
        kineticEnergyNew += diagnostic.energy[sp.name()];
        sp.lostEnergyZ = sp.lostEnergyXY = 0;
        sp.lostParticlesZ = sp.lostParticlesXY = 0;

        energyJe_ex +=
            dot_product_sum(fieldE_external, sp.currentOnGrid, irange);
        energyJe += dot_product_sum(fieldEn12, sp.currentOnGrid, irange);
    }

    diagnostic.addEnergy("energyFieldE", calc_energy_field(fieldEn, irange));
    diagnostic.addEnergy("energyFieldB", calc_energy_field(fieldBn, irange));
    fieldBFull.data() = fieldBn.data() + fieldBInit.data();
    diagnostic.addEnergy("energyFieldBFull",
                         calc_energy_field(fieldBFull, irange));
    double energyFieldEold = calc_energy_field(fieldE, irange);
    double energyFieldBold = calc_energy_field(fieldB, irange);

    double energyFieldDifference = diagnostic.energy["energyFieldB"] +
                                   diagnostic.energy["energyFieldE"] -
                                   energyFieldBold - energyFieldEold;

    const double dt = get_checked<double>(system_config, "Dt");

    std::cout << "Energy " << kineticEnergyNew - kineticEnergy << " "
              << energyFieldDifference << " " << dt * energyJe << " "
              << dt * energyJe_ex << "\n";

    diagnostic.addEnergy("energyConserve",
                         std::abs(kineticEnergyNew - kineticEnergy +
                                  energyFieldDifference - dt * energyJe_ex));
}
