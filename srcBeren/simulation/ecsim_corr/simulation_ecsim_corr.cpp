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
    const double n0 = parameters.get_double("n0");

    // ColliderWithNeutrals obj(n0);
    // std::pair<double, double> f = find_maximum_universal(
    //     [&obj](double x) { return obj.Sigma_e(x); }, 1e-16, 50, 1000000);
    
    // std::cout << "sigma_e " << f.first << " " << f.second << "\n";
    // f = find_maximum_universal(
    //     [&obj](double x) { return obj.Sigma_p(x); }, 1e-16, 50, 1000000);

    // std::cout << "sigma_p " << f.first << " " << f.second << "\n";
    // f = find_maximum_universal(
    //     [&obj](double x) { return obj.Sigma_cx(x); }, 1e-16, 50, 1000000);

    // std::cout << "sigma_cx " << f.first << " " << f.second << "\n";
    // std::cout << "sigma_p " << obj.Sigma_p(50*(0.001)*(0.001)) << "\n";
    const double dt = parameters.get_double("Dt");
    globalTimer.start("Total");
    std::cout << "timestep CORRECTION"<< "\n";
    globalTimer.start("densityCalc");
    for (auto &sp : species) {
        sp->density_on_grid_update(SHAPE_CH);   // calculate dendity field
    }
    globalTimer.finish("densityCalc");

    collect_charge_density(mesh.chargeDensityOld);
    mesh.apply_density_boundaries(mesh.chargeDensityOld, domain);

    globalTimer.start("particles1");
    fieldBFull.data() = fieldB.data() + fieldBInit.data();

    for (auto &sp_ref : charged_species) {
        auto &sp = sp_ref.get();
        sp.move_and_calc_current(0.5 * dt, sp.currentOnGrid, SHAPE_CH);
        //  +++ x_n -> x_{n+1/2}
        double t1 = omp_get_wtime();
        sp.update_cells(domain);
        std::cout << "time update cells " << omp_get_wtime() - t1 << "\n";
        // +++ get J(x_{n+1/2},v_n)_predict

        sp.predict_current(fieldBFull, fieldJp, domain, dt, SHAPE);
    }
    globalTimer.finish("particles1");

    globalTimer.start("particlesLmat2");


    // for (auto &sp : species) {
    //     sp->fill_matrixL(mesh, fieldBFull, domain, dt, SHAPE);
    // }

    prepare_block_matrix(SHAPE);

    for (auto &sp_ref : charged_species) {
        auto& sp = sp_ref.get();
        sp.fill_matrixL2(mesh, fieldBFull, domain, dt, SHAPE);
    }
    // todo: zeros Lmat + current
    globalTimer.finish("particlesLmat2");

   mesh.apply_boundaries(fieldJp, domain);

    globalTimer.start("bound1");
   // mesh.apply_boundaries(mesh.LmatX, domain);
    globalTimer.finish("bound1");

    globalTimer.start("stencilLmat2");

    mesh.stencil_Lmat2(mesh.Lmat2, domain);
   // convert_block_matrix(SHAPE);


    globalTimer.finish("stencilLmat2");

    //mesh.print_Lmat(Lmat2);
    globalTimer.start("bound2");
    mesh.apply_boundaries(mesh.Lmat2, domain);
    globalTimer.finish("bound2");
    //mesh.apply_boundaries(mesh.LmatX, domain);

   // Operator Lmat_compare(domain.total_size() * 3, domain.total_size() * 3); 
   // mesh.stencil_Lmat(Lmat_compare, domain);
  //  std::cout << "norm block convert_matrix " << (Lmat_compare - mesh.Lmat2).norm() << std::endl;

    globalTimer.start("FieldsPredict");
    // --- solve A*E'_{n+1}=f(E_n, B_n, J(x_{n+1/2})). mesh consist En,
    // En+1_predict
    // solve system of linear equations using Lmat for find fieldE
    mesh.predictE(fieldEp, fieldE, fieldB, fieldJp, dt);

    globalTimer.finish("FieldsPredict");

    globalTimer.start("particles2");

    fieldBFull.data() = fieldB.data() + fieldBInit.data();
    for (auto &sp_ref : charged_species) {
        auto& sp = sp_ref.get();
        // +++ get v'_{n+1} from v_{n} and E'_{n+1}
        sp.predict_velocity(fieldE, fieldEp, fieldBFull, domain, dt, SHAPE);
        // calc new particles velocity using new fieldE
        // +++ x_{n+1/2} -> x_{n+1}
        sp.move_and_calc_current(0.5 * dt, sp.currentOnGrid, SHAPE_CH);
        sp.update_cells(domain);
        sp.currentOnGrid.data() *= 0.5;
        mesh.apply_boundaries(sp.currentOnGrid, domain);
    }
    globalTimer.finish("particles2");

    collect_current(fieldJe);

    std::cout << "Current " << fieldJe.data().norm() << "\n";
        // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final

        globalTimer.start("FieldsCorr");
        // solve simple systeomof linear equations for correct fieldE
        mesh.correctE(fieldEn, fieldE, fieldB, fieldJe, dt);
        globalTimer.finish("FieldsCorr");

    fieldJp_full.data() =
        fieldJp.data() + mesh.Lmat2 * (fieldE.data() + fieldEn.data()) / dt;

    globalTimer.start("particles3");

        for (auto &sp_ref : charged_species) {
            auto& sp = sp_ref.get();
            // correct particles velocity for save energy
            correctv(sp, fieldJp_full, dt);
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
    mesh.apply_density_boundaries(mesh.chargeDensity, domain);

    std::cout << mesh.chargeDensity.data().norm()
              << " norm mesh.chargeDensity \n";

    auto divJ = mesh.divE * fieldJe.data();

    auto delta =
        (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (dt) +
        divJ;
    std::cout << delta.norm() << " norm drho / Dt - divJ \n";
    globalTimer.finish("Total");
}
