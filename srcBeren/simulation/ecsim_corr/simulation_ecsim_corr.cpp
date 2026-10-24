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
#include "algorithms_ecsim_corr.h"
#include "collision.h"
#include "collisions_with_neutrals.h"
#include "containers.h"
#include "cross_section.h"
#include "log_macros.h"
#include "recovery.h"
#include "timer.h"

void SimulationEcsimCorr::deposit_density(ParticlesArray &sp) {
    sp.density_on_grid_update<SHAPE_CH>();
    bc_handler.apply_to_fields(sp.densityOnGrid, FieldType::DENSITY, domain);
}

void SimulationEcsimCorr::second_push() {
    RECORD_TIMER;
    const double dt = dt_;

    globalTimer.start("particles2");

    fieldBFull.data() = fieldB.data() + fieldBInit.data();
    Field3d fieldE_full = fieldEp + fieldE_external;

    boundary_events_.clear();
    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        double pred_w = 0;
        fused_push_and_deposit<SHAPE_CH>(sp, fieldE_full, fieldBFull, dt, 0.5 * dt, sp.currentOnGrid, pred_w,
                                         bc_handler.path_policy(sp, domain, &boundary_events_[sp.name()]));
        pred_work_[sp.name()] = pred_w;
        sp.currentOnGrid.data() *= 0.5;
        bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
    }

    globalTimer.finish("particles2");
}

// Particles have ccordinates and velocities. Mesh have 3D fields in nodes (each
// field stored in 1D array with 4d index x,y,z,d)
void SimulationEcsimCorr::make_step([[maybe_unused]] const int timestep) {
    RECORD_TIMER;

    const double dt = dt_;
    globalTimer.start("Total");
    globalTimer.start("densityCalc");
    for (auto &kv : species) {
        auto &sp = *kv.second;
        deposit_density(sp);
    }
    globalTimer.finish("densityCalc");

    collect_charge_density(mesh.chargeDensityOld);
    advance_neutral_particles(dt);

    globalTimer.start("particles1");
    fieldBFull.data() = fieldB.data() + fieldBInit.data();

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        ParticleBoundaryEvents events;
        move_and_calc_current<SHAPE_CH>(sp, 0.5 * dt, sp.currentOnGrid, bc_handler.path_policy(sp, domain, &events));
        // x_n -> x_{n+1/2}: consume the same wall event as the current deposit,
        // before particle addresses change during cell redistribution.
        bc_handler.apply_to_particles(sp, species, domain, &events);
        sp.update_cells(domain);

        // +++ get J(x_{n+1/2},v_n)_predict

        algorithmsECSIM::predict_current<SHAPE>(sp, fieldBFull, fieldJp, dt);
    }
    globalTimer.finish("particles1");
    // <E0, J_raw> до apply_to_fields(CURRENT): учётная работа отфильтрованного
    // граничного тока (не мощность порта/Poynting)
    er0Work_ = dot_product_sum(fieldE0, fieldJp, bc_handler.active_range(domain.grid));
    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    globalTimer.start("bound1");
    globalTimer.finish("bound1");

    globalTimer.start("stencilLmat2");
    assembleLmat2(dt);
    globalTimer.finish("stencilLmat2");

    globalTimer.start("bound2");
    bc_handler.apply_to_operator(mesh.Lmat2, domain);
    globalTimer.finish("bound2");

    globalTimer.start("FieldsPredict");
    // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})).
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);

    globalTimer.finish("FieldsPredict");

    globalTimer.start("particles2");
    // Keep outgoing particles through the velocity correction; their boundary
    // decisions are saved by the deposit and applied afterwards.
    second_push();
    globalTimer.finish("particles2");

    collect_current(fieldJe);

    // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final

    globalTimer.start("FieldsCorr");
    // solve simple systeomof linear equations for correct fieldE
    correctE(fieldEn, fieldE, fieldB, fieldJe, dt);
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    globalTimer.finish("FieldsCorr");

    globalTimer.start("particles3");

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        // correct particles velocity for save energy
        // Freeze the fate at crossing velocity: correctv must not change the
        // absorption/reflection decision after the matching current deposit.
        correctv(sp, dt);
        bc_handler.apply_to_particles(sp, species, domain, &boundary_events_.at(sp.name()));
        sp.update_cells(domain);
    }

    for (auto &kv : species) {
        auto &sp = *kv.second;
        deposit_density(sp);
    }
    globalTimer.finish("particles3");

    globalTimer.start("computeB");
    // calculate fieldB
    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn, dt_);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
    globalTimer.finish("computeB");

    // check charge conservation: drho/dt + divJ = 0
    collect_charge_density(mesh.chargeDensity);

    auto divJ = mesh.divE * fieldJe.data();
    auto delta = (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (dt) + divJ;
    LOG_STEP("  |drho/dt + divJ| = " << delta.norm() << "\n");

    // Secondaries are flushed after this step's deposits (J, L-matrix, rho):
    // they enter the simulation consistently starting from the next step.
    bc_handler.flush_species(species);

    globalTimer.finish("Total");
}

void SimulationEcsimCorr::diagnostic_energy(Diagnostics &diagnostic) {
    IndexRange irange = bc_handler.active_range(domain.grid);
    Field3d fieldEn12 = 0.5 * (fieldE + fieldEn);
    double kineticEnergy = 0;
    double kineticEnergyNew = 0;
    double energyJe_ex = 0;
    double energyJe = 0;
    double totalLostEnergy = 0;

    collect_per_species_diagnostics(diagnostic, kineticEnergy, kineticEnergyNew, totalLostEnergy);

    for (auto &kv : species) {
        auto &sp = *kv.second;
        energyJe_ex += dot_product_sum(fieldE_external, sp.currentOnGrid, irange);
        energyJe += dot_product_sum(fieldEn12, sp.currentOnGrid, irange);
    }

    const double dt = dt_;
    compute_field_energy_and_conservation(diagnostic, irange, dt, kineticEnergy, kineticEnergyNew, totalLostEnergy,
                                          diagnostic.energy["totalInjectEnergy"], energyJe_ex, dampingEnergy_,
                                          er0Work_);

    double energyFieldDifference = diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] -
                                   calc_energy_field(fieldE, irange) - calc_energy_field(fieldB, irange);

    diagnostic.addEnergy("deltaKinetic", kineticEnergyNew - kineticEnergy);
    diagnostic.addEnergy("deltaField", energyFieldDifference);
    diagnostic.addEnergy("JeCorrWork", dt * energyJe);
    diagnostic.addEnergy("JeExWork", dt * energyJe_ex);
    diagnostic.addEnergy("energyRemovedJWork", dt * er0Work_);
}
