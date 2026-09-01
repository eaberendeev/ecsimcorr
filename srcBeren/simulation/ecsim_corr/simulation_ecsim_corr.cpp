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

void SimulationEcsimCorr::init_operators() {
    SimulationEcsim::init_operators();

    // M3: предобуславливатель фиксированного оператора corrector (IMmat).
    // Выбор через PRECONDITIONER: jacobi | lu | lu32 | none
    // (алиасы "cholesky"/"cholesky32" сохранены для совместимости)
    const char *pc = getenv("PRECONDITIONER");
    const std::string name = pc ? pc : "jacobi";
    if (name == "lu" || name == "cholesky") {
        precond_ = std::make_unique<SparseLUPreconditioner>(mesh.IMmat);
    } else if (name == "lu32" || name == "cholesky32") {
        precond_ = std::make_unique<SparseLUPreconditionerFp32>(mesh.IMmat);
    } else if (name == "jacobi") {
        precond_ = std::make_unique<JacobiPreconditioner>(mesh.IMmatDiag);
    } else {
        precond_ = nullptr;
    }
    std::cout << "Preconditioner: " << name << "\n";

    // Предобуславливатель predictor (независимый выбор; по умолчанию пусто —
    // тогда predictor идёт через M2 Jacobi + mixed precision)
    const char *pcp = getenv("PRECONDITIONER_PRED");
    const std::string namep = pcp ? pcp : "none";
    if (namep == "lu" || namep == "cholesky") {
        precond_pred_ = std::make_unique<SparseLUPreconditioner>(mesh.IMmat);
    } else if (namep == "lu32" || namep == "cholesky32") {
        precond_pred_ = std::make_unique<SparseLUPreconditionerFp32>(mesh.IMmat);
    } else if (namep == "jacobi") {
        precond_pred_ = std::make_unique<JacobiPreconditioner>(mesh.IMmatDiag);
    } else {
        precond_pred_ = nullptr;
    }
    std::cout << "Preconditioner pred: " << namep << "\n";

    // Кэш фиксированного оператора corrector (thread-partitioned, NUMA-friendly)
    IMmatTP_ = std::make_unique<ThreadPartitionedSparseMatrix<double>>(mesh.IMmat);
}

void SimulationEcsimCorr::second_push() {
    RECORD_TIMER;
    const double dt = get_checked<double>(system_config, "Dt");

    globalTimer.start("particles2");

    fieldBFull.data() = fieldB.data() + fieldBInit.data();
    Field3d fieldE_full = fieldEp + fieldE_external;

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        double pred_w = 0;
        fused_push_and_deposit(sp, fieldE_full, fieldBFull, dt, 0.5 * dt, sp.currentOnGrid, pred_w, SHAPE_CH);
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

    const double dt = get_checked<double>(system_config, "Dt");
    globalTimer.start("Total");
    globalTimer.start("densityCalc");
    for (auto &kv : species) {
        auto &sp = *kv.second;
        sp.density_on_grid_update(SHAPE_CH);   // calculate dendity field
        bc_handler.apply_to_fields(sp.densityOnGrid, FieldType::DENSITY, domain);
    }
    globalTimer.finish("densityCalc");

    collect_charge_density(mesh.chargeDensityOld);

    globalTimer.start("particles1");
    fieldBFull.data() = fieldB.data() + fieldBInit.data();

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        move_and_calc_current(sp, 0.5 * dt, sp.currentOnGrid, SHAPE_CH);
        //  +++ x_n -> x_{n+1/2}
        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);

        // +++ get J(x_{n+1/2},v_n)_predict

        algorithmsECSIM::predict_current(sp, fieldBFull, fieldJp, dt, SHAPE);
    }
    globalTimer.finish("particles1");
    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    globalTimer.start("particlesLmat2");

    // for (auto& kv : charged_species) {
    //     auto& sp = kv.second.get();
    //     fill_matrixL(sp, LmatX, fieldBFull, domain, dt, SHAPE);
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

    mesh.stencil_Lmat2(mesh.Lmat2, domain, mesh.workspacePtr);
    // convert_block_matrix(SHAPE);

    globalTimer.finish("stencilLmat2");

    // mesh.print_Lmat(Lmat2);
    globalTimer.start("bound2");
    bc_handler.apply_to_operator(mesh.Lmat2, domain);
    globalTimer.finish("bound2");
    // mesh.apply_boundaries(mesh.LmatX, domain);

    // Operator Lmat_compare(domain.total_size() * 3, domain.total_size() * 3);
    // mesh.stencil_Lmat(Lmat_compare, domain);
    //  std::cout << "norm block convert_matrix " << (Lmat_compare -
    //  mesh.Lmat2).norm() << std::endl;

    globalTimer.start("FieldsPredict");
    // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})).
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);

    globalTimer.finish("FieldsPredict");

    globalTimer.start("particles2");
    // Hack!!!
    // We postponed updating the current positions and working with boundary conditions (particle removal) until the
    // velocities was corrected and the energy balance was balanced.
    second_push();
    globalTimer.finish("particles2");

    collect_current(fieldJe);

    // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final

    globalTimer.start("FieldsCorr");
    // ---- get E_{n+1} from E_n and J_e. mesh En changed to En+1_final

    // P1: дефектный пропуск corrector (следствие теоремы (T)):
    // если E_{n+1} := 2E' − E_n (т.е. δ=0), то истинная невязка corrector-системы
    // равна дефекту rhs_def = dt(J'−J_e) + 2L(E'+E_ex), и дрейф энергии за шаг
    // = −⟨E_avg, rhs_def⟩ ≤ |E_avg|·|rhs_def|. Поэтому corrector можно пропускать,
    // когда |rhs_def| ≤ τ·|rhs| (τ = SKIP_CORR_TOL, 0 — выключено).
    const double tau_skip = getenvParsed<double>("SKIP_CORR_TOL", 0.0);
    bool corr_skipped = false;
    if (tau_skip > 0) {
        Field3d rhs_corr = fieldE + dt * (mesh.curlB * fieldB - fieldJe) + mesh.Mmat * fieldE;
        Field3d rhs_def = dt * (fieldJp - fieldJe) + 2.0 * mesh.Lmat2 * (fieldEp + fieldE_external);
        const double rel_defect = rhs_def.norm() / rhs_corr.norm();
        if (rel_defect <= tau_skip) {
            fieldEn = 2.0 * fieldEp - fieldE;
            corr_skipped = true;
            LOG_STEP("  corr SKIPPED: |defect|/|rhs| = " << rel_defect << "\n");
        }
    }
    if (!corr_skipped) {
        correctE(fieldEn, fieldE, fieldEp, fieldB, fieldJe, dt);
    }
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    globalTimer.finish("FieldsCorr");

    globalTimer.start("particles3");

    for (auto &kv : charged_species) {
        auto &sp = kv.second.get();
        // correct particles velocity for save energy
        correctv(sp, dt);
        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);
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
    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn, get_checked<double>(system_config, "Dt"));
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
    globalTimer.finish("computeB");

    // check charge conservation: drho/dt + divJ = 0
    collect_charge_density(mesh.chargeDensity);

    auto divJ = mesh.divE * fieldJe.data();
    auto delta = (mesh.chargeDensity.data() - mesh.chargeDensityOld.data()) / (dt) + divJ;
    LOG_STEP("  |drho/dt + divJ| = " << delta.norm() << "\n");
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

    const double dt = get_checked<double>(system_config, "Dt");
    compute_field_energy_and_conservation(diagnostic, irange, dt, kineticEnergy, kineticEnergyNew, totalLostEnergy,
                                          diagnostic.energy["totalInjectEnergy"], energyJe_ex, dampingEnergy_);

    double energyFieldDifference = diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] -
                                   calc_energy_field(fieldE, irange) - calc_energy_field(fieldB, irange);

    diagnostic.addEnergy("deltaKinetic", kineticEnergyNew - kineticEnergy);
    diagnostic.addEnergy("deltaField", energyFieldDifference);
    diagnostic.addEnergy("JeCorrWork", dt * energyJe);
    diagnostic.addEnergy("JeExWork", dt * energyJe_ex);
}
