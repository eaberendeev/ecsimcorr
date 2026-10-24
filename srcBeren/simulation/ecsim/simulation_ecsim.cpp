// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_ecsim.h"

#include <iomanip>
#include <iostream>
#include <map>
#include <source_location>
#include <sstream>
#include <string>

#include "Coil.h"
#include "Damping.h"
#include "Diagnostic.h"
#include "DiagnosticOutput.h"
#include "Mesh.h"
#include "ParticlesArray.h"
#include "ParticlesDiagnostic.h"
#include "Read.h"
#include "World.h"
#include "collision.h"
#include "containers.h"
#include "external_fieldsB.h"
#include "external_fieldsE.h"
#include "log_macros.h"
#include "recovery.h"
#include "row_block.h"
#include "solverSLE.h"
#include "timer.h"

void SimulationEcsim::assembleLmat2(double dt) {

    Operator tmpMat;

    if (assembly_check_counter_ % envOptions::validationPeriodicity() == 0) {
        timer::commonTimer timerCopy("copy matrix");
        tmpMat = mesh.Lmat2;
    }
    timer::commonTimer timerTestAssemble("new optimized assemble");
    for (int i = 0; i < omp_get_max_threads(); ++i) {
        rowBlocksGlobal[i].resizeAndReset(0);
    }

    for (auto &kv : species) {
        ParticlesArray &sp = *kv.second;
        sp.fill_matrixL2<SHAPE>(mesh, fieldBFull, domain, dt, rowBlocksGlobal);
    }
    mesh.stencil_Lmat2(mesh.Lmat2, rowBlocksGlobal, mesh.workspacePtr);
    timerTestAssemble.finish();

    if (assembly_check_counter_ % envOptions::validationPeriodicity() == 0) {
        timer::commonTimer timerRefAssemble("old assemble");
        prepare_block_matrix<SHAPE>();

        for (auto &kv : species) {
            ParticlesArray &sp = *kv.second;
            sp.fill_matrixL2_Reference<SHAPE>(mesh, fieldBFull, domain, dt);
        }
        mesh.stencil_Lmat2_Reference(tmpMat, domain);
        timerRefAssemble.finish();

        const bool isSameShape = checkMatrixPortraitCoincidence(mesh.Lmat2, tmpMat);
        if (!isSameShape) {
            const std::source_location location = std::source_location::current();
            std::cerr << location.file_name() << ":" << location.line()
                      << " Warning: optimized and reference assembly produced different matrix portraits" << std::endl;
        }
        timer::commonTimer timerNorms("compute norms");
        const double diffNorm = (mesh.Lmat2 - tmpMat).norm();
        const double refNorm = tmpMat.norm();
        timerNorms.finish();
        const double normalizedErr = refNorm == 0.0 && diffNorm == 0.0 ? 0.0 : diffNorm / refNorm;
        // Round-off tolerance: the optimized and reference assembly differ only
        // by the order of summation, so relative differences up to a few 1e-13
        // are expected and must not be flagged.
        constexpr double kAssemblyRelTol = 1e-10;
        if (normalizedErr >= kAssemblyRelTol) {
            const std::source_location location = std::source_location::current();
            std::cerr << location.file_name() << ":" << location.line()
                      << " Error between optimized and reference assembly is too large: normalized error = "
                      << normalizedErr << " >= " << kAssemblyRelTol << std::endl;
        }
    }
    assembly_check_counter_ += 1;
}

void SimulationEcsim::first_push() {
    RECORD_TIMER;

    const double dt = dt_;

    globalTimer.start("particles1");
    blas::sum(1.0, fieldB.data(), 1.0, fieldBInit.data(), fieldBFull.data());

    for (auto &kv : species) {
        auto &sp = *kv.second;
        //  +++ x_{n-1/2} -> x_{n+1/2}
        sp.move(dt);

        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);

        // +++ get J(x_{n+1/2},v_n)_predict
        algorithmsECSIM::predict_current<SHAPE>(sp, fieldBFull, fieldJp, dt);
    }
    globalTimer.finish("particles1");

    // <E0, J_raw> до apply_to_fields(CURRENT): учётная работа отфильтрованного
    // граничного тока (не мощность порта/Poynting)
    er0Work_ = dot_product_sum(fieldE0, fieldJp, bc_handler.active_range(domain.grid));
    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    assembleLmat2(dt);

    globalTimer.start("bound2");
    bc_handler.apply_to_operator(mesh.Lmat2, domain);
    globalTimer.finish("bound2");
}

void SimulationEcsim::second_push() {
    RECORD_TIMER;

    const double dt = dt_;

    globalTimer.start("particles2");

    blas::sum(1.0, fieldB.data(), 1.0, fieldBInit.data(), fieldBFull.data());
    Field3d fieldE_full = fieldEp + fieldE_external;
    for (auto &kv : species) {
        auto &sp = *kv.second;
        // +++ get v'_{n+1} from v_{n} and E'_{n+1/2}
        algorithmsECSIM::predict_velocity<SHAPE>(sp, fieldE_full, fieldBFull, dt);
    }
    globalTimer.finish("particles2");
}

// Particles have ccordinates and velocities. Mesh have 3D fields in nodes (each
// field stored in 1D array with 4d index x,y,z,d)
void SimulationEcsim::make_step([[maybe_unused]] const int timestep) {
    RECORD_TIMER;
    globalTimer.start("Total");

    first_push();

    globalTimer.start("FieldsPredict");
    // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})).
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);

    globalTimer.finish("FieldsPredict");

    second_push();

    for (auto &kv : species) {
        auto &sp = *kv.second;
        deposit_density(sp);
    }

    globalTimer.start("computeB");
    // calculate fieldB
    blas::sum(2.0, fieldEp.data(), -1.0, fieldE.data(), fieldEn.data());
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);

    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn, dt_);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
    globalTimer.finish("computeB");

    // Secondaries are flushed after this step's deposits (J, L-matrix, rho):
    // they enter the simulation consistently starting from the next step.
    bc_handler.flush_species(species);

    globalTimer.finish("Total");
}

template <ShapeType ST>
void SimulationEcsim::prepare_block_matrix() {
    RECORD_TIMER;
    Array3D<int> countInCell(domain.size());
    countInCell.setZero();
    for (auto &kv : species) {
        auto &sp = *kv.second;
#pragma omp parallel for schedule(static, 32)
        for (int i = 0; i < sp.particlesData.capacity(); i++) {
            countInCell(i) += sp.particlesData(i).size();
        }
    }

    static_assert(ST == ShapeType::Linear, "prepare_block_matrix supports Linear shape only");
    mesh.LmatX2.prepare(countInCell);
}

template void SimulationEcsim::prepare_block_matrix<ShapeType::Linear>();

void SimulationEcsim::predict_electric_field(Field3d &Ep, const Field3d &E, const Field3d &B, Field3d &J) {
    const double dt = dt_;

    Operator A = mesh.IMmat + mesh.Lmat2;
    mesh.Lmat2.makeCompressed();

    Field3d rhs = E - 0.5 * dt * J + 0.5 * dt * (mesh.curlB * B);
    impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);

    // E(n+1/2) = (M-L) * E(n+1/2) + E - 0.5*dt*(J + rotB)
    solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);
    LOG_STEP("  solver error=" << (A * Ep - rhs).norm() << "\n");
}

void SimulationEcsim::predict_electric_field(Field3d &Ep, const Field3d &E, const Field3d &E_ex, const Field3d &B,
                                             Field3d &J) {
    RECORD_TIMER;

    const double dt = dt_;

    timer::flatTimer timerDestructors(timer::NoStart{});
    {
        timer::commonTimer timerA("construct A");
        Operator A = parallelSparseSum(mesh.IMmat, mesh.Lmat2);
        timerA.finish();

        timer::commonTimer compressTimer("compress Lmat2");
        mesh.Lmat2.makeCompressed();
        compressTimer.finish();

        timer::commonTimer timerRhs("make rhs");
        Field3d rhs = E + 0.5 * dt * (mesh.curlB * B - J) - mesh.Lmat2 * E_ex;
        // Единая Дирихле-маска er0 (mode "solv") в правой части
        impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);
        timerRhs.finish();

        // E(n+1/2) = (M-L) * E(n+1/2)  - L*E_ex + E - 0.5*dt*(J + rotB)
        // (M*Ex = 0)
        const double err = solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);
        LOG_STEP("  solver error=" << err << "\n");

        // A и rhs уничтожаются при выходе из этого scope — замеряем их деструкторы
        timerDestructors.start("destructor operator A and rhs");
    }
}

void SimulationEcsim::init_operators() {
    Simulation::init_operators();
}
void SimulationEcsim::init_fields() {
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
    fieldE_external.resize(domain.size(), 3);

    fieldJp.setZero();
    fieldJe.setZero();

    fieldE.setZero();
    fieldB.setZero();

    if (get_checked<int>(system_config, "StartFromTime") > 0) {
        read_fields_from_recovery(fieldE, fieldB);
    }

    bc_handler.init_electric_field(fieldE, domain);

    // fieldEn должен получить уже установленный Дирихле-профиль fieldE
    fieldEn = fieldE;
    fieldBn = fieldB;

    // Поле-шаблон предписанных er0-значений (нули вне кольца) для учёта работы
    fieldE0.resize(domain.size(), 3);
    fieldE0.setZero();
    bc_handler.init_electric_field(fieldE0, domain);

    fieldE_external.setZero();
    if (auto e_cfg = create_electric_field_config(system_config, "ExternalFieldE")) {
        e_cfg->apply(fieldE_external, domain);
        std::cout << "Electric field config: " << fieldE_external.norm() << std::endl;
    }

    fieldBInit.setZero();
    if (auto b_cfg = create_magnetic_field_config(system_config, "ExternalFieldB")) {
        b_cfg->apply(fieldBInit, domain);
    }
    LmatX.resize(domain.total_size() * 3);
}

void SimulationEcsim::prepare_step(const int timestep) {
    RECORD_TIMER;
    const double dt = dt_;
    for (auto &kv : species) {
        auto &sp = *kv.second;
        sp.prepare();
        sp.diag.injection_energy = sp.inject_particles_step(sp.get_injection_distributions(), timestep, domain, dt,
                                                            system_config.value("co_locate_species", true));
    }

    dampingEnergy_ = damping_fields(fieldEn, fieldBn, domain, system_config);
    fieldE = fieldEn;
    fieldB = fieldBn;
    fieldJp.setZero();
    fieldJe.setZero();

}

SimulationEcsim::~SimulationEcsim() = default;

void SimulationEcsim::make_diagnostic(const int timestep) {
    RECORD_TIMER;

    if (timestep == 0) {
        fieldEp.setZero();

        for (auto &kv : species) {
            auto &sp = *kv.second;
            sp.currentOnGrid.setZero();
        }
    }

    if (!diagnostic_ptr_) {
        nlohmann::json diagnostic_config =
            system_config.contains("diagnostics") ? system_config["diagnostics"] : nlohmann::json::object();
        diagnostic_ptr_ = std::make_unique<Diagnostics>(diagnostic_config, domain, species);
        outputs_ = OutputFactory::create(diagnostic_config, domain, species, system_config, fieldEn, fieldBn,
                                         fieldE_external, fieldBInit);
    }

    diagnostic_energy(*diagnostic_ptr_);

    for (auto &out : outputs_) {
        out->output(timestep, *diagnostic_ptr_);
    }
}

void SimulationEcsim::finalize() {
    for (auto &out : outputs_) {
        out->finalize();
    }
    Simulation::finalize();
}

void SimulationEcsim::collect_per_species_diagnostics(Diagnostics &diagnostic, double &kineticEnergy,
                                                      double &kineticEnergyNew, double &totalLostEnergy) {
    static const Face all_faces[] = {Face::XMIN, Face::XMAX, Face::YMIN,    Face::YMAX,
                                     Face::ZMIN, Face::ZMAX, Face::CYLINDER};
    static const char *face_names[] = {"XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX", "CYLINDER"};
    static const PerFaceStats zero_stats;
    double totalInjectEnergy = 0;
    double totalEmitEnergy = 0;

    for (auto &kv : species) {
        auto &sp = *kv.second;
        const auto moments = sp.diagnostic_moments();
        diagnostic.addEnergy(sp.name(), moments.energy());
        diagnostic.addEnergy(sp.name() + "Particles", static_cast<double>(moments.count));
        diagnostic.addEnergy(sp.name() + "Inject", sp.diag.injection_energy);
        diagnostic.addEnergy(sp.name() + "Z", moments.kinetic.z());
        diagnostic.addEnergy(sp.name() + "XY", moments.kinetic.x() + moments.kinetic.y());
        kineticEnergy += sp.kineticEnergy;
        kineticEnergyNew += diagnostic.energy[sp.name()];
        totalInjectEnergy += sp.diag.injection_energy;

        diagnostic.addBoundary(sp.name() + "_InjectE", sp.diag.injection_energy);
        diagnostic.addBoundary(sp.name() + "_InjectN", static_cast<double>(sp.diag.injection_count));

        double sp_lost_energy = 0;
        for (int fi = 0; fi < 7; ++fi) {
            Face f = all_faces[fi];
            auto it = sp.diag.boundary.find(f);
            const auto &s = (it != sp.diag.boundary.end()) ? it->second : zero_stats;
            std::string prefix = sp.name() + "_" + face_names[fi];
            diagnostic.addBoundary(prefix + "LostE", s.lost_energy);
            diagnostic.addBoundary(prefix + "LostN", static_cast<double>(s.lost_count));
            diagnostic.addBoundary(prefix + "ReflN", static_cast<double>(s.reflected_count));
            diagnostic.addBoundary(prefix + "EmitE", s.emitted_energy);
            diagnostic.addBoundary(prefix + "EmitN", static_cast<double>(s.emitted_count));
            sp_lost_energy += s.lost_energy;
            totalEmitEnergy += s.emitted_energy;
        }
        totalLostEnergy += sp_lost_energy;


    }
    diagnostic.addEnergy("totalInjectEnergy", totalInjectEnergy);
    diagnostic.addEnergy("totalEmitEnergy", totalEmitEnergy);
}

void SimulationEcsim::compute_field_energy_and_conservation(Diagnostics &diagnostic, const IndexRange &irange,
                                                            double dt, double kineticEnergy, double kineticEnergyNew,
                                                            double totalLostEnergy, double totalInjectEnergy,
                                                            double energyJe_ex, double dampingEnergy, double er0Work) {
    diagnostic.addEnergy("totalLostEnergy", totalLostEnergy);
    diagnostic.addEnergy("energyDamping", dampingEnergy);
    diagnostic.addEnergy("energyFieldE", calc_energy_field(fieldEn, irange));
    diagnostic.addEnergy("energyFieldB", calc_energy_field(fieldBn, irange));
    fieldBFull.data() = fieldBn.data() + fieldBInit.data();
    diagnostic.addEnergy("energyFieldBFull", calc_energy_field(fieldBFull, irange));
    double energyFieldEold = calc_energy_field(fieldE, irange);
    double energyFieldBold = calc_energy_field(fieldB, irange);

    double energyFieldDifference =
        diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] - energyFieldBold - energyFieldEold;
    // emitted (secondary) particles are an external source, like injection
    double totalEmitEnergy = diagnostic.energy["totalEmitEnergy"];

    // er0Work_ = <E0, J_raw>, где J_raw — fieldJp ДО apply_to_fields(CURRENT).
    // Ebar равно E0 только на кольце и нулю на остальном замаскированном слое,
    // а J_solver там равен нулю, поэтому это дискретная работа отфильтрованного
    // граничного тока dt*<Ebar, J_raw - J_solver>. Это НЕ физическая мощность
    // порта/Poynting и не параметр для подгонки баланса.
    diagnostic.addEnergy("energyConserve",
                         std::abs(kineticEnergyNew - kineticEnergy - totalInjectEnergy - totalEmitEnergy +
                                  energyFieldDifference - dt * energyJe_ex - dt * er0Work + totalLostEnergy));
}

void SimulationEcsim::diagnostic_energy(Diagnostics &diagnostic) {
    double kineticEnergy = 0;
    double kineticEnergyNew = 0;
    double energyJe_ex = 0;
    double energyJe = 0;
    double totalLostEnergy = 0;
    IndexRange irange = bc_handler.active_range(domain.grid);

    collect_per_species_diagnostics(diagnostic, kineticEnergy, kineticEnergyNew, totalLostEnergy);

    for (auto &kv : species) {
        auto &sp = *kv.second;
        algorithmsECSIM::calculate_current(sp, sp.currentOnGrid);
        bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
        energyJe_ex += dot_product_sum(fieldE_external, sp.currentOnGrid, irange);
        energyJe += dot_product_sum(fieldEp, sp.currentOnGrid, irange);
    }

    const double dt = dt_;
    compute_field_energy_and_conservation(diagnostic, irange, dt, kineticEnergy, kineticEnergyNew, totalLostEnergy,
                                          diagnostic.energy["totalInjectEnergy"], energyJe_ex, dampingEnergy_, er0Work_);

    fieldJp_full.data() = fieldJp.data() + mesh.Lmat2 * (fieldE.data() + fieldEn.data()) / dt;
    double energyJe2 = dot_product_sum(fieldEp, fieldJp_full, irange);
    double energyFieldDifference = diagnostic.energy["energyFieldB"] + diagnostic.energy["energyFieldE"] -
                                   calc_energy_field(fieldE, irange) - calc_energy_field(fieldB, irange);

    diagnostic.addEnergy("deltaKinetic", kineticEnergyNew - kineticEnergy);
    diagnostic.addEnergy("deltaField", energyFieldDifference);
    diagnostic.addEnergy("JePredWork", dt * energyJe2);
    diagnostic.addEnergy("JeRealWork", dt * energyJe);
    diagnostic.addEnergy("JeExWork", dt * energyJe_ex);
    diagnostic.addEnergy("energyRemovedJWork", dt * er0Work_);
}
