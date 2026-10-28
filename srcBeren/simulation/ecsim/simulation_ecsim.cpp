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
#include "DiagnosticMeasurements.h"
#include "EnergyMeasurements.h"
#include "Mesh.h"
#include "OutputFactory.h"
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

    for (auto& kv : species) {
        ParticlesArray& sp = *kv.second;
        sp.fill_matrixL2<SHAPE>(mesh, fieldBFull, domain, dt, rowBlocksGlobal);
    }
    mesh.stencil_Lmat2(mesh.Lmat2, rowBlocksGlobal, mesh.workspacePtr);
    timerTestAssemble.finish();

    if (assembly_check_counter_ % envOptions::validationPeriodicity() == 0) {
        timer::commonTimer timerRefAssemble("old assemble");
        prepare_block_matrix<SHAPE>();

        for (auto& kv : species) {
            ParticlesArray& sp = *kv.second;
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
        // Optimized/reference sums differ by round-off. This warning measures
        // Frobenius-relative deviation without changing the solver tolerance.
        constexpr double expectedTol = 1e-14;
        if (normalizedErr >= expectedTol) {
            const std::source_location location = std::source_location::current();
            std::cerr << location.file_name() << ":" << location.line() << " Warning: assembly deviation "
                      << normalizedErr << " exceeds expectedTol=" << expectedTol << std::endl;
        }
    }
    assembly_check_counter_ += 1;
}

void SimulationEcsim::drift_particles_and_assemble_response() {
    RECORD_TIMER;

    const double dt = dt_;

    globalTimer.start("drift_and_response");
    blas::sum(1.0, fieldB.data(), 1.0, fieldBInit.data(), fieldBFull.data());

    for (auto& kv : species) {
        auto& sp = *kv.second;
        //  +++ x_{n-1/2} -> x_{n+1/2}
        sp.move(dt);

        sp.update_cells(domain);
        bc_handler.apply_to_particles(sp, species, domain);

        // +++ get J(x_{n+1/2},v_n)_predict
        algorithmsECSIM::predict_current<SHAPE>(sp, fieldBFull, fieldJp, dt);
    }
    globalTimer.finish("drift_and_response");

    // Discrete filtered-current power <E0,Jp_raw>, measured before the mask.
    filtered_current_power_ = dot_product_sum(fieldE0, fieldJp, bc_handler.active_range(domain.grid));
    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    assembleLmat2(dt);

    globalTimer.start("response_boundary_constraints");
    bc_handler.apply_to_operator(mesh.Lmat2, domain);
    globalTimer.finish("response_boundary_constraints");
}

void SimulationEcsim::advance_particle_velocities() {
    RECORD_TIMER;

    const double dt = dt_;

    globalTimer.start("advance_velocities");

    blas::sum(1.0, fieldB.data(), 1.0, fieldBInit.data(), fieldBFull.data());
    Field3d fieldE_full = fieldEp + fieldE_external;
    for (auto& kv : species) {
        auto& sp = *kv.second;
        // +++ get v'_{n+1} from v_{n} and E'_{n+1/2}
        algorithmsECSIM::predict_velocity<SHAPE>(sp, fieldE_full, fieldBFull, dt);
    }
    globalTimer.finish("advance_velocities");
}

void SimulationEcsim::make_step([[maybe_unused]] const int timestep) {
    RECORD_TIMER;
    globalTimer.start("Total");

    drift_particles_and_assemble_response();

    globalTimer.start("predict_electric_field");
    // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})).
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);

    globalTimer.finish("predict_electric_field");

    advance_particle_velocities();

    globalTimer.start("advance_endpoint_fields");
    // calculate fieldB
    blas::sum(2.0, fieldEp.data(), -1.0, fieldE.data(), fieldEn.data());
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);

    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn, dt_);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
    globalTimer.finish("advance_endpoint_fields");

    // Secondaries did not contribute to this step's current/response.
    // Make them visible to endpoint diagnostics and the next step.
    bc_handler.flush_species(species);

    globalTimer.finish("Total");
}

template <ShapeType ST>
void SimulationEcsim::prepare_block_matrix() {
    RECORD_TIMER;
    Array3D<int> countInCell(domain.size());
    countInCell.setZero();
    for (auto& kv : species) {
        auto& sp = *kv.second;
#pragma omp parallel for schedule(static, 32)
        for (int i = 0; i < sp.particlesData.capacity(); i++) {
            countInCell(i) += sp.particlesData(i).size();
        }
    }

    static_assert(ST == ShapeType::Linear, "prepare_block_matrix supports Linear shape only");
    mesh.LmatX2.prepare(countInCell);
}

template void SimulationEcsim::prepare_block_matrix<ShapeType::Linear>();

void SimulationEcsim::predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& B, Field3d& J) {
    const double dt = dt_;

    Operator A = mesh.IMmat + mesh.Lmat2;
    mesh.Lmat2.makeCompressed();

    Field3d rhs = E - 0.5 * dt * J + 0.5 * dt * (mesh.curlB * B);
    impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);

    // E(n+1/2) = (M-L) * E(n+1/2) + E - 0.5*dt*(J + rotB)
    solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);
    LOG_STEP("  solver error=" << (A * Ep - rhs).norm() << "\n");
}

void SimulationEcsim::predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& E_ex, const Field3d& B,
                                             Field3d& J) {
    RECORD_TIMER;

    const double dt = dt_;

    timer::flatTimer timerDestructors(timer::NoStart{});
    {
        timer::commonTimer timerRhs("make rhs");
        Field3d rhs = E + 0.5 * dt * (mesh.curlB * B - J) - mesh.Lmat2 * E_ex;
        // Единая Дирихле-маска er0 (mode "solv") в правой части
        impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);
        timerRhs.finish();

        // E(n+1/2) = (M-L) * E(n+1/2)  - L*E_ex + E - 0.5*dt*(J + rotB)
        // (M*Ex = 0)
        timer::commonTimer timerA("construct A");
        const Operator A = parallelSparseSum(mesh.IMmat, mesh.Lmat2);
        timerA.finish();

        const double err = solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);
        LOG_STEP("  solver error = " << err << "\n");

        // A и rhs уничтожаются при выходе из этого scope — замеряем их деструкторы
        timerDestructors.start("destructor operator A and rhs");
    }
}

void SimulationEcsim::init_fields() {
    fieldJp.resize(domain.size(), 3);
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
}

void SimulationEcsim::prepare_step(const int timestep) {
    RECORD_TIMER;
    const double dt = dt_;
    for (auto& kv : species) {
        auto& sp = *kv.second;
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
        timer::commonTimer timerPrepare("zero-step diagnostic");
        fieldEp.setZero();

        for (auto& kv : species) {
            auto& sp = *kv.second;
            sp.currentOnGrid.setZero();
        }
    }

    initialize_diagnostic_outputs();

    prepare_diagnostic_currents();
    diagnostic_energy(*diagnostic_ptr_);
    write_diagnostics(
        timestep, *diagnostic_ptr_,
        {domain, species, fieldEn, fieldBn, fieldE_external, fieldBInit, bc_handler, diagnostic_density_shape()},
        outputs_);
}

void SimulationEcsim::initialize_diagnostic_outputs() {
    if (!diagnostic_ptr_) {
        nlohmann::json diagnostic_config =
            system_config.contains("diagnostics") ? system_config["diagnostics"] : nlohmann::json::object();
        outputs_ = OutputFactory::create(diagnostic_config, domain, species, system_config, fieldEn, fieldBn,
                                         fieldE_external, fieldBInit);
        diagnostic_ptr_ = std::make_unique<Diagnostics>(species);
    }
}

bool SimulationEcsim::charge_continuity_requested(int timestep) {
    initialize_diagnostic_outputs();
    for (const auto& output : outputs_)
        if (output->checks_charge_continuity(timestep))
            return true;
    return false;
}

void SimulationEcsim::finalize() {
    for (auto& out : outputs_) {
        out->finalize();
    }
    Simulation::finalize();
}

EnergyFields SimulationEcsim::energy_fields() const {
    return {.E = fieldE,
            .En = fieldEn,
            .Ep = fieldEp,
            .B = fieldB,
            .Bn = fieldBn,
            .B_prescribed = fieldBInit,
            .E_external = fieldE_external,
            .Jp = fieldJp};
}

void SimulationEcsim::prepare_diagnostic_currents() {
    for (auto& kv : species) {
        auto& sp = *kv.second;
        algorithmsECSIM::calculate_current(sp, sp.currentOnGrid);
        bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
    }
}

void SimulationEcsim::diagnostic_energy(Diagnostics& diagnostic) const {
    measure_ecsim_energy(diagnostic, species, energy_fields(), mesh.Lmat2, bc_handler.active_range(domain.grid), dt_,
                         dampingEnergy_, filtered_current_power_);
}
