// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_ecsim_base.h"

#include <iostream>
#include <source_location>

#include "Damping.h"
#include "Diagnostic.h"
#include "DiagnosticMeasurements.h"
#include "EnergyMeasurements.h"
#include "OutputFactory.h"
#include "config/field_config.h"
#include "config/json_settings.h"
#include "runtime/log_macros.h"
#include "io/recovery.h"
#include "solverSLE.h"
#include "runtime/profiler.h"

SimulationEcsimBase::SimulationEcsimBase(const nlohmann::json& system_config, const nlohmann::json& particles_config,
                                         int argc, char** argv)
    : Simulation(system_config, particles_config, argc, argv), response_row_blocks_(omp_get_max_threads()) {
    constexpr int reserve_bytes = 32 * 1024 * 1024;
    for (SimpleArrayBuffer<RowBlock<36>>& row_blocks : response_row_blocks_)
        row_blocks.reserve(reserve_bytes / sizeof(RowBlock<36>));
}

void SimulationEcsimBase::assembleLmat2(double dt) {
    Operator tmpMat;

    if (assembly_check_counter_ % envOptions::validationPeriodicity() == 0) {
        timer::commonTimer timerCopy("copy matrix");
        tmpMat = field_system.L;
    }
    timer::commonTimer timerAssemble("assemble L");
    ecsim_response::clear_row_blocks(response_row_blocks_);

    for (auto& kv : species) {
        ParticlesArray& sp = *kv.second;
        ecsim_response::accumulate_species<SHAPE>(sp, fieldBFull, domain, dt, response_row_blocks_);
    }
    ecsim_response::assemble_matrix(field_system.L, response_row_blocks_, response_workspace_);
    timerAssemble.finish();

    if (assembly_check_counter_ % envOptions::validationPeriodicity() == 0) {
        timer::commonTimer timerRefAssemble("reference assemble L");
        prepare_block_matrix<SHAPE>(reference_response_blocks_);

        for (auto& kv : species) {
            ParticlesArray& sp = *kv.second;
            ecsim_response::accumulate_species_reference<SHAPE>(sp, fieldBFull, domain, dt, reference_response_blocks_);
        }
        ecsim_response::assemble_reference(tmpMat, domain, reference_response_blocks_);
        timerRefAssemble.finish();

        const bool isSameShape = checkMatrixPortraitCoincidence(field_system.L, tmpMat);
        if (!isSameShape) {
            const std::source_location location = std::source_location::current();
            std::cerr << location.file_name() << ":" << location.line()
                      << " Warning: optimized and reference assembly produced different matrix portraits" << std::endl;
        }
        timer::commonTimer timerNorms("compute norms");
        const double diffNorm = (field_system.L - tmpMat).norm();
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

template <ShapeType ST>
void SimulationEcsimBase::prepare_block_matrix(BlockMatrix& reference_blocks) {
    RECORD_TIMER;
    Array3D<int> countInCell(domain.size());
    countInCell.setZero();
    for (auto& kv : species) {
        auto& sp = *kv.second;
#pragma omp parallel for schedule(static, 32)
        for (int i = 0; i < sp.cells().capacity(); i++) {
            countInCell(i) += sp(i).size();
        }
    }

    static_assert(ST == ShapeType::Linear, "prepare_block_matrix supports Linear shape only");
    reference_blocks.resize(domain.total_size());
    reference_blocks.prepare(countInCell);
}

template void SimulationEcsimBase::prepare_block_matrix<ShapeType::Linear>(BlockMatrix&);

void SimulationEcsimBase::predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& B, Field3d& J) {
    RECORD_TIMER;
    const double dt = dt_;

    Operator A = field_system.IM + field_system.L;
    field_system.L.makeCompressed();

    Field3d rhs = E - 0.5 * dt * J + 0.5 * dt * (maxwell.curlB * B);
    impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);

    // (I-M+L) Ep = E + dt/2 (curlB B - J).
    solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);
    LOG_STEP("  solver error=" << (A * Ep - rhs).norm() << "\n");
}

void SimulationEcsimBase::predict_electric_field(Field3d& Ep, const Field3d& E, const Field3d& E_ex, const Field3d& B,
                                                 Field3d& J) {
    RECORD_TIMER;

    const double dt = dt_;

    timer::flatTimer timerDestructors(timer::NoStart{});
    {
        timer::commonTimer timerRhs("make rhs");
        Field3d rhs = E + 0.5 * dt * (maxwell.curlB * B - J) - field_system.L * E_ex;
        // Единая Дирихле-маска er0 (mode "solv") в правой части
        impose_dirichlet_rhs(rhs, bc_handler.electric_dirichlet(domain), domain);
        timerRhs.finish();

        // (I-M+L) Ep = E + dt/2 (curlB B - J) - L E_external.
        timer::commonTimer timerA("construct A");
        const Operator A = parallelSparseSum(field_system.IM, field_system.L);
        timerA.finish();

        const double err = solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, Ep, E);
        LOG_STEP("  solver error = " << err << "\n");

        // A и rhs уничтожаются при выходе из этого scope — замеряем их деструкторы
        timerDestructors.start("destructor operator A and rhs");
    }
}

void SimulationEcsimBase::init_fields() {
    RECORD_TIMER;
    damping_parameters_ = damping_parameters_from_config(system_config);
    {
        MaxwellStencils stencils = assemble_maxwell_stencils(domain, bc_handler);
        field_system.init(stencils, domain, dt_, bc_handler);
        maxwell.init(std::move(stencils));
    }
    fieldJp.resize(domain.size(), 3);

    fieldE.resize(domain.size(), 3);
    fieldEn.resize(domain.size(), 3);
    fieldEp.resize(domain.size(), 3);
    fieldB.resize(domain.size(), 3);
    fieldBn.resize(domain.size(), 3);
    fieldBInit.resize(domain.size(), 3);
    fieldBFull.resize(domain.size(), 3);
    fieldE_external.resize(domain.size(), 3);

    fieldJp.setZero();

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

void SimulationEcsimBase::prepare_step(const int timestep) {
    RECORD_TIMER;
    const double dt = dt_;
    for (auto& kv : species) {
        ParticlesArray& sp = *kv.second;
        sp.begin_step();
        sp.step_state.injection_energy =
            inject_particles_step(sp, particle_sources.at(sp.name()).injection, timestep, domain, dt,
                                  system_config.value("co_locate_species", true));
    }

    dampingEnergy_ = damping_fields(fieldEn, fieldBn, domain, damping_parameters_);
    fieldE = fieldEn;
    fieldB = fieldBn;
    fieldJp.setZero();
}

SimulationEcsimBase::~SimulationEcsimBase() = default;

void SimulationEcsimBase::make_diagnostic(const int timestep) {
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

void SimulationEcsimBase::initialize_diagnostic_outputs() {
    RECORD_TIMER;
    if (!diagnostic_ptr_) {
        nlohmann::json diagnostic_config =
            system_config.contains("diagnostics") ? system_config["diagnostics"] : nlohmann::json::object();
        outputs_ = OutputFactory::create(diagnostic_config, domain, species, system_config, fieldEn, fieldBn,
                                         fieldE_external, fieldBInit);
        diagnostic_ptr_ = std::make_unique<Diagnostics>(species);
    }
}

bool SimulationEcsimBase::charge_continuity_requested(int timestep) {
    initialize_diagnostic_outputs();
    for (const auto& output : outputs_)
        if (output->checks_charge_continuity(timestep))
            return true;
    return false;
}

void SimulationEcsimBase::finalize() {
    for (auto& out : outputs_) {
        out->finalize();
    }
    Simulation::finalize();
}

EnergyFields SimulationEcsimBase::energy_fields() const {
    return {.E = fieldE,
            .En = fieldEn,
            .Ep = fieldEp,
            .B = fieldB,
            .Bn = fieldBn,
            .B_prescribed = fieldBInit,
            .E_external = fieldE_external,
            .Jp = fieldJp};
}
