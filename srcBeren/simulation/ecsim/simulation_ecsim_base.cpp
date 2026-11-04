// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_ecsim_base.h"

#include <iostream>

#include "Damping.h"
#include "Diagnostic.h"
#include "DiagnosticMeasurements.h"
#include "EnergyMeasurements.h"
#include "OutputFactory.h"
#include "config/field_config.h"
#include "config/json_settings.h"
#include "log_macros.h"
#include "recovery.h"
#include "solverSLE.h"
#include "timer.h"

SimulationEcsimBase::SimulationEcsimBase(const nlohmann::json& system_config, const nlohmann::json& particles_config,
                                         int argc, char** argv)
    : Simulation(system_config, particles_config, argc, argv) {
}

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
