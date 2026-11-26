// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_ecsim_corr.h"

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "algorithms_ecsim_corr.h"
#include "coupling/current.h"
#include "diagnostics/DensityMeasurements.h"
#include "diagnostics/Diagnostic.h"
#include "diagnostics/EnergyMeasurements.h"
#include "fields/Coil.h"
#include "fields/Damping.h"
#include "fields/solverSLE.h"
#include "grid/domain.h"
#include "io/recovery.h"
#include "io/text_input.h"
#include "particles/particle_array.h"
#include "runtime/log_macros.h"
#include "runtime/profiler.h"
#include "simulation/ecsim/algorithms_ecsim.h"
#include "utils/containers.h"

void SimulationEcsimCorr::init_fields() {
    SimulationEcsimBase::init_fields();
    corrector_passes_ = system_config.value("CorrectorPasses", 0);
    if (corrector_passes_ < 0)
        throw std::invalid_argument("CorrectorPasses must not be negative");
    fieldJe.resize(domain.size(), 3);
    fieldJe.setZero();
    if (corrector_passes_ > 0)
        trial_current_.resize(domain.size(), 3);
}

void SimulationEcsimCorr::prepare_step(int timestep) {
    SimulationEcsimBase::prepare_step(timestep);
    fieldJe.setZero();
}

SimulationEcsimCorr::SecondHalfStepResults SimulationEcsimCorr::advance_velocities_and_deposit_second_half_current() {
    RECORD_TIMER;
    const double dt = dt_;

    fieldBFull = fieldB + fieldBInit;
    Field3d fieldE_full = fieldEp + fieldE_external;

    SecondHalfStepResults second_half_step;
    for (auto& kv : charged_species) {
        auto& sp = kv.second.get();
        auto& result = second_half_step[sp.name()];
        fused_push_and_deposit<SHAPE_CH>(sp, fieldE_full, fieldBFull, dt, 0.5 * dt, sp.currentOnGrid,
                                         result.predictor_power,
                                         bc_handler.path_policy(sp, domain, &result.boundary_events));
        sp.currentOnGrid *= 0.5;
        bc_handler.apply_to_fields(sp.currentOnGrid, FieldType::CURRENT, domain);
    }

    return second_half_step;
}

// A corrector pass deposits the second half-paths in E' like
// advance_velocities_and_deposit_second_half_current, but leaves the particles
// at x^{n+1/2}, v^n and their first half-path current in currentOnGrid.
// fieldJe receives the averaged current of both half-paths.
void SimulationEcsimCorr::deposit_trial_current() {
    RECORD_TIMER;
    fieldBFull = fieldB + fieldBInit;
    const Field3d fieldE_full = fieldEp + fieldE_external;

    fieldJe.setZero();
    for (const std::pair<const std::string, std::reference_wrapper<ParticlesArray>>& kv : charged_species) {
        const ParticlesArray& sp = kv.second.get();
        trial_current_ = sp.currentOnGrid;
        trial_push_and_deposit<SHAPE_CH>(sp, fieldE_full, fieldBFull, dt_, 0.5 * dt_, trial_current_,
                                         bc_handler.path_policy(sp, domain));
        trial_current_ *= 0.5;
        bc_handler.apply_to_fields(trial_current_, FieldType::CURRENT, domain);
        fieldJe += trial_current_;
    }
}

// Particles are pushed by the predicted field E', while E^{n+1} follows from
// their Esirkepov current, so E' differs from (E^n+E^{n+1})/2. One quasi-Newton
// step towards equality, with the ECSIM response L standing for the derivative
// of that current: (I-M+L) dE = (I-M) ((E^n+E^{n+1})/2 - E'). E^{n+1} itself is
// not needed: (I-M) E^{n+1} = b on the field unknowns, so the right-hand side is
// (I-M) (E^n/2 - E') + b/2.
void SimulationEcsimCorr::refine_predicted_field(const Operator& A) {
    RECORD_TIMER;
    const Field3d endpoint = field_system.restrict_rhs(endpoint_rhs(fieldE, fieldB, fieldJe, dt_));
    const Field3d rhs = field_system.IM * field_system.restrict(0.5 * fieldE - fieldEp) + 0.5 * endpoint;
    if (rhs.norm() == 0.0)
        return;
    Field3d solution(field_system.unknowns());
    solution.setZero();
    const double err = solve_linear_system<BicgstabSolver<Field3d>>(A, rhs, solution, solution);
    LOG_STEP("  corrector pass: |(I-M)(E_mid-E')|=" << rhs.norm() << " solver error=" << err << "\n");
    Field3d correction(domain.size(), 3);
    field_system.prolong(solution, correction);
    fieldEp = fieldEp + correction;
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);
}

// First half: x^n -> x^{n+1/2}; v^n is unchanged except at walls.
void SimulationEcsimCorr::drift_half_step_and_assemble_response() {
    RECORD_TIMER;
    const double dt = dt_;

    fieldBFull = fieldB + fieldBInit;

    for (auto& kv : charged_species) {
        auto& sp = kv.second.get();
        ParticleBoundaryEvents events;
        coupling::drift_and_deposit_current<SHAPE_CH>(sp, 0.5 * dt, sp.currentOnGrid,
                                                      bc_handler.path_policy(sp, domain, &events));
        // Apply the x^n -> x^{n+1/2} decisions before redistribution changes
        // particle addresses. The snapshots are not used after this call.
        bc_handler.apply_to_particles(sp, species, domain, &events);
        sp.update_cells(domain);

        // +++ get J(x_{n+1/2},v_n)_predict

        algorithmsECSIM::predict_current<SHAPE>(sp, fieldBFull, fieldJp, dt);
    }
    // Discrete filtered-current power <E0,Jp_raw>, measured before the mask.
    filtered_current_power_ = dot_product_sum(fieldE0, fieldJp, bc_handler.active_range(domain.grid));
    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    assembleLmat2(dt);

    bc_handler.apply_to_operator(field_system.L, domain, field_system.storage_index);
}

// The deposited paths already fixed absorption/reflection. Rescale first,
// then apply those saved decisions before redistributing particles.
void SimulationEcsimCorr::correct_particle_energy_and_apply_boundaries(SecondHalfStepResults& second_half_step) {
    RECORD_TIMER;
    const double dt = dt_;

    for (auto& kv : charged_species) {
        auto& sp = kv.second.get();
        auto& result = second_half_step.at(sp.name());
        // Only velocities may change here. Apply the decisions made at the
        // crossing velocity before discarding the address snapshots.
        rescale_velocities_for_energy_balance(sp, dt, result.predictor_power);
        bc_handler.apply_to_particles(sp, species, domain, &result.boundary_events);
        second_half_step.erase(sp.name());
        sp.update_cells(domain);
    }
}

void SimulationEcsimCorr::make_step(const int timestep) {
    RECORD_TIMER_PARAMS(timestep);

    const double dt = dt_;
    // rho^n is sampled after this step's injection, before either half-drift.
    // No density grid is allocated on unchecked steps.
    continuity_.reset();
    const bool check_continuity = charge_continuity_requested(timestep);
    std::optional<Field3d> charge_before;
    if (check_continuity)
        charge_before = measure_charge_density(species, domain, bc_handler, diagnostic_density_shape());
    advance_neutral_particles(dt);

    drift_half_step_and_assemble_response();

    {
        // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})); A lives through the passes.
        const Operator A = predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
        bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);
        for (int pass = 0; pass < corrector_passes_; ++pass) {
            deposit_trial_current();
            refine_predicted_field(A);
        }
    }

    // Keep outgoing particles in the same storage until correction and
    // application of these saved decisions. Field solves do not move particles.
    SecondHalfStepResults second_half_step = advance_velocities_and_deposit_second_half_current();
    collect_current(fieldJe);
    solve_electric_field_from_deposited_current(fieldEn, fieldE, fieldB, fieldJe, dt_);
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);

    correct_particle_energy_and_apply_boundaries(second_half_step);

    // calculate fieldB
    maxwell.advance_magnetic_field(fieldBn, fieldB, fieldE, fieldEn, dt_);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);

    if (check_continuity) {
        const Field3d charge_after = measure_charge_density(species, domain, bc_handler, diagnostic_density_shape());
        continuity_ = measure_charge_continuity(*charge_before, charge_after, fieldJe, maxwell.divE, dt);
    }

    // Secondaries did not contribute to this step's current/response.
    // Make them visible to endpoint diagnostics and the next step.
    bc_handler.flush_species(species);
}

void SimulationEcsimCorr::prepare_diagnostic_currents() {
    // Corr retains the deposited path current. Recomputing point current after
    // velocity rescaling/collisions would change its temporal level and work.
}

void SimulationEcsimCorr::diagnostic_values(Diagnostics& diagnostic) const {
    diagnostic.charge_continuity = continuity_;
    measure_corr_energy(diagnostic, species, energy_fields(), bc_handler.active_range(domain.grid), dt_, dampingEnergy_,
                        filtered_current_power_);
}
