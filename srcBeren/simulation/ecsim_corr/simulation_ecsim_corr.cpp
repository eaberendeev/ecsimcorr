// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "simulation_ecsim_corr.h"

#include <iostream>
#include <map>
#include <string>

#include "Coil.h"
#include "Damping.h"
#include "DensityMeasurements.h"
#include "Diagnostic.h"
#include "EnergyMeasurements.h"
#include "Mesh.h"
#include "ParticlesArray.h"
#include "ParticlesDiagnostic.h"
#include "Read.h"
#include "algorithms_ecsim_corr.h"
#include "collision.h"
#include "collisions_with_neutrals.h"
#include "containers.h"
#include "cross_section.h"
#include "grid/Domain.h"
#include "log_macros.h"
#include "recovery.h"
#include "timer.h"

SimulationEcsimCorr::SecondHalfStepResults SimulationEcsimCorr::advance_velocities_and_deposit_second_half_current() {
    RECORD_TIMER;
    const double dt = dt_;

    globalTimer.start("advance_velocities_and_second_half_current");

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

    globalTimer.finish("advance_velocities_and_second_half_current");
    return second_half_step;
}

// First half: x^n -> x^{n+1/2}; v^n is unchanged except at walls.
void SimulationEcsimCorr::drift_half_step_and_assemble_response() {
    RECORD_TIMER;
    const double dt = dt_;

    globalTimer.start("first_half_drift_and_response");
    fieldBFull = fieldB + fieldBInit;

    for (auto& kv : charged_species) {
        auto& sp = kv.second.get();
        ParticleBoundaryEvents events;
        move_and_calc_current<SHAPE_CH>(sp, 0.5 * dt, sp.currentOnGrid, bc_handler.path_policy(sp, domain, &events));
        // Apply the x^n -> x^{n+1/2} decisions before redistribution changes
        // particle addresses. The snapshots are not used after this call.
        bc_handler.apply_to_particles(sp, species, domain, &events);
        sp.update_cells(domain);

        // +++ get J(x_{n+1/2},v_n)_predict

        algorithmsECSIM::predict_current<SHAPE>(sp, fieldBFull, fieldJp, dt);
    }
    globalTimer.finish("first_half_drift_and_response");
    // Discrete filtered-current power <E0,Jp_raw>, measured before the mask.
    filtered_current_power_ = dot_product_sum(fieldE0, fieldJp, bc_handler.active_range(domain.grid));
    bc_handler.apply_to_fields(fieldJp, FieldType::CURRENT, domain);

    globalTimer.start("assemble_particle_response");
    assembleLmat2(dt);
    globalTimer.finish("assemble_particle_response");

    globalTimer.start("response_boundary_constraints");
    bc_handler.apply_to_operator(mesh.Lmat2, domain);
    globalTimer.finish("response_boundary_constraints");
}

// The deposited paths already fixed absorption/reflection. Rescale first,
// then apply those saved decisions before redistributing particles.
void SimulationEcsimCorr::correct_particle_energy_and_apply_boundaries(SecondHalfStepResults& second_half_step) {
    RECORD_TIMER;
    const double dt = dt_;

    globalTimer.start("correct_particle_energy_and_boundaries");

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

    globalTimer.finish("correct_particle_energy_and_boundaries");
}

void SimulationEcsimCorr::make_step(const int timestep) {
    RECORD_TIMER;

    const double dt = dt_;
    globalTimer.start("Total");
    // rho^n is sampled after this step's injection, before either half-drift.
    // No density grid is allocated on unchecked steps.
    continuity_.reset();
    const bool check_continuity = charge_continuity_requested(timestep);
    std::optional<Field3d> charge_before;
    if (check_continuity)
        charge_before = measure_charge_density(species, domain, bc_handler, diagnostic_density_shape());
    advance_neutral_particles(dt);

    drift_half_step_and_assemble_response();

    globalTimer.start("predict_electric_field");
    // --- solve A*E'_{n+1/2}=f(E_n, B_n, J(x_{n+1/2})).
    predict_electric_field(fieldEp, fieldE, fieldE_external, fieldB, fieldJp);
    bc_handler.apply_to_fields(fieldEp, FieldType::ELECTRIC, domain);

    globalTimer.finish("predict_electric_field");

    // Keep outgoing particles in the same storage until correction and
    // application of these saved decisions. Field solves do not move particles.
    auto second_half_step = advance_velocities_and_deposit_second_half_current();

    collect_current(fieldJe);

    globalTimer.start("solve_electric_field_from_current");
    solve_electric_field_from_deposited_current(fieldEn, fieldE, fieldB, fieldJe, dt);
    bc_handler.apply_to_fields(fieldEn, FieldType::ELECTRIC, domain);
    globalTimer.finish("solve_electric_field_from_current");

    correct_particle_energy_and_apply_boundaries(second_half_step);

    globalTimer.start("advance_endpoint_fields");
    // calculate fieldB
    mesh.compute_fieldB(fieldBn, fieldB, fieldE, fieldEn, dt_);
    bc_handler.apply_to_fields(fieldBn, FieldType::MAGNETIC, domain);
    globalTimer.finish("advance_endpoint_fields");

    if (check_continuity) {
        const Field3d charge_after = measure_charge_density(species, domain, bc_handler, diagnostic_density_shape());
        continuity_ = measure_charge_continuity(*charge_before, charge_after, fieldJe, mesh.divE, dt);
    }

    // Secondaries did not contribute to this step's current/response.
    // Make them visible to endpoint diagnostics and the next step.
    bc_handler.flush_species(species);

    globalTimer.finish("Total");
}

void SimulationEcsimCorr::prepare_diagnostic_currents() {
    // Corr retains the deposited path current. Recomputing point current after
    // velocity rescaling/collisions would change its temporal level and work.
}

void SimulationEcsimCorr::diagnostic_energy(Diagnostics& diagnostic) const {
    diagnostic.charge_continuity = continuity_;
    measure_corr_energy(diagnostic, species, energy_fields(), bc_handler.active_range(domain.grid), dt_, dampingEnergy_,
                        filtered_current_power_);
}
