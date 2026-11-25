// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#ifndef SIMULATION_ECSIM_CORR_H
#define SIMULATION_ECSIM_CORR_H

#include <map>
#include <string>

#include "grid/domain.h"
#include "particles/particle_array.h"
#include "simulation/ecsim/simulation_ecsim_base.h"
#include "simulation/simulation.h"
#include "utils/containers.h"

// Corr uses the ECSIM predictor, but owns both half-paths and the correction.
class SimulationEcsimCorr : public SimulationEcsimBase {
   public:
    using SimulationEcsimBase::SimulationEcsimBase;
    void init_fields() override;
    void prepare_step(int timestep) override;
    // Current averaged over both half-paths; reset before each Corr step.
    Field3d fieldJe;
    void make_step(const int timestep) override;

    ShapeType diagnostic_density_shape() const override {
        return SHAPE_CH;
    }

    void rescale_velocities_for_energy_balance(ParticlesArray& sort, const double dt, double predictor_power);
    void diagnostic_values(Diagnostics& diagnostic) const override;
    // Storage right-hand side b of (I-M) E^{n+1} = b for the deposited current J.
    Field3d endpoint_rhs(const Field3d& E, const Field3d& B, const Field3d& J, const double dt);
    void solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J,
                                                     const double dt);

   private:
    std::optional<ChargeContinuityResidual> continuity_;
    void prepare_diagnostic_currents() override;

    // One species' Boris power and wall decisions belong to the same half-step.
    // Particle addresses in boundary_events remain valid until boundary application;
    // do not reorder, inject or redistribute particles between these two stages.
    struct SpeciesSecondHalfStep {
        double predictor_power = 0;   // sum(q*mpw*v_mid . E_predictor), in code energy/time units
        ParticleBoundaryEvents boundary_events;
    };

    using SecondHalfStepResults = std::map<std::string, SpeciesSecondHalfStep>;

    // The state every corrector pass starts from: v^n of the particles in cell
    // order and their first half-path current. The second half only moves a
    // particle by dt/2 v^{n+1} and leaves it in the same cell storage, so
    // x^{n+1/2} = x^{n+1} - dt/2 v^{n+1} need not be stored.
    struct MidpointState {
        std::vector<int64_t> offsets;   // first particle of each cell in velocities
        std::vector<Vector3R> velocities;
        Field3d current;
    };

    std::map<std::string, MidpointState> midpoint_;

    // "CorrectorPasses": extra velocity pushes per step, each with the predicted
    // field moved towards (E^n+E^{n+1})/2. Zero keeps the single-pass scheme.
    int corrector_passes_ = 0;

    void drift_half_step_and_assemble_response();
    void save_midpoint_states();
    void restore_midpoint_states();
    void refine_predicted_field();
    SecondHalfStepResults advance_velocities_and_deposit_second_half_current();
    void correct_particle_energy_and_apply_boundaries(SecondHalfStepResults& second_half_step);
};

#endif
