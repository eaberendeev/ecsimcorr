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
    void init() override;
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
    void solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J,
                                                     const double dt);

   private:
    std::optional<ChargeContinuityResidual> continuity_;
    // "EnergyCorrection": "global" scales each species by one factor; "local"
    // gives every particle the corrected work on its own path (local_correction.h)
    // and scales only the remainder.
    bool local_correction_ = false;

    struct LocalCorrectionStatistics {
        double target = 0;
        double local = 0;
        double remainder = 0;
        int64_t pooled = 0;
    } local_statistics_;

    // `corrected` is (E^n+E^{n+1})/2 + E_ext, `predictor` E' + E_ext.
    void correct_energy_locally(ParticlesArray& sort, double dt, double predictor_power, const Field3d& corrected,
                                const Field3d& predictor);

    // Particles absorbed in the first half-drift (x^{n+1/2} beyond the wall,
    // v^n): the corrected field's work on their path leaves with them.
    struct AbsorbedPath {
        Particle particle;
        Face face;
    };

    std::map<std::string, std::vector<AbsorbedPath>> absorbed_first_half_;
    void prepare_diagnostic_currents() override;

    // One species' Boris power and wall decisions belong to the same half-step.
    // Particle addresses in boundary_events remain valid until boundary application;
    // do not reorder, inject or redistribute particles between these two stages.
    struct SpeciesSecondHalfStep {
        double predictor_power = 0;   // sum(q*mpw*v_mid . E_predictor), in code energy/time units
        ParticleBoundaryEvents boundary_events;
    };

    using SecondHalfStepResults = std::map<std::string, SpeciesSecondHalfStep>;

    // Particles at x^{n+1/2} with v^n and their first half-path current: the
    // state every corrector pass starts from.
    struct MidpointState {
        Array3D<std::vector<Particle>> cells;
        Field3d current;
    };

    using MidpointStates = std::map<std::string, MidpointState>;

    // "CorrectorPasses": extra velocity pushes per step, each with the predicted
    // field moved towards (E^n+E^{n+1})/2. Zero keeps the single-pass scheme.
    int corrector_passes_ = 0;

    void drift_half_step_and_assemble_response();
    MidpointStates save_midpoint_states() const;
    void restore_midpoint_states(const MidpointStates& states);
    void refine_predicted_field();
    SecondHalfStepResults push_velocities_and_solve_endpoint_field();
    SecondHalfStepResults advance_velocities_and_deposit_second_half_current();
    void correct_particle_energy_and_apply_boundaries(SecondHalfStepResults& second_half_step);
};

#endif
