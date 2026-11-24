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
    void solve_electric_field_from_deposited_current(Field3d& En, const Field3d& E, const Field3d& B, Field3d& J,
                                                     const double dt);

   protected:
    // One species' Boris power and wall decisions belong to the same half-step.
    // Particle addresses in boundary_events remain valid until boundary application;
    // do not reorder, inject or redistribute particles between these two stages.
    struct SpeciesSecondHalfStep {
        double predictor_power = 0;   // sum(q*mpw*v_mid . E_predictor), in code energy/time units
        ParticleBoundaryEvents boundary_events;
    };

    using SecondHalfStepResults = std::map<std::string, SpeciesSecondHalfStep>;

    // "CorrectorPasses": extra velocity pushes per step, each with the predicted
    // field moved towards (E^n+E^{n+1})/2. Zero keeps the single-pass scheme.
    int corrector_passes_ = 0;

    // Steps of the scheme that a variant replaces (ecsim_path). The defaults
    // are the ECSIM-Corr ones.
    // Predictor current of one species at x^{n+1/2} with v^n into fieldJp.
    virtual void predict_species_current(ParticlesArray& sp);
    // After the first half-step of all species (particles at x^{n+1/2}, v^n).
    virtual void after_first_half_step() {
    }
    // Current that is not deposited by the particle paths (added to the
    // predictor current and to the endpoint current).
    virtual void add_extra_current(Field3d&) const {
    }
    // v^n -> v^{n+1}, x^{n+1/2} -> x^{n+1} and the second half-path current.
    virtual void push_species_and_deposit(ParticlesArray& sp, const Field3d& electric, const Field3d& magnetic,
                                          SpeciesSecondHalfStep& result);
    // After B^{n+1}, before the energy correction and the wall decisions.
    virtual void after_magnetic_field_advance() {
    }
    // Energy correction of one species (velocities only).
    virtual void correct_species_energy(ParticlesArray& sp, double predictor_power);

   private:
    std::optional<ChargeContinuityResidual> continuity_;
    void prepare_diagnostic_currents() override;

    // Particles at x^{n+1/2} with v^n and their first half-path current: the
    // state every corrector pass starts from.
    struct MidpointState {
        Array3D<std::vector<Particle>> cells;
        Field3d current;
    };

    using MidpointStates = std::map<std::string, MidpointState>;

    void drift_half_step_and_assemble_response();
    MidpointStates save_midpoint_states() const;
    void restore_midpoint_states(const MidpointStates& states);
    void refine_predicted_field();
    SecondHalfStepResults push_velocities_and_solve_endpoint_field();
    SecondHalfStepResults advance_velocities_and_deposit_second_half_current();
    void correct_particle_energy_and_apply_boundaries(SecondHalfStepResults& second_half_step);
};

#endif
