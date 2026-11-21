#pragma once

#include <memory>

#include "quadratic_point_pattern.h"
#include "simulation/ecsim/simulation_ecsim.h"
#include "whitney_shape.h"

// ECSIM time staggering and Maxwell solve with compatible spatial forms.
// No orbit iteration or Gauss correction is introduced.
class SimulationEcsimWhitney : public SimulationEcsim {
   public:
    using SimulationEcsim::SimulationEcsim;
    void init() override;
    ShapeType diagnostic_density_shape() const override {
        return kernel_ == ecsimWhitney::KernelOrder::Quadratic ? ShapeType::Quadratic : ShapeType::Linear;
    }
    void drift_particles_and_assemble_response() override;
    void advance_particle_velocities() override;

   protected:
    // Assemble Jp and the orbit-averaged response L at the current particle
    // positions (stored in field_system.L). The caller owns the drift: full
    // step for ECSIM, first half-step for ECSIM-Corr.
    void assemble_point_response();
    // Direct CSR scatter for the fully periodic quadratic point response.
    void assemble_point_response_pattern();
    // field_system.L on the field unknowns from a response on storage indices.
    // sum_images adds rows and columns of periodic images to their owners, as
    // the ECSIM deposit does; otherwise only owner rows are read, because they
    // already hold the whole response. Wall rows are constrained afterwards.
    void set_storage_response(const Operator& storage, bool sum_images);
    // Non-periodic axes require a wall-constrained Maxwell response. Point
    // response drops crossing particles before assembly; the frozen ballistic
    // path instead clips their orbits and keeps their truncated current.
    bool wall_boundaries_ = false;
    virtual bool supports_wall_boundaries() const {
        return true;
    }
    bool parallel_point_ = false;   // sequential reference remains available
    ecsimWhitney::KernelOrder kernel_ = ecsimWhitney::KernelOrder::Linear;
    // Built only for the fully periodic quadratic point response.
    std::unique_ptr<ecsimWhitney::QuadraticPointPattern> pattern_;
    // Dense-pattern storage; field_system.L receives a pruned copy per step.
    Operator pattern_matrix_;
    // currentOnGrid already contains q w W (v^n+v^{n+1})/2 from the push.
    void prepare_diagnostic_currents() override {
    }
};
