#pragma once

#include <array>
#include <memory>
#include <string>
#include <unordered_map>

#include "drift_kinetic_grid.h"
#include "drift_kinetic_particles.h"
#include "drift_kinetic_push.h"
#include "simulation/eccapfim/simulation_eccapfim.h"

// Self-consistent drift-kinetic scheme on the Yee grid.
//
// The electric field is the ECCAPFIM midpoint field E^{n+1/2}; the magnetic
// field is the ECCAPFIM self-consistent B^{n+1/2} plus the external B.  The
// nonlinear source is J_gc + curl(M), with J_gc and M deposited by the
// guiding-centre particles.  Yee and cell-centred fields are bridged by
// averaging/scattering in this file; the guiding-centre model itself remains
// cell-centred.
class SimulationDriftKinetic : public SimulationEccapfim {
   public:
    using SimulationEccapfim::SimulationEccapfim;
    void init() override;
    void init_fields() override;
    void make_diagnostic(int timestep) override;
    void finalize() override;

   protected:
    void begin_nonlinear_step() override;
    void evaluate_current(const Field3d& electric, const Field3d& magnetic) override;
    void accept_nonlinear_step() override;

   private:
    void convert_initial_particles();
    void sync_grid_from_fields(const Field3d& electric, const Field3d& magnetic);
    void deposit_cell_field_to_yee(const driftKinetic::VectorField& cell_field, Field3d& yee_field);
    Vector3R drift_velocity(const driftKinetic::GuidingCenter& state, const driftKinetic::Fields& fields) const;

    std::unique_ptr<driftKinetic::DriftKineticGrid> grid_;
    std::unordered_map<std::string, std::unique_ptr<driftKinetic::DriftKineticParticles>> particles_;
    std::unordered_map<std::string, std::unique_ptr<driftKinetic::DriftKineticParticles>> initial_particles_;
    std::unordered_map<std::string, std::unique_ptr<driftKinetic::DriftKineticParticles>> trial_particles_;
    driftKinetic::GuidingCenterPush push_;
};
