#pragma once

#include <array>
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

#include "drift_kinetic_adaptive.h"
#include "drift_kinetic_grid.h"
#include "drift_kinetic_particles.h"
#include "drift_kinetic_push.h"
#include "simulation/eccapfim/simulation_eccapfim.h"

// Self-consistent drift-kinetic/hybrid scheme on the Yee grid.
//
// The electric field is the ECCAPFIM midpoint field E^{n+1/2}; the magnetic
// field is the ECCAPFIM self-consistent B^{n+1/2} plus the external B.  The
// nonlinear source is the kinetic current plus J_gc + curl(M), with J_gc and
// M deposited by the guiding-centre particles.  `Scheme: "drift_kinetic"`
// converts every charged species; `Scheme: "hybrid_kinetic"` converts only
// species whose particle config has `model: "drift_kinetic"` and keeps the
// rest on the kinetic ECCAPFIM backend.
class SimulationDriftKinetic : public SimulationEccapfim {
   public:
    using SimulationEccapfim::SimulationEccapfim;
    void init() override;
    void init_fields() override;
    void prepare_step(int timestep) override;
    void make_diagnostic(int timestep) override;
    double extra_kinetic_energy() const override;
    double extra_previous_kinetic_energy() const override;
    void add_extra_diagnostics(Diagnostics& diag) override;
    void finalize() override;

   protected:
    void begin_nonlinear_step() override;
    void evaluate_current(const Field3d& electric, const Field3d& magnetic) override;
    void accept_nonlinear_step() override;
    void assemble_mass_matrix(std::vector<Eigen::Triplet<double>>& triplets) override;

   private:
    void convert_initial_particles();
    void sync_grid_from_fields(const Field3d& electric, const Field3d& magnetic);
    void add_cell_field_to_yee(const driftKinetic::VectorField& cell_field, Field3d& yee_field);
    Vector3R drift_velocity(const driftKinetic::GuidingCenter& state, const driftKinetic::Fields& fields) const;

    std::unique_ptr<driftKinetic::DriftKineticGrid> grid_;
    std::unordered_map<std::string, std::unique_ptr<driftKinetic::DriftKineticParticles>> particles_;
    std::unordered_map<std::string, std::unique_ptr<driftKinetic::DriftKineticParticles>> initial_particles_;
    std::unordered_map<std::string, std::unique_ptr<driftKinetic::DriftKineticParticles>> trial_particles_;
    driftKinetic::GuidingCenterPush push_;
    driftKinetic::AdaptiveSubstepping adaptive_;
    double previous_energy_ = 0.0;
    double current_energy() const;
};
