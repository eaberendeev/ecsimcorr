#pragma once

#include <map>
#include <memory>
#include <string>

#include "drift_kinetic_adaptive.h"
#include "drift_kinetic_particles.h"
#include "drift_kinetic_push.h"
#include "drift_kinetic_yee.h"
#include "simulation/eccapfim/simulation_eccapfim.h"

// Self-consistent drift-kinetic and hybrid PIC on the Yee grid (xpic
// `drift_kinetic` and `hybrid_drift_kinetic`, without MPI and PETSc).
//
// Unknowns of a step are the midpoint fields E^{n+1/2}, B^{n+1/2}. For every
// field iterate the guiding centres are advanced implicitly with
//   E = E^{n+1/2} along the orbit, B^n at r^n, B^{n+1} = 2 B^{n+1/2} - B^n at r^{n+1},
// and deposit the free current J_gc = q w (r^{n+1} - r^n) / dt and the
// magnetization M = -mu w b. Maxwell's equations are
//   (E^{n+1/2} - E^n) / (dt/2) = curl_B (B^{n+1/2} - M) - J,
//   (B^{n+1/2} - B^n) / (dt/2) = -curl_E E^{n+1/2},
// with J = J_kinetic + J_gc. Since curl_B = curl_E^T, the discrete energy
// theorem is dK_gc = dt E . J_gc - M . (B^{n+1} - B^n), and the total
// particle + field energy is conserved to the nonlinear tolerance.
//
// `Scheme: "drift_kinetic"` converts every charged species; `"hybrid_kinetic"`
// converts the species with `model: "drift_kinetic"` and keeps the rest on the
// kinetic ECCAPFIM backend.
class SimulationDriftKinetic : public SimulationEccapfim {
   public:
    using SimulationEccapfim::SimulationEccapfim;
    void init() override;
    void init_fields() override;
    void prepare_step(int timestep) override;
    void make_step(int timestep) override;
    void make_diagnostic(int timestep) override;
    double extra_kinetic_energy() const override;
    double extra_previous_kinetic_energy() const override;
    double extra_external_work() const override;
    void add_extra_diagnostics(Diagnostics& diag) override;

    // Guiding-centre energy sum w (m v_par^2 / 2 + mu |B(r)|) for the full field `magnetic`.
    double guiding_center_energy(const Field3d& magnetic) const;
    driftKinetic::DriftKineticParticles& guiding_centers(const std::string& name);
    const driftKinetic::PushStatistics& push_statistics(const std::string& name) const;

    // Sources of the last field evaluation (the accepted one after a step).
    Field3d fieldJ_guiding;   // free guiding-centre current on E positions
    Field3d fieldM;           // magnetization on B positions

   protected:
    void begin_nonlinear_step() override;
    void evaluate_current(const Field3d& electric, const Field3d& magnetic) override;
    void accept_nonlinear_step() override;
    void assemble_mass_matrix(const Field3d& magnetic) override;

   private:
    using Species = std::map<std::string, std::unique_ptr<driftKinetic::DriftKineticParticles>>;
    void convert_initial_particles();

    std::unique_ptr<driftKinetic::YeeCoupling> coupling_;   // created in init(), after the boundaries
    Species particles_, initial_particles_, trial_particles_;
    std::map<std::string, driftKinetic::PushStatistics> push_stats_;
    driftKinetic::GuidingCenterPush push_;
    driftKinetic::AdaptiveSubstepping adaptive_;
    Field3d magnetic_begin_, magnetic_end_;   // full B^n and B^{n+1} of the current iterate
    double previous_energy_ = 0.0;
    bool use_aitken_ = false;
    bool coord_is_gc_ = false;
    bool evolve_fields_ = true;
    bool fail_on_terminal_nonconvergence_ = false;
    bool fail_on_final_nonconvergence_ = false;
    bool warned_this_step_ = false;
};
