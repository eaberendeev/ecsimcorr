#pragma once

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "drift_kinetic_adaptive.h"
#include "drift_kinetic_moments.h"
#include "drift_kinetic_particles.h"
#include "drift_kinetic_push.h"
#include "drift_kinetic_switch.h"
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
//
// A species with `model: "adaptive"` is held in both descriptions: the kinetic
// markers stay in the species `<name>`, the guiding centres form the species
// `<name>_gc`. Before a step every marker is moved to the model fitting its
// magnetization (`switch_models`); see drift_kinetic_switch.h for the rules.
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

    double extra_lost_energy() const override {
        return lost_energy_;
    }

    // sum m w v_par b(r) with B^{n+1}; the drift momentum is omitted, like the drift energy.
    Vector3R extra_momentum() const override;
    void add_extra_diagnostics(Diagnostics& diag) override;

    const MomentSources* extra_moment_sources() const override {
        return &moment_sources_;
    }

    // Guiding-centre energy sum w (m v_par^2 / 2 + mu |B(r)|) for the full field `magnetic`.
    double guiding_center_energy(const Field3d& magnetic) const;
    driftKinetic::DriftKineticParticles& guiding_centers(const std::string& name);
    const driftKinetic::PushStatistics& push_statistics(const std::string& name) const;

    // Model change of the markers of the adaptive species at the time level n
    // (fields E^n, B^n = fieldBn). With w = omega_c dt at the marker and rho
    // its Larmor radius in units of the smallest cell size:
    //   kinetic -> guiding centre: w > omega_dt_to_guiding_center and rho < larmor_radius_to_guiding_center,
    //   guiding centre -> kinetic: w < omega_dt_to_kinetic or rho > larmor_radius_to_kinetic.
    // The gaps between the thresholds keep a marker from flipping every step.
    // Charge density and energy are unchanged. Called from prepare_step before
    // the energies of the step are recorded.
    void switch_models(int timestep);

    struct SwitchStatistics {
        long to_guiding_center = 0, to_kinetic = 0;   // markers moved by the last call
        double energy_defect = 0.0;                   // energy after - before the last call (rounding)
    };

    // Statistics of the last switch for the guiding-centre part `<name>_gc`.
    const SwitchStatistics& switch_statistics(const std::string& guiding_name) const {
        return switching_.at(guiding_name).statistics;
    }

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
    // Open walls: the orbit of a guiding centre that leaves the box is cut at
    // the wall, at the step fraction theta found by bisection.
    bool outside_box(const Vector3R& position) const;
    template <typename FieldFn>
    bool cut_at_wall(const driftKinetic::GuidingCenter& start, const FieldFn& fields,
                     std::vector<driftKinetic::Leaf>& leaves, double& theta) const;
    bool walls_ = false;   // some axis is nonperiodic (open faces)

    // Adaptive species, by the name of the guiding-centre part.
    struct SwitchingSpecies {
        std::string kinetic;   // name of the kinetic part in `species`
        SwitchStatistics statistics;
    };

    std::map<std::string, SwitchingSpecies> switching_;
    double switch_omega_to_guiding_ = 1.0, switch_omega_to_kinetic_ = 0.5;     // thresholds of omega_c dt
    double switch_radius_to_guiding_ = 1.0, switch_radius_to_kinetic_ = 1.5;   // thresholds of rho / h
    int switch_interval_ = 1;                                                  // steps between model changes

    // Guiding centres absorbed by the walls in the last evaluation: energy at
    // the wall hit (< 0: kept) and the face, per cell and marker.
    struct Absorption {
        double energy = -1.0;
        int face = -1;
    };

    std::map<std::string, std::vector<std::vector<Absorption>>> absorbed_;
    double lost_energy_ = 0.0;                                                     // of the accepted step
    std::map<std::string, std::array<std::pair<double, long>, 6>> lost_by_face_;   // energy, count
    // Recovery file of a guiding-centre species under `root` (".", or ".." on restart).
    static std::string recovery_file(const std::string& root, const std::string& name);
    int recovery_interval_ = 0;   // of the `recovery` output; 0: none

    std::unique_ptr<driftKinetic::YeeCoupling> coupling_;   // created in init(), after the boundaries
    Species particles_, initial_particles_, trial_particles_;
    std::map<std::string, driftKinetic::PushStatistics> push_stats_;
    driftKinetic::GuidingCenterPush push_;
    driftKinetic::AdaptiveSubstepping adaptive_;
    Field3d magnetic_begin_, magnetic_end_;                        // full B^n and B^{n+1} of the current iterate
    std::vector<Field3d> thread_current_, thread_magnetization_;   // per-thread deposits of one evaluation
    // Per-species sources of the last field evaluation; the current
    // J_gc + curl_B M of the accepted step feeds the diagnostics.
    std::map<std::string, Field3d> species_free_current_, species_magnetization_, species_current_;
    std::map<std::string, std::unique_ptr<driftKinetic::GuidingCenterMoments>> moment_adapters_;
    MomentSources moment_sources_;
    double previous_energy_ = 0.0;
    bool use_aitken_ = false;
    bool coord_is_gc_ = false;
    bool evolve_fields_ = true;
    bool fail_on_terminal_nonconvergence_ = false;
    bool fail_on_final_nonconvergence_ = false;
    bool warned_this_step_ = false;
};
