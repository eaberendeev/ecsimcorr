#pragma once

#include "simulation_ecsim_whitney.h"
#include "whitney_orbit.h"

class SimulationEcsimWhitneyFrozen : public SimulationEcsimWhitney {
   public:
    using SimulationEcsimWhitney::SimulationEcsimWhitney;
    void init() override;
    void drift_particles_and_assemble_response() override;
    void advance_particle_velocities() override;

   protected:
    bool supports_wall_boundaries() const override {
        return predictor_ == ecsimWhitney::OrbitPredictor::Ballistic;
    }
    ecsimWhitney::OrbitPredictor predictor_ = ecsimWhitney::OrbitPredictor::Ballistic;
    bool cell_assembly_ = false;   // reference particle triplets remain the default
    struct ParticleResponse {
        std::vector<ecsimWhitney::EdgeWeight> electric;
        Vector3R magnetic;
        // Absorbed at the first wall: its truncated current/L drive this step.
        // The wall velocity update accounts for field work in loss energy, but
        // no separate push current is deposited; the particle is removed later.
        bool removed = false;
    };
    // Rebuilt after drift/rebin every step, consumed before any further rebin
    // or collisions. Indexed in species/bucket/particle traversal order; owns
    // no particle pointers. Values are invalid after the push; allocations are
    // retained and overwritten at the next assembly. Not part of restart.
    std::map<std::string, std::vector<ParticleResponse>> response_;
    bool response_ready_ = false;
    // Assemble response at the current coordinates. The caller owns the drift;
    // ECSIM-Corr calls this after its half-step, the full scheme after dt.
    void assemble_frozen_response();
    void append_response_matrix(std::vector<Triplet>& entries, const std::vector<ParticleResponse>& responses,
                                size_t begin, size_t end, double qm, double qw, std::vector<int>& local_index) const;
};
