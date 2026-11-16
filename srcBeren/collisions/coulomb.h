#pragma once

#include <cstdint>
#include <random>

#include "collisions/collision_manager.h"
#include "runtime/random.h"

// Per-thread counters collected only when COLLISION_STATS=1.
struct CoulombCollisionStats {
    int64_t pairs = 0;
    int64_t variance_ge1 = 0;
    int64_t modu_small = 0;
    double sum_var = 0.0;
    double max_var = 0.0;
    int64_t max_pairs_per_cell = 0;
    int64_t cells_skipped = 0;
};

class CoulombCollisionOperator : public CollisionOperator {
   public:
    // energy_correction: largest fraction of a pair's centre-of-mass energy changed when
    // unequal species weights restore exact energy conservation (Angus et al. 2024).
    CoulombCollisionOperator(const std::string& species1, const std::string& species2, double n0,
                             double coulomb_log = 15.0, double energy_correction = 0.05);
    void apply(Species& species, double dt) override;
    std::string info() const override;

    void collide_same_type(ParticlesArray& sp, double dt);
    void collide_diff_type(ParticlesArray& sp1, ParticlesArray& sp2, double dt);

   private:
    // Change of the relative velocity v1 - v2 in one binary scattering.
    Vector3R scatter(const Vector3R& v1, const Vector3R& v2, double variance_numerator, double plasma_factor,
                     double denominator_const, double variance_factor,
                     std::uniform_real_distribution<double>& uniformDistr,
                     std::normal_distribution<double>& normalDistr, LehmerEngine& rndEng, CoulombCollisionStats* stats);

    std::string species1_name_;
    std::string species2_name_;
    bool is_same_type_;
    double n0_;
    double coulomb_log_;
    double energy_correction_;
    // Deterministic per-channel seed derived from the ordered species pair.
    uint64_t stream_seed_;
    // Advanced once per collide_*() call so successive steps use fresh streams.
    uint64_t apply_index_ = 0;
};
