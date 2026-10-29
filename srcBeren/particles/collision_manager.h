#pragma once

#ifndef COLLISION_MANAGER_H
#define COLLISION_MANAGER_H

#include <cstdint>
#include <memory>
#include <nlohmann/json.hpp>
#include <random>
#include <string>
#include <vector>

#include "ParticlesArray.h"
#include "grid/Domain.h"
#include "random.h"
#include "sgs.h"

class CollisionOperator {
   public:
    virtual ~CollisionOperator() = default;
    virtual void apply(Species& species, const Domain& domain, double dt) = 0;
    virtual std::string info() const = 0;
};

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
    CoulombCollisionOperator(const std::string& species1, const std::string& species2, double n0,
                             double coulomb_log = 15.0);
    void apply(Species& species, const Domain& domain, double dt) override;
    std::string info() const override;

    void collide_same_type(ParticlesArray& sp, double dt);
    void collide_diff_type(ParticlesArray& sp1, ParticlesArray& sp2, double dt);

   private:
    void bin_collide(Vector3R& v1, Vector3R& v2, double variance_numerator, double plasma_factor,
                     double denominator_const, double variance_factor, double mu_over_m1, double mu_over_m2,
                     std::uniform_real_distribution<double>& uniformDistr,
                     std::normal_distribution<double>& normalDistr, LehmerEngine& rndEng, CoulombCollisionStats* stats);

    std::string species1_name_;
    std::string species2_name_;
    bool is_same_type_;
    double n0_;
    double coulomb_log_;
    // Deterministic per-channel seed derived from the ordered species pair.
    uint64_t stream_seed_;
    // Advanced once per collide_*() call so successive steps use fresh streams.
    uint64_t apply_index_ = 0;
};

struct NeutralCollisionConfig {
    std::string charged_species;
    std::string neutral_species;
    std::string electron_product;
    std::string ion_product;
    std::string scheme = "physical_only";
    bool electron_ionization = true;
    bool proton_ionization = true;
    bool proton_charge_exchange = true;
};

class NeutralCollisionOperator : public CollisionOperator {
   public:
    NeutralCollisionOperator(const NeutralCollisionConfig& config, double n0);
    void apply(Species& species, const Domain& domain, double dt) override;
    std::string info() const override;

   private:
    NeutralCollisionConfig config_;
    double n0_;
    LehmerEngine baseRndEng_;
};

class CollisionManager {
   public:
    CollisionManager() = default;
    void init_from_json(const nlohmann::json& config, double n0);
    void apply(Species& species, const Domain& domain, double dt);
    bool has_collisions() const {
        return !operators_.empty();
    }

   private:
    std::vector<std::unique_ptr<CollisionOperator>> operators_;
};

#endif
