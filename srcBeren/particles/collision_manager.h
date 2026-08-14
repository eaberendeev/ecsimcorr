#pragma once

#ifndef COLLISION_MANAGER_H
#define COLLISION_MANAGER_H

#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "random.h"
#include "sgs.h"

class CollisionOperator {
   public:
    virtual ~CollisionOperator() = default;
    virtual void apply(Species &species, const Domain &domain, double dt) = 0;
    virtual std::string info() const = 0;
};

class CoulombCollisionOperator : public CollisionOperator {
   public:
    CoulombCollisionOperator(const std::string &species1, const std::string &species2, double n0,
                             double coulomb_log = 15.0);
    void apply(Species &species, const Domain &domain, double dt) override;
    std::string info() const override;

    void collide_same_type(ParticlesArray &sp, double dt);
    void collide_diff_type(ParticlesArray &sp1, ParticlesArray &sp2, double dt);

   private:
    void bin_collide(Vector3R &v1, Vector3R &v2, double q1, double q2, double n1, double n2, double m1, double m2,
                     double dt, double variance_factor, LehmerEngine &rndGen);
    double get_variance(double u, double q1, double q2, double n, double m, double dt);

    std::string species1_name_;
    std::string species2_name_;
    bool is_same_type_;
    double n0_;
    double coulomb_log_;
    LehmerEngine baseEng_;
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
    NeutralCollisionOperator(const NeutralCollisionConfig &config, double n0);
    void apply(Species &species, const Domain &domain, double dt) override;
    std::string info() const override;

   private:
    NeutralCollisionConfig config_;
    double n0_;
    LehmerEngine baseEng_;
};

class CollisionManager {
   public:
    CollisionManager() = default;
    void init_from_json(const nlohmann::json &config, double n0);
    void apply(Species &species, const Domain &domain, double dt);
    bool has_collisions() const {
        return !operators_.empty();
    }

   private:
    std::vector<std::unique_ptr<CollisionOperator>> operators_;
};

#endif
